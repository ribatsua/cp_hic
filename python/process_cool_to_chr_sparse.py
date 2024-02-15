#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@File    :   process_cool_to_chr_sparse.py
@Time    :   2024/01/04 16:05:11
@Author  :   cp 
@Version :   1.0
@Desc    :   None
'''
'''
#从高分辨cool文件获取不同分辨率稀疏矩阵按染色体保存所有细胞的稀疏矩阵文件
#输入：.cool文件所在文件夹路径 分辨率参数 输出路径
#通过直接对行号和列号除以分辨率与原始分辨率的比值，转换成pandas数据帧后 使用group by 合并行列索引相同的元素值，从而得到不同分辨率稀疏矩阵
#稀疏矩阵使用来自scipy.sparse的csr_matrix函数保存，保存类型为float32 以适应higashi 计算A/B区室部分代码
#运行结果为按染色体保存的整个数据集的稀疏矩阵文件
'''

from tqdm import tqdm
from scipy.sparse import csr_matrix
import subprocess
import pickle
import os
import math
import cooler
import numpy as np
from collections import Counter
import pandas as pd



def find_file(folder_path,data_dir):
    file_path=[]
    for root, dirs, files in os.walk(folder_path):
            for file in files:
                if file.endswith('.cool'): #添加.split("_")[-1]后通用于以cell_1_chr1.txt命名文件
                # if file.endswith('.txt'):
                    file_path.append(os.path.join(root, file))
    # sorted_paths = sorted(file_path, key=sort_key)
    cell_count = []
    for path in file_path:
        cell_type=path.split('/')[-1].split('_')[-3]
        cell_count.append(cell_type)

    result = Counter(cell_count)
    temp_list = {}
    cell_info = []
    start_cell_id = 0
    for key in result.keys():
         temp_list[key]=[n for n in file_path if n.split('/')[-1].split('_')[-3] == key ]
         cell_info.extend(list([key]*len(temp_list[key])))
         start_cell_id+=len(temp_list[key])
    file_path = []
    for key in temp_list.keys():
        file_path.extend(temp_list[key])
    label_info = {}
    label_info['cell_type'] = cell_info
    outpath =f"{data_dir}/label_info.pickle"
    if not os.path.exists(outpath):
        print("generate label_info",outpath)
        with open(outpath, "wb") as f:
            pickle.dump(label_info, f)

                    
    return file_path

def create_bed(res,genome_reference_path):
    file_path=genome_reference_path
    # file_path = '/home/python/higashi/cellcycle/250k/config/mm9.chrom.sizes.txt'  
    with open(file_path, 'r') as file:
        lines = file.readlines()
    chromosome_bins = {}
    # 特定分辨率
    resolution = res 
    # 遍历文件中的每一行
    for line in lines:
        # 分割每一行的染色体和对应的长度
        chromosome, length = line.strip().split('\t')    
        # 计算染色体对应的bin数量
        # bin_count = int(length) // resolution
        bin_count = math.ceil(int(length) / resolution)
        
        # 存储结果到字典中
        chromosome_bins[chromosome] = bin_count
    return chromosome_bins


def get_dimension(lst):
    if type(lst) != type([]):
        return 0
    else:
        return 1 + get_dimension(lst[0])
    


def download_file_with_curl(url, destination):
    subprocess.run(['curl', '-O', url, '-L'], check=True, cwd=destination)


#将高分辨率数据行列编号直接除以与需要的分辨率的倍数。向下取整 ，然后转换成pandas数据帧使用使用 Pandas 的 groupby 进行分组求和得到低分辨率数据
def process_sparse_matrix(sparse_matrix, divisor,bin_count):
    # 将行和列索引除以指定数字并转换为整数
    # sparse_matrix.row=np.ceil(sparse_matrix.row / divisor).astype(int)
    # sparse_matrix.col=np.ceil(sparse_matrix.col / divisor).astype(int)
    # //向下取整
    sparse_matrix.row=sparse_matrix.row // divisor
    sparse_matrix.col=sparse_matrix.col // divisor
    # 将稀疏矩阵转换为 Pandas DataFrame
    df = pd.DataFrame({'row': sparse_matrix.row, 'col': sparse_matrix.col, 'value': sparse_matrix.data})

    # 使用 Pandas 的 groupby 进行分组求和
    grouped_df = df.groupby(['row', 'col']).sum().reset_index()

    # 将结果重新转换为稀疏矩阵
    grouped_sparse_matrix = csr_matrix((grouped_df['value'], (grouped_df['row'], grouped_df['col'])),
                                          shape=(bin_count, bin_count), dtype='float32')

    return grouped_sparse_matrix
    
#从原始cool文件

def process_chrom(chrom,res,sorted_paths,bin_count,output_dir,win_ins):
    #初始化空列表用于存放该染色体的稀疏矩阵
    sparse_matrices_list=[]
    # chrNum = chrom.split("chr")[1]
 
    for path in tqdm(sorted_paths,desc="Processing %s"%chrom):
        clr = cooler.Cooler(path)
        selector = clr.matrix(balance=False,sparse=True)
        data = selector.fetch(chrom,chrom)  # same as fetch(chrm, chrm)
        if win_ins==1:
            sparse_matrix_csr = csr_matrix((data.data, (data.row, data.col)),
                                          shape=(bin_count, bin_count), dtype='float32')
        else:       
            sparse_matrix_csr = process_sparse_matrix(data, win_ins,bin_count)
        sparse_matrices_list.append(sparse_matrix_csr)

    np.save(os.path.join(output_dir,'%s_sparse_adj.npy'% chrom), sparse_matrices_list)
    del sparse_matrices_list
    print (chrom, "finished")


def main():
    folder_path = '/home/dataset/snm3c/10k_cool' 
    data_dir = '/home/dataset/snm3c/250k'
    output_dir = os.path.join(data_dir,'raw')
    if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    # res = 10000
    genome = 'hg19'
    genome_reference_path=f"{data_dir}/{genome}.chrom.sizes"
    if not os.path.exists(f"{data_dir}/{genome}.chrom.sizes"):
        print(f"download {genome}.chrom.sizes")
        download_file_with_curl(f"http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.chrom.sizes", data_dir)

    if not os.path.exists(f"{data_dir}/cytoBand.txt.gz"):
        print("download cytoBand.txt.gz")
        download_file_with_curl(f"http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/database/cytoBand.txt.gz", data_dir)
        
    
    
   

    res = 250000
    win_ins = res/10000
    chrom_list = ["chr6","chr7","chr8","chr9","chr10",
    "chr11","chr12","chr13","chr14","chr15",
    "chr16","chr17","chr18","chr19","chrX"]
    # chrom_list = ["chr1"]
    chrom_list = np.array(chrom_list)
    print(chrom_list)
    sorted_paths=find_file(folder_path,data_dir)
    chromosome_bins = create_bed(res,genome_reference_path)
    for chrom in chrom_list:
        bin_count=chromosome_bins[chrom]
        print(bin_count)
        process_chrom(chrom,res,sorted_paths,bin_count,output_dir,win_ins)

if __name__ == '__main__':
    main()

