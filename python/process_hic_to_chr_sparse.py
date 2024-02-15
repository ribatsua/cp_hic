#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@File    :   process_hic_to_chr_sparse.py
@Time    :   2024/01/04 16:19:04
@Author  :   cp 
@Version :   1.0
@Desc    :   None
'''


#由hic文件提取稀疏矩阵按染色体保存整个数据集
import hicstraw
import os
import numpy as np
from scipy.sparse import csr_matrix
from tqdm import tqdm
import math

#计算每条染色体在分辨率res下划分的bin的数量
def create_bed(res,genome_reference_path):
    file_path=genome_reference_path
    # file_path = '/home/python/higashi/cellcycle/250k/config/mm9.chrom.sizes.txt'  
    with open(file_path, 'r') as file:
        lines = file.readlines()
    chromosome_bins = {}
    # 特定分辨率
    resolution = res # 请替换成你的分辨率
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
#按文件名对文件列表排序
def sort_key(path):
    # 先使用斜杠分割路径，然后使用下划线分割每个部分，获取数字作为排序关键字
    parts = []
    for part in path.split('/'):
        parts.extend(part.split('_'))
    return [int(part) if part.isdigit() else part for part in parts]
#获取数据集根目录下所有包含特定命名的hic文件并返回排序后的文件列表
def find_file(folder_path):
    file_path=[]
    for root, dirs, files in os.walk(folder_path):
            for file in files:
                # if file.endswith('.txt') and chrom == file.split(".")[0].split("_")[-1]: #添加.split("_")[-1]后通用于以cell_1_chr1.txt命名文件
                if file.endswith('.hic') and 'cortex' == file.split(".")[0].split("-")[0].split("_")[-1]:
                # if file.endswith('.txt'):
                    file_path.append(os.path.join(root, file))
    sorted_paths = sorted(file_path, key=sort_key)
    return sorted_paths          

def process_chrom(chrom,res,sorted_paths,bin_count,output_dir):
    sparse_matrices_list=[]
    chrNum = chrom.split("chr")[1]
    for path in tqdm(sorted_paths,desc="Processing %s"%chrom):
        result = hicstraw.straw("observed", 'NONE', path, chrNum,chrNum, 'BP', res)
        sorted_result = sorted(result, key=lambda x: (x.binX, x.binY))
        rows = []
        cols = []
        values = []
        for elem in sorted_result:
            row, col, value = (int(elem.binX/res),int(elem.binY/res),elem.counts)
            rows.append(row)
            cols.append(col)
            values.append(value)
        sparse_matrix_csr = csr_matrix((values, (rows, cols)), shape=(bin_count, bin_count), dtype=np.float32)#统一稀疏矩阵形状为bin数量
        sparse_matrices_list.append(sparse_matrix_csr)

    np.save(os.path.join(output_dir,'%s_sparse_adj.npy'% chrom), sparse_matrices_list)
    print (chrom, "finished")

def main():
    genome_reference_path='/home/annotation/mm10/mm10.chrom.sizes'
    folder_path = '/home/dataset/tan/hic_data'
    output_dir = '/home/dataset/tan/cortex/100k/raw'

    res = 100000
    chrom_list = ["chr1","chr2","chr3","chr4","chr5",
    "chr6","chr7","chr8","chr9","chr10",
    "chr11","chr12","chr13","chr14","chr15",
    "chr16","chr17","chr18","chr19","chrX"]
    chrom_list = np.array(chrom_list)
    print(chrom_list)
    sorted_paths=find_file(folder_path)
    chromosome_bins = create_bed(res,genome_reference_path)
    for chrom in chrom_list:
        bin_count=chromosome_bins[chrom]
        process_chrom(chrom,res,sorted_paths,bin_count,output_dir)

if __name__ == '__main__':
    main()