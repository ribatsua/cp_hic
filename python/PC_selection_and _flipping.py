#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@File    :   PC_selection_and _flipping.py
@Time    :   2024/01/05 16:05:47
@Author  :   cp 
@Version :   1.0
@Desc    :   pc挑选和翻转
'''

import os
import subprocess
import h5py 
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from itertools import chain
import re


import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr


# 激活自动转换为pandas数据框
pandas2ri.activate()

# 安装并导入必要的R包
base = importr('base')
stats = importr('stats')





def download_file_with_curl(url, destination):
    subprocess.run(['curl', '-O', url, '-L'], check=True, cwd=destination)


def download_and_extract(genome, folder, resolution):
    if folder is None or folder.strip() == "":
        folder = f"{genome}_{int(resolution)}_goldenpathData"
        if not os.path.exists(folder):
            os.makedirs(folder)
            current_path = os.getcwd()
            os.chdir(current_path)
        folder = os.path.normpath(folder)
    else:
        folder = os.path.normpath(folder)

    genome_fa_file = f"{folder}/{genome}.fa"
    if not os.path.exists(f"{folder}/cytoBand.txt.gz"):
        print("download cytoBand.txt.gz")
        download_file_with_curl(f"http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/database/cytoBand.txt.gz", folder)
        

    if not os.path.exists(f"{folder}/{genome}.chrom.sizes"):
        print(f"download {genome}.chrom.sizes")
        download_file_with_curl(f"http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.chrom.sizes", folder)
    if not os.path.exists(f"{folder}/{genome}.refGene.gtf.gz"):
        print(f"download {genome}.refGene.gtf.gz")
        download_file_with_curl(f"http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/genes/{genome}.refGene.gtf.gz", folder)
        

    if not os.path.exists(f"{folder}/{genome}.fa.gz"):
        genome_fa_url = f"http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.fa.gz"
        print(f"download {genome}.fa.gz")
        download_file_with_curl(genome_fa_url, folder)
    if not os.path.exists(f"{folder}/{genome}.fa"):  
        print(f"Unzipping {genome}.fa")
        #解压文件输出定向到{folder}文件夹下 {genome}.fa文件
        subprocess.run(f"gunzip -c {folder}/{genome}.fa.gz > {folder}/{genome}.fa", shell=True, check=True)


    tss_bed_file = f"{folder}/{genome}.tss.bed"
    if not os.path.exists(tss_bed_file):
        cmd = f"gunzip -c {folder}/{genome}.refGene.gtf.gz | awk -v OFS='\\t' '{{if($3==\"transcript\"){{if($7==\"+\"){{print $1,$4,$4+1}}else{{print $1,$5-1,$5}}}}}}' | grep -v 'alt' | grep -v 'random' | sort |uniq |sort -k 1,1 -k2,2n > {folder}/{genome}.tss.bed"
    
        print(f"Running {cmd}")
        subprocess.run(cmd, shell=True, check=True)

    binned_bed_file = f"{folder}/{genome}.binned.bed"
    if not os.path.exists(binned_bed_file):
        cmd = f"bedtools makewindows -g {folder}/{genome}.chrom.sizes -w {int(resolution)} > {folder}/{genome}.binned.bed"
        print(f"Running {cmd}")
        subprocess.run(cmd, shell=True, check=True)

    gcpt_bedgraph_file = f"{folder}/{genome}.GCpt.bedGraph"
    if not os.path.exists(gcpt_bedgraph_file):
        cmd = f"bedtools nuc -fi {folder}/{genome}.fa -bed  {folder}/{genome}.binned.bed | grep -v '#' | awk -v OFS='\\t' '{{print $1,$2,$3,$5}}' | grep -v 'alt' | grep -v 'random' | sort -k 1,1 -k2,2n > {folder}/{genome}.GCpt.bedGraph"
        print(f"Running {cmd}")
        subprocess.run(cmd, shell=True, check=True)

    gcpt_tss_bedgraph_file = f"{folder}/{genome}.GCpt.tss.bedGraph"
    if not os.path.exists(gcpt_tss_bedgraph_file):
        cmd = f"bedtools map -a {folder}/{genome}.GCpt.bedGraph -b {folder}/{genome}.tss.bed -c 1 -o count -null 0 > {folder}/{genome}.GCpt.tss.bedGraph"
        print(f"Running {cmd}")
        subprocess.run(cmd, shell=True, check=True)
    print("GCpt.tss.bedGraph already exists")
    goldenpath = f"{folder}/{genome}.GCpt.tss.bedGraph"
    return goldenpath

def process_and_transpose(gcc_values):
    # 将字典中的每个值转换为 DataFrame 
    dfs = {key: pd.DataFrame(value) for key, value in gcc_values.items()}

    # 使用 pd.concat 在轴 1 上合并 DataFrame
    result_df = pd.concat(dfs.values(), axis=1)

    # 转置 DataFrame
    result_df = result_df.transpose()
    # 在列名后添加 .cor 后缀
    result_df.columns = result_df.columns + '.cor'

    # 添加 'name' 列，该列包含原始 DataFrame 的索引
    result_df['name'] = result_df.index

    return result_df
vals_path = 'Python_vals.txt'

if os.path.exists(vals_path):
    os.remove(vals_path)
def run_r_function(i,pc_mat, sam, gcc_values, tss_values, len_values, chr_list):
    # 将Python对象转换为R对象
    r_pc_mat=pandas2ri.py2rpy(pc_mat)
    # r_pc_mat = robjects.r['as.matrix'](r_pc_mat)
    sam_length = len(sam)
    gcc_values_r = pandas2ri.py2rpy(gcc_values)
    tss_values_r = pandas2ri.py2rpy(tss_values)
    len_values_r = pandas2ri.py2rpy(len_values)

    # 在R环境中定义r_pc_mat
    robjects.globalenv['i'] = i
    robjects.globalenv['r_pc_mat'] = r_pc_mat
    robjects.globalenv['sam_length'] = sam_length
    # 在设置 chr_list 到 R 环境之前，将其转换为字符向量
    robjects.globalenv['chr_list'] = robjects.StrVector(chr_list)

    # robjects.globalenv['chr_list'] = chr_list
    robjects.globalenv['gcc_values_r'] = gcc_values_r
    robjects.globalenv['tss_values_r'] = tss_values_r
    robjects.globalenv['len_values_r'] = len_values_r
    

    # 执行R代码
    robjects.r('''
        # print(head(r_pc_mat))
        # print(chr_list)
        # print(gcc_values_r)
    
        hcl <- hclust(as.dist(round(1-cor(r_pc_mat),4)))
		cl<- list()
		k <- 1
		h <- 0.05   
		while (h < 1) {
			g <- data.frame(group=cutree(hcl, h=h))	
			g[,"name"]  <- rownames(g)
			g[,"name2"] <- paste0(g$group,"-",do.call(rbind,strsplit(g$name,"[.]"))[,1])
			v <- data.frame(member=table(g$name2))
			v$member.Var1 <- do.call(rbind,strsplit(as.character(v$member.Var1),"-"))[,1]
			v <- data.frame(member=table(v$member.Var1))
            # print(g)
			s <- v[v$member.Freq == sam_length,]$member.Var1
			if (length(s) > 0) {

				for(u in 1:length(s)) {
					if (length(unique(g[g$group==s[u],]$name)) == sam_length) {
						cl[[k]] <- g[g$group==s[u],]
						k <- k + 1
					}
				}
			}
			h <- h + 0.025
            # print(paste("Current h value:", h))
            # print(paste("Number of groups:", length(unique(g$group))))
            # print(paste("Number of unique s values:", length(unique(s))))
		}
        if (length(cl) > 0) {
			cl <- as.data.frame(do.call(rbind, cl))
			cl <- unique(cl)
			cl[,"chr"] <- chr_list[i+1]
			cl[,"gcc.cor"] <- 0
			cl[,"tss.cor"] <- 0

			for(v in 1:nrow(gcc_values_r)) {
				cl[cl$name == gcc_values_r$name[v],"gcc.cor"] <- gcc_values_r$gcc.cor[v]
				cl[cl$name == gcc_values_r$name[v],"tss.cor"] <- tss_values_r$tss.cor[v]
				cl[cl$name == gcc_values_r$name[v],"len.cor"] <- abs(len_values_r$len.cor[v])
			}
			data$clus[[d]] <- cl
			data$cor[[d]] <- as.data.frame(aggregate(gcc.cor ~ group + chr, mean, data=cl))
			data$cor[[d]][,"tss.cor"] <- aggregate(tss.cor ~ group + chr, mean, data=cl)$tss.cor
			data$cor[[d]][,"len.cor"] <- aggregate(len.cor ~ group + chr, mean, data=cl)$len.cor
			data$cor[[d]][,"score"] <- apply(data$cor[[d]][,c("tss.cor","gcc.cor")],1,sum)
			data$cor[[d]]$gcc.cor <- round(data$cor[[d]]$gcc.cor, 4)
			data$cor[[d]]$tss.cor <- round(data$cor[[d]]$tss.cor, 4)
			data$cor[[d]]$len.cor <- round(data$cor[[d]]$len.cor, 4)
			data$cor[[d]]$score   <- round(data$cor[[d]]$score, 4)
			data$chrom[[d]] <- chr_list[i+1]
			d <- d + 1
		}
        length_cl<-length(cl)''')

    # 从R中检索结果
    data= robjects.globalenv['data']
    length_cl = robjects.globalenv['length_cl']
    # data_clus = robjects.globalenv['data'].rx2('clus')
    # data_cor = robjects.globalenv['data'].rx2('cor')
    # data_chrom = robjects.globalenv['data'].rx2('chrom')
    # print(data_clus)
    # 将R数据框转换为Pandas数据框
    # data_clus = pandas2ri.ri2py(data_clus)
    # data_cor = pandas2ri.ri2py(data_cor)
    # data_chrom = pandas2ri.ri2py(data_chrom)
    # 使用提取数据的方法
    # data_clus_df = pandas2ri.rpy2py_dataframe(data_clus)
    # data_cor_df = pandas2ri.rpy2py_dataframe(data_cor)
    # data_chrom_df = pandas2ri.rpy2py_dataframe(data_chrom)

    return data,length_cl
robjects.r('''
        d <- 1
        data <- list(clus=list(),cor=list(),chrom=list())
        ''')


genome_name = "mm10"
resolution_value = 250000
pc_k = 2
data_folder = None
#下载参考基因组文件并生成 GCpt.tss.bedGraph文件（用于后续挑选pc）
goldenpath = download_and_extract(genome_name, data_folder, resolution_value)
group_id = 'compartment_raw'
df = '/home/python/higashi/dataset_hic/dataset2/cortex250k/compartment_raw.h5'
file = h5py.File(df, 'r')
bin_chr = np.array([chrom.decode('utf-8') for chrom in file['compartment_raw']['bin']['chrom']])
bin_start = file[group_id + '/bin/start'][:]
bin_end = file[group_id + '/bin/end'][:]

chr_list = list(np.unique(bin_chr))
chr_list.sort(key=lambda l: int(re.findall('\d+', l)[0]))
# chr_list = ['chr1','chr2']
print(chr_list)
vals = {}

cell_id = [1,2,3]
sam = ['cell_' + str(cell) for cell in cell_id]

for i in range(len(chr_list)):
    chr_index = np.where(bin_chr == chr_list[i])[0]
    gcc_values = {}
    tss_values = {}
    len_values = {}
    chrom_list = {}
    count_vect = {}

    for j in range(len(sam)):
        print(f"在 {sam[j]} 样本中运行 {chr_list[i]}")
        cell_data = file[group_id + f'/cell_{j}/'][:]
        #按染色体名取对应位置数据
        pca = pd.DataFrame({'chr': bin_chr[chr_index], 'start': bin_start[chr_index], 'end': bin_end[chr_index],'index': bin_end[chr_index]/resolution_value})
        lhs = pd.DataFrame(data=cell_data)
        # print(pca.shape)
        pca[['pc1', 'pc2', 'pc3']] = lhs.iloc[chr_index, :].values
        pca =pca.iloc[:,0:(4+pc_k)]
        
        #添加行名如：chr_250000
        pca.index = pca["chr"].astype(str) + "_" + pca["start"].astype(str)
        
        gcc = pd.read_table(goldenpath, header=None, names=["chr", "start", "end", "gcc", "tss"])
        # 为GCC数据添加行名
        gcc.index = gcc["chr"].astype(str) + "_" + gcc["start"].astype(str)

        # print(pca)

        # 将GCC和TSS信息添加到PCA中
        pca[["gcc", "tss"]] = gcc.loc[pca.index, ["gcc", "tss"]]
        pca['len'] = range(1, pca.shape[0] + 1)
        # print(pca.iloc[:, 4])
        # 生成bedGraph文件
        for k in range(4, pca.shape[1] - 3):

            # 计算pc值和gcc的相关性 决定pc值是否乘以-1
            pca.iloc[:, k] = np.sign(pearsonr(pca.iloc[:, k], pca["gcc"])[0])*pca.iloc[:, k]

        # 存储染色体的pc值(正负校正后)
        chrom_list[j] = pca.iloc[:, 4:(pca.shape[1] - 3)]
        
        colnames_chrom = [f"{sam[j]}.PC{i}" for i in range(1, chrom_list[j].shape[1] + 1)]
        # print(colnames_chrom)
        chrom_list[j].columns = colnames_chrom

        # 存储行名
        count_vect[j] = chrom_list[j].index
        # print(count_vect[j])
        # 计算pc值和GCC、TSS、行数的相关性
        gcc_values[j]=[pd.concat([chrom_list[j], pca["gcc"]], axis=1).corr().iloc[:-1, -1].round(4).transpose()]
        # print(gcc_values)
        tss_values[j]=[pd.concat([chrom_list[j], pca["tss"]], axis=1).corr().iloc[:-1, -1].round(4).transpose()]
        len_values[j]=[pd.concat([chrom_list[j], pca["len"]], axis=1).corr().iloc[:-1, -1].round(4).transpose()]
        # print(gcc_values[j])

    gcc_df=process_and_transpose(gcc_values)
    tss_df = process_and_transpose(tss_values)
    len_df = process_and_transpose(len_values)
    #合并相关性结果
    vals[i] = pd.DataFrame()
    vals[i]['gcc.cor'] = gcc_df['gcc.cor']
    vals[i]['tss.cor'] = tss_df['tss.cor']
    vals[i]['len.cor'] = len_df['len.cor']
    vals[i]['chr'] = [chr_list[i]]*len(gcc_df['gcc.cor'])
    vals[i]['sample'] = vals[i].index.str.split('.').str[0]
    vals[i]['pc']= vals[i].index.str.split('.').str[1]

    vals[i].to_csv(vals_path, sep="\t", header = False,index=False,mode='a')
    # print(vals[i])

    # print(count_vect[0]+count_vect[1])
    #字典中的列表转换成pd.Series并合并成一列
    # values, counts = np.unique( pd.concat([pd.Series(count_vect[0]), pd.Series(count_vect[1])], ignore_index=True), return_counts=True)
    #改成适用于一个或多个列表
    values, counts = np.unique(pd.concat([pd.Series(lst) for lst in count_vect.values()], ignore_index=True), return_counts=True)
    # 将结果转换为数据框
    count_vect_df = pd.DataFrame({'value': values, 'freq': counts})
    # print(count_vect_df)
    # print(count_vect_df)
    pc_mat=[]
    # 遍历 sam 列表的每个元素
  
    for j in range(len(sam)):
        # 选择 chrom_list 中的第 j 个元素的子集
        subset = chrom_list[j].loc[count_vect_df['value'], :]

        # 将子集添加到 pc_mat 列表中
        pc_mat.append(subset)
    # 使用 numpy.c_ 将 pc_mat 列表中的矩阵按列合并
    pc_mat_combined = np.c_[tuple(pc_mat)]

    # 将结果转换为 DataFrame 行索引使用count_vect_df['value']列 列名使用vals[i].index
    pc_mat_df = pd.DataFrame(pc_mat_combined,index=count_vect_df['value'],columns=vals[i].index)


    data,length_cl = run_r_function(i,pc_mat_df, sam, gcc_df, tss_df, len_df, chr_list)

robjects.r('''
    clus.df <- as.data.frame(unique(do.call(rbind,data$clus)))
	cor.df  <- as.data.frame(do.call(rbind,data$cor))
	chr.vec <- unlist(data$chrom)
        ''')
clus_df = pandas2ri.rpy2py_dataframe(robjects.globalenv['clus.df'])
cor_df = pandas2ri.rpy2py_dataframe(robjects.globalenv['cor.df'])

#############################处理完所有染色体后
chr_max = []
for i in range(len(chr_list)):
    cor_df_chrom = cor_df[cor_df['chr']==chr_list[i]]
    clus_df_chrom = clus_df[clus_df['chr']==chr_list[i]]
    if len(cor_df_chrom)>0:
        selected_rows = cor_df_chrom.loc[cor_df_chrom['score'].idxmax()]
        name = clus_df_chrom[clus_df_chrom['group']==selected_rows['group']]['name'].str.split('[.]', expand=True)

    #    name_split = selected_rows['name'].str.split('[.]', expand=True)


        chr_max.append(pd.DataFrame({
            'group': selected_rows['group'],
            'chr': selected_rows['chr'],
            'gcc.cor': selected_rows['gcc.cor'],
            'tss.cor': selected_rows['tss.cor'],
            'len.cor': selected_rows['len.cor'],
            'sample': [','.join(name[0])],
            'pcs': [','.join(map(str, name[1]))]
        }))
    # print(chr_list[i])
    else:
        vals_chrom = vals[i]
        vals_chrom['score'] = vals_chrom['gcc.cor']+vals_chrom['tss.cor']
        # print(vals_chrom)
        sample_df = {}
        #选择每个sample中score最大的pc所在的行
        for j in range(len(sam)):
            sample_df[j] = vals_chrom[vals_chrom['sample']==sam[j]]
            
            sample_df[j] = sample_df[j].loc[sample_df[j]['score'].idxmax()]
            
            # print(type(sample_df[j]))
        #按行合并并转换成数据帧
        sample_df = pd.DataFrame(pd.concat(sample_df,axis=1).T)
        # print(sample_df)
        # 计算每列的均值
        gcc_cor_mean = round(sample_df['gcc.cor'].mean(), 4)
        tss_cor_mean = round(sample_df['tss.cor'].mean(), 4)
        len_cor_mean = round(sample_df['len.cor'].mean(), 4)
        
        # 将结果添加到 chr_max 列表中
        chr_max.append(pd.DataFrame({
            'group': [1],
            'chr': [chr_list[i]],
            'gcc.cor': [gcc_cor_mean],
            'tss.cor': [tss_cor_mean],
            'len.cor': [len_cor_mean],
            'sample': [','.join(sample_df['sample'])],
            'pcs': [','.join(map(str, sample_df['pc']))]
        }))
# 使用 pd.concat 将列表中的数据框按行合并
chr_max_df = pd.concat(chr_max, ignore_index=True)

# 将数据框写入文本文件
chr_max_df.to_csv("Python_test_chr_pc_selected1.txt", sep="\t", index=False)