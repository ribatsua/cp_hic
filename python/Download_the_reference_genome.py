#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@File    :   Download_the_reference_genome.py
@Time    :   2024/01/05 15:11:50
@Author  :   cp 
@Version :   1.0
@Desc    :   None
'''
#从UCSC下载参考基因组文件

import subprocess
import os


def download_file_with_curl(url, destination):
    subprocess.run(['curl', '-O', url, '-L'], check=True, cwd=destination)

def download_and_extract(genome, folder=None):
    if folder is None or folder.strip() == "":
        folder = f"/home/annotation/{genome}"
        if not os.path.exists(folder):
            os.makedirs(folder)
    if not os.path.exists(f"{folder}/{genome}.chrom.sizes"):
        print(f"download {genome}.chrom.sizes")
        download_file_with_curl(f"http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.chrom.sizes", folder)
    if not os.path.exists(f"{folder}/cytoBand.txt.gz"):
        print("download cytoBand.txt.gz")
        download_file_with_curl(f"http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/database/cytoBand.txt.gz", folder)
        print(f"Unzipping {genome}.cytoBand.txt")
        subprocess.run(f"gunzip -c {folder}/cytoBand.txt.gz > {folder}/cytoBand.txt", shell=True, check=True)

genome = 'mm9'
download_and_extract(genome)
