#!/bin/bash

# 指定文件列表文件路径
file_list="new.txt"

# 检查文件列表文件是否存在
if [ ! -f "$file_list" ]; then
  echo "File list '$file_list' not found."
  exit 1
fi


#/usr/local/anaconda/envs/HiC-Pro_v3.0.0/etc/asperaweb_id_dsa.openssh 


# 读取文件列表文件的每一行并下载文件
while IFS= read -r file_path; do
  echo "Downloading $file_path"
  ascp -i /usr/local/anaconda/envs/HiC-Pro_v3.0.0/etc/asperaweb_id_dsa.openssh -l 300M -k 1 -T anonftp@ftp.ncbi.nlm.nih.gov:$file_path
  echo "Downloaded $file_path"
done < "$file_list"

echo "All downloads completed."