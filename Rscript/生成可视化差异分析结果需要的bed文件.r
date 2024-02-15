library(data.table)
library(rhdf5)

data_dir='/home/dataset/snm3c/250k/raw/compartment_raw/'
bed_path<-paste0(data_dir, "padjust_snm3c_250k.bedGraph")
dist_path<-paste0(data_dir, "dist_snm3c_250k.bedGraph")

data_path<-paste0("/home/dataset/snm3c/250k/raw/compartment_raw/shap_chord.csv")
# label_path<-paste0(data_dir,"label_bed.csv")
label_path<-'/home/dataset/snm3c/250k/raw/test_chr1_2_compartment.hdf5'

data<-fread(data_path)

label_file<-H5Fopen(label_path)

group_id <- 'compartment_raw'

# 读该组数据
group <- h5read(label_file, group_id)
# 关闭h5文件
H5Fclose(label_file)

bin<-group$bin
# bin$chr
label <- data.frame(matrix(nrow = length(bin$chr), ncol = 3))
colnames(label)<-c('chr','start','end')
label[,c('chr','end','start')]<-bin
data1<-label
data1[,'pvalue']<-log10(data$p_adjust+10^(-200))
data1[,'effect_size']<-data$effect_size
# str(label)
# label<-fread(label_path)

lines <- c(
  "# locus chr19:49302001-49304701",
  "# refGene encodeRegions",
  "# zero-based, half-open coords",
  'track name="-log10(shap_Padjust)" description="BedGraph format" visibility=full color=0,153,0 priority=20 plotType="points"'
)
lines1 <- c(
  "# locus chr19:49302001-49304701",
  "# refGene encodeRegions",
  "# zero-based, half-open coords",
  'track name="shap_dist" description="BedGraph format" visibility=full color=0,153,0 priority=20 plotType="points"'
)

# 将特定文本写入文件
writeLines(lines, bed_path)
writeLines(lines1, dist_path)

############################################################################################################################
# data<-fread("/home/python/higashi/notebook/cortex250k_compartment_shap_chord_dmanova_w2_default_b150.csv")
# data<-fread("/home/python/higashi/notebook/cortex250k_compartment_shap_chord_dmanova_w2_default_b150_remove_chrX_ncb.csv")
# data<-fread('cortex250k_compartment_shap_chord_dmanova_w2_default_b150_remove_chrX_ncb.csv')
# data<-fread('/home/compartment_raw_sim_shap_manhattan.csv')
# data<-fread('shap_hic.csv')
# data<-fread('/home/compartment_zscore_shap_hic.csv')
############################################################################################################################

data1 <- subset(data1, effect_size>0.16)[,1:4]
str(data1)

# filtered_data <- subset(data, effect_size>0.16)
# ############################################################################################################################
# # filtered_data <- subset(data, p_value<10^(-6))
# ############################################################################################################################


# data1=data.frame(matrix(nrow = length(filtered_data$effect_size)))
# data1 <- data.frame(p_value = -log10(filtered_data$p_adjust+10^(-200)), index = as.numeric(sub("bin_(\\d+)", "\\1", filtered_data$bin)) - 1)
# ############################################################################################################################
# # data1 <- data.frame(p_value = -log10(filtered_data$p_value), index = as.numeric(sub("bin_(\\d+)", "\\1", filtered_data$bin)) - 1)
# ############################################################################################################################
# # str(data1)

# merged_data <- merge(label, data1, by.x = "index", by.y = "index")


# data2 <- data.frame(
#   chrom = merged_data$chr,
#   start = merged_data$start,
#   end = merged_data$end,
#   p_value = merged_data$p_value
# )

# str(data2)
# data3=data.frame(matrix(nrow = length(filtered_data$effect_size)))
# data3 <- data.frame(dist = filtered_data$F, index = as.numeric(sub("bin_(\\d+)", "\\1", filtered_data$bin)) - 1)

# merged_data2<-merge(label, data3, by.x = "index", by.y = "index")
# str(merged_data)
# data4<- data.frame(
#   chrom = merged_data2$chr,
#   start = merged_data2$start,
#   end = merged_data2$end,
#   dist = merged_data2$dist
# )
# str(data4)
# merged_data2$index
# write.table(data2, bed_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE,append = TRUE)

# write.table(data4, dist_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE,append = TRUE)