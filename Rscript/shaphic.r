library(rhdf5)
library(vegan)
library(GUniFrac)

group_id<-'compartment_raw'

data_dir<-'/home/dataset/snm3c/250k/raw'
df<-paste(data_dir,group_id,'compartment_250k_shap.h5',sep = '/')
file_path <- paste(data_dir, group_id, "shap_chord.csv",sep = '/')

file <- H5Fopen(df)
group <- H5Gopen(file, group_id)
datasets_info <- h5ls(group)
bin_datasets <- datasets_info$name[grep("^bin_", datasets_info$name)]
results_table <- data.frame(matrix(ncol = 5, nrow = length(bin_datasets)))
colnames(results_table) <- c("bin","p_value",'effect_size','p_adjust','F')
#获取细胞id
# label <- group$label$cell_id
#获取细胞分类标签
rhs <- group$label$cell_type

#计算程序执行剩余时间
calculateRemainingTime <- function(start_time, n, i) {
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  time_per_iteration <- elapsed_time / i
  remaining_iterations <- n - i
  remaining_time <- remaining_iterations * time_per_iteration
  
  # 将时间转换为易读格式，保留一位小数，单位为分钟
  remaining_time_str <- as.difftime(round(remaining_time, 1), units = "mins")
  
  # 返回计算出的时间
  return(remaining_time_str)
}
#使用R基础包设置进度条风格3
pb <- txtProgressBar(style=3)
start_time <- Sys.time() ## 记录开始时间
for (i in 1:length(bin_datasets)) {
  
  bin_name <- paste0("bin_", i)
  lhs <- h5read(file, paste0(group_id,"/", bin_name))
  
  if (all(lhs == 0)) {
    print(paste("Bin", bin_name, "contains all zeros. Skipping adonis2 calculation and result writing..."))
    results_table[i, "bin"] <- bin_name
    next  # 跳出循环
    
  } else {
    # 如果 lhs 不全为 0，进行 adonis2 或者（dmanova）计算和结果写入
    # table <- data.frame(data = t(lhs), row.names = label)
    table <- data.frame(data = t(lhs))
    distance_matrix <- vegdist(t(lhs), method = "chord")
    # distance_matrix <- vegdist(t(lhs), method = "manhattan")
    adonis_result <- dmanova(distance_matrix ~ rhs, data = table)
   
    N <- adonis_result$aov.tab[rownames(adonis_result$aov.tab) == "Total", "Df"] + 1
    MS_res <- adonis_result$aov.tab[pmatch("Residual", rownames(adonis_result$aov.tab)), "MeanSqs"]
    SS_tot <- adonis_result$aov.tab[rownames(adonis_result$aov.tab) == "Total", "SumsOfSqs"]
    
    omega <- apply(adonis_result$aov.tab, 1, function(x) (x["Df"]*(x["MeanSqs"]-MS_res))/(x["Df"]*x["MeanSqs"]+(N-x["Df"])*MS_res))
    #omega <- apply(adonis_result$aov.tab, 1, function(x) (x["SumsOfSqs"]-x["Df"]*MS_res)/(SS_tot+MS_res))#伪w2
    
    # 将结果写入数据帧第i行对应列
    results_table[i, "bin"] <- bin_name
    results_table[i, "p_value"] <- adonis_result$aov.tab$'Pr(>F)'[1]
    results_table[i, "effect_size"] <- omega[1]
    results_table[i, "F"] <- adonis_result$aov.tab$F.Model[1]
    # print(results_table[i,"p_value"])
  }
  ## 第二个位置：实时反映进度
  setTxtProgressBar(pb, i/length(bin_datasets))
  remaining_time_str <- calculateRemainingTime(start_time, length(bin_datasets), i)
  # 打印剩余时间
  cat(sprintf("\r剩余时间: %s min", remaining_time_str))
}
cat("\n")
pvalue_adjust <- p.adjust(results_table$p_value, method="BH")
results_table$p_adjust <- pvalue_adjust

#关闭进度条
close(pb)
end_time <- Sys.time()  ## 记录结束时间
run_time <- end_time - start_time  ## 计算程序运行耗时
run_time
# file_path <- paste0(data_dir, group_id, "_ori_pc1_shap_chord_no_softmax.csv")

write.csv(results_table, file = file_path, row.names = FALSE)
print(paste("csv save to:", file_path))