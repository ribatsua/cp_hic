library(data.table)
library(ggplot2)
library(png)
args <- commandArgs(trailingOnly = TRUE)
# print(args)

if (length(args) != 1) {
  stop("Usage: Rscript script.R <parameter>")
}

parameter <- args
# parameter<-'mean4'
print(paste0("处理参数:", parameter, "\n"))


data_dir<-'/home/dataset/sim_data/cellcycle/1m/png/'
A<-fread('/home/dataset/sim_data/cellcycle/1m/compartment_raw_0_1171.csv_mean1+var1.csv')
B<-fread(paste0('/home/dataset/sim_data/cellcycle/1m/compartment_raw_0_1171.csv_',parameter,'.csv'))

# max(A$Varianceraw)
test_res<-fread(paste0('/home/dataset/sim_data/cellcycle/1m/res/',parameter,'_manhattan.csv'))

mean_list<-list(c(0,0.5),c(0.5,1),c(1,1.5),c(1.5,2),c(2,2.5),c(2.5,3))
# mean_list<-list(c(0,0.5),c(0.5,1),c(1,1.5))

mode<-c('Meanraw')


test_res[[mode]]<-A[[mode]]
# A[[mode]]
# test_res[[mode]]
# min(A[[mode]])
# # A[[mode]]


# mean_list<-list(c(0,0.5),c(0.5,1),c(1,1.5))
length(mean_list)
search_term <- "wd|vd"
threshold<-0.16
recall_list <- numeric(length = length(mean_list))
precision_list<-numeric(length = length(mean_list))
acc_list<-numeric(length = length(mean_list))
f1_list<-numeric(length = length(mean_list))



for (i in 1:length(mean_list)){
    # print(i[1])
    # print(i[2])
    a<-A[abs(A[[mode]])>mean_list[[i]][1]&abs(A[[mode]])<=mean_list[[i]][2]]
    # print(length(a))
    b<-B[abs(B[[mode]])>mean_list[[i]][1]&abs(B[[mode]])<=mean_list[[i]][2]]
    c<-test_res[abs(test_res[[mode]])>mean_list[[i]][1]&abs(test_res[[mode]])<=mean_list[[i]][2]]  
    # str(a)
    # print(paste0('b_bin:',b$bin))
    resultA <- a[grep(search_term, a$datalabel)]$bin
    resultB <- b[grep(search_term, b$datalabel)]$bin
    # print(paste('result_B',resultB))
    combined_result <- union(resultA, resultB) #实际有显著差异
    N<-length(c$effect_size)
    print(paste("N",N))
    # data<-fread('/home/compartment_raw_sim_var_shap_manhattan.csv')
    filtered_data <- subset(c, effect_size>threshold)
    filtered_data$bin_number <- as.numeric(sub("bin_", "", filtered_data$bin))-1#预测为有显著差异
    c$bin_number <- as.numeric(sub("bin_", "", c$bin))-1
    # print(paste("c$bin_number:", c$bin_number))

    # print(paste("combined_result:", combined_result))

    combined_nd<-setdiff(c$bin_number,combined_result) #实际无显著差异
    # print(paste("combined_nd:", combined_nd))
    N_number<-setdiff(c$bin_number,filtered_data$bin_number) #预测为无显著差异
    length(N_number)
    length(filtered_data$bin_number)
    length(combined_nd)
    length(combined_result)
    # print(paste("length(N_number):", length(N_number)))
    # print(paste("length(filtered_data$bin_number)", length(filtered_data$bin_number)))
    # print(paste("length(combined_nd)",length(combined_nd)))
    # print(paste("length(combined_result)", length(combined_result)))
    # 合并resultA和resultB并保留唯一值
    # combined_result <- union(resultA, resultB)
    TP=length(intersect(filtered_data$bin_number, combined_result))
    TN=length(intersect(combined_nd, N_number))
    FN=length(intersect(combined_result, N_number))
    FP=length(intersect(filtered_data$bin_number,combined_nd))

    print(paste("TP:", TP))
    print(paste("TN:", TN))
    print(paste("FP:",FP))
    print(paste("FN:", FN))
    recall <- ifelse(TP == 0 | (TP + FN)==0, 0, TP / (TP + FN))
    precision<-ifelse(TP == 0 | (TP + FP)==0, 0, TP / (TP + FP))
    # precision=TP/(TP+FP)
    # recall=TP/(TP+FN)
    acc=(TP+TN)/N
    f1=(2*precision*recall)/(precision+recall)
    FPR=FP/(TN+FP)
    recall_list[i]<-recall
    precision_list[i]<-precision
    acc_list[i]<-acc
    f1_list[i]<-f1
    # print(acc)
    print(paste('precision',precision))
    print(paste('recall',recall))
    print(paste('combined_result',length(combined_result)))

   

}

upper = sapply(mean_list, "[[", 2)
results_df <- data.frame(
  mean_count = sapply(mean_list, "[[", 2),
  Precision = precision_list,
  Recall = recall_list,
  Accuracy = acc_list,
  f1 = f1_list
)
png(paste0(data_dir,parameter,"_manhattan_f1_1d_std0.png"), width = 800, height = 600, units = "px", res = 100)

# recalls <- results_df$Recall
recalls<-results_df$f1
midpoints <- sapply(mean_list, function(x) mean(x))
print(midpoints)
# 计算条形的宽度
bar_width <- 0.5  # 适当调整宽度

# 创建垂直条形图，每个条形图放在一个区间的中点
barplot(recalls, 
        names.arg = rep("", length(midpoints)),  # 空白标签
        col = "skyblue", 
        main = parameter, 
        xlab = mode, 
        ylab = "f1_score",
        ylim = c(0, 1.2),
        space = 0,
        width = rep(bar_width, length(midpoints)),  # 设置条形宽度
        axisnames = FALSE  # 不显示默认的坐标轴标签
)

# 在区间中点绘制折线图的点
points(
  midpoints,
  recalls,
  col = "red",
  pch = 16
)

# 绘制折线图
lines(
  midpoints,
  recalls,
  col = "red",
  type = "b",
  lwd = 2
)

# 添加标签
text(
  midpoints,
  recalls,
  labels = round(recalls, 2),
  pos = 3
)

# 手动设置 x 轴刻度
axis(1, at = midpoints, labels = sapply(mean_list, function(x) paste(x[1], "-", x[2])))

dev.off()  # 保存图形
###############################################################################################
# plot(
#   results_df$mean_count, results_df$Recall,
#   type = "l", col = "red", xlab = "mean_count", ylab = "sensitivity",
#   xlim = c(0,3.5), ylim = c(0, 1.2),
#   main = parameter
# )

###############################################################################################
# png('test.png', width = 800, height = 600, units = "px", res = 100)
# plot(
#   results_df$mean_count, results_df$Precision,
#   type = "l", col = "red", xlab = "mean_count", ylab = "sensitivity",
#   xlim = c(0,3.5), ylim = c(0, 1.2),
#   main = parameter
# )

# Adding points with shape icons
# points(
#   results_df$mean_count, results_df$Precision,
#   col = "red", pch = 16  # 16 represents a solid circle, you can choose a different value
# )

# lines(results_df$mean_count, results_df$Recall, col = "blue")
######################################################################################################
# points(results_df$mean_count, results_df$Recall, col = "blue", pch = 17)  # 17 represents a solid triangle
# text(
#   results_df$mean_count, results_df$Recall,
#   labels = round(results_df$Recall, 2),
#   pos = 3
# )
###############################################################################################
# lines(results_df$mean_count, results_df$Accuracy, col = "green")
# points(results_df$mean_count, results_df$Accuracy, col = "green", pch = 18)  # 18 represents a solid diamond

# Adding legend
# legend("bottomleft", legend = c("Precision", "Recall", "Accuracy"), col = c("red", "blue", "green"), lty = 1, pch = c(16, 17, 18))


# dev.off()  # Close the PNG device