library(data.table)
library(ggplot2)
library(png)
library(rhdf5)
args <- commandArgs(trailingOnly = TRUE)
# print(args)

if (length(args) != 1) {
  stop("Usage: Rscript script.R <parameter>")
}

parameter <- args

print(paste0("处理参数:", parameter, "\n"))

if (grepl("mean", parameter)) {
  mode<-c('Meanraw')
} else {
  mode<-c('Varianceraw')
}

data_dir<-'/home/dataset/sim_data/cellcycle/250k/png/'
threshold<-0.10
png_path<-paste0(data_dir,parameter,'_',threshold,"_manhattan_f1.png")

test_res<-fread(paste0('/home/dataset/sim_data/cellcycle/250k/res/',parameter,'_manhattan.csv'))

label_file<-H5Fopen('/home/dataset/sim_data/cellcycle/250k/compartment250k_simulate.hdf5')

test_res[[mode]]<-h5read(label_file, paste0('compartment250k_',parameter,"/", mode))
max_value<-ceiling(max(abs(test_res[[mode]])))
test_list<-seq(0,max_value,0.5)
#取前n-1个元素和后n-1个元素，两两合并返回列表
mean_list <- mapply(function(x, y) c(x, y), head(test_list, -1), tail(test_list, -1), SIMPLIFY = FALSE)


length(mean_list)
search_term <- "wd|vd"

recall_list <- numeric(length = length(mean_list))
precision_list<-numeric(length = length(mean_list))
acc_list<-numeric(length = length(mean_list))
f1_list<-numeric(length = length(mean_list))



for (i in 1:length(mean_list)){
    ID_res<-test_res[abs(test_res[[mode]])>mean_list[[i]][1]&abs(test_res[[mode]])<=mean_list[[i]][2]]
    N<-length(ID_res$effect_size)  #均值区间内bin总数
    print(paste("bin num:",N))
    T_dc<-ID_res[ID_res$is_dc==1]$bin #实际加差异的bin
    T_ndc<-ID_res[ID_res$is_dc==0]$bin#实际无差异的Bin
    filtered_data <- subset(ID_res, effect_size>threshold)
    ID_dc<-filtered_data$bin #预测有差异的bin
    ID_ndc<-setdiff(ID_res$bin,ID_dc) #预测无差异的Bin
    #TP 预测有差异 实际也有差异  intersect取交集
    TP=length(intersect(ID_dc, T_dc))
    # 预测无差异，实际也无差异
    TN=length(intersect(ID_ndc, T_ndc))
    #FP 预测有差异 实际无差异
    FP=length(intersect(ID_dc,T_ndc))
    #fn 预测无差异 实际有差异  未识别出的差异bin的数量
    FN=length(intersect(ID_ndc, T_dc))

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
    print(paste('T_dc',length(T_dc)))
}

upper = sapply(mean_list, "[[", 2)
results_df <- data.frame(
  mean_count = sapply(mean_list, "[[", 2),
  Precision = precision_list,
  Recall = recall_list,
  Accuracy = acc_list,
  f1 = f1_list
)
png(png_path, width = 800, height = 600, units = "px", res = 100)


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