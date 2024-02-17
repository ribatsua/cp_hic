library(data.table)
moc <- function(P, Q) {
  NP <- length(P)
  NQ <- length(Q)

  if (NP == NQ & NP == 1) {
    moc <- 1
  } else {
    moc <- 0
    for (i in 1:NP) {
      for (j in 1:NQ) {
        #注意语言从1开始编号
        Pi <- seq(P[[i]][1], P[[i]][2]-1)
        Qj <- seq(Q[[j]][1], Q[[j]][2]-1)
        Fij <- intersect(Pi, Qj)
        moc <- moc + (length(Fij)^2 / (length(Pi) * length(Qj)))
      }
    }
    moc <- (moc - 1) / ((NP * NQ)^(1/2) - 1)
  }
  return(moc)
}

#获取满足条件的连续区间的起始和终止坐标以计算moc
get_continuous_intervals <- function(data, condition_column,comparison_operator,target_value) {
  # 获取所有满足条件的行号
  if (comparison_operator == "==") {
    selected_rows <- which(data[[condition_column]] == target_value)
  } else if (comparison_operator == ">") {
    selected_rows <- which(data[[condition_column]] > target_value)
  }

  # 如果没有满足条件的行，返回空列表
  if (length(selected_rows) == 0) {
    return(list())
  }

  # 初始化变量
  intervals <- list()

  # 获取连续区间的起始和终止行号
  start_row <- selected_rows[1]

  for (i in 2:length(selected_rows)) {
    if (selected_rows[i] != selected_rows[i - 1] + 1) {
      # 当不连续时，存储区间的起始和终止行号
      intervals <- append(intervals, list(c(start_row, selected_rows[i - 1])))
      # 更新起始行号
      start_row <- selected_rows[i]
    }
  }

  # 处理最后一个区间
  intervals <- append(intervals, list(c(start_row, selected_rows[length(selected_rows)])))

  return(intervals)
}

#测试
# parameter <- 'mean0.5'
# file_path <- paste0('/home/dataset/sim_data/cellcycle/1m/res/', parameter, '_manhattan.csv')

# # 读取CSV文件
# test_res <- fread(file_path)
# sim_data<-fread('/home/dataset/sim_data/cellcycle/1m/compartment_raw_0_1171.csv_mean4.csv')
# #预测的差异区间
# pre_dc<-get_continuous_intervals(test_res, "effect_size",">",0.16)
# #真实差异区间
# sim_data$is_dc <- sapply(sim_data$datalabel, function(x) ifelse(x %in% c('wd', 'vd'), 1, 0))

# real_dc<-get_continuous_intervals(sim_data, "is_dc","==",1)
# # print(sim_data$is_dc)
# result <- moc(real_dc,pre_dc)
# print(result)
# a<-real_dc

# b<-pre_dc


# print(real_dc)
# print(pre_dc)

# 测试
parameter <- 'mean4'
file_path <- paste0('/home/dataset/sim_data/cellcycle/250k/res/', parameter, '_manhattan.csv')

# 读取CSV文件
test_res <- fread(file_path)
result_intervals <- get_continuous_intervals(test_res, "is_dc","==",1)
pre_dc<-get_continuous_intervals(test_res, "effect_size",">",0.16)
# 输出结果
result <- moc(result_intervals,pre_dc)
print(result)
a<-result_intervals

b<-pre_dc
#########区间可视化


#  提取所有区间的起始点和终止点
all_points <- unlist(c(a, b))

# 获取图形的坐标范围
x_min <- min(all_points)
x_max <- max(all_points)
y_min <- 0
y_max <- 1

# 设置PNG文件路径
png_file_path <- "./overlap_plot.png"

# 创建一个空白的绘图区域，并保存为PNG
png(png_file_path, width = 800, height = 400)
plot(0, xlim = c(x_min, x_max), ylim = c(y_min, y_max), type = "n", xlab = "Position", ylab = "", main = "Overlap of Intervals")

# 添加第一组区间
for (i in 1:length(a)) {
  rect(a[[i]][1], 0.5, a[[i]][2], 0.6, col = "blue", border = NA)
}

# 添加第二组区间
for (i in 1:length(b)) {
  rect(b[[i]][1], 0.4, b[[i]][2], 0.5, col = "red", border = NA)
}

# 添加图例
legend("topright", legend = c("real_dc", "pre_dc"), fill = c("blue", "red"))

text(mean(c(max(sapply(a, max)), min(sapply(b, min)))), 0.8, labels = paste("MOC: ", result), col = "black", cex = 1.2)

# 关闭PNG设备
dev.off()

# 测试
# parameter <- 'mean0.5'
# file_path <- paste0('/home/dataset/sim_data/cellcycle/250k/res/', parameter, '_manhattan.csv')

# # 读取CSV文件
# test_res <- fread(file_path)
# result_intervals <- get_continuous_intervals(test_res, "is_dc","==",1)
# pre_dc<-get_continuous_intervals(test_res, "effect_size",">",0.16)
# # 输出结果
# result <- moc(result_intervals,pre_dc)
# print(result)





# # 测试
# a <- list(c(1, 5), c(6, 7),c(8,10))
# b <- list(c(1, 2), c(3, 7),c(8,10))
# result <- moc(a, b)
# print(result)
