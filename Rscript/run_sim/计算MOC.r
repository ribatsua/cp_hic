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



# 测试
parameter <- 'mean0.5'
file_path <- paste0('/home/dataset/sim_data/cellcycle/250k/res/', parameter, '_manhattan.csv')

# 读取CSV文件
test_res <- fread(file_path)
result_intervals <- get_continuous_intervals(test_res, "is_dc","==",1)
pre_dc<-get_continuous_intervals(test_res, "effect_size",">",0.16)
# 输出结果
# print(result_intervals)
# print(pre_dc)
result <- moc(result_intervals,pre_dc)
print(result)





# # 测试
# a <- list(c(1, 5), c(6, 7),c(8,10))
# b <- list(c(1, 2), c(3, 7),c(8,10))
# result <- moc(a, b)
# print(result)
