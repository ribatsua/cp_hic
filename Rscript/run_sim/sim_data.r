# The function simulate_data is designed to simulate data based on a given CSV file. It takes four parameters: file_path, r_mean, r_var, and sample_N.
#
# file_path is the path to the input CSV file.
# r_mean (default value = 1) represents the multiple by which the mean will be changed during simulation.
# r_var (default value = 1) represents the multiple by which the variance will be changed during simulation.
# sample_N (default value = 1171) represents the number of samples to be generated for each row of data.
# The main steps of the function are as follows:
#
# Read the data from the CSV file using read.csv.
# Perform parameter estimation for each row of data by calculating the mean and variance using mean and sd functions.
# Simulate changes in mean and variance based on the given r_mean and r_var values.
# Create labels to indicate data partitions and whether mean or variance has been changed.
# Generate a new data frame with simulated mean and variance values.
# Set a random seed and sample data from a normal distribution using rnorm.
# Save the sampled data to a new CSV file.
# The output of the function is the path to the output file.
# It's important to note that the function uses specific labels such as 'nd', 'wd', and 'vd' to mark data partitions and indicate whether mean or variance has been changed.

Simulate <- function(file_path, r_mean = 1, r_var=1, sample_N = 1171) {
  if (missing(file_path)) {
    stop("Error: file path is missing.")
  }
  if(is.null(r_mean)){
    r_mean <- 1
  }
  if(is.null(r_var)){
    r_var <- 1
  }
  if(is.null(sample_N)){
    sample_N <- 1171
  }
  # parameter estimation
  data <- read.csv(file_path)
  param_estimates <- data.frame(matrix(ncol = 3, nrow = nrow(data)))
  colnames(param_estimates) <- c("bin", "Meanraw", "Varianceraw")
  
  param_estimates$bin <- seq(0, nrow(param_estimates) - 1)
  
  for (i in 1:nrow(data)) {
    row_data <- as.numeric(unlist(data[i, -1]))
    param_estimates[i, "Meanraw"] <- mean(row_data)
    param_estimates[i, "Varianceraw"] <- sd(row_data)
  }
  
  # parameter simulaton
  mean_data <- param_estimates[, "Meanraw"]
  var_data <- param_estimates[, "Varianceraw"]
  
  # Create empty labels vector and partition the data
  labelsAB <- rep("", length(mean_data))
  label_block <- rep("", length(mean_data))
  labels <- rep("nd", length(mean_data))
  labelsvar <- rep("", length(var_data))
  
  labelsAB <- ifelse(mean_data > 0, "1", "0")
  label_block[1] <- "0"
  mean_n <- 0
  
  for (i in 2:length(mean_data)) {
    product <- mean_data[i] * mean_data[i-1]
    if (product > 0) {
      label_block[i] <- label_block[i-1]
    } else {
      mean_n <- mean_n + 1
      label_block[i] <- floor(mean_n)
    }
  }
  
  for (i in 1:length(labelsAB)) {
    if (labelsAB[i] == "1") {
      max_value_A <- label_block[i]
    }
    if (labelsAB[i] == "0") {
      max_value_B <- label_block[i]
    }
  }
  
  # mean change
  set.seed(123)
  # 按labelsAB和label_block标签值相同的分组计算最大值和最小值
  grouped_max <- tapply(mean_data[labelsAB==1], list(labelsAB[labelsAB==1], label_block[labelsAB==1]), max)
  grouped_min <- tapply(mean_data[labelsAB==0], list(labelsAB[labelsAB==0], label_block[labelsAB==0]), min)
  
  # 计算labels3的值
  labels3 <- rep(1, length(mean_data))
  num_LA <- 0
  num_LB <- 0
  max_val <- max(mean_data)
  min_val <- min(mean_data)
  # 创建一个新的向量进行标记强弱区,1-S强，0-W弱
  for (i in 1:length(mean_data)) {
    curr_label1 <- labelsAB[i]
    curr_label2 <- label_block[i]
    if (curr_label2 != 0) {
      if (curr_label1 == 1) {
        if (mean_data[i] == grouped_max[curr_label1, curr_label2] & mean_data[i] > max_val*2/3) {
          labels3[labelsAB==curr_label1 & label_block==curr_label2] <- 2
          num_LA = num_LA + 1
        }
        if (mean_data[i] == grouped_max[curr_label1, curr_label2] & mean_data[i] < max_val/3) {
          labels3[labelsAB==curr_label1 & label_block==curr_label2] <- 0
        }
      } else {
        if (abs(mean_data[i]) == abs(grouped_min[curr_label1, curr_label2]) & abs(mean_data[i]) > abs(min_val)*2/3) {
          labels3[labelsAB==curr_label1 & label_block==curr_label2] <- 2
          num_LB = num_LB + 1
        }
        if (abs(mean_data[i]) == abs(grouped_min[curr_label1, curr_label2]) & abs(mean_data[i]) < abs(min_val)/3) {
          labels3[labelsAB==curr_label1 & label_block==curr_label2] <- 0
        }
      }
    }
  }
  
  
  # 循环遍历label1和label_block
  for (i in 1:length(labelsAB)) {
    if (labelsAB[i] == "0") {
      # 如果label1[i]==0，则更新max_value_B为label_block[i]
      max_value_B <-label_block[i]
    }
  }
  # 循环遍历label1和label_block
  for (i in 1:length(labelsAB)) {
    if (labelsAB[i] == "1") {
      # 如果label1[i]==1，则更新max_value_A为label_block[i]
      max_value_A <- label_block[i]
    }
    if (labelsAB[i] == "0") {
      # 如果label1[i]==0，则更新max_value_B为label_block[i]
      max_value_B <- label_block[i]
    }
  }
  
  # 设置随机种子
  set.seed(123) 
  
  # 获取满足条件的label_block取值
  possible_values_LA <- label_block[labels3 == "2"&labelsAB == "1" ]
  possible_values_LB <- label_block[labels3 == "2"&labelsAB == "0" ]
  
  # 随机挑选数值
  random_values_LA <- sample(possible_values_LA, ceiling(as.numeric(num_LA) /2))
  random_values_LB <- sample(possible_values_LB, ceiling(as.numeric(num_LB) /2))

  
  # 从（0~N）中随机抽取M个数
  N_A <- max_value_A  # 上限值
  M_A <- ceiling(as.numeric(N_A) %/%8) # 需要抽取的个数
  
  # 从（0~N）中随机抽取M个数
  N_B <- max_value_B  # 上限值
  M_B<- ceiling(as.numeric(N_B) %/%8) # 需要抽取的个数
  

  extra_nums <- sample(random_values_LA,ceiling(as.numeric(num_LA) %/%4))
  extra_nums1 <- sample(random_values_LB,ceiling(as.numeric(num_LB) %/%4))

  
  sampled_nums <- sample(0:N_A, M_A)
  sampled_nums1 <- sample(0:N_B,M_B)
  extra_nums2 <- sample(sampled_nums,ceiling(as.numeric(M_A/2)))
  extra_nums3 <- sample(sampled_nums1,ceiling(as.numeric(M_B/2)))
  mean_data_new <- mean_data  # 创建新的列数据
  var_data_new <- var_data 
  if(r_mean!=1&&r_var==1){
    for (i in 1:length(mean_data)) {
    if (label_block[i] %in% random_values_LA||label_block[i] %in% sampled_nums&&labelsAB[i]==1) {
      if(!(label_block[i] %in% extra_nums||label_block[i] %in% extra_nums2))
      {
        mean_data_new[i] <- mean_data[i] * r_mean*-1
        labels[i] <- 'wd'
      }
      else{
        mean_data_new[i] <- mean_data[i] * r_mean
        labels[i] <- 'wd'
      }
    }
    if (label_block[i] %in% random_values_LB||label_block[i] %in% sampled_nums1&&labelsAB[i]==0){
      if(!(label_block[i] %in% extra_nums1||label_block[i] %in% extra_nums3))
      {
        mean_data_new[i] <- mean_data[i] * r_mean*-1
        labels[i] <- 'wd'
      }
      else{
        mean_data_new[i] <- mean_data[i] * r_mean
        labels[i] <- 'wd'
      }
    }
  }
  }
  if(r_var!=1&&r_mean==1){
   for (i in 1:length(mean_data)) {
    if (label_block[i] %in% random_values_LA||label_block[i] %in% sampled_nums&&labelsAB[i]==1) {
        var_data_new[i] <- var_data[i] * r_var
        labels[i] <- 'vd'
    }
    if (label_block[i] %in% random_values_LB||label_block[i] %in% sampled_nums1&&labelsAB[i]==0){
        var_data_new[i] <- var_data[i] * r_var
        labels[i] <- 'vd'
    }
  }  
  }
 
  data <- cbind(data, Meannew = mean_data_new,Variancenew=var_data_new,datalabel = labels)
  # sampling
  sampled_data_df <- data.frame()
  for (i in 1:nrow(data)) {
    mu_param <-data[i, "Meannew"]
    sigma_param <- data[i, "Variancenew"]
    set.seed(1407)
    sampled_values <- rnorm(sample_N, mean = mu_param, sd = sigma_param)
    sampled_data_df <- rbind(sampled_data_df, sampled_values)
  }
  
  colnames(sampled_data_df) <- paste("cell_", 0:(ncol(sampled_data_df)-1), sep = "")
  sampled_data_df$bin <- 0:(nrow(sampled_data_df)-1)
  sampled_data_df <- sampled_data_df[, c("bin", setdiff(names(sampled_data_df), "bin"))]
  sampled_data_df <- cbind(sampled_data_df,Meanraw = mean_data,Varianceraw=var_data,Meannew = mean_data_new,Variancenew=var_data_new,datalabel = labels)
  # if(r_mean==1&&r_var==1){
  #   output_file <- paste0(file_path, "_mean", r_mean, "+var", r_var, ".csv")
  #   }
  if(r_mean==1){
    if(r_var!=1){
      output_file <- paste0(file_path, "_var", r_var, ".csv")
    }
    else{
      output_file <- paste0(file_path, "_mean", r_mean, "+var", r_var, ".csv")
    }
  } 
  else{
    if(r_var==1){
      output_file <- paste0(file_path, "_mean", r_mean, ".csv")
    }
    else{
    output_file <- paste0(file_path, "_mean", r_mean, "+var", r_var, ".csv")
    }
  }
   
  write.csv(sampled_data_df, file = output_file, row.names = FALSE)
  return(output_file)
}
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",1,1)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",0.25,1)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",0.3,1)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",0.5,1)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",2,1)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",3,1)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",4,1)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",5,1)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",6,1)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",7,1)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",8,1)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",9,1)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",10,1)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",1,0.25)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",1,0.3)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",1,0.5)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",1,2)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",1,3)
# Simulate("D:/Desktop/data12.15/data12.15/compartment_raw_0_1171.csv",1,4)


