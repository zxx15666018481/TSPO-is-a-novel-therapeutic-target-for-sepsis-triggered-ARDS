# ====================== 1. 设置文件路径 ======================
count_file <- "C:/Users/Administrator/Desktop/1/count.txt"   # count数据文件（第一列是基因ID，列是样本）
length_file <- "C:/Users/Administrator/Desktop/1/length.txt" # 基因长度文件（id列是基因ID，eff_length是有效长度）
output_tpm <- "C:/Users/Administrator/Desktop/1/tpm.txt"     # 输出TPM结果
output_log2tpm <- "C:/Users/Administrator/Desktop/1/log2tpm.txt" # 输出log2(TPM+1)结果


# ====================== 2. 读取数据 ======================
# 读取count数据（假设制表符分隔，表头存在）
count_df <- read.table(
  count_file, 
  header = TRUE, 
  sep = "\t",        # 如果是逗号分隔，改为sep=","
  stringsAsFactors = FALSE
)

# 读取基因长度数据（假设制表符分隔，表头存在）
length_df <- read.table(
  length_file, 
  header = TRUE, 
  sep = "\t",        # 如果是逗号分隔，改为sep=","
  stringsAsFactors = FALSE
)


# ====================== 3. 数据检查与合并 ======================
# 检查基因长度是否合法（避免除以0）
if (any(length_df$eff_length <= 0)) {
  stop("发现有效长度为0或负数的基因，请检查length文件：\n", 
       paste(length_df$id[length_df$eff_length <= 0], collapse = ", "))
}

# 合并count和长度数据（按基因ID匹配，假设count的基因ID列是"id"，length的基因ID列是"id"；若列名不同，需修改by参数）
merged_df <- merge(
  count_df, 
  length_df, 
  by = "id"  # 若count的基因列是"gene_id"，length是"id"，则改为 by.x="gene_id", by.y="id"
)


# ====================== 4. 计算TPM ======================
# 提取样本列（排除基因ID和长度列）
sample_cols <- setdiff(colnames(merged_df), c("id", "eff_length"))

# 检查是否有样本列
if (length(sample_cols) == 0) {
  stop("未找到样本列，请检查count文件的列名！")
}

# 初始化TPM结果表
tpm_df <- data.frame(
  id = merged_df$id, 
  eff_length = merged_df$eff_length, 
  stringsAsFactors = FALSE
)

# 逐样本计算TPM
for (col in sample_cols) {
  # 1. 提取该样本的count值
  count_vals <- merged_df[[col]]
  
  # 2. 转换基因长度为kb（eff_length单位是bp，除以1000）
  length_kb <- merged_df$eff_length / 1000
  
  # 3. 计算“每千碱基的count”（count / 长度(kb)）
  per_kb <- count_vals / length_kb
  
  # 4. 计算该样本的总“每千碱基count”（用于标准化）
  sum_per_kb <- sum(per_kb, na.rm = TRUE)  # 忽略可能的NA（若有）
  
  # 5. 计算TPM：(per_kb / 总per_kb) * 1e6
  tpm_vals <- (per_kb / sum_per_kb) * 1e6
  
  # 6. 避免极端值（理论上sum_per_kb不会为0，但若全0则TPM全为0）
  if (sum_per_kb == 0) {
    warning(paste("样本", col, "的count全为0，TPM结果全为0！"))
  }
  
  # 7. 保存到结果表
  tpm_df[[col]] <- tpm_vals
}


# ====================== 5. 计算log2(TPM + 1) ======================
log2tpm_df <- data.frame(id = tpm_df$id, stringsAsFactors = FALSE)
for (col in sample_cols) {
  # TPM + 1 避免log2(0)，然后取log2
  log2tpm_df[[col]] <- log2(tpm_df[[col]] + 1)
}


# ====================== 6. 保存结果 ======================
# 保存TPM
write.table(
  tpm_df, 
  output_tpm, 
  sep = "\t",      # 制表符分隔，可改为","（csv格式）
  quote = FALSE,   # 不添加引号
  row.names = FALSE# 不保存行名
)

# 保存log2(TPM + 1)
write.table(
  log2tpm_df, 
  output_log2tpm, 
  sep = "\t",      
  quote = FALSE,   
  row.names = FALSE
)

# 提示完成
cat("转换成功！\n", 
    "TPM结果保存到：", output_tpm, "\n", 
    "log2(TPM+1)结果保存到：", output_log2tpm, "\n")