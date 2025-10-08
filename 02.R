#!/usr/bin/env Rscript

# FerrDb 基因整合脚本
# 目标：读取 Driver / Suppressor / Marker 三类基因，并合并输出去重表

suppressMessages({
  # 本脚本仅依赖 base R，无需额外包
})

# -------- 基础配置 --------
roles <- c("Driver", "Suppressor", "Marker")
local_files <- file.path("data", paste0("ferrdb_", tolower(roles), ".csv"))
names(local_files) <- roles

# FerrDb v1 legacy 下载链接（若本地缺失）
base_urls <- c(
  Driver = "https://www.zhounan.org/ferrdb/assets/fordownload/1_info_driver.csv",
  Suppressor = "https://www.zhounan.org/ferrdb/assets/fordownload/2_info_suppressor.csv",
  Marker = "https://www.zhounan.org/ferrdb/assets/fordownload/3_info_marker.csv"
)

# 若 data 目录不存在则创建（避免后续写入失败）
if (!dir.exists("data")) {
  dir.create("data", recursive = TRUE, showWarnings = FALSE)
}

read_role_table <- function(role) {
  local_path <- local_files[[role]]
  url <- base_urls[[role]]

  message(sprintf("处理 %s 基因表...", role))

  if (!file.exists(local_path)) {
    message(sprintf("未找到本地文件 %s，尝试从网络下载。", local_path))
    ok <- FALSE
    try({
      utils::download.file(url, destfile = local_path, mode = "wb", quiet = TRUE)
      ok <- file.exists(local_path)
    }, silent = TRUE)

    if (!ok) {
      stop(sprintf(
        paste0(
          "无法获取 %s 数据。请将 FerrDb legacy 页面的 CSV 手动下载为 %s 后重试。\n",
          "下载链接：%s"
        ),
        role, local_path, url
      ), call. = FALSE)
    }
  }

  df <- tryCatch(
    utils::read.csv(local_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      stop(sprintf("读取 %s 失败：%s", local_path, e$message), call. = FALSE)
    }
  )

  # 寻找基因符号列，优先匹配标准命名
  candidate_cols <- c("Symbol", "symbol", "Gene", "Gene.Symbol", "Gene_symbol", "Gene Symbol")
  symbol_col <- candidate_cols[candidate_cols %in% names(df)]
  if (length(symbol_col) == 0) {
    stop(sprintf("文件 %s 中未找到基因符号列，请检查格式。", local_path), call. = FALSE)
  }

  symbol_values <- df[[symbol_col[1]]]
  symbol_values <- trimws(as.character(symbol_values))
  symbol_values <- symbol_values[symbol_values != ""]

  out <- data.frame(Symbol = symbol_values, Role = role, stringsAsFactors = FALSE)
  return(out)
}

# 逐类读取并合并
merged_list <- lapply(roles, read_role_table)
merged_df <- do.call(rbind, merged_list)

# 去重并保留原始列顺序
merged_unique <- unique(merged_df[, c("Symbol", "Role")])

# 输出统计信息
counts <- table(merged_unique$Role)
for (role in roles) {
  role_count <- if (role %in% names(counts)) counts[[role]] else 0
  message(sprintf("%s 基因数：%d", role, role_count))
}
message(sprintf("去重后总数：%d", nrow(merged_unique)))

# 验证数量合理
stopifnot(nrow(merged_unique) > 200)

# 保存结果
output_path <- "data/ferrdb_genes_merged.csv"
utils::write.csv(merged_unique, output_path, row.names = FALSE)
message(sprintf("已保存合并结果至 %s", output_path))
