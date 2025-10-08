#!/usr/bin/env Rscript

# ---------------------------------------------------------------
# 01.R - GSE109169 数据下载与标准化处理脚本
# 本脚本自动下载 GEO 系列 GSE109169，构建分组信息并输出标准化表达矩阵。
# ---------------------------------------------------------------

#==============================
# 工具函数：安装缺失的 R 包
#==============================
install_if_missing <- function(pkgs) {
  pkgs <- unique(pkgs)
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("[INFO] 检测到缺失包 %s，正在安装...", pkg))
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      tryCatch({
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      }, error = function(e) {
        stop(sprintf("[ERROR] 安装包 %s 失败：%s", pkg, e$message))
      })
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

#==============================
# 必要的包
#==============================
required_pkgs <- c("GEOquery", "limma", "annotate", "Biobase", "readr")
install_if_missing(required_pkgs)

#==============================
# 创建目录结构
#==============================
create_dir_if_missing <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("[INFO] 已创建目录：%s", path))
  }
}

create_dir_if_missing("data")
create_dir_if_missing("results")
create_dir_if_missing("figures")

#==============================
# 下载 GEO 数据（带重试）
#==============================
download_gse <- function(gse_id, retries = 3, sleep_sec = 5) {
  attempt <- 1
  while (attempt <= retries) {
    message(sprintf("[INFO] 正在尝试第 %d/%d 次下载 %s...", attempt, retries, gse_id))
    result <- tryCatch({
      GEOquery::getGEO(gse_id, GSEMatrix = TRUE)
    }, error = function(e) {
      message(sprintf("[WARN] 下载失败：%s", e$message))
      NULL
    })
    if (!is.null(result)) {
      message("[INFO] 下载成功。")
      return(result)
    }
    if (attempt < retries) {
      message(sprintf("[INFO] %d 秒后重试...", sleep_sec))
      Sys.sleep(sleep_sec)
    }
    attempt <- attempt + 1
  }
  stop(sprintf("[ERROR] 多次尝试后仍无法下载 %s，请检查网络或稍后重试。", gse_id))
}

#==============================
# 主流程
#==============================
message("[STEP] 1. 下载并读取 GSE109169 数据...")
geo_list <- download_gse("GSE109169")

if (length(geo_list) > 1) {
  message("[INFO] 检测到多个 ExpressionSet，默认使用第一个对象。")
}
eset <- geo_list[[1]]

message("[STEP] 2. 提取表达矩阵与样本信息...")
expr_raw <- Biobase::exprs(eset)
pheno <- Biobase::pData(eset)

if (is.null(expr_raw) || nrow(expr_raw) == 0) {
  stop("[ERROR] 表达矩阵为空，请检查 GEO 数据是否完整。")
}

#==============================
# 构建分组信息
#==============================
message("[STEP] 3. 构建 Tumor 与 Normal 分组...")
sample_ids <- colnames(expr_raw)
if (is.null(sample_ids)) {
  stop("[ERROR] 无法从表达矩阵提取样本名。")
}

group <- ifelse(grepl("-Tumor$", sample_ids, ignore.case = TRUE), "Tumor",
                ifelse(grepl("-Normal$", sample_ids, ignore.case = TRUE), "Normal", NA))

if (any(is.na(group))) {
  warning("[WARN] 存在无法判定分组的样本，请手动检查。")
}

pheno_table <- data.frame(
  sample_id = sample_ids,
  group = factor(group, levels = c("Tumor", "Normal")),
  pheno,
  check.names = FALSE
)

#==============================
# 数据标准化
#==============================
message("[STEP] 4. 对表达矩阵进行 log2 转换和分位数标准化...")
if (max(expr_raw, na.rm = TRUE) > 100) {
  message("[INFO] 检测到原始值可能未经 log2 处理，执行 log2 转换 (log2(x + 1))。")
  expr_log <- log2(expr_raw + 1)
} else {
  expr_log <- expr_raw
}

expr_norm <- limma::normalizeBetweenArrays(expr_log, method = "quantile")

#==============================
# 基础 QC 图形
#==============================
message("[STEP] 5. 生成基础 QC 图形...")

qc_boxplot_path <- file.path("figures", "boxplot_expression.png")
png(filename = qc_boxplot_path, width = 2400, height = 1600, res = 300)
boxplot(expr_norm, outline = FALSE, las = 2, main = "Normalized Expression Distribution",
        ylab = "Expression", col = "lightblue")
dev.off()
message(sprintf("[INFO] 箱线图已保存：%s", qc_boxplot_path))

qc_cor_path <- file.path("figures", "correlation_heatmap.png")
cor_mat <- stats::cor(expr_norm)
png(filename = qc_cor_path, width = 2400, height = 2000, res = 300)
stats::heatmap(cor_mat, symm = TRUE, main = "Sample Correlation Heatmap")
dev.off()
message(sprintf("[INFO] 样本相关性热图已保存：%s", qc_cor_path))

#==============================
# 保存结果
#==============================
message("[STEP] 6. 保存标准化表达矩阵与样本信息...")
expr_out_path <- file.path("data", "expr_matrix.rds")
pheno_out_path <- file.path("data", "pheno_table.csv")

saveRDS(expr_norm, file = expr_out_path)
readr::write_csv(pheno_table, file = pheno_out_path)

message(sprintf("[INFO] 表达矩阵已保存至：%s", expr_out_path))
message(sprintf("[INFO] 样本信息已保存至：%s", pheno_out_path))

#==============================
# 控制台输出关键日志
#==============================
message("[STEP] 7. 打印关键日志...")
print_dim <- dim(expr_norm)
message(sprintf("[LOG] 表达矩阵维度：%d 行 x %d 列", print_dim[1], print_dim[2]))
message("[LOG] 分组计数：")
print(table(pheno_table$group, useNA = "ifany"))

message("[INFO] 脚本执行完毕。")
