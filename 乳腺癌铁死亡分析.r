# =============================================================
# Project: Breast Cancer Ferroptosis (Bulk RNA-seq + WGCNA + Cox)
# Author: (Your Name)
# Date: 2025-09-03 (Asia/Tokyo)
# R >= 4.3 recommended
# =============================================================

# ---------------------------
# 0) Setup & Utilities
# ---------------------------
options(stringsAsFactors = FALSE)
set.seed(20250903)

# Check and install required packages. When possible, catch errors and
# continue gracefully. The purrr package is required for map/map_dfr.
required_packages <- c(
  "dplyr", "tibble", "stringr", "readr", "data.table", "ggplot2",
  "BiocParallel", "limma", "edgeR", "DESeq2", "WGCNA",
  "clusterProfiler", "org.Hs.eg.db", "msigdbr", "GSVA",
  "SummarizedExperiment", "TCGAbiolinks", "survival", "survminer",
  "glmnet", "timeROC", "pROC", "rms", "GEOquery", "hgu133plus2.db",
  "AnnotationDbi", "janitor", "broom", "purrr"
)

optional_packages <- c("STRINGdb", "igraph")

install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      try(BiocManager::install(pkg), silent = TRUE)
    }
  }
}

install_if_missing(required_packages)

# Load required libraries; suppress messages for clarity
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(data.table)
  library(ggplot2)
  library(BiocParallel)
  library(limma)
  library(edgeR)
  library(DESeq2)
  library(WGCNA)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(GSVA)
  library(SummarizedExperiment)
  library(TCGAbiolinks)
  library(survival)
  library(survminer)
  library(glmnet)
  library(timeROC)
  library(pROC)
  library(rms)
  library(GEOquery)
  library(AnnotationDbi)
  library(janitor)
  library(broom)
  library(purrr)
})

# Optional packages; load quietly
suppressWarnings({
  try(library(STRINGdb), silent = TRUE)
  try(library(igraph), silent = TRUE)
  try(library(hgu133plus2.db), silent = TRUE)
})

# Allow parallel threads for WGCNA if available
allowWGCNAThreads()

# Set up project directories
proj_dir <- "~/BrCa_Ferroptosis_Project"
dir.create(proj_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path(proj_dir, "figures")
res_dir <- file.path(proj_dir, "results")
cache_dir <- file.path(proj_dir, "cache")
for (d in c(fig_dir, res_dir, cache_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Helper functions for saving files
safe_saveRDS <- function(obj, path) {
  saveRDS(obj, path)
  message("Saved: ", path)
}

safe_write <- function(df, path) {
  data.table::fwrite(df, path)
  message("Saved: ", path)
}

# Helper: convert counts to TPM, given gene lengths in kb
counts_to_tpm <- function(counts, gene_length_kb) {
  rpk <- sweep(counts, 1, gene_length_kb, FUN = "/")
  tpm <- sweep(rpk, 2, colSums(rpk), FUN = "/") * 1e6
  return(tpm)
}

# Helper: collapse duplicated gene symbols by keeping the row with highest mean expression
collapse_by_symbol <- function(expr, symbols) {
  stopifnot(nrow(expr) == length(symbols))
  valid_idx <- !is.na(symbols) & symbols != ""
  expr <- expr[valid_idx, , drop = FALSE]
  symbols <- symbols[valid_idx]
  dt <- as.data.table(expr)
  dt[, symbol := symbols]
  expr_cols <- setdiff(colnames(dt), "symbol")
  dt[, meanExpr := rowMeans(.SD, na.rm = TRUE), .SDcols = expr_cols]
  setorder(dt, symbol, -meanExpr)
  dt_unique <- dt[!duplicated(symbol)]
  mat <- as.matrix(dt_unique[, ..expr_cols])
  rownames(mat) <- dt_unique$symbol
  return(mat)
}

# Helper: detect OS time/status in GEO phenotype data
extract_surv_from_pdata <- function(pdat) {
  nm <- tolower(colnames(pdat))
  time_patterns <- c(
    "os_months", "os_time", "ostime", "futime", "time",
    "months_to_last_followup", "followup_months", "t.rfs",
    "rfs_time", "dfs_time", "overall_survival_months"
  )
  status_patterns <- c(
    "os_event", "os_status", "vital_status", "status",
    "event", "rfs_event", "dfs_event", "death_event"
  )
  time_idx <- which(nm %in% time_patterns)
  if (length(time_idx) == 0) {
    time_idx <- grep("(os|overall|survival|futime|time|follow).*(month|day|time)", nm)
  }
  status_idx <- which(nm %in% status_patterns)
  if (length(status_idx) == 0) {
    status_idx <- grep("(os|overall|survival|vital|status|event)", nm)
  }
  if (length(time_idx) > 0) {
    time_vec <- suppressWarnings(as.numeric(as.character(pdat[[time_idx[1]]])))
  } else {
    time_vec <- rep(NA_real_, nrow(pdat))
  }
  if (length(status_idx) > 0) {
    status_raw <- as.character(pdat[[status_idx[1]]])
    status_num <- rep(NA_real_, length(status_raw))
    status_num[grepl("dead|deceased|1|event|yes|death", tolower(status_raw))] <- 1
    status_num[grepl("alive|0|censored|no|living", tolower(status_raw))] <- 0
    if (all(is.na(status_num))) {
      status_num <- suppressWarnings(as.numeric(status_raw))
    }
  } else {
    status_num <- rep(NA_real_, nrow(pdat))
  }
  list(
    time = time_vec,
    status = status_num,
    status_note = if (length(status_idx) > 0) colnames(pdat)[status_idx[1]] else "not_found"
  )
}

# ---------------------------
# 1) TCGA-BRCA: download & prep
# ---------------------------
message("Step 1: Downloading TCGA-BRCA data...")

# Cached file paths
tcga_counts_rds <- file.path(cache_dir, "tcga_brca_counts.rds")
tcga_clin_rds   <- file.path(cache_dir, "tcga_brca_clinical.rds")

if (!file.exists(tcga_counts_rds) || !file.exists(tcga_clin_rds)) {
  tryCatch({
    query <- GDCquery(
      project = "TCGA-BRCA",
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    GDCdownload(query, method = "api")
    se <- GDCprepare(query)
    safe_saveRDS(se, tcga_counts_rds)
    clin <- colData(se) %>% as.data.frame()
    safe_saveRDS(clin, tcga_clin_rds)
  }, error = function(e) {
    message("Error downloading TCGA data: ", e$message)
    stop("Please check internet connection and try again")
  })
} else {
  message("Loading cached TCGA data...")
  se <- readRDS(tcga_counts_rds)
  clin <- readRDS(tcga_clin_rds)
}

# Extract counts matrix and gene info
count_mat <- assay(se, "unstranded")
gene_info <- rowData(se) %>% as.data.frame()
ens_ids <- gene_info$gene_id
if (is.null(ens_ids)) {
  ens_ids <- rownames(count_mat)
}
ens_clean <- gsub("\\..*$", "", ens_ids)

# Map Ensembl IDs to gene symbols
symbols <- mapIds(org.Hs.eg.db,
                  keys = ens_clean,
                  keytype = "ENSEMBL",
                  column = "SYMBOL",
                  multiVals = "first")

valid_idx <- !is.na(symbols) & symbols != ""
count_mat_filt <- count_mat[valid_idx, ]
symbols_filt <- symbols[valid_idx]
count_mat_collapsed <- collapse_by_symbol(count_mat_filt, symbols_filt)

# Sample information
sample_info <- colData(se) %>% as.data.frame()
sample_type <- sample_info$shortLetterCode
names(sample_type) <- sample_info$barcode

# Match sample names
common_samples <- intersect(names(sample_type), colnames(count_mat_collapsed))
count_mat_final <- count_mat_collapsed[, common_samples]
sample_type_final <- sample_type[common_samples]
clin_sub <- sample_info[rownames(sample_info) %in% common_samples, ]

# Keep tumor (TP) and normal (NT) samples
keep_samples <- sample_type_final %in% c("TP", "NT")
count_mat_final <- count_mat_final[, keep_samples]
sample_type_final <- sample_type_final[keep_samples]
clin_sub <- clin_sub[rownames(clin_sub) %in% names(sample_type_final), ]

message("Final matrix: ", nrow(count_mat_final), " genes x ", ncol(count_mat_final), " samples")

# ---------------------------
# 2) Differential Expression Analysis
# ---------------------------
message("Step 2: Differential expression analysis...")

# Construct coldata
coldata <- data.frame(
  sample = colnames(count_mat_final),
  group = factor(ifelse(sample_type_final == "TP", "Tumor", "Normal")),
  row.names = colnames(count_mat_final)
)

# Filter low-expressed genes
dge <- DGEList(counts = count_mat_final)
dge <- calcNormFactors(dge)
keep_genes <- filterByExpr(dge, group = coldata$group, min.count = 10, min.total.count = 15)
count_filt <- count_mat_final[keep_genes, ]
message("After filtering: ", nrow(count_filt), " genes retained")

# Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = count_filt,
  colData = coldata,
  design = ~ group
)
# Pre-filter
dds <- dds[rowSums(counts(dds)) >= 10, ]
# DE analysis
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "Tumor", "Normal"))

# Organize results
res_df <- as.data.frame(res) %>%
  rownames_to_column("SYMBOL") %>%
  arrange(padj) %>%
  mutate(
    sig = !is.na(padj) & abs(log2FoldChange) >= 1 & padj < 0.05,
    direction = case_when(
      sig & log2FoldChange > 0 ~ "Up",
      sig & log2FoldChange < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )

safe_write(res_df, file.path(res_dir, "TCGA_BRCA_DE_Tumor_vs_Normal.csv"))

# VST transformation
vst_mat <- assay(vst(dds, blind = FALSE))

# Select top variable genes for WGCNA
gene_vars <- apply(vst_mat, 1, var, na.rm = TRUE)
top_var_genes <- head(order(gene_vars, decreasing = TRUE), 15000)
vst_top <- vst_mat[top_var_genes, ]

message("DE analysis completed: ", sum(res_df$sig, na.rm = TRUE), " significant genes")

# ---------------------------
# 3) Ferroptosis gene sets
# ---------------------------
message("Step 3: Loading ferroptosis gene sets...")

fer_v1_path <- file.path(proj_dir, "FerrDb_v1_259.csv")
fer_v2_path <- file.path(proj_dir, "FerrDb_v2_621.csv")

# Attempt to load FerrDb v1
if (file.exists(fer_v1_path)) {
  fer_symbols <- tryCatch({
    fer_v1 <- readr::read_csv(fer_v1_path, show_col_types = FALSE) %>%
      janitor::clean_names()
    if ("symbol" %in% colnames(fer_v1)) {
      unique(na.omit(fer_v1$symbol))
    } else {
      stop("Column 'symbol' not found in FerrDb v1 file")
    }
  }, error = function(e) {
    message("Error loading FerrDb v1: ", e$message)
    NULL
  })
} else {
  message("FerrDb v1 file not found at: ", fer_v1_path)
  fer_symbols <- NULL
}

# If FerrDb v1 load failed or has too few genes, fallback to MSigDB
if (is.null(fer_symbols) || length(fer_symbols) < 50) {
  message("Using MSigDB ferroptosis gene sets...")
  fer_symbols <- tryCatch({
    msig <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
    fer_sets <- msig %>% filter(grepl("FERROPTOSIS", gs_name, ignore.case = TRUE))
    unique(fer_sets$gene_symbol)
  }, error = function(e) {
    message("Error loading from MSigDB: ", e$message)
    c("ACSL4", "GPX4", "SLC7A11", "TFRC", "FTH1", "PTGS2",
      "ALOX15", "NOX1", "LPCAT3", "AIFM2")
  })
}

# Optionally load FerrDb v2 (unused but available)
if (file.exists(fer_v2_path)) {
  tryCatch({
    fer2 <- readr::read_csv(fer_v2_path, show_col_types = FALSE) %>%
      janitor::clean_names()
    message("Loaded FerrDb v2: ", nrow(fer2), " entries")
  }, error = function(e) {
    message("Error loading FerrDb v2: ", e$message)
  })
}

# ---------------------------
# 4) Ferroptosis-related DEGs & enrichment
# ---------------------------
message("Step 4: Analyzing ferroptosis-related DEGs...")

# Intersect DE results with ferroptosis genes
fe_deg <- res_df %>%
  filter(SYMBOL %in% fer_symbols) %>%
  arrange(padj)

safe_write(fe_deg, file.path(res_dir, "Fe_related_DEGs_TCGA_BRCA.csv"))
message("Ferroptosis DEGs found: ", nrow(fe_deg))

# Identify up- and down-regulated ferroptosis genes
fe_up <- fe_deg %>% filter(sig & log2FoldChange > 0) %>% pull(SYMBOL) %>% unique()
fe_down <- fe_deg %>% filter(sig & log2FoldChange < 0) %>% pull(SYMBOL) %>% unique()
message("Ferroptosis genes - Up: ", length(fe_up), ", Down: ", length(fe_down))

# Perform enrichment analysis if enough genes
if (length(c(fe_up, fe_down)) >= 5) {
  tryCatch({
    all_fe_genes <- unique(c(fe_up, fe_down))
    sym2ent <- bitr(all_fe_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    if (length(fe_up) >= 5) {
      up_entrez <- sym2ent$ENTREZID[sym2ent$SYMBOL %in% fe_up]
      ego_up <- enrichGO(
        gene = up_entrez,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID", ont = "BP",
        pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE
      )
      if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
        safe_write(as.data.frame(ego_up), file.path(res_dir, "GO_BP_up_in_Ferro_DEGs.csv"))
      }
    }
    if (length(fe_down) >= 5) {
      down_entrez <- sym2ent$ENTREZID[sym2ent$SYMBOL %in% fe_down]
      ego_down <- enrichGO(
        gene = down_entrez,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID", ont = "BP",
        pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE
      )
      if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
        safe_write(as.data.frame(ego_down), file.path(res_dir, "GO_BP_down_in_Ferro_DEGs.csv"))
      }
    }
  }, error = function(e) {
    message("Error in enrichment analysis: ", e$message)
  })
}

# ---------------------------
# 5) Ferroptosis score & WGCNA
# ---------------------------
message("Step 5: Calculating ferroptosis score and running WGCNA...")

# Identify ferroptosis genes present in expression data
fer_in_vst <- intersect(fer_symbols, rownames(vst_top))
message("Ferroptosis genes in expression data: ", length(fer_in_vst))

if (length(fer_in_vst) < 5) {
  stop("Too few ferroptosis genes present in expression data. Please check gene sets.")
}

# Compute ssGSEA-based ferroptosis score
fer_score <- tryCatch({
  fer_ssgsea <- gsva(
    vst_top,
    list(FERROPTOSIS = fer_in_vst),
    method = "ssgsea", kcdf = "Gaussian", parallel.sz = 1
  )
  as.numeric(fer_ssgsea[1, ])
}, error = function(e) {
  message("Error calculating ferroptosis score via GSVA: ", e$message)
  colMeans(vst_top[fer_in_vst, ], na.rm = TRUE)
})
names(fer_score) <- colnames(vst_top)

# Prepare WGCNA data (samples x genes)
wgcna_dat <- t(vst_top)
message("WGCNA input: ", nrow(wgcna_dat), " samples x ", ncol(wgcna_dat), " genes")

# Check sample and gene quality
gsg <- goodSamplesGenes(wgcna_dat, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    wgcna_dat <- wgcna_dat[, gsg$goodGenes]
  }
  if (sum(!gsg$goodSamples) > 0) {
    wgcna_dat <- wgcna_dat[gsg$goodSamples, ]
    fer_score <- fer_score[gsg$goodSamples]
  }
}

# Sample clustering to detect outliers
sampleTree <- hclust(dist(wgcna_dat), method = "average")
pdf(file.path(fig_dir, "WGCNA_SampleClustering.pdf"), width = 10, height = 6)
par(mfrow = c(1, 1))
plot(sampleTree, main = "Sample clustering to detect outliers",
     xlab = "", sub = "", cex = 0.7)
dev.off()

# Choose soft-thresholding power
message("Choosing soft-thresholding power...")
powers <- c(1:20)
sft <- pickSoftThreshold(wgcna_dat, powerVector = powers,
                         networkType = "signed", verbose = 0)

# Plot soft-threshold selection
pdf(file.path(fig_dir, "WGCNA_SoftThreshold.pdf"), width = 12, height = 5)
par(mfrow = c(1, 2))
cex1 <- 0.9
# Scale-free topology fit index
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.80, col = "red", lty = 2)
# Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")
dev.off()

# Determine soft threshold
sft_power <- sft$powerEstimate
if (is.na(sft_power)) {
  good_powers <- which(-sign(sft$fitIndices[,3]) * sft$fitIndices[,2] >= 0.80)
  if (length(good_powers) > 0) {
    sft_power <- min(sft$fitIndices[good_powers, 1])
  } else {
    sft_power <- 6
  }
}
message("Selected soft-thresholding power: ", sft_power)

# Build co-expression network
message("Building co-expression network...")
net <- blockwiseModules(
  wgcna_dat,
  power = sft_power,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  deepSplit = 2,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3,
  maxBlockSize = 5000
)

MEs <- net$MEs
moduleColors <- labels2colors(net$colors)
MEs_col <- orderMEs(MEs)
message("Identified ", ncol(MEs_col), " modules")

# Module–trait correlations
trait_df <- data.frame(FerroScore = fer_score[rownames(wgcna_dat)])
modTraitCor <- cor(MEs_col, trait_df$FerroScore, use = "pairwise.complete.obs")
modTraitP <- corPvalueStudent(modTraitCor, nrow(wgcna_dat))

mtc_df <- data.frame(
  module = rownames(modTraitCor),
  cor = as.numeric(modTraitCor),
  p = as.numeric(modTraitP)
) %>% arrange(p)

safe_write(mtc_df, file.path(res_dir, "WGCNA_Module_vs_FerroScore.csv"))

# Plot module–trait relationships
pdf(file.path(fig_dir, "WGCNA_ModuleTrait_FerroScore.pdf"), width = 6, height = 8)
labeledHeatmap(
  Matrix = modTraitCor,
  xLabels = "FerroScore",
  yLabels = rownames(modTraitCor),
  ySymbols = rownames(modTraitCor),
  colors = blueWhiteRed(50),
  textMatrix = paste0(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")"),
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)
dev.off()

# Identify key module
key_mod <- mtc_df %>% filter(p < 0.05) %>% arrange(desc(abs(cor))) %>% slice_head(n = 1)
if (nrow(key_mod) == 0) {
  key_mod <- mtc_df %>% arrange(desc(abs(cor))) %>% slice_head(n = 1)
}
key_color <- labels2colors(as.numeric(gsub("ME", "", key_mod$module)))
message("Key module: ", key_mod$module, " (", key_color, ") - cor=", round(key_mod$cor, 3), ", p=", signif(key_mod$p, 3))

# Gene-module relationships
kME <- signedKME(wgcna_dat, MEs_col)
colnames(kME) <- paste0("kME", colnames(MEs_col))

GS_ferro <- as.numeric(cor(wgcna_dat, trait_df$FerroScore, use = "pairwise.complete.obs"))
GS_ferro_p <- as.numeric(corPvalueStudent(GS_ferro, nrow(wgcna_dat)))

wgcna_annot <- data.frame(
  SYMBOL = colnames(wgcna_dat),
  module = moduleColors,
  GS_Ferro = GS_ferro,
  GSp_Ferro = GS_ferro_p,
  kME_key = kME[, paste0("kME", key_mod$module)],
  stringsAsFactors = FALSE
)

safe_write(wgcna_annot, file.path(res_dir, "WGCNA_GeneAnnotations.csv"))

# Candidate hub genes
cand_hubs <- wgcna_annot %>%
  filter(module == key_color) %>%
  mutate(isFerro = SYMBOL %in% fer_symbols, abs_GS = abs(GS_Ferro))
if (nrow(cand_hubs) > 10) {
  kME_threshold <- quantile(cand_hubs$kME_key, 0.7, na.rm = TRUE)
  GS_threshold <- quantile(cand_hubs$abs_GS, 0.7, na.rm = TRUE)
} else {
  kME_threshold <- quantile(cand_hubs$kME_key, 0.5, na.rm = TRUE)
  GS_threshold <- quantile(cand_hubs$abs_GS, 0.5, na.rm = TRUE)
}
cand_hubs_filtered <- cand_hubs %>%
  filter(kME_key >= kME_threshold & abs_GS >= GS_threshold) %>%
  arrange(desc(kME_key))

safe_write(cand_hubs_filtered, file.path(res_dir, "Candidate_Hubs_inKeyModule.csv"))
message("Candidate hub genes identified: ", nrow(cand_hubs_filtered))

# ---------------------------
# 6) Survival analysis and modeling
# ---------------------------
message("Step 6: Survival analysis and modeling...")

# Extract survival data
os_days <- pmax(as.numeric(clin_sub$days_to_death), as.numeric(clin_sub$days_to_last_follow_up), na.rm = TRUE)
os_status <- ifelse(!is.na(as.numeric(clin_sub$days_to_death)) & as.numeric(clin_sub$days_to_death) > 0, 1, 0)

valid_surv <- !is.na(os_days) & !is.na(os_status) & os_days > 0
os_days <- os_days[valid_surv]
os_status <- os_status[valid_surv]
clin_surv <- clin_sub[valid_surv, ]
expr_tcga <- vst_top[, rownames(clin_surv)]

# Candidate genes for Cox regression
cox_candidates <- unique(c(cand_hubs_filtered$SYMBOL, fe_deg$SYMBOL[fe_deg$sig]))
cox_genes <- intersect(cox_candidates, rownames(expr_tcga))
message("Cox regression candidate genes: ", length(cox_genes))
if (length(cox_genes) < 5) {
  cox_genes <- fe_deg %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(50) %>% pull(SYMBOL)
  cox_genes <- intersect(cox_genes, rownames(expr_tcga))
}

expr_cox <- t(expr_tcga[cox_genes, , drop = FALSE])
cox_data <- data.frame(
  sample = rownames(expr_cox),
  time = os_days,
  status = os_status,
  expr_cox,
  stringsAsFactors = FALSE
)

# Univariate Cox regression
message("Running univariate Cox regression...")
univ_results <- map_dfr(cox_genes, function(gene) {
  tryCatch({
    form <- as.formula(paste("Surv(time, status) ~", gene))
    fit <- coxph(form, data = cox_data)
    summary_fit <- summary(fit)
    data.frame(
      gene = gene,
      hr = summary_fit$coefficients[1, "exp(coef)"],
      hr_lower = summary_fit$conf.int[1, "lower .95"],
      hr_upper = summary_fit$conf.int[1, "upper .95"],
      z_score = summary_fit$coefficients[1, "z"],
      p_value = summary_fit$coefficients[1, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(gene = gene, hr = NA, hr_lower = NA, hr_upper = NA,
               z_score = NA, p_value = NA, stringsAsFactors = FALSE)
  })
}) %>% filter(!is.na(p_value)) %>% arrange(p_value) %>% mutate(fdr = p.adjust(p_value, method = "BH"))

safe_write(univ_results, file.path(res_dir, "Univariate_Cox_TCGA.csv"))

# LASSO Cox regression
message("Running LASSO Cox regression...")
set.seed(20250903)
x_lasso <- as.matrix(expr_cox)
y_lasso <- Surv(cox_data$time, cox_data$status)
# Remove genes with zero variance
gene_vars <- apply(x_lasso, 2, var, na.rm = TRUE)
x_lasso <- x_lasso[, gene_vars > 0]
cv_fit <- cv.glmnet(x_lasso, y_lasso, family = "cox", alpha = 1, nfolds = 10, standardize = TRUE, type.measure = "C")
# Save CV curve
pdf(file.path(fig_dir, "LASSO_CV_Curve.pdf"), width = 8, height = 6)
plot(cv_fit)
title("LASSO Cross-Validation")
dev.off()

lambda_opt <- cv_fit$lambda.1se
lasso_fit <- glmnet(x_lasso, y_lasso, family = "cox", alpha = 1, lambda = lambda_opt, standardize = TRUE)
coef_lasso <- coef(lasso_fit, s = lambda_opt)
selected_genes <- rownames(coef_lasso)[which(coef_lasso != 0)]
selected_coefs <- as.numeric(coef_lasso[coef_lasso != 0])
if (length(selected_genes) == 0) {
  lambda_opt <- cv_fit$lambda.min
  lasso_fit <- glmnet(x_lasso, y_lasso, family = "cox", alpha = 1, lambda = lambda_opt, standardize = TRUE)
  coef_lasso <- coef(lasso_fit, s = lambda_opt)
  selected_genes <- rownames(coef_lasso)[which(coef_lasso != 0)]
  selected_coefs <- as.numeric(coef_lasso[coef_lasso != 0])
}

lasso_results <- data.frame(gene = selected_genes, coefficient = selected_coefs, stringsAsFactors = FALSE) %>% arrange(desc(abs(coefficient)))
safe_write(lasso_results, file.path(res_dir, "LASSO_Selected_Genes_TCGA.csv"))
message("LASSO selected genes: ", nrow(lasso_results))

# Compute risk score and perform survival analysis if any genes were selected
if (nrow(lasso_results) > 0) {
  risk_scores <- as.numeric(x_lasso[, selected_genes, drop = FALSE] %*% selected_coefs)
  risk_df <- data.frame(sample = rownames(x_lasso), risk_score = risk_scores, time = cox_data$time, status = cox_data$status, stringsAsFactors = FALSE)
  safe_write(risk_df, file.path(res_dir, "RiskScore_TCGA.csv"))
  risk_median <- median(risk_scores, na.rm = TRUE)
  risk_df$risk_group <- ifelse(risk_df$risk_score >= risk_median, "High", "Low")
  surv_fit <- survfit(Surv(time, status) ~ risk_group, data = risk_df)
  pdf(file.path(fig_dir, "KM_Survival_TCGA.pdf"), width = 8, height = 6)
  ggsurvplot(
    surv_fit, data = risk_df, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
    xlab = "Time (days)", ylab = "Overall survival probability",
    title = "Kaplan-Meier Survival Curves (TCGA-BRCA)", risk.table.height = 0.3
  )
  dev.off()
  # Time-dependent ROC
  time_points <- c(365, 1095, 1825)
  tryCatch({
    roc_obj <- timeROC(T = risk_df$time, delta = risk_df$status, marker = risk_df$risk_score,
                       cause = 1, times = time_points, iid = TRUE)
    roc_results <- data.frame(
      time_point = paste0(time_points/365, "_years"),
      time_days = time_points,
      auc = roc_obj$AUC,
      stringsAsFactors = FALSE
    )
    safe_write(roc_results, file.path(res_dir, "TimeROC_TCGA.csv"))
    pdf(file.path(fig_dir, "TimeROC_TCGA.pdf"), width = 8, height = 6)
    plot(roc_obj, time = 1095, col = "red", lwd = 2, title = "Time-dependent ROC (3 years)")
    dev.off()
    message("Time-dependent ROC AUCs: ", paste(paste0(time_points/365, "y: ", round(roc_obj$AUC, 3)), collapse = ", "))
  }, error = function(e) {
    message("Error in time-dependent ROC: ", e$message)
  })
} else {
  message("No genes selected by LASSO. Skipping risk score calculation.")
}

# ---------------------------
# 7) External validation
# ---------------------------
message("Step 7: External validation using GEO data...")

graph validation_results <- list()
geo_datasets <- c("GSE20685", "GSE45255", "GSE58812")

for (geo_id in geo_datasets) {
  message("Attempting validation with ", geo_id, "...")
  geo_cache_file <- file.path(cache_dir, paste0(geo_id, "_eset.rds"))
  tryCatch({
    if (!file.exists(geo_cache_file)) {
      gset <- getGEO(geo_id, GSEMatrix = TRUE, getGPL = FALSE)
      if (length(gset) > 0) {
        eset <- gset[[1]]
        saveRDS(eset, geo_cache_file)
      } else {
        next
      }
    } else {
      eset <- readRDS(geo_cache_file)
    }
    expr_geo <- exprs(eset)
    pdata_geo <- pData(eset)
    if (geo_id %in% c("GSE20685", "GSE45255")) {
      if (requireNamespace("hgu133plus2.db", quietly = TRUE)) {
        annot <- AnnotationDbi::select(
          hgu133plus2.db::hgu133plus2.db,
          keys = rownames(expr_geo), keytype = "PROBEID", columns = "SYMBOL"
        )
        annot <- annot[!is.na(annot$SYMBOL) & annot$SYMBOL != "", ]
        expr_geo_filt <- expr_geo[annot$PROBEID, ]
        expr_geo_symbol <- collapse_by_symbol(expr_geo_filt, annot$SYMBOL)
      } else {
        message("Package hgu133plus2.db not available, skipping ", geo_id)
        next
      }
    } else {
      expr_geo_symbol <- expr_geo
    }
    expr_geo_norm <- t(scale(t(expr_geo_symbol), center = TRUE, scale = TRUE))
    surv_info <- extract_surv_from_pdata(pdata_geo)
    geo_time <- surv_info$time
    geo_status <- surv_info$status
    valid_surv_geo <- !is.na(geo_time) & !is.na(geo_status) & geo_time > 0
    if (sum(valid_surv_geo) < 50) {
      message("Insufficient survival data in ", geo_id, ", skipping...")
      next
    }
    if (exists("selected_genes") && length(selected_genes) > 0) {
      common_genes <- intersect(selected_genes, rownames(expr_geo_norm))
      if (length(common_genes) >= max(3, length(selected_genes) * 0.5)) {
        risk_geo <- as.numeric(t(expr_geo_norm[common_genes, valid_surv_geo, drop = FALSE]) %*%
                                 selected_coefs[match(common_genes, selected_genes)])
        validation_df <- data.frame(
          sample = colnames(expr_geo_norm)[valid_surv_geo], risk_score = risk_geo,
          time = geo_time[valid_surv_geo], status = geo_status[valid_surv_geo],
          dataset = geo_id, stringsAsFactors = FALSE
        )
        validation_df$risk_group <- ifelse(validation_df$risk_score >= median(validation_df$risk_score, na.rm = TRUE), "High", "Low")
        surv_fit_val <- survfit(Surv(time, status) ~ risk_group, data = validation_df)
        surv_test <- survdiff(Surv(time, status) ~ risk_group, data = validation_df)
        p_val <- 1 - pchisq(surv_test$chisq, length(surv_test$n) - 1)
        pdf(file.path(fig_dir, paste0("KM_Validation_", geo_id, ".pdf")), width = 8, height = 6)
        print(ggsurvplot(
          surv_fit_val, data = validation_df, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
          xlab = "Time", ylab = "Survival probability", title = paste("Validation in", geo_id),
          risk.table.height = 0.3
        ))
        dev.off()
        validation_results[[geo_id]] <- list(
          data = validation_df, p_value = p_val, n_patients = nrow(validation_df),
          common_genes = length(common_genes)
        )
        message("Validation in ", geo_id, " completed. P-value: ", signif(p_val, 3), ", N=", nrow(validation_df))
      } else {
        message("Insufficient gene overlap with ", geo_id, " (", length(common_genes), "/", length(selected_genes), ")")
      }
    }
  }, error = function(e) {
    message("Error processing ", geo_id, ": ", e$message)
  })
}

# Summarize validation results
if (length(validation_results) > 0) {
  validation_summary <- map_dfr(names(validation_results), function(dataset) {
    res <- validation_results[[dataset]]
    data.frame(dataset = dataset, n_patients = res$n_patients, p_value = res$p_value,
               significant = res$p_value < 0.05, common_genes = res$common_genes,
               stringsAsFactors = FALSE)
  })
  safe_write(validation_summary, file.path(res_dir, "Validation_Summary.csv"))
  message("Validation completed in ", nrow(validation_summary), " datasets")
} else {
  message("No successful validation datasets")
}

# ---------------------------
# 8) Multivariable Cox regression and nomogram
# ---------------------------
message("Step 8: Multivariable Cox regression and nomogram...")

if (exists("risk_df") && nrow(risk_df) > 0) {
  clinical_vars <- data.frame(
    sample = rownames(expr_cox),
    age = as.numeric(clin_surv$age_at_diagnosis) / 365.25,
    stage = clin_surv$ajcc_pathologic_stage,
    grade = clin_surv$neoplasm_histologic_grade,
    er_status = clin_surv$er_status_by_ihc,
    pr_status = clin_surv$pr_status_by_ihc,
    her2_status = clin_surv$her2_status_by_ihc,
    stringsAsFactors = FALSE
  )
  clinical_vars$stage_clean <- case_when(
    grepl("Stage I[^IV]", clinical_vars$stage) ~ "Stage_I",
    grepl("Stage II", clinical_vars$stage) ~ "Stage_II",
    grepl("Stage III", clinical_vars$stage) ~ "Stage_III",
    grepl("Stage IV", clinical_vars$stage) ~ "Stage_IV",
    TRUE ~ NA_character_
  )
  clinical_vars$grade_clean <- case_when(
    grepl("G1|Grade 1", clinical_vars$grade) ~ "Grade_1",
    grepl("G2|Grade 2", clinical_vars$grade) ~ "Grade_2",
    grepl("G3|Grade 3", clinical_vars$grade) ~ "Grade_3",
    TRUE ~ NA_character_
  )
  clinical_vars$er_clean <- case_when(
    grepl("Positive", clinical_vars$er_status) ~ "Positive",
    grepl("Negative", clinical_vars$er_status) ~ "Negative",
    TRUE ~ NA_character_
  )
  
  multiv_df <- risk_df %>% left_join(clinical_vars, by = "sample") %>%
    filter(!is.na(time) & !is.na(status) & time > 0 & !is.na(risk_score) & !is.na(age) & age > 0)
  if (nrow(multiv_df) >= 100) {
    cox_formula <- "Surv(time, status) ~ risk_score + age"
    if (sum(!is.na(multiv_df$stage_clean)) >= 50) {
      cox_formula <- paste(cox_formula, "+ stage_clean")
      multiv_df$stage_clean <- factor(multiv_df$stage_clean)
    }
    if (sum(!is.na(multiv_df$er_clean)) >= 50) {
      cox_formula <- paste(cox_formula, "+ er_clean")
      multiv_df$er_clean <- factor(multiv_df$er_clean)
    }
    tryCatch({
      cox_multiv <- coxph(as.formula(cox_formula), data = multiv_df)
      cox_summary <- broom::tidy(cox_multiv, exponentiate = TRUE, conf.int = TRUE)
      safe_write(cox_summary, file.path(res_dir, "Multivariable_Cox_TCGA.csv"))
      cox_zph <- cox.zph(cox_multiv)
      safe_write(broom::tidy(cox_zph), file.path(res_dir, "Cox_Proportional_Hazards_Test.csv"))
      if (requireNamespace("rms", quietly = TRUE)) {
        tryCatch({
          dd <- datadist(multiv_df)
          options(datadist = "dd")
          cox_rms <- cph(as.formula(cox_formula), data = multiv_df, x = TRUE, y = TRUE, surv = TRUE)
          surv_times <- c(365, 1095, 1825)
          surv_funcs <- map(surv_times, function(t) {
            function(x) survest(cox_rms, newdata = multiv_df, times = t)$surv
          })
          nom <- nomogram(cox_rms, fun = surv_funcs,
                          funlabel = c("1-year Survival", "3-year Survival", "5-year Survival"),
                          maxscale = 100)
          pdf(file.path(fig_dir, "Nomogram_TCGA.pdf"), width = 12, height = 8)
          plot(nom, xfrac = 0.3)
          dev.off()
          tryCatch({
            cal_3y <- calibrate(cox_rms, u = 1095, method = "boot", B = 100)
            pdf(file.path(fig_dir, "Calibration_3year_TCGA.pdf"), width = 8, height = 6)
            plot(cal_3y, main = "3-year Survival Calibration")
            dev.off()
          }, error = function(e) {
            message("Error in calibration plot: ", e$message)
          })
          message("Nomogram created successfully")
        }, error = function(e) {
          message("Error creating nomogram: ", e$message)
        })
      }
    }, error = function(e) {
      message("Error in multivariable Cox regression: ", e$message)
    })
  } else {
    message("Insufficient samples for multivariable analysis")
  }
}

# ---------------------------
# 9) Protein-protein interaction network analysis (optional)
# ---------------------------
message("Step 9: Protein-protein interaction analysis (optional)...")

if (exists("selected_genes") && length(selected_genes) >= 5) {
  tryCatch({
    if (requireNamespace("STRINGdb", quietly = TRUE) && requireNamespace("igraph", quietly = TRUE)) {
      string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400, input_directory = cache_dir)
      gene_df <- data.frame(gene = selected_genes, stringsAsFactors = FALSE)
      mapped_genes <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)
      if (nrow(mapped_genes) >= 3) {
        interactions <- string_db$get_interactions(mapped_genes$STRING_id)
        if (nrow(interactions) > 0) {
          g <- igraph::graph_from_data_frame(
            interactions[, c("from", "to", "combined_score")], directed = FALSE
          )
          degree_centrality <- igraph::degree(g)
          betweenness_centrality <- igraph::betweenness(g)
          closeness_centrality <- igraph::closeness(g)
          network_stats <- data.frame(
            STRING_id = names(degree_centrality),
            degree = degree_centrality,
            betweenness = betweenness_centrality,
            closeness = closeness_centrality,
            stringsAsFactors = FALSE
          ) %>% left_join(mapped_genes[, c("STRING_id", "gene")], by = "STRING_id") %>% arrange(desc(degree))
          safe_write(network_stats, file.path(res_dir, "PPI_Network_Statistics.csv"))
          pdf(file.path(fig_dir, "PPI_Network.pdf"), width = 10, height = 10)
          plot(g,
               vertex.size = sqrt(degree_centrality) * 3,
               vertex.label = mapped_genes$gene[match(V(g)$name, mapped_genes$STRING_id)],
               vertex.label.cex = 0.8,
               edge.width = E(g)$combined_score / 200,
               layout = layout_with_fr,
               main = "Protein-Protein Interaction Network")
          dev.off()
          message("PPI network analysis completed with ", nrow(mapped_genes), " genes and ", nrow(interactions), " interactions")
        }
      }
    }
  }, error = function(e) {
    message("Error in PPI analysis: ", e$message)
  })
}

# ---------------------------
# 10) Generate summary and cleanup
# ---------------------------
message("Step 10: Generating final reports...")

analysis_summary <- list(
  project_info = list(
    date = Sys.Date(),
    r_version = R.version.string,
    working_directory = proj_dir
  ),
  data_summary = list(
    tcga_samples = ncol(count_mat_final),
    tcga_genes = nrow(count_mat_final),
    ferroptosis_genes = length(fer_symbols),
    de_genes = sum(res_df$sig, na.rm = TRUE),
    ferroptosis_deg = nrow(fe_deg)
  ),
  wgcna_summary = list(
    input_genes = ncol(vst_top),
    modules_identified = if (exists("net")) ncol(net$MEs) else 0,
    key_module = if (exists("key_mod")) key_mod$module else "none",
    hub_genes = if (exists("cand_hubs_filtered")) nrow(cand_hubs_filtered) else 0
  ),
  survival_summary = list(
    lasso_genes = if (exists("selected_genes")) length(selected_genes) else 0,
    validation_datasets = if (exists("validation_summary")) nrow(validation_summary) else 0
  )
)

writeLines(capture.output(str(analysis_summary)), file.path(res_dir, "Analysis_Summary.txt"))

if (exists("selected_genes") && length(selected_genes) > 0) {
  key_results <- data.frame(
    gene = selected_genes,
    coefficient = selected_coefs,
    is_ferroptosis = selected_genes %in% fer_symbols,
    stringsAsFactors = FALSE
  ) %>% arrange(desc(abs(coefficient)))
  safe_write(key_results, file.path(res_dir, "Key_Prognostic_Genes.csv"))
}

sink(file.path(res_dir, "SessionInfo.txt"))
cat("Analysis completed on:", as.character(Sys.time()), "\n\n")
sessionInfo()
sink()

# Clean up large objects to free memory
rm(list = c("count_mat", "count_mat_final", "vst_mat", "expr_tcga"))

# Display completion messages
message("\n", strrep("=", 60))
message("ANALYSIS COMPLETED SUCCESSFULLY!")
message(strrep("=", 60))
message("Results saved in: ", res_dir)
message("Figures saved in: ", fig_dir)
message("Cache files in: ", cache_dir)
message(strrep("=", 60))
