############################################################
# SU2C-MARK
# Macro_FCGBP ssGSEA + Response boxplot + Logistic regression
############################################################

library(GSVA)
library(dplyr)
library(ggplot2)

############################################################
# 1. 路径
############################################################

root <- "/data3/wangmeiheng/pan-cancer/FCGBP_public_data_20260924"

expr_file <- file.path(
  root,
  "processed/SU2C_MARK/expression_source_scale.tsv.gz"
)

clinical_file <- file.path(
  root,
  "processed/SU2C_MARK/clinical_aligned.tsv"
)

signature_file <- "/data3/wangmeiheng/pan-cancer/Macro_FCGBP_signature.csv"

outdir <- file.path(
  root,
  "raw/SU2C_MARK/results/Macro_FCGBP_ssGSEA_simple"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# 2. 读取表达矩阵
# 行 = gene
# 列 = sample
############################################################

con <- gzfile(expr_file, "rt")

tpm <- read.delim(
  con,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

close(con)

tpm <- as.matrix(tpm)
storage.mode(tpm) <- "numeric"

cat(
  "Expression matrix:",
  nrow(tpm), "genes x",
  ncol(tpm), "samples\n"
)

############################################################
# 3. 读取临床数据
############################################################

clinical <- read.delim(
  clinical_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

cat(
  "Clinical:",
  nrow(clinical),
  "samples\n"
)

# 检查response_binary是否存在
if (!"response_binary" %in% colnames(clinical)) {
  stop("clinical中不存在 response_binary 列")
}

############################################################
# 4. 读取Macro_FCGBP signature
############################################################

signature <- read.csv(
  signature_file,
  check.names = FALSE
)

if (!"gene" %in% colnames(signature)) {
  stop("signature文件中不存在 gene 列")
}

markerGenes <- unique(
  trimws(
    as.character(signature$gene)
  )
)

markerGenes <- markerGenes[
  !is.na(markerGenes) &
    markerGenes != ""
]

############################################################
# 5. signature与表达矩阵取交集
############################################################

genes_use <- intersect(
  markerGenes,
  rownames(tpm)
)

cat(
  "Signature genes:", length(markerGenes), "\n",
  "Matched genes:", length(genes_use), "\n",
  "Coverage:",
  round(
    length(genes_use) / length(markerGenes) * 100,
    2
  ),
  "%\n"
)

if (length(genes_use) < 2) {
  stop("匹配到的signature基因太少，无法进行ssGSEA")
}

write.table(
  data.frame(gene = genes_use),
  file = file.path(
    outdir,
    "Macro_FCGBP_genes_used.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# 6. log2(TPM + 1)
############################################################

expr <- log2(tpm + 1)

gene_sets <- list(
  Macro_FCGBP = genes_use
)

############################################################
# 7. ssGSEA
# 自动兼容新版/旧版GSVA
############################################################

if (packageVersion("GSVA") >= "1.50.0") {
  
  param <- ssgseaParam(
    expr,
    gene_sets,
    alpha = 0.25,
    normalize = TRUE
  )
  
  ssgsea_score <- gsva(
    param,
    verbose = FALSE
  )
  
} else {
  
  ssgsea_score <- gsva(
    expr,
    gene_sets,
    method = "ssgsea",
    kcdf = "Gaussian",
    abs.ranking = FALSE,
    ssgsea.norm = TRUE,
    parallel.sz = 1,
    verbose = FALSE
  )
}

############################################################
# 8. 整理ssGSEA结果
############################################################
############################################################
# 8. 整理ssGSEA结果 + PDGFB + ADGRV1表达
############################################################

# 检查两个基因是否存在
genes_single <- c("PDGFB", "ADGRV1")

missing_genes <- setdiff(
  genes_single,
  rownames(expr)
)

if (length(missing_genes) > 0) {
  stop(
    paste0(
      "表达矩阵中缺少基因: ",
      paste(missing_genes, collapse = ", ")
    )
  )
}

# sample × score/expression
score_df <- data.frame(
  sample = colnames(ssgsea_score),
  
  Macro_FCGBP_score = as.numeric(
    ssgsea_score[
      "Macro_FCGBP",
      colnames(ssgsea_score)
    ]
  ),
  
  PDGFB = as.numeric(
    expr[
      "PDGFB",
      colnames(ssgsea_score)
    ]
  ),
  
  ADGRV1 = as.numeric(
    expr[
      "ADGRV1",
      colnames(ssgsea_score)
    ]
  ),
  
  stringsAsFactors = FALSE
)

head(score_df)

# 确认
colnames(score_df)

summary(score_df$PDGFB)
summary(score_df$ADGRV1)

write.table(
  score_df,
  file = file.path(
    outdir,
    "Macro_FCGBP_PDGFB_ADGRV1_scores.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
############################################################
# 9. 合并临床数据
############################################################

clinical$sample <- rownames(clinical)

dat <- inner_join(
  score_df,
  clinical,
  by = "sample"
)

cat(
  "Matched expression + clinical:",
  nrow(dat),
  "samples\n"
)

############################################################
# 10. 可选：只保留Pre-treatment RNA QC = Keep
# 如果你想全部样本分析，把这一段注释掉
############################################################

if ("Pre-treatment_RNA_Sample_QC" %in% colnames(dat)) {
  
  dat_keep <- dat %>%
    filter(
      !is.na(`Pre-treatment_RNA_Sample_QC`),
      `Pre-treatment_RNA_Sample_QC` == "Keep"
    )
  
} else {
  
  dat_keep <- dat
}

cat(
  "Samples used:",
  nrow(dat_keep),
  "\n"
)

############################################################
# 11. Response数据
############################################################

dat_response <- dat_keep %>%
  filter(
    !is.na(Macro_FCGBP_score),
    !is.na(response_binary)
  )

# 转为numeric
dat_response$response_binary <- as.numeric(
  as.character(
    dat_response$response_binary
  )
)

# 检查0/1
dat_response <- dat_response %>%
  filter(
    response_binary %in% c(0, 1)
  )

if (nrow(dat_response) == 0) {
  stop("没有可用于response分析的样本")
}

############################################################
# 12. 定义Responder / Non-responder
############################################################

dat_response$response_group <- factor(
  dat_response$response_binary,
  levels = c(0, 1),
  labels = c(
    "Non-responder",
    "Responder"
  )
)

cat("\nResponse group:\n")
print(
  table(
    dat_response$response_group
  )
)

############################################################
# 13. Wilcoxon检验
############################################################

wilcox_res <- wilcox.test(
  Macro_FCGBP_score ~ response_group,
  data = dat_response,
  exact = FALSE
)

wilcox_p <- wilcox_res$p.value

cat(
  "\nWilcoxon P =",
  wilcox_p,
  "\n"
)

############################################################
# 14. 箱线图
############################################################
setwd("/data3/wangmeiheng/pan-cancer/FCGBP_public_data_20260924/raw/SU2C_MARK/results/Macro_FCGBP_ssGSEA_simple")
y_max <- max(
  dat_response$Macro_FCGBP_score,
  na.rm = TRUE
)

y_min <- min(
  dat_response$Macro_FCGBP_score,
  na.rm = TRUE
)

y_range <- y_max - y_min

if (y_range == 0) {
  y_range <- 1
}

p_box <- ggplot(
  dat_response,
  aes(
    x = response_group,
    y = Macro_FCGBP_score,
    fill = response_group
  )
) +
  
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    linewidth = 0.8
  ) +
  
  geom_jitter(
    aes(
      color = response_group
    ),
    width = 0.12,
    size = 2.2,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(
    values = c(
      "Non-responder" = "#386B9B",
      "Responder" = "#CB713A"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "Non-responder" = "#386B9B",
      "Responder" = "#CB713A"
    )
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = y_max + 0.08 * y_range,
    label = sprintf(
      "Wilcoxon P = %.3g",
      wilcox_p
    ),
    size = 4.5
  ) +
  
  expand_limits(
    y = y_max + 0.15 * y_range
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 12
    ),
    axis.title.y = element_text(
      color = "black",
      size = 14
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.8
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.8
    )
  ) +
  
  labs(
    y = "Macro_FCGBP ssGSEA score"
  )

ggsave(
  filename = file.path(
    outdir,
    "Macro_FCGBP_response_boxplot.pdf"
  ),
  plot = p_box,
  width = 5,
  height = 5
)


############################################################
# 15. Logistic regression
#
# 这里定义：
# response_binary = 1 = Responder
# response_binary = 0 = Non-responder
#
# 因此：
# OR < 1 表示score越高越倾向Non-responder
############################################################

dat_response$score_Z <- as.numeric(
  scale(
    dat_response$Macro_FCGBP_score
  )
)

fit_logistic <- glm(
  response_binary ~ score_Z,
  data = dat_response,
  family = binomial()
)

summary_logistic <- summary(
  fit_logistic
)

coef_tab <- summary_logistic$coefficients

beta <- coef_tab[
  "score_Z",
  "Estimate"
]

se <- coef_tab[
  "score_Z",
  "Std. Error"
]

logistic_p <- coef_tab[
  "score_Z",
  "Pr(>|z|)"
]

OR <- exp(beta)

CI_low <- exp(
  beta - 1.96 * se
)

CI_high <- exp(
  beta + 1.96 * se
)

############################################################
# 16. Logistic结果
############################################################

logistic_result <- data.frame(
  n = nrow(dat_response),
  
  Responder = sum(
    dat_response$response_binary == 1
  ),
  
  Non_responder = sum(
    dat_response$response_binary == 0
  ),
  
  OR_per_1SD_for_response = OR,
  
  CI_low = CI_low,
  
  CI_high = CI_high,
  
  Logistic_P = logistic_p,
  
  Wilcoxon_P = wilcox_p
)

print(
  logistic_result
)

write.table(
  logistic_result,
  file = file.path(
    outdir,
    "Macro_FCGBP_response_logistic.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# 17. 保存用于分析的数据
############################################################

write.table(
  dat_response,
  file = file.path(
    outdir,
    "Macro_FCGBP_response_analysis_data.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# 18. 输出结果
############################################################

cat("\n==============================\n")

cat(
  "Wilcoxon P =",
  signif(
    wilcox_p,
    4
  ),
  "\n"
)

cat(
  "Logistic OR per 1 SD =",
  round(
    OR,
    3
  ),
  "\n"
)

cat(
  "95% CI =",
  round(
    CI_low,
    3
  ),
  "-",
  round(
    CI_high,
    3
  ),
  "\n"
)

cat(
  "Logistic P =",
  signif(
    logistic_p,
    4
  ),
  "\n"
)

cat(
  "Results saved in:",
  outdir,
  "\n"
)

cat(
  "==============================\n"
)
###############################################################################
################################三指标#########################################
###############################################################################
############################################################
# 8. 整理 Macro_FCGBP + PDGFB + ADGRV1
############################################################

# 检查基因是否存在
genes_single <- c("PDGFB", "ADGRV1")

missing_genes <- setdiff(
  genes_single,
  rownames(expr)
)

if (length(missing_genes) > 0) {
  stop(
    paste0(
      "表达矩阵中缺少基因: ",
      paste(missing_genes, collapse = ", ")
    )
  )
}

feature_df <- data.frame(
  sample = colnames(expr),
  
  Macro_FCGBP_score = as.numeric(
    ssgsea_score[
      "Macro_FCGBP",
      colnames(expr)
    ]
  ),
  
  PDGFB = as.numeric(
    expr[
      "PDGFB",
      colnames(expr)
    ]
  ),
  
  ADGRV1 = as.numeric(
    expr[
      "ADGRV1",
      colnames(expr)
    ]
  ),
  
  stringsAsFactors = FALSE
)

head(feature_df)

write.table(
  feature_df,
  file = file.path(
    outdir,
    "Macro_FCGBP_PDGFB_ADGRV1_scores.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# 9. 合并临床数据
############################################################

clinical$sample <- rownames(clinical)

dat <- inner_join(
  feature_df,
  clinical,
  by = "sample"
)

cat(
  "Matched expression + clinical:",
  nrow(dat),
  "samples\n"
)

############################################################
# 10. 只保留 Pre-treatment RNA QC = Keep
############################################################

if ("Pre-treatment_RNA_Sample_QC" %in% colnames(dat)) {
  
  dat_keep <- dat %>%
    filter(
      `Pre-treatment_RNA_Sample_QC` == "Keep"
    )
  
} else {
  
  dat_keep <- dat
}

cat(
  "Samples used:",
  nrow(dat_keep),
  "\n"
)

############################################################
# 11. Response数据
############################################################

dat_response <- dat_keep %>%
  filter(
    !is.na(response_binary)
  )

dat_response$response_binary <- as.numeric(
  as.character(
    dat_response$response_binary
  )
)

dat_response <- dat_response %>%
  filter(
    response_binary %in% c(0, 1)
  )

if (nrow(dat_response) == 0) {
  stop("没有可用于response分析的样本")
}

############################################################
# 12. 定义Responder / Non-responder
############################################################

dat_response$response_group <- factor(
  dat_response$response_binary,
  levels = c(0, 1),
  labels = c(
    "Non-responder",
    "Responder"
  )
)

cat("\nResponse group:\n")

print(
  table(
    dat_response$response_group
  )
)

############################################################
# 13. 设置需要分析的3个指标
############################################################

features <- c(
  "Macro_FCGBP_score",
  "PDGFB",
  "ADGRV1"
)

feature_labels <- c(
  Macro_FCGBP_score = "Macro_FCGBP ssGSEA score",
  PDGFB = "PDGFB expression",
  ADGRV1 = "ADGRV1 expression"
)

############################################################
# 14. 创建结果对象
############################################################

result_list <- list()

############################################################
# 15. 循环：箱线图 + Wilcoxon + Logistic
############################################################

for (feature in features) {
  
  cat(
    "\n=============================\n",
    "Analyzing:", feature, "\n",
    "=============================\n"
  )
  
  # 当前指标的数据
  tmp <- dat_response %>%
    filter(
      !is.na(.data[[feature]])
    )
  
  if (nrow(tmp) < 10) {
    cat("Too few samples, skip:", feature, "\n")
    next
  }
  
  ##########################################################
  # Wilcoxon
  ##########################################################
  
  wilcox_res <- wilcox.test(
    tmp[[feature]] ~ tmp$response_group,
    exact = FALSE
  )
  
  wilcox_p <- wilcox_res$p.value
  
  ##########################################################
  # Z-score
  ##########################################################
  
  tmp$feature_Z <- as.numeric(
    scale(
      tmp[[feature]]
    )
  )
  
  ##########################################################
  # Logistic regression
  #
  # response_binary:
  # 1 = Responder
  # 0 = Non-responder
  #
  # OR < 1:
  # 表达/评分越高，越倾向Non-responder
  ##########################################################
  
  fit_logistic <- glm(
    response_binary ~ feature_Z,
    data = tmp,
    family = binomial()
  )
  
  coef_tab <- summary(
    fit_logistic
  )$coefficients
  
  beta <- coef_tab[
    "feature_Z",
    "Estimate"
  ]
  
  se <- coef_tab[
    "feature_Z",
    "Std. Error"
  ]
  
  logistic_p <- coef_tab[
    "feature_Z",
    "Pr(>|z|)"
  ]
  
  OR <- exp(beta)
  
  CI_low <- exp(
    beta - 1.96 * se
  )
  
  CI_high <- exp(
    beta + 1.96 * se
  )
  
  ##########################################################
  # 保存统计结果
  ##########################################################
  
  result_list[[feature]] <- data.frame(
    
    Feature = feature,
    
    n = nrow(tmp),
    
    Responder = sum(
      tmp$response_binary == 1
    ),
    
    Non_responder = sum(
      tmp$response_binary == 0
    ),
    
    Non_responder_median = median(
      tmp[[feature]][
        tmp$response_binary == 0
      ],
      na.rm = TRUE
    ),
    
    Responder_median = median(
      tmp[[feature]][
        tmp$response_binary == 1
      ],
      na.rm = TRUE
    ),
    
    OR_per_1SD_for_response = OR,
    
    CI_low = CI_low,
    
    CI_high = CI_high,
    
    Logistic_P = logistic_p,
    
    Wilcoxon_P = wilcox_p
  )
  
  ##########################################################
  # 箱线图
  ##########################################################
  
  y_max <- max(
    tmp[[feature]],
    na.rm = TRUE
  )
  
  y_min <- min(
    tmp[[feature]],
    na.rm = TRUE
  )
  
  y_range <- y_max - y_min
  
  if (
    !is.finite(y_range) ||
    y_range == 0
  ) {
    y_range <- 1
  }
  
  p_box <- ggplot(
    tmp,
    aes(
      x = response_group,
      y = .data[[feature]],
      fill = response_group
    )
  ) +
    
    geom_boxplot(
      width = 0.55,
      outlier.shape = NA,
      linewidth = 0.8
    ) +
    
    geom_jitter(
      aes(
        color = response_group
      ),
      width = 0.12,
      size = 2.2,
      alpha = 0.7,
      show.legend = FALSE
    ) +
    
    scale_fill_manual(
      values = c(
        "Non-responder" = "#386B9B",
        "Responder" = "#CB713A"
      )
    ) +
    
    scale_color_manual(
      values = c(
        "Non-responder" = "#386B9B",
        "Responder" = "#CB713A"
      )
    ) +
    
    annotate(
      "text",
      x = 1.5,
      y = y_max + 0.08 * y_range,
      label = sprintf(
        "Wilcoxon P = %.3g",
        wilcox_p
      ),
      size = 4.5
    ) +
    
    expand_limits(
      y = y_max + 0.15 * y_range
    ) +
    
    theme_classic(
      base_size = 14
    ) +
    
    theme(
      aspect.ratio = 1,
      
      legend.position = "none",
      
      axis.title.x = element_blank(),
      
      axis.text = element_text(
        color = "black",
        size = 12
      ),
      
      axis.title.y = element_text(
        color = "black",
        size = 14
      ),
      
      axis.line = element_line(
        color = "black",
        linewidth = 0.8
      ),
      
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.8
      )
    ) +
    
    labs(
      y = feature_labels[[feature]]
    )
  
  ##########################################################
  # 保存PDF
  ##########################################################
  
  ggsave(
    filename = file.path(
      outdir,
      paste0(
        feature,
        "_response_boxplot.pdf"
      )
    ),
    plot = p_box,
    width = 5,
    height = 5
  )
  
  ##########################################################
  
  
  ##########################################################
  # 输出当前结果
  ##########################################################
  
  cat(
    feature,
    "\nWilcoxon P =",
    signif(
      wilcox_p,
      4
    ),
    "\n"
  )
  
  cat(
    "Logistic OR per 1 SD =",
    round(
      OR,
      3
    ),
    "\n"
  )
  
  cat(
    "95% CI =",
    round(
      CI_low,
      3
    ),
    "-",
    round(
      CI_high,
      3
    ),
    "\n"
  )
  
  cat(
    "Logistic P =",
    signif(
      logistic_p,
      4
    ),
    "\n"
  )
}


############################################################
# 16. 合并统计结果
############################################################

logistic_results <- bind_rows(
  result_list
)

print(
  logistic_results
)

############################################################
# 17. 多重检验FDR
############################################################

logistic_results$Logistic_FDR <- p.adjust(
  logistic_results$Logistic_P,
  method = "BH"
)

logistic_results$Wilcoxon_FDR <- p.adjust(
  logistic_results$Wilcoxon_P,
  method = "BH"
)

############################################################
# 18. 保存结果
############################################################

write.table(
  logistic_results,
  file = file.path(
    outdir,
    "Macro_FCGBP_PDGFB_ADGRV1_response_results.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  dat_response,
  file = file.path(
    outdir,
    "Macro_FCGBP_PDGFB_ADGRV1_response_analysis_data.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# 19. 最终结果
############################################################

logistic_results
###############################################################################
########################## 额外分析1：PDGFB + ADGRV1 四分组与ICB响应 ################
###############################################################################

library(dplyr)
library(ggplot2)

############################################################
# 1. 保留PDGFB、ADGRV1和response完整样本
############################################################

dat_joint <- dat_response %>%
  filter(
    !is.na(PDGFB),
    !is.na(ADGRV1),
    !is.na(response_binary)
  )

############################################################
# 2. PDGFB和ADGRV1分别按中位数分High/Low
############################################################

PDGFB_cutoff <- median(
  dat_joint$PDGFB,
  na.rm = TRUE
)

ADGRV1_cutoff <- median(
  dat_joint$ADGRV1,
  na.rm = TRUE
)

cat(
  "PDGFB cutoff =", PDGFB_cutoff, "\n",
  "ADGRV1 cutoff =", ADGRV1_cutoff, "\n"
)

dat_joint$PDGFB_group <- ifelse(
  dat_joint$PDGFB >= PDGFB_cutoff,
  "High",
  "Low"
)

dat_joint$ADGRV1_group <- ifelse(
  dat_joint$ADGRV1 >= ADGRV1_cutoff,
  "High",
  "Low"
)

############################################################
# 3. 四分组
############################################################

dat_joint$Joint_group <- paste0(
  dat_joint$PDGFB_group,
  "_",
  dat_joint$ADGRV1_group
)

dat_joint$Joint_group <- factor(
  dat_joint$Joint_group,
  levels = c(
    "Low_Low",
    "High_Low",
    "Low_High",
    "High_High"
  )
)

cat("\nGroup number:\n")

print(
  table(
    dat_joint$Joint_group
  )
)

############################################################
# 4. 四组Responder / Non-responder数量
############################################################

cat("\nResponse table:\n")

response_table <- table(
  dat_joint$Joint_group,
  dat_joint$response_group
)

print(response_table)

############################################################
# 5. 计算各组响应率和不响应率
############################################################

response_stat <- dat_joint %>%
  group_by(Joint_group) %>%
  summarise(
    n = n(),
    
    Responder = sum(
      response_binary == 1
    ),
    
    Non_responder = sum(
      response_binary == 0
    ),
    
    Response_rate =
      Responder / n,
    
    Non_response_rate =
      Non_responder / n,
    
    .groups = "drop"
  )

print(response_stat)

############################################################
# 6. Fisher检验
# 检验四组与ICB response是否总体相关
############################################################

tab <- table(
  dat_joint$Joint_group,
  dat_joint$response_binary
)

fisher_result <- fisher.test(
  tab
)

fisher_p <- fisher_result$p.value

cat(
  "\nFisher exact test P =",
  fisher_p,
  "\n"
)

############################################################
# 7. Logistic regression
#
# 这里重新定义：
# NR = 1：Non-responder
# NR = 0：Responder
#
# 因此：
# OR > 1 表示更倾向ICB non-response
############################################################

dat_joint$NR <- ifelse(
  dat_joint$response_binary == 0,
  1,
  0
)

# 双低作为参考组
dat_joint$Joint_group <- relevel(
  dat_joint$Joint_group,
  ref = "Low_Low"
)

fit_joint <- glm(
  NR ~ Joint_group,
  data = dat_joint,
  family = binomial()
)

summary(fit_joint)

############################################################
# 8. 提取OR、95%CI、P
############################################################

coef_joint <- summary(
  fit_joint
)$coefficients

# 去掉截距
coef_joint <- coef_joint[
  rownames(coef_joint) != "(Intercept)",
  ,
  drop = FALSE
]

joint_logistic <- data.frame(
  Comparison = rownames(
    coef_joint
  ),
  
  beta = coef_joint[
    ,
    "Estimate"
  ],
  
  SE = coef_joint[
    ,
    "Std. Error"
  ],
  
  P = coef_joint[
    ,
    "Pr(>|z|)"
  ],
  
  stringsAsFactors = FALSE
)

joint_logistic$OR <- exp(
  joint_logistic$beta
)

joint_logistic$CI_low <- exp(
  joint_logistic$beta -
    1.96 * joint_logistic$SE
)

joint_logistic$CI_high <- exp(
  joint_logistic$beta +
    1.96 * joint_logistic$SE
)

joint_logistic <- joint_logistic %>%
  select(
    Comparison,
    OR,
    CI_low,
    CI_high,
    P
  )

print(
  joint_logistic
)

############################################################
# 9. 四组Non-response rate柱状图
############################################################

p_joint <- ggplot(
  response_stat,
  aes(
    x = Joint_group,
    y = Non_response_rate,
    fill = Joint_group
  )
) +
  
  geom_col(
    width = 0.7
  ) +
  
  geom_text(
    aes(
      label = paste0(
        Non_responder,
        "/",
        n
      )
    ),
    vjust = -0.4,
    size = 4
  ) +
  
  scale_fill_manual(
    values = c(
      "Low_Low" = "#386B9B",
      "High_Low" = "#70A1D7",
      "Low_High" = "#E6A15C",
      "High_High" = "#CB713A"
    )
  ) +
  
  scale_y_continuous(
    limits = c(
      0,
      max(
        response_stat$Non_response_rate
      ) * 1.2
    ),
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    
    axis.text.y = element_text(
      color = "black"
    )
  ) +
  
  labs(
    x = NULL,
    y = "ICB non-response rate",
    title = sprintf(
      "Fisher P = %.3g",
      fisher_p
    )
  )

p_joint

############################################################
# 10. 保存图片
############################################################

ggsave(
  filename = file.path(
    outdir,
    "PDGFB_ADGRV1_4group_nonresponse.pdf"
  ),
  plot = p_joint,
  width = 5,
  height = 5
)

############################################################
# 11. 保存统计结果
############################################################

write.table(
  response_stat,
  file = file.path(
    outdir,
    "PDGFB_ADGRV1_4group_response_rate.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  joint_logistic,
  file = file.path(
    outdir,
    "PDGFB_ADGRV1_4group_logistic.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  dat_joint,
  file = file.path(
    outdir,
    "PDGFB_ADGRV1_4group_response_data.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# 12. 输出
############################################################

cat(
  "\n====================================\n"
)

cat(
  "Fisher P =",
  signif(
    fisher_p,
    4
  ),
  "\n"
)

print(
  response_stat
)

print(
  joint_logistic
)

cat(
  "====================================\n"
)