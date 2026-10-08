dir.create("/data3/wangmeiheng/pan-cancer/FCGBP_public_data_20260924/raw/TCGA/TCGA_Macro_FCGBP_OS_KM_PDGFB", showWarnings = FALSE)
setwd("/data3/wangmeiheng/pan-cancer/FCGBP_public_data_20260924/raw/TCGA/TCGA_Macro_FCGBP_OS_KM_PDGFB")
dat1 <- tcga[c("PDGFB", "ADGRV1"), , drop = FALSE] %>%
  t() %>%
  as.data.frame()

dat1$sample <- rownames(dat1)

dat1 <- inner_join(
  dat1,
  tcga_clinical,
  by = "sample"
)
result <- list()
six_cancer<-c("BRCA","KIRC","COAD","READ","GBM","HNSC","LUAD","LUSC")
for (ca in intersect(six_cancer,unique(dat1$cancer.type.abbreviation))) {
  #ca=six_cancer[1]
  cat("Processing:", ca, "\n")
  
  tmp <- dat1 %>%
    filter(
      cancer.type.abbreviation == ca,
      !is.na(OS.time),
      !is.na(OS)
    )
  
  # 样本和死亡事件太少则跳过
  if (nrow(tmp) < 30 || sum(tmp$OS == 1) < 10)
    next
  
  # ① KM：每个癌种内部按中位数分组
  cutoff <- median(tmp$PDGFB)
  
  tmp$Score_group <- ifelse(
    tmp$PDGFB >= cutoff,
    "High", "Low"
  )
  
  tmp$Score_group <- factor(
    tmp$Score_group,
    levels = c("Low", "High")
  )
  
  fit <- survfit(
    Surv(OS.time, OS) ~ Score_group,
    data = tmp
  )
  
  # ② log-rank P
  logrank <- survdiff(
    Surv(OS.time, OS) ~ Score_group,
    data = tmp
  )
  
  logrank_p <- 1 - pchisq(
    logrank$chisq,
    df = length(logrank$n) - 1
  )
  
  # ③ 连续评分Cox：每增加1 SD
  tmp$score_Z <- as.numeric(scale(tmp$PDGFB))
  
  cox <- coxph(
    Surv(OS.time, OS) ~ score_Z,
    data = tmp
  )
  
  cs <- summary(cox)
  
  HR      <- cs$conf.int[1, "exp(coef)"]
  CI_low  <- cs$conf.int[1, "lower .95"]
  CI_high <- cs$conf.int[1, "upper .95"]
  cox_p   <- cs$coefficients[1, "Pr(>|z|)"]
  
  # ④ 提取KM数据
  ss <- summary(fit)
  
  plot_dat <- data.frame(
    time = ss$time,
    surv = ss$surv,
    strata = ss$strata
  )
  
  plot_dat$strata <- sub(
    "Score_group=", "",
    plot_dat$strata
  )
  
  # 补充t=0，使曲线从1开始
  start_dat <- data.frame(
    time = 0,
    surv = 1,
    strata = c("Low", "High")
  )
  
  plot_dat <- bind_rows(start_dat, plot_dat)
  
  # ⑤ ggplot画KM
  p <- ggplot(
    plot_dat,
    aes(
      x = time,
      y = surv,
      group = strata,
      color = strata
    )
  ) +
    geom_step(linewidth = 1.1) +
    
    # 红蓝配色
    scale_color_manual(
      values = c(
        "Low"  = "#3B6FB6",   # 蓝色
        "High" = "#D94A4A"    # 红色
      )
    ) +
    
    # Y轴固定0-1
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = c(0, 0)
    ) +
    
    # X轴稍微留白
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) +
    
    theme_classic(base_size = 14) +
    
    # 保持方形
    theme(
      aspect.ratio = 1,
      
      axis.title = element_text(
        size = 14,
        color = "black"
      ),
      
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      
      legend.title = element_text(
        size = 12
      ),
      
      legend.text = element_text(
        size = 11
      ),
      
      legend.position = c(0.80, 0.78),
      legend.background = element_blank(),
      
      plot.margin = margin(
        10, 15, 10, 10
      )
    ) +
    
    labs(
      title = ca,
      x = "Time (days)",
      y = "Overall survival probability",
      color = "Macro_FCGBP"
    ) +
    
    # 统计结果放右上角
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.4,
      size = 4.0,
      color = "black",
      label = sprintf(
        "Log-rank P = %.3g",
        logrank_p
        
      )
    )
  
  ggsave(
    paste0(
      "TCGA_Macro_FCGBP_OS_KM/",
      ca,
      "_Macro_FCGBP_OS.pdf"
    ),
    p,
    width = 5,
    height = 5
  )
  
  result[[ca]] <- data.frame(
    cancer = ca,
    n = nrow(tmp),
    death = sum(tmp$OS == 1),
    cutoff = cutoff,
    logrank_P = logrank_p,
    HR = HR,
    CI_low = CI_low,
    CI_high = CI_high,
    Cox_P = cox_p
  )
}

# 汇总所有癌种
result <- bind_rows(result) %>%
  mutate(
    FDR = p.adjust(Cox_P, method = "BH")
  ) %>%
  arrange(Cox_P)

write.csv(
  result,
  "TCGA_Macro_FCGBP_OS_result.csv",
  row.names = FALSE
)

######################################
dir.create(
  "TCGA_Macro_FCGBP_OS_KM_3year",
  showWarnings = FALSE
)

result_3y <- list()

for (ca in intersect(
  six_cancer,
  unique(dat1$cancer.type.abbreviation)
)) {
  
  cat("Processing 3-year OS:", ca, "\n")
  
  tmp <- dat1 %>%
    filter(
      cancer.type.abbreviation == ca,
      !is.na(PDGFB),
      !is.na(OS.time),
      !is.na(OS)
    )
  
  # 天 -> 月
  tmp$OS.month <- tmp$OS.time / 30.44
  
  # ---------------------------------
  # 3年截尾：36个月
  # ---------------------------------
  
  # 随访时间超过36个月的统一截到36个月
  tmp$time_3y <- pmin(
    tmp$OS.month,
    36
  )
  
  # 只有36个月以内发生的死亡算event
  # 36个月之后死亡，按36个月删失
  tmp$event_3y <- ifelse(
    tmp$OS == 1 &
      tmp$OS.month <= 36,
    1,
    0
  )
  
  # 样本数或3年内死亡事件太少则跳过
  if (
    nrow(tmp) < 30 ||
    sum(tmp$event_3y == 1) < 10
  ) {
    cat(
      "Skip:",
      ca,
      "n =", nrow(tmp),
      "3y deaths =", sum(tmp$event_3y),
      "\n"
    )
    next
  }
  
  # ---------------------------------
  # ① 每个癌种内部按中位数分High/Low
  # ---------------------------------
  
  cutoff <- median(
    tmp$PDGFB,
    na.rm = TRUE
  )
  
  tmp$Score_group <- ifelse(
    tmp$PDGFB >= cutoff,
    "High",
    "Low"
  )
  
  tmp$Score_group <- factor(
    tmp$Score_group,
    levels = c("Low", "High")
  )
  
  # ---------------------------------
  # ② 3年KM
  # ---------------------------------
  
  fit <- survfit(
    Surv(time_3y, event_3y) ~ Score_group,
    data = tmp
  )
  
  # log-rank
  logrank <- survdiff(
    Surv(time_3y, event_3y) ~ Score_group,
    data = tmp
  )
  
  logrank_p <- 1 - pchisq(
    logrank$chisq,
    df = length(logrank$n) - 1
  )
  
  # ---------------------------------
  # ③ 连续评分Cox：每增加1 SD
  # 同样使用3年截尾结局
  # ---------------------------------
  
  tmp$score_Z <- as.numeric(
    scale(tmp$PDGFB)
  )
  
  cox <- coxph(
    Surv(time_3y, event_3y) ~ score_Z,
    data = tmp
  )
  
  cs <- summary(cox)
  
  HR <- cs$conf.int[
    1,
    "exp(coef)"
  ]
  
  CI_low <- cs$conf.int[
    1,
    "lower .95"
  ]
  
  CI_high <- cs$conf.int[
    1,
    "upper .95"
  ]
  
  cox_p <- cs$coefficients[
    1,
    "Pr(>|z|)"
  ]
  
  # ---------------------------------
  # ④ 提取KM曲线
  # ---------------------------------
  
  ss <- summary(fit)
  
  plot_dat <- data.frame(
    time = ss$time,
    surv = ss$surv,
    strata = ss$strata
  )
  
  plot_dat$strata <- sub(
    "Score_group=",
    "",
    plot_dat$strata
  )
  
  # 从t=0、survival=1开始
  start_dat <- data.frame(
    time = 0,
    surv = 1,
    strata = c(
      "Low",
      "High"
    )
  )
  
  plot_dat <- bind_rows(
    start_dat,
    plot_dat
  )
  
  plot_dat$strata <- factor(
    plot_dat$strata,
    levels = c(
      "Low",
      "High"
    )
  )
  
  # ---------------------------------
  # ⑤ 绘制3年KM
  # ---------------------------------
  
  p <- ggplot(
    plot_dat,
    aes(
      x = time,
      y = surv,
      group = strata,
      color = strata
    )
  ) +
    geom_step(
      linewidth = 1.1
    ) +
    
    scale_color_manual(
      values = c(
        "Low" = "#3B6FB6",
        "High" = "#D94A4A"
      )
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        0.2
      ),
      expand = c(0, 0)
    ) +
    
    scale_x_continuous(
      limits = c(0, 36),
      breaks = c(
        0,
        12,
        24,
        36
      ),
      expand = c(0, 0)
    ) +
    
    theme_classic(
      base_size = 14
    ) +
    
    theme(
      aspect.ratio = 1,
      
      axis.title = element_text(
        size = 14,
        color = "black"
      ),
      
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      
      legend.title = element_text(
        size = 12
      ),
      
      legend.text = element_text(
        size = 11
      ),
      
      legend.position = c(
        0.80,
        0.78
      ),
      
      legend.background =
        element_blank(),
      
      plot.margin = margin(
        10,
        15,
        10,
        10
      )
    ) +
    
    labs(
      title = ca,
      x = "Time (months)",
      y = "Overall survival probability",
      color = "Macro_FCGBP"
    ) +
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.4,
      size = 4,
      color = "black",
      label = sprintf(
        paste0(
          "Log-rank P = %.3g\n",
          "HR per 1 SD = %.2f ",
          "(95%% CI %.2f–%.2f)\n",
          "Cox P = %.3g"
        ),
        logrank_p,
        HR,
        CI_low,
        CI_high,
        cox_p
      )
    )
  
  # 保存图
  ggsave(
    filename = paste0(
      "TCGA_Macro_FCGBP_OS_KM_3year/",
      ca,
      "_Macro_FCGBP_OS_3year.pdf"
    ),
    plot = p,
    width = 5,
    height = 5
  )
  
  # 保存统计结果
  result_3y[[ca]] <- data.frame(
    cancer = ca,
    n = nrow(tmp),
    death_3y = sum(
      tmp$event_3y == 1
    ),
    cutoff = cutoff,
    logrank_P = logrank_p,
    HR = HR,
    CI_low = CI_low,
    CI_high = CI_high,
    Cox_P = cox_p
  )
}

result_3y <- bind_rows(
  result_3y
)

result_3y$FDR <- p.adjust(
  result_3y$Cox_P,
  method = "BH"
)

write.csv(
  result_3y,
  "TCGA_Macro_FCGBP_OS_3year_results.csv",
  row.names = FALSE
)

result_3y
####################################
dir.create(
  "TCGA_Macro_FCGBP_OS_KM_5year",
  showWarnings = FALSE
)

result_5y <- list()

for (ca in intersect(
  six_cancer,
  unique(dat1$cancer.type.abbreviation)
)) {
  
  cat("Processing 5-year OS:", ca, "\n")
  
  tmp <- dat1 %>%
    filter(
      cancer.type.abbreviation == ca,
      !is.na(PDGFB),
      !is.na(OS.time),
      !is.na(OS)
    )
  
  # 天 -> 月
  tmp$OS.month <- tmp$OS.time / 30.44
  
  # ---------------------------------
  # 5年截尾：60个月
  # ---------------------------------
  
  tmp$time_5y <- pmin(
    tmp$OS.month,
    60
  )
  
  # 只有60个月以内死亡算event
  tmp$event_5y <- ifelse(
    tmp$OS == 1 &
      tmp$OS.month <= 60,
    1,
    0
  )
  
  # 样本数或5年内死亡事件太少则跳过
  if (
    nrow(tmp) < 30 ||
    sum(tmp$event_5y == 1) < 10
  ) {
    cat(
      "Skip:",
      ca,
      "n =", nrow(tmp),
      "5y deaths =", sum(tmp$event_5y),
      "\n"
    )
    next
  }
  
  # ---------------------------------
  # ① 中位数分High/Low
  # ---------------------------------
  
  cutoff <- median(
    tmp$PDGFB,
    na.rm = TRUE
  )
  
  tmp$Score_group <- ifelse(
    tmp$PDGFB >= cutoff,
    "High",
    "Low"
  )
  
  tmp$Score_group <- factor(
    tmp$Score_group,
    levels = c(
      "Low",
      "High"
    )
  )
  
  # ---------------------------------
  # ② 5年KM
  # ---------------------------------
  
  fit <- survfit(
    Surv(time_5y, event_5y) ~ Score_group,
    data = tmp
  )
  
  logrank <- survdiff(
    Surv(time_5y, event_5y) ~ Score_group,
    data = tmp
  )
  
  logrank_p <- 1 - pchisq(
    logrank$chisq,
    df = length(logrank$n) - 1
  )
  
  # ---------------------------------
  # ③ 连续评分Cox：每增加1 SD
  # ---------------------------------
  
  tmp$score_Z <- as.numeric(
    scale(tmp$PDGFB)
  )
  
  cox <- coxph(
    Surv(time_5y, event_5y) ~ score_Z,
    data = tmp
  )
  
  cs <- summary(cox)
  
  HR <- cs$conf.int[
    1,
    "exp(coef)"
  ]
  
  CI_low <- cs$conf.int[
    1,
    "lower .95"
  ]
  
  CI_high <- cs$conf.int[
    1,
    "upper .95"
  ]
  
  cox_p <- cs$coefficients[
    1,
    "Pr(>|z|)"
  ]
  
  # ---------------------------------
  # ④ 提取KM曲线
  # ---------------------------------
  
  ss <- summary(fit)
  
  plot_dat <- data.frame(
    time = ss$time,
    surv = ss$surv,
    strata = ss$strata
  )
  
  plot_dat$strata <- sub(
    "Score_group=",
    "",
    plot_dat$strata
  )
  
  start_dat <- data.frame(
    time = 0,
    surv = 1,
    strata = c(
      "Low",
      "High"
    )
  )
  
  plot_dat <- bind_rows(
    start_dat,
    plot_dat
  )
  
  plot_dat$strata <- factor(
    plot_dat$strata,
    levels = c(
      "Low",
      "High"
    )
  )
  
  # ---------------------------------
  # ⑤ 绘制5年KM
  # ---------------------------------
  
  p <- ggplot(
    plot_dat,
    aes(
      x = time,
      y = surv,
      group = strata,
      color = strata
    )
  ) +
    geom_step(
      linewidth = 1.1
    ) +
    
    scale_color_manual(
      values = c(
        "Low" = "#3B6FB6",
        "High" = "#D94A4A"
      )
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        0.2
      ),
      expand = c(0, 0)
    ) +
    
    scale_x_continuous(
      limits = c(0, 60),
      breaks = c(
        0,
        12,
        24,
        36,
        48,
        60
      ),
      expand = c(0, 0)
    ) +
    
    theme_classic(
      base_size = 14
    ) +
    
    theme(
      aspect.ratio = 1,
      
      axis.title = element_text(
        size = 14,
        color = "black"
      ),
      
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      
      legend.title = element_text(
        size = 12
      ),
      
      legend.text = element_text(
        size = 11
      ),
      
      legend.position = c(
        0.80,
        0.78
      ),
      
      legend.background =
        element_blank(),
      
      plot.margin = margin(
        10,
        15,
        10,
        10
      )
    ) +
    
    labs(
      title = ca,
      x = "Time (months)",
      y = "Overall survival probability",
      color = "Macro_FCGBP"
    ) +
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.4,
      size = 4,
      color = "black",
      label = sprintf(
        paste0(
          "Log-rank P = %.3g\n",
          "HR per 1 SD = %.2f ",
          "(95%% CI %.2f–%.2f)\n",
          "Cox P = %.3g"
        ),
        logrank_p,
        HR,
        CI_low,
        CI_high,
        cox_p
      )
    )
  
  # 保存图
  ggsave(
    filename = paste0(
      "TCGA_Macro_FCGBP_OS_KM_5year/",
      ca,
      "_Macro_FCGBP_OS_5year.pdf"
    ),
    plot = p,
    width = 5,
    height = 5
  )
  
  # 保存结果
  result_5y[[ca]] <- data.frame(
    cancer = ca,
    n = nrow(tmp),
    death_5y = sum(
      tmp$event_5y == 1
    ),
    cutoff = cutoff,
    logrank_P = logrank_p,
    HR = HR,
    CI_low = CI_low,
    CI_high = CI_high,
    Cox_P = cox_p
  )
}

result_5y <- bind_rows(
  result_5y
)

result_5y$FDR <- p.adjust(
  result_5y$Cox_P,
  method = "BH"
)

write.csv(
  result_5y,
  "TCGA_Macro_FCGBP_OS_5year_results.csv",
  row.names = FALSE
)

result_5y
#############################################################DFI
###############################################################
result <- list()
six_cancer<-c("BRCA","KIRC","COAD","READ","GBM","HNSC","LUAD","LUSC")
for (ca in intersect(six_cancer,unique(dat1$cancer.type.abbreviation))) {
  #ca=six_cancer[1]
  cat("Processing:", ca, "\n")
  
  tmp <- dat1 %>%
    filter(
      cancer.type.abbreviation == ca,
      !is.na(PDGFB),
      !is.na(DFI.time),
      !is.na(DFI)
    )
  
  # 样本和死亡事件太少则跳过
  if (nrow(tmp) < 30 || sum(tmp$DFI == 1) < 10)
    next
  
  # ① KM：每个癌种内部按中位数分组
  cutoff <- median(tmp$PDGFB)
  
  tmp$Score_group <- ifelse(
    tmp$PDGFB >= cutoff,
    "High", "Low"
  )
  
  tmp$Score_group <- factor(
    tmp$Score_group,
    levels = c("Low", "High")
  )
  
  fit <- survfit(
    Surv(DFI.time, DFI) ~ Score_group,
    data = tmp
  )
  
  # ② log-rank P
  logrank <- survdiff(
    Surv(DFI.time, DFI) ~ Score_group,
    data = tmp
  )
  
  logrank_p <- 1 - pchisq(
    logrank$chisq,
    df = length(logrank$n) - 1
  )
  
  # ③ 连续评分Cox：每增加1 SD
  tmp$score_Z <- as.numeric(scale(tmp$PDGFB))
  
  cox <- coxph(
    Surv(DFI.time, DFI) ~ score_Z,
    data = tmp
  )
  
  cs <- summary(cox)
  
  HR      <- cs$conf.int[1, "exp(coef)"]
  CI_low  <- cs$conf.int[1, "lower .95"]
  CI_high <- cs$conf.int[1, "upper .95"]
  cox_p   <- cs$coefficients[1, "Pr(>|z|)"]
  
  # ④ 提取KM数据
  ss <- summary(fit)
  
  plot_dat <- data.frame(
    time = ss$time,
    surv = ss$surv,
    strata = ss$strata
  )
  
  plot_dat$strata <- sub(
    "Score_group=", "",
    plot_dat$strata
  )
  
  # 补充t=0，使曲线从1开始
  start_dat <- data.frame(
    time = 0,
    surv = 1,
    strata = c("Low", "High")
  )
  
  plot_dat <- bind_rows(start_dat, plot_dat)
  
  # ⑤ ggplot画KM
  p <- ggplot(
    plot_dat,
    aes(
      x = time,
      y = surv,
      group = strata,
      color = strata
    )
  ) +
    geom_step(linewidth = 1.1) +
    
    # 红蓝配色
    scale_color_manual(
      values = c(
        "Low"  = "#3B6FB6",   # 蓝色
        "High" = "#D94A4A"    # 红色
      )
    ) +
    
    # Y轴固定0-1
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = c(0, 0)
    ) +
    
    # X轴稍微留白
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) +
    
    theme_classic(base_size = 14) +
    
    # 保持方形
    theme(
      aspect.ratio = 1,
      
      axis.title = element_text(
        size = 14,
        color = "black"
      ),
      
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      
      legend.title = element_text(
        size = 12
      ),
      
      legend.text = element_text(
        size = 11
      ),
      
      legend.position = c(0.80, 0.78),
      legend.background = element_blank(),
      
      plot.margin = margin(
        10, 15, 10, 10
      )
    ) +
    
    labs(
      title = ca,
      x = "Time (days)",
      y = "Overall survival probability",
      color = "Macro_FCGBP"
    ) +
    
    # 统计结果放右上角
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.4,
      size = 4.0,
      color = "black",
      label = sprintf(
        "Log-rank P = %.3g",
        logrank_p
        
      )
    )
  
  ggsave(
    paste0(
      "TCGA_Macro_FCGBP_DFI_KM/",
      ca,
      "_Macro_FCGBP_DFI.pdf"
    ),
    p,
    width = 5,
    height = 5
  )
  
  result[[ca]] <- data.frame(
    cancer = ca,
    n = nrow(tmp),
    death = sum(tmp$DFI == 1),
    cutoff = cutoff,
    logrank_P = logrank_p,
    HR = HR,
    CI_low = CI_low,
    CI_high = CI_high,
    Cox_P = cox_p
  )
}

# 汇总所有癌种
result <- bind_rows(result) %>%
  mutate(
    FDR = p.adjust(Cox_P, method = "BH")
  ) %>%
  arrange(Cox_P)

write.csv(
  result,
  "TCGA_Macro_FCGBP_DFI_result.csv",
  row.names = FALSE
)

######################################
dir.create(
  "TCGA_Macro_FCGBP_DFI_KM_3year",
  showWarnings = FALSE
)

result_3y <- list()

for (ca in intersect(
  six_cancer,
  unique(dat1$cancer.type.abbreviation)
)) {
  
  cat("Processing 3-year DFI:", ca, "\n")
  
  tmp <- dat1 %>%
    filter(
      cancer.type.abbreviation == ca,
      !is.na(PDGFB),
      !is.na(DFI.time),
      !is.na(DFI)
    )
  
  # 天 -> 月
  tmp$DFI.month <- tmp$DFI.time / 30.44
  
  # ---------------------------------
  # 3年截尾：36个月
  # ---------------------------------
  
  # 随访时间超过36个月的统一截到36个月
  tmp$time_3y <- pmin(
    tmp$DFI.month,
    36
  )
  
  # 只有36个月以内发生的死亡算event
  # 36个月之后死亡，按36个月删失
  tmp$event_3y <- ifelse(
    tmp$DFI == 1 &
      tmp$DFI.month <= 36,
    1,
    0
  )
  
  # 样本数或3年内死亡事件太少则跳过
  if (
    nrow(tmp) < 30 ||
    sum(tmp$event_3y == 1) < 10
  ) {
    cat(
      "Skip:",
      ca,
      "n =", nrow(tmp),
      "3y deaths =", sum(tmp$event_3y),
      "\n"
    )
    next
  }
  
  # ---------------------------------
  # ① 每个癌种内部按中位数分High/Low
  # ---------------------------------
  
  cutoff <- median(
    tmp$PDGFB,
    na.rm = TRUE
  )
  
  tmp$Score_group <- ifelse(
    tmp$PDGFB >= cutoff,
    "High",
    "Low"
  )
  
  tmp$Score_group <- factor(
    tmp$Score_group,
    levels = c("Low", "High")
  )
  
  # ---------------------------------
  # ② 3年KM
  # ---------------------------------
  
  fit <- survfit(
    Surv(time_3y, event_3y) ~ Score_group,
    data = tmp
  )
  
  # log-rank
  logrank <- survdiff(
    Surv(time_3y, event_3y) ~ Score_group,
    data = tmp
  )
  
  logrank_p <- 1 - pchisq(
    logrank$chisq,
    df = length(logrank$n) - 1
  )
  
  # ---------------------------------
  # ③ 连续评分Cox：每增加1 SD
  # 同样使用3年截尾结局
  # ---------------------------------
  
  tmp$score_Z <- as.numeric(
    scale(tmp$PDGFB)
  )
  
  cox <- coxph(
    Surv(time_3y, event_3y) ~ score_Z,
    data = tmp
  )
  
  cs <- summary(cox)
  
  HR <- cs$conf.int[
    1,
    "exp(coef)"
  ]
  
  CI_low <- cs$conf.int[
    1,
    "lower .95"
  ]
  
  CI_high <- cs$conf.int[
    1,
    "upper .95"
  ]
  
  cox_p <- cs$coefficients[
    1,
    "Pr(>|z|)"
  ]
  
  # ---------------------------------
  # ④ 提取KM曲线
  # ---------------------------------
  
  ss <- summary(fit)
  
  plot_dat <- data.frame(
    time = ss$time,
    surv = ss$surv,
    strata = ss$strata
  )
  
  plot_dat$strata <- sub(
    "Score_group=",
    "",
    plot_dat$strata
  )
  
  # 从t=0、survival=1开始
  start_dat <- data.frame(
    time = 0,
    surv = 1,
    strata = c(
      "Low",
      "High"
    )
  )
  
  plot_dat <- bind_rows(
    start_dat,
    plot_dat
  )
  
  plot_dat$strata <- factor(
    plot_dat$strata,
    levels = c(
      "Low",
      "High"
    )
  )
  
  # ---------------------------------
  # ⑤ 绘制3年KM
  # ---------------------------------
  
  p <- ggplot(
    plot_dat,
    aes(
      x = time,
      y = surv,
      group = strata,
      color = strata
    )
  ) +
    geom_step(
      linewidth = 1.1
    ) +
    
    scale_color_manual(
      values = c(
        "Low" = "#3B6FB6",
        "High" = "#D94A4A"
      )
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        0.2
      ),
      expand = c(0, 0)
    ) +
    
    scale_x_continuous(
      limits = c(0, 36),
      breaks = c(
        0,
        12,
        24,
        36
      ),
      expand = c(0, 0)
    ) +
    
    theme_classic(
      base_size = 14
    ) +
    
    theme(
      aspect.ratio = 1,
      
      axis.title = element_text(
        size = 14,
        color = "black"
      ),
      
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      
      legend.title = element_text(
        size = 12
      ),
      
      legend.text = element_text(
        size = 11
      ),
      
      legend.position = c(
        0.80,
        0.78
      ),
      
      legend.background =
        element_blank(),
      
      plot.margin = margin(
        10,
        15,
        10,
        10
      )
    ) +
    
    labs(
      title = ca,
      x = "Time (months)",
      y = "Overall survival probability",
      color = "Macro_FCGBP"
    ) +
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.4,
      size = 4,
      color = "black",
      label = sprintf(
        paste0(
          "Log-rank P = %.3g\n",
          "HR per 1 SD = %.2f ",
          "(95%% CI %.2f–%.2f)\n",
          "Cox P = %.3g"
        ),
        logrank_p,
        HR,
        CI_low,
        CI_high,
        cox_p
      )
    )
  
  # 保存图
  ggsave(
    filename = paste0(
      "TCGA_Macro_FCGBP_DFI_KM_3year/",
      ca,
      "_Macro_FCGBP_DFI_3year.pdf"
    ),
    plot = p,
    width = 5,
    height = 5
  )
  
  # 保存统计结果
  result_3y[[ca]] <- data.frame(
    cancer = ca,
    n = nrow(tmp),
    death_3y = sum(
      tmp$event_3y == 1
    ),
    cutoff = cutoff,
    logrank_P = logrank_p,
    HR = HR,
    CI_low = CI_low,
    CI_high = CI_high,
    Cox_P = cox_p
  )
}

result_3y <- bind_rows(
  result_3y
)

result_3y$FDR <- p.adjust(
  result_3y$Cox_P,
  method = "BH"
)

write.csv(
  result_3y,
  "TCGA_Macro_FCGBP_DFI_3year_results.csv",
  row.names = FALSE
)

result_3y
####################################
dir.create(
  "TCGA_Macro_FCGBP_DFI_KM_5year",
  showWarnings = FALSE
)

result_5y <- list()

for (ca in intersect(
  six_cancer,
  unique(dat1$cancer.type.abbreviation)
)) {
  
  cat("Processing 5-year DFI:", ca, "\n")
  
  tmp <- dat1 %>%
    filter(
      cancer.type.abbreviation == ca,
      !is.na(PDGFB),
      !is.na(DFI.time),
      !is.na(DFI)
    )
  
  # 天 -> 月
  tmp$DFI.month <- tmp$DFI.time / 30.44
  
  # ---------------------------------
  # 5年截尾：60个月
  # ---------------------------------
  
  tmp$time_5y <- pmin(
    tmp$DFI.month,
    60
  )
  
  # 只有60个月以内死亡算event
  tmp$event_5y <- ifelse(
    tmp$DFI == 1 &
      tmp$DFI.month <= 60,
    1,
    0
  )
  
  # 样本数或5年内死亡事件太少则跳过
  if (
    nrow(tmp) < 30 ||
    sum(tmp$event_5y == 1) < 10
  ) {
    cat(
      "Skip:",
      ca,
      "n =", nrow(tmp),
      "5y deaths =", sum(tmp$event_5y),
      "\n"
    )
    next
  }
  
  # ---------------------------------
  # ① 中位数分High/Low
  # ---------------------------------
  
  cutoff <- median(
    tmp$PDGFB,
    na.rm = TRUE
  )
  
  tmp$Score_group <- ifelse(
    tmp$PDGFB >= cutoff,
    "High",
    "Low"
  )
  
  tmp$Score_group <- factor(
    tmp$Score_group,
    levels = c(
      "Low",
      "High"
    )
  )
  
  # ---------------------------------
  # ② 5年KM
  # ---------------------------------
  
  fit <- survfit(
    Surv(time_5y, event_5y) ~ Score_group,
    data = tmp
  )
  
  logrank <- survdiff(
    Surv(time_5y, event_5y) ~ Score_group,
    data = tmp
  )
  
  logrank_p <- 1 - pchisq(
    logrank$chisq,
    df = length(logrank$n) - 1
  )
  
  # ---------------------------------
  # ③ 连续评分Cox：每增加1 SD
  # ---------------------------------
  
  tmp$score_Z <- as.numeric(
    scale(tmp$PDGFB)
  )
  
  cox <- coxph(
    Surv(time_5y, event_5y) ~ score_Z,
    data = tmp
  )
  
  cs <- summary(cox)
  
  HR <- cs$conf.int[
    1,
    "exp(coef)"
  ]
  
  CI_low <- cs$conf.int[
    1,
    "lower .95"
  ]
  
  CI_high <- cs$conf.int[
    1,
    "upper .95"
  ]
  
  cox_p <- cs$coefficients[
    1,
    "Pr(>|z|)"
  ]
  
  # ---------------------------------
  # ④ 提取KM曲线
  # ---------------------------------
  
  ss <- summary(fit)
  
  plot_dat <- data.frame(
    time = ss$time,
    surv = ss$surv,
    strata = ss$strata
  )
  
  plot_dat$strata <- sub(
    "Score_group=",
    "",
    plot_dat$strata
  )
  
  start_dat <- data.frame(
    time = 0,
    surv = 1,
    strata = c(
      "Low",
      "High"
    )
  )
  
  plot_dat <- bind_rows(
    start_dat,
    plot_dat
  )
  
  plot_dat$strata <- factor(
    plot_dat$strata,
    levels = c(
      "Low",
      "High"
    )
  )
  
  # ---------------------------------
  # ⑤ 绘制5年KM
  # ---------------------------------
  
  p <- ggplot(
    plot_dat,
    aes(
      x = time,
      y = surv,
      group = strata,
      color = strata
    )
  ) +
    geom_step(
      linewidth = 1.1
    ) +
    
    scale_color_manual(
      values = c(
        "Low" = "#3B6FB6",
        "High" = "#D94A4A"
      )
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        0.2
      ),
      expand = c(0, 0)
    ) +
    
    scale_x_continuous(
      limits = c(0, 60),
      breaks = c(
        0,
        12,
        24,
        36,
        48,
        60
      ),
      expand = c(0, 0)
    ) +
    
    theme_classic(
      base_size = 14
    ) +
    
    theme(
      aspect.ratio = 1,
      
      axis.title = element_text(
        size = 14,
        color = "black"
      ),
      
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      
      legend.title = element_text(
        size = 12
      ),
      
      legend.text = element_text(
        size = 11
      ),
      
      legend.position = c(
        0.80,
        0.78
      ),
      
      legend.background =
        element_blank(),
      
      plot.margin = margin(
        10,
        15,
        10,
        10
      )
    ) +
    
    labs(
      title = ca,
      x = "Time (months)",
      y = "Overall survival probability",
      color = "Macro_FCGBP"
    ) +
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.4,
      size = 4,
      color = "black",
      label = sprintf(
        paste0(
          "Log-rank P = %.3g\n",
          "HR per 1 SD = %.2f ",
          "(95%% CI %.2f–%.2f)\n",
          "Cox P = %.3g"
        ),
        logrank_p,
        HR,
        CI_low,
        CI_high,
        cox_p
      )
    )
  
  # 保存图
  ggsave(
    filename = paste0(
      "TCGA_Macro_FCGBP_DFI_KM_5year/",
      ca,
      "_Macro_FCGBP_DFI_5year.pdf"
    ),
    plot = p,
    width = 5,
    height = 5
  )
  
  # 保存结果
  result_5y[[ca]] <- data.frame(
    cancer = ca,
    n = nrow(tmp),
    death_5y = sum(
      tmp$event_5y == 1
    ),
    cutoff = cutoff,
    logrank_P = logrank_p,
    HR = HR,
    CI_low = CI_low,
    CI_high = CI_high,
    Cox_P = cox_p
  )
}

result_5y <- bind_rows(
  result_5y
)

result_5y$FDR <- p.adjust(
  result_5y$Cox_P,
  method = "BH"
)

write.csv(
  result_5y,
  "TCGA_Macro_FCGBP_DFI_5year_results.csv",
  row.names = FALSE
)

result_5y
#############################################################PFI
###############################################################
result <- list()
six_cancer<-c("BRCA","KIRC","COAD","READ","GBM","HNSC","LUAD","LUSC")
for (ca in intersect(six_cancer,unique(dat1$cancer.type.abbreviation))) {
  #ca=six_cancer[1]
  cat("Processing:", ca, "\n")
  
  tmp <- dat1 %>%
    filter(
      cancer.type.abbreviation == ca,
      !is.na(PDGFB),
      !is.na(PFI.time),
      !is.na(PFI)
    )
  
  # 样本和死亡事件太少则跳过
  if (nrow(tmp) < 30 || sum(tmp$PFI == 1) < 10)
    next
  
  # ① KM：每个癌种内部按中位数分组
  cutoff <- median(tmp$PDGFB)
  
  tmp$Score_group <- ifelse(
    tmp$PDGFB >= cutoff,
    "High", "Low"
  )
  
  tmp$Score_group <- factor(
    tmp$Score_group,
    levels = c("Low", "High")
  )
  
  fit <- survfit(
    Surv(PFI.time, PFI) ~ Score_group,
    data = tmp
  )
  
  # ② log-rank P
  logrank <- survdiff(
    Surv(PFI.time, PFI) ~ Score_group,
    data = tmp
  )
  
  logrank_p <- 1 - pchisq(
    logrank$chisq,
    df = length(logrank$n) - 1
  )
  
  # ③ 连续评分Cox：每增加1 SD
  tmp$score_Z <- as.numeric(scale(tmp$PDGFB))
  
  cox <- coxph(
    Surv(PFI.time, PFI) ~ score_Z,
    data = tmp
  )
  
  cs <- summary(cox)
  
  HR      <- cs$conf.int[1, "exp(coef)"]
  CI_low  <- cs$conf.int[1, "lower .95"]
  CI_high <- cs$conf.int[1, "upper .95"]
  cox_p   <- cs$coefficients[1, "Pr(>|z|)"]
  
  # ④ 提取KM数据
  ss <- summary(fit)
  
  plot_dat <- data.frame(
    time = ss$time,
    surv = ss$surv,
    strata = ss$strata
  )
  
  plot_dat$strata <- sub(
    "Score_group=", "",
    plot_dat$strata
  )
  
  # 补充t=0，使曲线从1开始
  start_dat <- data.frame(
    time = 0,
    surv = 1,
    strata = c("Low", "High")
  )
  
  plot_dat <- bind_rows(start_dat, plot_dat)
  
  # ⑤ ggplot画KM
  p <- ggplot(
    plot_dat,
    aes(
      x = time,
      y = surv,
      group = strata,
      color = strata
    )
  ) +
    geom_step(linewidth = 1.1) +
    
    # 红蓝配色
    scale_color_manual(
      values = c(
        "Low"  = "#3B6FB6",   # 蓝色
        "High" = "#D94A4A"    # 红色
      )
    ) +
    
    # Y轴固定0-1
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = c(0, 0)
    ) +
    
    # X轴稍微留白
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) +
    
    theme_classic(base_size = 14) +
    
    # 保持方形
    theme(
      aspect.ratio = 1,
      
      axis.title = element_text(
        size = 14,
        color = "black"
      ),
      
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      
      legend.title = element_text(
        size = 12
      ),
      
      legend.text = element_text(
        size = 11
      ),
      
      legend.position = c(0.80, 0.78),
      legend.background = element_blank(),
      
      plot.margin = margin(
        10, 15, 10, 10
      )
    ) +
    
    labs(
      title = ca,
      x = "Time (days)",
      y = "Overall survival probability",
      color = "Macro_FCGBP"
    ) +
    
    # 统计结果放右上角
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.4,
      size = 4.0,
      color = "black",
      label = sprintf(
        "Log-rank P = %.3g",
        logrank_p
        
      )
    )
  
  ggsave(
    paste0(
      "TCGA_Macro_FCGBP_PFI_KM/",
      ca,
      "_Macro_FCGBP_PFI.pdf"
    ),
    p,
    width = 5,
    height = 5
  )
  
  result[[ca]] <- data.frame(
    cancer = ca,
    n = nrow(tmp),
    death = sum(tmp$PFI == 1),
    cutoff = cutoff,
    logrank_P = logrank_p,
    HR = HR,
    CI_low = CI_low,
    CI_high = CI_high,
    Cox_P = cox_p
  )
}

# 汇总所有癌种
result <- bind_rows(result) %>%
  mutate(
    FDR = p.adjust(Cox_P, method = "BH")
  ) %>%
  arrange(Cox_P)

write.csv(
  result,
  "TCGA_Macro_FCGBP_PFI_result.csv",
  row.names = FALSE
)

######################################
dir.create(
  "TCGA_Macro_FCGBP_PFI_KM_3year",
  showWarnings = FALSE
)

result_3y <- list()

for (ca in intersect(
  six_cancer,
  unique(dat1$cancer.type.abbreviation)
)) {
  
  cat("Processing 3-year PFI:", ca, "\n")
  
  tmp <- dat1 %>%
    filter(
      cancer.type.abbreviation == ca,
      !is.na(PDGFB),
      !is.na(PFI.time),
      !is.na(PFI)
    )
  
  # 天 -> 月
  tmp$PFI.month <- tmp$PFI.time / 30.44
  
  # ---------------------------------
  # 3年截尾：36个月
  # ---------------------------------
  
  # 随访时间超过36个月的统一截到36个月
  tmp$time_3y <- pmin(
    tmp$PFI.month,
    36
  )
  
  # 只有36个月以内发生的死亡算event
  # 36个月之后死亡，按36个月删失
  tmp$event_3y <- ifelse(
    tmp$PFI == 1 &
      tmp$PFI.month <= 36,
    1,
    0
  )
  
  # 样本数或3年内死亡事件太少则跳过
  if (
    nrow(tmp) < 30 ||
    sum(tmp$event_3y == 1) < 10
  ) {
    cat(
      "Skip:",
      ca,
      "n =", nrow(tmp),
      "3y deaths =", sum(tmp$event_3y),
      "\n"
    )
    next
  }
  
  # ---------------------------------
  # ① 每个癌种内部按中位数分High/Low
  # ---------------------------------
  
  cutoff <- median(
    tmp$PDGFB,
    na.rm = TRUE
  )
  
  tmp$Score_group <- ifelse(
    tmp$PDGFB >= cutoff,
    "High",
    "Low"
  )
  
  tmp$Score_group <- factor(
    tmp$Score_group,
    levels = c("Low", "High")
  )
  
  # ---------------------------------
  # ② 3年KM
  # ---------------------------------
  
  fit <- survfit(
    Surv(time_3y, event_3y) ~ Score_group,
    data = tmp
  )
  
  # log-rank
  logrank <- survdiff(
    Surv(time_3y, event_3y) ~ Score_group,
    data = tmp
  )
  
  logrank_p <- 1 - pchisq(
    logrank$chisq,
    df = length(logrank$n) - 1
  )
  
  # ---------------------------------
  # ③ 连续评分Cox：每增加1 SD
  # 同样使用3年截尾结局
  # ---------------------------------
  
  tmp$score_Z <- as.numeric(
    scale(tmp$PDGFB)
  )
  
  cox <- coxph(
    Surv(time_3y, event_3y) ~ score_Z,
    data = tmp
  )
  
  cs <- summary(cox)
  
  HR <- cs$conf.int[
    1,
    "exp(coef)"
  ]
  
  CI_low <- cs$conf.int[
    1,
    "lower .95"
  ]
  
  CI_high <- cs$conf.int[
    1,
    "upper .95"
  ]
  
  cox_p <- cs$coefficients[
    1,
    "Pr(>|z|)"
  ]
  
  # ---------------------------------
  # ④ 提取KM曲线
  # ---------------------------------
  
  ss <- summary(fit)
  
  plot_dat <- data.frame(
    time = ss$time,
    surv = ss$surv,
    strata = ss$strata
  )
  
  plot_dat$strata <- sub(
    "Score_group=",
    "",
    plot_dat$strata
  )
  
  # 从t=0、survival=1开始
  start_dat <- data.frame(
    time = 0,
    surv = 1,
    strata = c(
      "Low",
      "High"
    )
  )
  
  plot_dat <- bind_rows(
    start_dat,
    plot_dat
  )
  
  plot_dat$strata <- factor(
    plot_dat$strata,
    levels = c(
      "Low",
      "High"
    )
  )
  
  # ---------------------------------
  # ⑤ 绘制3年KM
  # ---------------------------------
  
  p <- ggplot(
    plot_dat,
    aes(
      x = time,
      y = surv,
      group = strata,
      color = strata
    )
  ) +
    geom_step(
      linewidth = 1.1
    ) +
    
    scale_color_manual(
      values = c(
        "Low" = "#3B6FB6",
        "High" = "#D94A4A"
      )
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        0.2
      ),
      expand = c(0, 0)
    ) +
    
    scale_x_continuous(
      limits = c(0, 36),
      breaks = c(
        0,
        12,
        24,
        36
      ),
      expand = c(0, 0)
    ) +
    
    theme_classic(
      base_size = 14
    ) +
    
    theme(
      aspect.ratio = 1,
      
      axis.title = element_text(
        size = 14,
        color = "black"
      ),
      
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      
      legend.title = element_text(
        size = 12
      ),
      
      legend.text = element_text(
        size = 11
      ),
      
      legend.position = c(
        0.80,
        0.78
      ),
      
      legend.background =
        element_blank(),
      
      plot.margin = margin(
        10,
        15,
        10,
        10
      )
    ) +
    
    labs(
      title = ca,
      x = "Time (months)",
      y = "Overall survival probability",
      color = "Macro_FCGBP"
    ) +
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.4,
      size = 4,
      color = "black",
      label = sprintf(
        paste0(
          "Log-rank P = %.3g\n",
          "HR per 1 SD = %.2f ",
          "(95%% CI %.2f–%.2f)\n",
          "Cox P = %.3g"
        ),
        logrank_p,
        HR,
        CI_low,
        CI_high,
        cox_p
      )
    )
  
  # 保存图
  ggsave(
    filename = paste0(
      "TCGA_Macro_FCGBP_PFI_KM_3year/",
      ca,
      "_Macro_FCGBP_PFI_3year.pdf"
    ),
    plot = p,
    width = 5,
    height = 5
  )
  
  # 保存统计结果
  result_3y[[ca]] <- data.frame(
    cancer = ca,
    n = nrow(tmp),
    death_3y = sum(
      tmp$event_3y == 1
    ),
    cutoff = cutoff,
    logrank_P = logrank_p,
    HR = HR,
    CI_low = CI_low,
    CI_high = CI_high,
    Cox_P = cox_p
  )
}

result_3y <- bind_rows(
  result_3y
)

result_3y$FDR <- p.adjust(
  result_3y$Cox_P,
  method = "BH"
)

write.csv(
  result_3y,
  "TCGA_Macro_FCGBP_PFI_3year_results.csv",
  row.names = FALSE
)

result_3y
####################################
dir.create(
  "TCGA_Macro_FCGBP_PFI_KM_5year",
  showWarnings = FALSE
)

result_5y <- list()

for (ca in intersect(
  six_cancer,
  unique(dat1$cancer.type.abbreviation)
)) {
  
  cat("Processing 5-year PFI:", ca, "\n")
  
  tmp <- dat1 %>%
    filter(
      cancer.type.abbreviation == ca,
      !is.na(PDGFB),
      !is.na(PFI.time),
      !is.na(PFI)
    )
  
  # 天 -> 月
  tmp$PFI.month <- tmp$PFI.time / 30.44
  
  # ---------------------------------
  # 5年截尾：60个月
  # ---------------------------------
  
  tmp$time_5y <- pmin(
    tmp$PFI.month,
    60
  )
  
  # 只有60个月以内死亡算event
  tmp$event_5y <- ifelse(
    tmp$PFI == 1 &
      tmp$PFI.month <= 60,
    1,
    0
  )
  
  # 样本数或5年内死亡事件太少则跳过
  if (
    nrow(tmp) < 30 ||
    sum(tmp$event_5y == 1) < 10
  ) {
    cat(
      "Skip:",
      ca,
      "n =", nrow(tmp),
      "5y deaths =", sum(tmp$event_5y),
      "\n"
    )
    next
  }
  
  # ---------------------------------
  # ① 中位数分High/Low
  # ---------------------------------
  
  cutoff <- median(
    tmp$PDGFB,
    na.rm = TRUE
  )
  
  tmp$Score_group <- ifelse(
    tmp$PDGFB >= cutoff,
    "High",
    "Low"
  )
  
  tmp$Score_group <- factor(
    tmp$Score_group,
    levels = c(
      "Low",
      "High"
    )
  )
  
  # ---------------------------------
  # ② 5年KM
  # ---------------------------------
  
  fit <- survfit(
    Surv(time_5y, event_5y) ~ Score_group,
    data = tmp
  )
  
  logrank <- survdiff(
    Surv(time_5y, event_5y) ~ Score_group,
    data = tmp
  )
  
  logrank_p <- 1 - pchisq(
    logrank$chisq,
    df = length(logrank$n) - 1
  )
  
  # ---------------------------------
  # ③ 连续评分Cox：每增加1 SD
  # ---------------------------------
  
  tmp$score_Z <- as.numeric(
    scale(tmp$PDGFB)
  )
  
  cox <- coxph(
    Surv(time_5y, event_5y) ~ score_Z,
    data = tmp
  )
  
  cs <- summary(cox)
  
  HR <- cs$conf.int[
    1,
    "exp(coef)"
  ]
  
  CI_low <- cs$conf.int[
    1,
    "lower .95"
  ]
  
  CI_high <- cs$conf.int[
    1,
    "upper .95"
  ]
  
  cox_p <- cs$coefficients[
    1,
    "Pr(>|z|)"
  ]
  
  # ---------------------------------
  # ④ 提取KM曲线
  # ---------------------------------
  
  ss <- summary(fit)
  
  plot_dat <- data.frame(
    time = ss$time,
    surv = ss$surv,
    strata = ss$strata
  )
  
  plot_dat$strata <- sub(
    "Score_group=",
    "",
    plot_dat$strata
  )
  
  start_dat <- data.frame(
    time = 0,
    surv = 1,
    strata = c(
      "Low",
      "High"
    )
  )
  
  plot_dat <- bind_rows(
    start_dat,
    plot_dat
  )
  
  plot_dat$strata <- factor(
    plot_dat$strata,
    levels = c(
      "Low",
      "High"
    )
  )
  
  # ---------------------------------
  # ⑤ 绘制5年KM
  # ---------------------------------
  
  p <- ggplot(
    plot_dat,
    aes(
      x = time,
      y = surv,
      group = strata,
      color = strata
    )
  ) +
    geom_step(
      linewidth = 1.1
    ) +
    
    scale_color_manual(
      values = c(
        "Low" = "#3B6FB6",
        "High" = "#D94A4A"
      )
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        0.2
      ),
      expand = c(0, 0)
    ) +
    
    scale_x_continuous(
      limits = c(0, 60),
      breaks = c(
        0,
        12,
        24,
        36,
        48,
        60
      ),
      expand = c(0, 0)
    ) +
    
    theme_classic(
      base_size = 14
    ) +
    
    theme(
      aspect.ratio = 1,
      
      axis.title = element_text(
        size = 14,
        color = "black"
      ),
      
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      
      legend.title = element_text(
        size = 12
      ),
      
      legend.text = element_text(
        size = 11
      ),
      
      legend.position = c(
        0.80,
        0.78
      ),
      
      legend.background =
        element_blank(),
      
      plot.margin = margin(
        10,
        15,
        10,
        10
      )
    ) +
    
    labs(
      title = ca,
      x = "Time (months)",
      y = "Overall survival probability",
      color = "Macro_FCGBP"
    ) +
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.4,
      size = 4,
      color = "black",
      label = sprintf(
        paste0(
          "Log-rank P = %.3g\n",
          "HR per 1 SD = %.2f ",
          "(95%% CI %.2f–%.2f)\n",
          "Cox P = %.3g"
        ),
        logrank_p,
        HR,
        CI_low,
        CI_high,
        cox_p
      )
    )
  
  # 保存图
  ggsave(
    filename = paste0(
      "TCGA_Macro_FCGBP_PFI_KM_5year/",
      ca,
      "_Macro_FCGBP_PFI_5year.pdf"
    ),
    plot = p,
    width = 5,
    height = 5
  )
  
  # 保存结果
  result_5y[[ca]] <- data.frame(
    cancer = ca,
    n = nrow(tmp),
    death_5y = sum(
      tmp$event_5y == 1
    ),
    cutoff = cutoff,
    logrank_P = logrank_p,
    HR = HR,
    CI_low = CI_low,
    CI_high = CI_high,
    Cox_P = cox_p
  )
}

result_5y <- bind_rows(
  result_5y
)

result_5y$FDR <- p.adjust(
  result_5y$Cox_P,
  method = "BH"
)

write.csv(
  result_5y,
  "TCGA_Macro_FCGBP_PFI_5year_results.csv",
  row.names = FALSE
)

result_5y
#############################################################PFI
###############################################################
result <- list()
six_cancer<-c("BRCA","KIRC","COAD","READ","GBM","HNSC","LUAD","LUSC")
for (ca in intersect(six_cancer,unique(dat1$cancer.type.abbreviation))) {
  #ca=six_cancer[1]
  cat("Processing:", ca, "\n")
  
  tmp <- dat1 %>%
    filter(
      cancer.type.abbreviation == ca,
      !is.na(PDGFB),
      !is.na(PFI.time),
      !is.na(PFI)
    )
  
  # 样本和死亡事件太少则跳过
  if (nrow(tmp) < 30 || sum(tmp$PFI == 1) < 10)
    next
  
  # ① KM：每个癌种内部按中位数分组
  cutoff <- median(tmp$PDGFB)
  
  tmp$Score_group <- ifelse(
    tmp$PDGFB >= cutoff,
    "High", "Low"
  )
  
  tmp$Score_group <- factor(
    tmp$Score_group,
    levels = c("Low", "High")
  )
  
  fit <- survfit(
    Surv(PFI.time, PFI) ~ Score_group,
    data = tmp
  )
  
  # ② log-rank P
  logrank <- survdiff(
    Surv(PFI.time, PFI) ~ Score_group,
    data = tmp
  )
  
  logrank_p <- 1 - pchisq(
    logrank$chisq,
    df = length(logrank$n) - 1
  )
  
  # ③ 连续评分Cox：每增加1 SD
  tmp$score_Z <- as.numeric(scale(tmp$PDGFB))
  
  cox <- coxph(
    Surv(PFI.time, PFI) ~ score_Z,
    data = tmp
  )
  
  cs <- summary(cox)
  
  HR      <- cs$conf.int[1, "exp(coef)"]
  CI_low  <- cs$conf.int[1, "lower .95"]
  CI_high <- cs$conf.int[1, "upper .95"]
  cox_p   <- cs$coefficients[1, "Pr(>|z|)"]
  
  # ④ 提取KM数据
  ss <- summary(fit)
  
  plot_dat <- data.frame(
    time = ss$time,
    surv = ss$surv,
    strata = ss$strata
  )
  
  plot_dat$strata <- sub(
    "Score_group=", "",
    plot_dat$strata
  )
  
  # 补充t=0，使曲线从1开始
  start_dat <- data.frame(
    time = 0,
    surv = 1,
    strata = c("Low", "High")
  )
  
  plot_dat <- bind_rows(start_dat, plot_dat)
  
  # ⑤ ggplot画KM
  p <- ggplot(
    plot_dat,
    aes(
      x = time,
      y = surv,
      group = strata,
      color = strata
    )
  ) +
    geom_step(linewidth = 1.1) +
    
    # 红蓝配色
    scale_color_manual(
      values = c(
        "Low"  = "#3B6FB6",   # 蓝色
        "High" = "#D94A4A"    # 红色
      )
    ) +
    
    # Y轴固定0-1
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = c(0, 0)
    ) +
    
    # X轴稍微留白
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) +
    
    theme_classic(base_size = 14) +
    
    # 保持方形
    theme(
      aspect.ratio = 1,
      
      axis.title = element_text(
        size = 14,
        color = "black"
      ),
      
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      
      legend.title = element_text(
        size = 12
      ),
      
      legend.text = element_text(
        size = 11
      ),
      
      legend.position = c(0.80, 0.78),
      legend.background = element_blank(),
      
      plot.margin = margin(
        10, 15, 10, 10
      )
    ) +
    
    labs(
      title = ca,
      x = "Time (days)",
      y = "Overall survival probability",
      color = "Macro_FCGBP"
    ) +
    
    # 统计结果放右上角
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.4,
      size = 4.0,
      color = "black",
      label = sprintf(
        "Log-rank P = %.3g",
        logrank_p
        
      )
    )
  
  ggsave(
    paste0(
      "TCGA_Macro_FCGBP_PFI_KM/",
      ca,
      "_Macro_FCGBP_PFI.pdf"
    ),
    p,
    width = 5,
    height = 5
  )
  
  result[[ca]] <- data.frame(
    cancer = ca,
    n = nrow(tmp),
    death = sum(tmp$PFI == 1),
    cutoff = cutoff,
    logrank_P = logrank_p,
    HR = HR,
    CI_low = CI_low,
    CI_high = CI_high,
    Cox_P = cox_p
  )
}

# 汇总所有癌种
result <- bind_rows(result) %>%
  mutate(
    FDR = p.adjust(Cox_P, method = "BH")
  ) %>%
  arrange(Cox_P)

write.csv(
  result,
  "TCGA_Macro_FCGBP_PFI_result.csv",
  row.names = FALSE
)

######################################
dir.create(
  "TCGA_Macro_FCGBP_PFI_KM_3year",
  showWarnings = FALSE
)

result_3y <- list()

for (ca in intersect(
  six_cancer,
  unique(dat1$cancer.type.abbreviation)
)) {
  
  cat("Processing 3-year PFI:", ca, "\n")
  
  tmp <- dat1 %>%
    filter(
      cancer.type.abbreviation == ca,
      !is.na(PDGFB),
      !is.na(PFI.time),
      !is.na(PFI)
    )
  
  # 天 -> 月
  tmp$PFI.month <- tmp$PFI.time / 30.44
  
  # ---------------------------------
  # 3年截尾：36个月
  # ---------------------------------
  
  # 随访时间超过36个月的统一截到36个月
  tmp$time_3y <- pmin(
    tmp$PFI.month,
    36
  )
  
  # 只有36个月以内发生的死亡算event
  # 36个月之后死亡，按36个月删失
  tmp$event_3y <- ifelse(
    tmp$PFI == 1 &
      tmp$PFI.month <= 36,
    1,
    0
  )
  
  # 样本数或3年内死亡事件太少则跳过
  if (
    nrow(tmp) < 30 ||
    sum(tmp$event_3y == 1) < 10
  ) {
    cat(
      "Skip:",
      ca,
      "n =", nrow(tmp),
      "3y deaths =", sum(tmp$event_3y),
      "\n"
    )
    next
  }
  
  # ---------------------------------
  # ① 每个癌种内部按中位数分High/Low
  # ---------------------------------
  
  cutoff <- median(
    tmp$PDGFB,
    na.rm = TRUE
  )
  
  tmp$Score_group <- ifelse(
    tmp$PDGFB >= cutoff,
    "High",
    "Low"
  )
  
  tmp$Score_group <- factor(
    tmp$Score_group,
    levels = c("Low", "High")
  )
  
  # ---------------------------------
  # ② 3年KM
  # ---------------------------------
  
  fit <- survfit(
    Surv(time_3y, event_3y) ~ Score_group,
    data = tmp
  )
  
  # log-rank
  logrank <- survdiff(
    Surv(time_3y, event_3y) ~ Score_group,
    data = tmp
  )
  
  logrank_p <- 1 - pchisq(
    logrank$chisq,
    df = length(logrank$n) - 1
  )
  
  # ---------------------------------
  # ③ 连续评分Cox：每增加1 SD
  # 同样使用3年截尾结局
  # ---------------------------------
  
  tmp$score_Z <- as.numeric(
    scale(tmp$PDGFB)
  )
  
  cox <- coxph(
    Surv(time_3y, event_3y) ~ score_Z,
    data = tmp
  )
  
  cs <- summary(cox)
  
  HR <- cs$conf.int[
    1,
    "exp(coef)"
  ]
  
  CI_low <- cs$conf.int[
    1,
    "lower .95"
  ]
  
  CI_high <- cs$conf.int[
    1,
    "upper .95"
  ]
  
  cox_p <- cs$coefficients[
    1,
    "Pr(>|z|)"
  ]
  
  # ---------------------------------
  # ④ 提取KM曲线
  # ---------------------------------
  
  ss <- summary(fit)
  
  plot_dat <- data.frame(
    time = ss$time,
    surv = ss$surv,
    strata = ss$strata
  )
  
  plot_dat$strata <- sub(
    "Score_group=",
    "",
    plot_dat$strata
  )
  
  # 从t=0、survival=1开始
  start_dat <- data.frame(
    time = 0,
    surv = 1,
    strata = c(
      "Low",
      "High"
    )
  )
  
  plot_dat <- bind_rows(
    start_dat,
    plot_dat
  )
  
  plot_dat$strata <- factor(
    plot_dat$strata,
    levels = c(
      "Low",
      "High"
    )
  )
  
  # ---------------------------------
  # ⑤ 绘制3年KM
  # ---------------------------------
  
  p <- ggplot(
    plot_dat,
    aes(
      x = time,
      y = surv,
      group = strata,
      color = strata
    )
  ) +
    geom_step(
      linewidth = 1.1
    ) +
    
    scale_color_manual(
      values = c(
        "Low" = "#3B6FB6",
        "High" = "#D94A4A"
      )
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        0.2
      ),
      expand = c(0, 0)
    ) +
    
    scale_x_continuous(
      limits = c(0, 36),
      breaks = c(
        0,
        12,
        24,
        36
      ),
      expand = c(0, 0)
    ) +
    
    theme_classic(
      base_size = 14
    ) +
    
    theme(
      aspect.ratio = 1,
      
      axis.title = element_text(
        size = 14,
        color = "black"
      ),
      
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      
      legend.title = element_text(
        size = 12
      ),
      
      legend.text = element_text(
        size = 11
      ),
      
      legend.position = c(
        0.80,
        0.78
      ),
      
      legend.background =
        element_blank(),
      
      plot.margin = margin(
        10,
        15,
        10,
        10
      )
    ) +
    
    labs(
      title = ca,
      x = "Time (months)",
      y = "Overall survival probability",
      color = "Macro_FCGBP"
    ) +
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.4,
      size = 4,
      color = "black",
      label = sprintf(
        paste0(
          "Log-rank P = %.3g\n",
          "HR per 1 SD = %.2f ",
          "(95%% CI %.2f–%.2f)\n",
          "Cox P = %.3g"
        ),
        logrank_p,
        HR,
        CI_low,
        CI_high,
        cox_p
      )
    )
  
  # 保存图
  ggsave(
    filename = paste0(
      "TCGA_Macro_FCGBP_PFI_KM_3year/",
      ca,
      "_Macro_FCGBP_PFI_3year.pdf"
    ),
    plot = p,
    width = 5,
    height = 5
  )
  
  # 保存统计结果
  result_3y[[ca]] <- data.frame(
    cancer = ca,
    n = nrow(tmp),
    death_3y = sum(
      tmp$event_3y == 1
    ),
    cutoff = cutoff,
    logrank_P = logrank_p,
    HR = HR,
    CI_low = CI_low,
    CI_high = CI_high,
    Cox_P = cox_p
  )
}

result_3y <- bind_rows(
  result_3y
)

result_3y$FDR <- p.adjust(
  result_3y$Cox_P,
  method = "BH"
)

write.csv(
  result_3y,
  "TCGA_Macro_FCGBP_PFI_3year_results.csv",
  row.names = FALSE
)

result_3y
####################################
dir.create(
  "TCGA_Macro_FCGBP_PFI_KM_5year",
  showWarnings = FALSE
)

result_5y <- list()

for (ca in intersect(
  six_cancer,
  unique(dat1$cancer.type.abbreviation)
)) {
  
  cat("Processing 5-year PFI:", ca, "\n")
  
  tmp <- dat1 %>%
    filter(
      cancer.type.abbreviation == ca,
      !is.na(PDGFB),
      !is.na(PFI.time),
      !is.na(PFI)
    )
  
  # 天 -> 月
  tmp$PFI.month <- tmp$PFI.time / 30.44
  
  # ---------------------------------
  # 5年截尾：60个月
  # ---------------------------------
  
  tmp$time_5y <- pmin(
    tmp$PFI.month,
    60
  )
  
  # 只有60个月以内死亡算event
  tmp$event_5y <- ifelse(
    tmp$PFI == 1 &
      tmp$PFI.month <= 60,
    1,
    0
  )
  
  # 样本数或5年内死亡事件太少则跳过
  if (
    nrow(tmp) < 30 ||
    sum(tmp$event_5y == 1) < 10
  ) {
    cat(
      "Skip:",
      ca,
      "n =", nrow(tmp),
      "5y deaths =", sum(tmp$event_5y),
      "\n"
    )
    next
  }
  
  # ---------------------------------
  # ① 中位数分High/Low
  # ---------------------------------
  
  cutoff <- median(
    tmp$PDGFB,
    na.rm = TRUE
  )
  
  tmp$Score_group <- ifelse(
    tmp$PDGFB >= cutoff,
    "High",
    "Low"
  )
  
  tmp$Score_group <- factor(
    tmp$Score_group,
    levels = c(
      "Low",
      "High"
    )
  )
  
  # ---------------------------------
  # ② 5年KM
  # ---------------------------------
  
  fit <- survfit(
    Surv(time_5y, event_5y) ~ Score_group,
    data = tmp
  )
  
  logrank <- survdiff(
    Surv(time_5y, event_5y) ~ Score_group,
    data = tmp
  )
  
  logrank_p <- 1 - pchisq(
    logrank$chisq,
    df = length(logrank$n) - 1
  )
  
  # ---------------------------------
  # ③ 连续评分Cox：每增加1 SD
  # ---------------------------------
  
  tmp$score_Z <- as.numeric(
    scale(tmp$PDGFB)
  )
  
  cox <- coxph(
    Surv(time_5y, event_5y) ~ score_Z,
    data = tmp
  )
  
  cs <- summary(cox)
  
  HR <- cs$conf.int[
    1,
    "exp(coef)"
  ]
  
  CI_low <- cs$conf.int[
    1,
    "lower .95"
  ]
  
  CI_high <- cs$conf.int[
    1,
    "upper .95"
  ]
  
  cox_p <- cs$coefficients[
    1,
    "Pr(>|z|)"
  ]
  
  # ---------------------------------
  # ④ 提取KM曲线
  # ---------------------------------
  
  ss <- summary(fit)
  
  plot_dat <- data.frame(
    time = ss$time,
    surv = ss$surv,
    strata = ss$strata
  )
  
  plot_dat$strata <- sub(
    "Score_group=",
    "",
    plot_dat$strata
  )
  
  start_dat <- data.frame(
    time = 0,
    surv = 1,
    strata = c(
      "Low",
      "High"
    )
  )
  
  plot_dat <- bind_rows(
    start_dat,
    plot_dat
  )
  
  plot_dat$strata <- factor(
    plot_dat$strata,
    levels = c(
      "Low",
      "High"
    )
  )
  
  # ---------------------------------
  # ⑤ 绘制5年KM
  # ---------------------------------
  
  p <- ggplot(
    plot_dat,
    aes(
      x = time,
      y = surv,
      group = strata,
      color = strata
    )
  ) +
    geom_step(
      linewidth = 1.1
    ) +
    
    scale_color_manual(
      values = c(
        "Low" = "#3B6FB6",
        "High" = "#D94A4A"
      )
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        0.2
      ),
      expand = c(0, 0)
    ) +
    
    scale_x_continuous(
      limits = c(0, 60),
      breaks = c(
        0,
        12,
        24,
        36,
        48,
        60
      ),
      expand = c(0, 0)
    ) +
    
    theme_classic(
      base_size = 14
    ) +
    
    theme(
      aspect.ratio = 1,
      
      axis.title = element_text(
        size = 14,
        color = "black"
      ),
      
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      
      legend.title = element_text(
        size = 12
      ),
      
      legend.text = element_text(
        size = 11
      ),
      
      legend.position = c(
        0.80,
        0.78
      ),
      
      legend.background =
        element_blank(),
      
      plot.margin = margin(
        10,
        15,
        10,
        10
      )
    ) +
    
    labs(
      title = ca,
      x = "Time (months)",
      y = "Overall survival probability",
      color = "Macro_FCGBP"
    ) +
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.4,
      size = 4,
      color = "black",
      label = sprintf(
        paste0(
          "Log-rank P = %.3g\n",
          "HR per 1 SD = %.2f ",
          "(95%% CI %.2f–%.2f)\n",
          "Cox P = %.3g"
        ),
        logrank_p,
        HR,
        CI_low,
        CI_high,
        cox_p
      )
    )
  
  # 保存图
  ggsave(
    filename = paste0(
      "TCGA_Macro_FCGBP_PFI_KM_5year/",
      ca,
      "_Macro_FCGBP_PFI_5year.pdf"
    ),
    plot = p,
    width = 5,
    height = 5
  )
  
  # 保存结果
  result_5y[[ca]] <- data.frame(
    cancer = ca,
    n = nrow(tmp),
    death_5y = sum(
      tmp$event_5y == 1
    ),
    cutoff = cutoff,
    logrank_P = logrank_p,
    HR = HR,
    CI_low = CI_low,
    CI_high = CI_high,
    Cox_P = cox_p
  )
}

result_5y <- bind_rows(
  result_5y
)

result_5y$FDR <- p.adjust(
  result_5y$Cox_P,
  method = "BH"
)

write.csv(
  result_5y,
  "TCGA_Macro_FCGBP_PFI_5year_results.csv",
  row.names = FALSE
)

result_5y