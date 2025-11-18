#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  matches <- grep(file_arg, cmd_args)
  if (length(matches) > 0) {
    return(normalizePath(sub(file_arg, "", cmd_args[matches[1]])))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
  stop("Cannot determine script location; please run with Rscript.")
}

script_path <- get_script_path()
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))

run_analysis_scripts <- function(project_root) {
  scripts <- c(
    "R/scripts/analyze_1all.R",
    "R/scripts/tests_1all.R",
    "R/scripts/analyze_2before_after.R",
    "R/scripts/tests_2before_after.R"
  )
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_root)
  for (script in scripts) {
    system2("Rscript", script, stdout = TRUE, stderr = TRUE)
  }
}

read_dataset <- function(file_path, skip_rows = 0) {
  read_csv(file_path, show_col_types = FALSE, skip = skip_rows) %>%
    transmute(
      CT = as.numeric(.data$CT),
      PET = as.numeric(.data$PET)
    ) %>%
    tidyr::drop_na(CT, PET)
}

format_p <- function(value) {
  if (is.null(value) || is.na(value)) return("NA")
  sprintf("%.3g", value)
}

datasets <- list(
  list(
    name = "1all",
    label = "1all データ",
    file = file.path(project_root, "data", "raw", "1all_data.csv"),
    skip = 0L,
    figure_prefix = "1all"
  ),
  list(
    name = "2before_after",
    label = "2before_after データ",
    file = file.path(project_root, "data", "raw", "2before_after.csv"),
    skip = 2L,
    figure_prefix = "2before_after"
  )
)

run_analysis_scripts(project_root)

results <- lapply(datasets, function(ds) {
  data <- read_dataset(ds$file, ds$skip)
  if (nrow(data) < 3) {
    stop("Dataset ", ds$name, " has fewer than 3 rows.")
  }
  lm_fit <- lm(CT ~ PET, data = data)
  cor_test <- cor.test(data$CT, data$PET, method = "pearson")
  diff_values <- data$PET - data$CT
  diff_nonzero <- diff_values[diff_values != 0]

  if (length(diff_nonzero) == 0) {
    sign_test <- NULL
  } else {
    sign_test <- binom.test(sum(diff_nonzero > 0), length(diff_nonzero), alternative = "greater")
  }

  if (length(unique(diff_values)) <= 1) {
    wilcox_res <- NULL
  } else {
    wilcox_res <- suppressWarnings(
      wilcox.test(data$PET, data$CT, alternative = "greater", paired = TRUE, exact = FALSE)
    )
  }

  list(
    dataset = ds,
    n = nrow(data),
    cor = cor_test$estimate,
    cor_p = cor_test$p.value,
    lm_intercept = coef(lm_fit)[1],
    lm_slope = coef(lm_fit)[2],
    adj_r2 = summary(lm_fit)$adj.r.squared,
    sign_p = if (!is.null(sign_test)) sign_test$p.value else NA_real_,
    wilcox_p = if (!is.null(wilcox_res)) wilcox_res$p.value else NA_real_
  )
})

report_path <- file.path(project_root, "output", "reports", "summary.md")
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

summary_table <- c(
  "| データセット | n | Pearson r | r の p値 | Sign test p値 | Wilcoxon p値 |",
  "| --- | --- | --- | --- | --- | --- |"
)
for (res in results) {
  summary_table <- c(
    summary_table,
    sprintf(
      "| %s | %d | %.3f | %.3g | %s | %s |",
      res$dataset$label,
      res$n,
      res$cor,
      res$cor_p,
      format_p(res$sign_p),
      format_p(res$wilcox_p)
    )
  )
}

report_lines <- c("# 1all / 2before_after 解析サマリー", "")
report_lines <- c(report_lines, "## 指標サマリー", summary_table, "")

for (res in results) {
  prefix <- res$dataset$figure_prefix
  report_lines <- c(
    report_lines,
    sprintf("## %s", res$dataset$label),
    "",
    sprintf("- サンプル数: %d", res$n),
    sprintf("- Pearson 相関係数: %.3f (p = %.3g)", res$cor, res$cor_p),
    sprintf("- 回帰式: CT = %.3f + %.3f * PET (調整済R² = %.3f)", res$lm_intercept, res$lm_slope, res$adj_r2),
    sprintf("- 符号検定 (PET > CT) p値: %s", format_p(res$sign_p)),
    sprintf("- Wilcoxon 符号付順位検定 (PET > CT) p値: %s", format_p(res$wilcox_p)),
    "",
    "### 図",
    sprintf("![Scatter %s](../figures/%s_scatter.png)", res$dataset$name, prefix),
    "",
    sprintf("![Jitter %s](../figures/%s_jitter.png)", res$dataset$name, prefix),
    "",
    sprintf("![Histogram %s](../figures/%s_difference_hist.png)", res$dataset$name, prefix),
    "",
    sprintf("![Case Bar %s](../figures/%s_difference_bar.png)", res$dataset$name, prefix),
    ""
  )
}

writeLines(report_lines, report_path)

message("Summary report written to ", report_path)
