#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tibble)
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

input_csv <- file.path(project_root, "data", "raw", "2before_after.csv")
log_path <- file.path(project_root, "output", "logs", "2before_after_nonparametric.txt")
table_path <- file.path(project_root, "output", "tables", "2before_after_nonparametric.csv")
figure_hist_base <- file.path(project_root, "output", "figures", "2before_after_difference_hist")
figure_bar_base <- file.path(project_root, "output", "figures", "2before_after_difference_bar")

if (!file.exists(input_csv)) {
  stop("Data file not found at ", input_csv)
}

data <- read_csv(input_csv, show_col_types = FALSE, skip = 2) %>%
  transmute(
    CT = as.numeric(.data$CT),
    PET = as.numeric(.data$PET),
    diff = PET - CT
  ) %>%
  tidyr::drop_na(CT, PET)

if (nrow(data) < 3) {
  stop("Need at least 3 paired observations for non-parametric tests.")
}

diff_nonzero <- data$diff[data$diff != 0]

if (length(diff_nonzero) == 0) {
  sign_test <- NULL
} else {
  sign_test <- binom.test(sum(diff_nonzero > 0), length(diff_nonzero), alternative = "greater")
}

if (length(unique(data$diff)) <= 1) {
  wilcox_res <- NULL
} else {
  wilcox_res <- suppressWarnings(
    wilcox.test(data$PET, data$CT, alternative = "greater", paired = TRUE, exact = FALSE)
  )
}

log_lines <- c(
  sprintf("Date: %s", Sys.time()),
  sprintf("Input: %s", input_csv),
  sprintf("Pairs: %d", nrow(data)),
  sprintf("Non-zero differences: %d", length(diff_nonzero))
)

if (!is.null(sign_test)) {
  log_lines <- c(
    log_lines,
    "Sign test (PET > CT):",
    sprintf("  positives: %d / %d", sum(diff_nonzero > 0), length(diff_nonzero)),
    sprintf("  p-value: %.6g", sign_test$p.value)
  )
} else {
  log_lines <- c(log_lines, "Sign test: skipped (no non-zero differences)")
}

if (!is.null(wilcox_res)) {
  log_lines <- c(
    log_lines,
    "Wilcoxon signed-rank test (PET > CT):",
    sprintf("  statistic V: %.6f", wilcox_res$statistic),
    sprintf("  p-value: %.6g", wilcox_res$p.value)
  )
} else {
  log_lines <- c(log_lines, "Wilcoxon signed-rank test: skipped (insufficient variation)")
}

diff_plot <- ggplot(data, aes(x = diff)) +
  geom_histogram(fill = "#2E86C1", color = "white", bins = 15, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#C0392B") +
  labs(
    title = "Distribution of PET - CT (2before_after)",
    x = "PET - CT",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)

bar_plot_data <- data %>%
  mutate(case = row_number())

bar_plot <- ggplot(bar_plot_data, aes(x = case, y = diff)) +
  geom_col(fill = "#E67E22", width = 0.8) +
  geom_hline(yintercept = 0, color = "#566573", linetype = "dashed", linewidth = 0.5) +
  labs(
    title = "Difference (PET - CT) by Case (2before_after)",
    x = "Case",
    y = "PET - CT"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank())

summary_tbl <- tibble(
  test = character(),
  statistic = character(),
  parameter = character(),
  p_value = numeric(),
  note = character()
)

if (!is.null(sign_test)) {
  summary_tbl <- add_row(
    summary_tbl,
    test = "Sign (PET > CT)",
    statistic = sprintf("k = %s", signif(sign_test$statistic, 6)),
    parameter = sprintf("n = %s", signif(sign_test$parameter, 6)),
    p_value = sign_test$p.value,
    note = ""
  )
} else {
  summary_tbl <- add_row(
    summary_tbl,
    test = "Sign (PET > CT)",
    statistic = NA_character_,
    parameter = NA_character_,
    p_value = NA_real_,
    note = "Skipped (no non-zero differences)"
  )
}

if (!is.null(wilcox_res)) {
  summary_tbl <- add_row(
    summary_tbl,
    test = "Wilcoxon signed-rank (PET > CT)",
    statistic = sprintf("V = %s", signif(unname(wilcox_res$statistic), 6)),
    parameter = sprintf("pairs = %d", nrow(data)),
    p_value = wilcox_res$p.value,
    note = ""
  )
} else {
  summary_tbl <- add_row(
    summary_tbl,
    test = "Wilcoxon signed-rank (PET > CT)",
    statistic = NA_character_,
    parameter = sprintf("pairs = %d", nrow(data)),
    p_value = NA_real_,
    note = "Skipped (insufficient variation)"
  )
}

dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(figure_hist_base), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(figure_bar_base), recursive = TRUE, showWarnings = FALSE)

writeLines(log_lines, log_path)
write_csv(summary_tbl, table_path)
formats <- c("png", "pdf", "svg")
for (ext in formats) {
  ggsave(sprintf("%s.%s", figure_hist_base, ext), diff_plot, width = 6, height = 4, dpi = 300)
  ggsave(sprintf("%s.%s", figure_bar_base, ext), bar_plot, width = 5.5, height = 4, dpi = 300)
}

message("2before_after non-parametric tests complete. Log: ", log_path)
message("Summary table saved to ", table_path)
message("Difference histogram base: ", figure_hist_base)
message("Difference bar plot base: ", figure_bar_base)
