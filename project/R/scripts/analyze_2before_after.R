#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
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
figure_count_base <- file.path(project_root, "output", "figures", "2before_after_scatter")
figure_jitter_base <- file.path(project_root, "output", "figures", "2before_after_jitter")
log_path <- file.path(project_root, "output", "logs", "2before_after_correlation.txt")

if (!file.exists(input_csv)) {
  stop("Data file not found at ", input_csv)
}

data <- read_csv(input_csv, show_col_types = FALSE, skip = 2) %>%
  transmute(
    CT = as.numeric(.data$CT),
    PET = as.numeric(.data$PET)
  ) %>%
  tidyr::drop_na(CT, PET)

if (nrow(data) < 3) {
  stop("Need at least 3 observations for regression and correlation test.")
}

lm_fit <- lm(CT ~ PET, data = data)
cor_test <- cor.test(data$CT, data$PET, method = "pearson")

plot_count <- ggplot(data, aes(x = PET, y = CT)) +
  geom_count(color = "#2C3E50", alpha = 0.8) +
  scale_size_area(max_size = 8, guide = guide_legend(title = "Count")) +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C", fill = "#F9EBEA") +
  labs(
    title = "PET vs CT (2before_after)",
    subtitle = sprintf("Pearson r = %.3f (p = %.3g)", cor_test$estimate, cor_test$p.value),
    x = "PET",
    y = "CT"
  ) +
  theme_minimal(base_size = 14)

plot_jitter <- ggplot(data, aes(x = PET, y = CT)) +
  geom_point(
    size = 2,
    alpha = 0.8,
    color = "#2C3E50",
    position = position_jitter(width = 0.15, height = 0.15)
  ) +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C", fill = "#F9EBEA") +
  labs(
    title = "PET vs CT (2before_after, Jittered)",
    subtitle = sprintf("Pearson r = %.3f (p = %.3g)", cor_test$estimate, cor_test$p.value),
    x = "PET",
    y = "CT"
  ) +
  theme_minimal(base_size = 14)

dir.create(dirname(figure_count_base), recursive = TRUE, showWarnings = FALSE)
formats <- c("png", "pdf", "svg")
for (ext in formats) {
  ggsave(sprintf("%s.%s", figure_count_base, ext), plot_count, width = 6, height = 4, dpi = 300)
  ggsave(sprintf("%s.%s", figure_jitter_base, ext), plot_jitter, width = 6, height = 4, dpi = 300)
}

log_lines <- c(
  sprintf("Date: %s", Sys.time()),
  sprintf("Input: %s", input_csv),
  sprintf("Pearson correlation: %.6f", cor_test$estimate),
  sprintf("p-value: %.6g", cor_test$p.value),
  sprintf("95%% CI: [%.6f, %.6f]", cor_test$conf.int[1], cor_test$conf.int[2]),
  sprintf("Linear model: CT = %.6f + %.6f * PET", coef(lm_fit)[1], coef(lm_fit)[2]),
  sprintf("Adjusted R-squared: %.6f", summary(lm_fit)$adj.r.squared)
)

dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
writeLines(log_lines, log_path)

message("2before_after analysis complete. Count figure base: ", figure_count_base)
message("2before_after analysis complete. Jitter figure base: ", figure_jitter_base)
message("Correlation stats saved to ", log_path)
