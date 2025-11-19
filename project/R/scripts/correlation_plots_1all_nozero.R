#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# Determine project root based on this script's location
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

input_csv <- file.path(project_root, "data", "raw", "1all_data.csv")
figure_count_base <- file.path(project_root, "output", "figures", "1all_scatter_nozero")
figure_jitter_base <- file.path(project_root, "output", "figures", "1all_jitter_nozero")
log_path <- file.path(project_root, "output", "logs", "1all_correlation_nozero.txt")

if (!file.exists(input_csv)) {
  stop("Data file not found at ", input_csv)
}

# Read data (trim BOM if present) and exclude CT=0 / PET=0 cases
data <- read_csv(input_csv, show_col_types = FALSE) %>%
  mutate(
    CT = as.numeric(CT),
    PET = as.numeric(PET)
  ) %>%
  tidyr::drop_na(CT, PET)

original_n <- nrow(data)
data <- data %>%
  filter(CT != 0, PET != 0)
filtered_n <- nrow(data)

if (nrow(data) < 3) {
  stop("Need at least 3 observations for regression and correlation test after filtering.")
}

lm_fit <- lm(CT ~ PET, data = data)
cor_test <- cor.test(data$CT, data$PET, method = "pearson")

plot_count <- ggplot(data, aes(x = PET, y = CT)) +
  geom_count(color = "#2C3E50", alpha = 0.8) +
  scale_size_area(max_size = 8, guide = guide_legend(title = "Count")) +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C", fill = "#F9EBEA") +
  labs(
    title = "PET vs CT (CT/PET != 0)",
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
    title = "PET vs CT (CT/PET != 0, Jittered)",
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
  sprintf("Rows used: %d (from %d)", filtered_n, original_n),
  "Filter: excluded rows where CT = 0 or PET = 0.",
  sprintf("Pearson correlation: %.6f", cor_test$estimate),
  sprintf("p-value: %.6g", cor_test$p.value),
  sprintf("95%% CI: [%.6f, %.6f]", cor_test$conf.int[1], cor_test$conf.int[2]),
  sprintf("Linear model: CT = %.6f + %.6f * PET", coef(lm_fit)[1], coef(lm_fit)[2]),
  sprintf("Adjusted R-squared: %.6f", summary(lm_fit)$adj.r.squared)
)

dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
writeLines(log_lines, log_path)

message("Excluded cases with CT = 0 or PET = 0 (", filtered_n, " of ", original_n, " rows used).")
message("Analysis complete. Count figure base: ", figure_count_base)
message("Analysis complete. Jitter figure base: ", figure_jitter_base)
message("Correlation stats saved to ", log_path)
