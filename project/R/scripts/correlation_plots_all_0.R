#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(tidyr)
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

input_csv <- file.path(project_root, "data", "raw", "1all_data_no00.csv")
figure_count_base <- file.path(project_root, "output", "figures", "1all_data_no00_scatter")
log_path <- file.path(project_root, "output", "logs", "1all_data_no00_correlation.txt")

if (!file.exists(input_csv)) {
  stop("Data file not found at ", input_csv)
}

data <- read_csv(input_csv, show_col_types = FALSE) %>%
  transmute(
    CT = as.numeric(.data$CT),
    PET = as.numeric(.data$PET)
  ) %>%
  tidyr::drop_na(CT, PET)

bowker_test <- function(tab) {
  if (nrow(tab) != ncol(tab)) {
    stop("Bowker test requires a square contingency table.")
  }
  stat <- 0
  df <- 0
  k <- nrow(tab)
  for (i in seq_len(k - 1)) {
    for (j in seq.int(i + 1, k)) {
      nij <- tab[i, j]
      nji <- tab[j, i]
      if (nij + nji > 0) {
        stat <- stat + (nij - nji)^2 / (nij + nji)
        df <- df + 1
      }
    }
  }
  p <- pchisq(stat, df, lower.tail = FALSE)
  list(statistic = stat, parameter = df, p.value = p)
}

if (nrow(data) < 3) {
  stop("Need at least 3 observations for regression and correlation test.")
}

lm_fit <- lm(CT ~ PET, data = data)
cor_test <- cor.test(data$CT, data$PET, method = "pearson")
p_label <- sprintf("%.3g", cor_test$p.value)
tab_ct_pet <- xtabs(~ CT + PET, data = data)
bowker_res <- tryCatch(bowker_test(tab_ct_pet), error = function(e) NULL)
bowker_subtitle <- if (!is.null(bowker_res)) {
  sprintf("Bowker: X2 = %.2f, df = %d, p = %.3g", bowker_res$statistic, bowker_res$parameter, bowker_res$p.value)
} else {
  "Bowker: not computed"
}

plot_count <- ggplot(data, aes(x = PET, y = CT)) +
  geom_count(color = "#2C3E50", alpha = 0.8) +
  scale_size_area(max_size = 12, guide = guide_legend(title = "Count")) +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C", fill = "#F9EBEA") +
  scale_x_continuous(labels = label_number(accuracy = 1), expand = expansion(mult = c(0.02, 0.1))) +
  scale_y_continuous(labels = label_number(accuracy = 1)) +
  labs(
    title = "PET vs CT (1all_data_no00)",
    subtitle = sprintf("Pearson r = %.3f (P = %s)", cor_test$estimate, p_label),
    x = "PET",
    y = "CT"
  ) +
  theme_minimal(base_size = 14)

dir.create(dirname(figure_count_base), recursive = TRUE, showWarnings = FALSE)
formats <- c("png", "pdf", "svg")
for (ext in formats) {
  ggsave(sprintf("%s.%s", figure_count_base, ext), plot_count, width = 6, height = 4, dpi = 300)
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
if (!is.null(bowker_res)) {
  log_lines <- c(
    log_lines,
    sprintf("Bowker test statistic: %.6f", unname(bowker_res$statistic)),
    sprintf("Bowker df: %s", bowker_res$parameter),
    sprintf("Bowker p-value: %.6g", bowker_res$p.value)
  )
} else {
  log_lines <- c(log_lines, "Bowker test: not computed (error)")
}

dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
writeLines(log_lines, log_path)

# Cross-tab (CT rows, PET columns) and heatmap
ct_levels <- sort(unique(data$CT))
pet_levels <- sort(unique(data$PET))
crosstab_long <- count(data, CT, PET, name = "count")
crosstab_wide <- crosstab_long %>%
  mutate(
    CT = factor(CT, levels = ct_levels),
    PET = factor(PET, levels = pet_levels)
  ) %>%
  pivot_wider(names_from = PET, values_from = count, values_fill = 0) %>%
  arrange(CT) %>%
  select(CT, all_of(as.character(pet_levels))) %>%
  mutate(CT = as.character(CT))

crosstab_path <- file.path(project_root, "output", "tables", "1all_data_no00_crosstab.csv")
dir.create(dirname(crosstab_path), recursive = TRUE, showWarnings = FALSE)
write_csv(crosstab_wide, crosstab_path)

heatmap_plot <- ggplot(
  crosstab_long %>%
    tidyr::complete(
      CT = ct_levels,
      PET = pet_levels,
      fill = list(count = 0)
    ) %>%
    mutate(
      CT = factor(CT, levels = ct_levels),
      PET = factor(PET, levels = pet_levels)
    ),
  aes(x = PET, y = CT, fill = count)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = count), color = "black", size = 4) +
  scale_fill_gradient(low = "#F5EEF8", high = "#5DADE2", name = "Count") +
  labs(
    title = "Heatmap: CT vs PET (1all_data_no00)",
    subtitle = bowker_subtitle,
    x = "PET score",
    y = "CT score"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())

heatmap_base <- file.path(project_root, "output", "figures", "1all_data_no00_crosstab_heatmap")
for (ext in formats) {
  ggsave(sprintf("%s.%s", heatmap_base, ext), heatmap_plot, width = 6, height = 5, dpi = 300)
}

message("all_0 analysis complete. Scatter figure base: ", figure_count_base)
message("Cross-tab saved to ", crosstab_path, "; heatmap saved to ", heatmap_base)
message("Correlation stats saved to ", log_path)
