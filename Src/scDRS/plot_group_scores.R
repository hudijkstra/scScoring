#' This script reads a CSV file containing FDR-filtered single-cell association results,
#' formats cell type names, reshapes the data, and creates two plots:
#'   1. A bar plot of the number of significant cells per cell type
#'   2. A bar plot of the percentage of significant cells per cell type
#'
#'
#' Usage:
#' ------
#' Rscript plot_significant_cells.R <input_csv> <output_dir> <trait_name>
#'
#' Arguments:
#' ----------
#' <input_csv>   : Path to the input CSV file with columns such as:
#'                 - group
#'                 - n_cell
#'                 - n_fdr_0.05 / n_fdr_0.1 / n_fdr_0.2
#'                 - assoc_mcp
#' 
#' <output_dir>  : Path to the directory where output PNGs will be saved.
#' <trait_name>  : Name of the trait to include in plot titles.
#'
#' Author: Hessel Dijkstra


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

csv_file <- args[1]
output_dir <- args[2]
trait <- args[3]

# Define the colors for the FDR thresholds
colors <- c("FDR < 0.05" = "#f2d7d5",  # Light red
            "FDR < 0.1" = "#c0392b",   # Medium red
            "FDR < 0.2" = "#7b241c")   # Dark red

# Read in the CSV file
df <- read.csv(csv_file)

# Format cell type names
df <- df %>%
  mutate(
    group = recode(group,
      "CD4T"           = "CD4+ T",
      "CD8T"           = "CD8+ T",
      "monocyte"       = "Monocyte",
      "B"              = "B",
      "CD14 Mono"      = "CD14+ Mono",
      "CD16 Mono"      = "CD16+ Mono",
      "CD4 Naive"      = "CD4+ Naive T",
      "CD4 TCM"        = "CD4+ TCM",
      "CD8 Naive"      = "CD8+ Naive T",
      "CD8 TEM_1"      = "CD8+ TEM",
      "Intermediate B" = "Intermediate B",
      "Memory B"       = "Memory B",
      "NK"             = "NK",
      "DC"             = "DC",
      "cDC"            = "cDC",
      "gdT"            = "γδ T",
      "pDC"            = "pDC"
    )
  )


df_long <- df %>%
  pivot_longer(cols = starts_with("n_fdr_"),
               names_to = "FDR_threshold",
               values_to = "count")

df_long <- df_long %>%
  mutate(percentage = round((count / df$n_cell[match(df_long$group, df$group)]) * 100, 2))

df_long <- df_long %>%
  mutate(signif_label = case_when(
    assoc_mcp < 0.001 & FDR_threshold == "n_fdr_0.2" ~ "**",
    assoc_mcp < 0.05  & FDR_threshold == "n_fdr_0.2" ~ "*",
    TRUE ~ NA_character_
  ))


df_long$FDR_threshold <- factor(df_long$FDR_threshold,
                                levels = c("n_fdr_0.05", "n_fdr_0.1", "n_fdr_0.2"),
                                labels = c("FDR < 0.05", "FDR < 0.1", "FDR < 0.2"))

########################
# Count plot
########################
df_cnt_long <- df_long %>%
  group_by(group) %>%
  mutate(cum_height_count = cumsum(count)) %>%  
  ungroup()  

count_plot <- ggplot(df_cnt_long, aes(x = group, y = count, fill = FDR_threshold)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = signif_label, y = cum_height_count),
            vjust = -0.5,  
            size = 5, 
            na.rm = TRUE) +
  scale_fill_manual(values = colors) +  
  labs(
    x = "Cell Type",
    y = "Number of Significant Cells",
    fill = "FDR Threshold",
    title = paste("Significant Cells and Disease Association for", trait)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  theme_minimal()

########################
# Percentage plot
########################
df_pct_long <- df_long %>%
  group_by(group) %>%
  mutate(cum_heigh_pct = cumsum(percentage)) %>%  
  ungroup()  

percentage_plot <- ggplot(df_pct_long, aes(x = group, y = percentage, fill = FDR_threshold)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = signif_label, y = cum_heigh_pct),  
            vjust = -0.5,  
            size = 5, 
            na.rm = TRUE) +
  scale_fill_manual(values = colors) +
  labs(
    x = "Cell Type",
    y = "Percentage of Significant Cells (%)",
    fill = "FDR Threshold",
    title = paste("Significant Cells and Disease Association for", trait)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  theme_minimal()

# Adjust x-axis text angle if there are more than 6 groups
if (length(unique(df_long$group)) > 6) {
  percentage_plot <- percentage_plot +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  count_plot <- count_plot +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# Save the plots
basefile_no_ext <- tools::file_path_sans_ext(basename(csv_file))
count_plot_file <- paste0(output_dir, basefile_no_ext, "_count.png")
percentage_plot_file <- paste0(output_dir, basefile_no_ext, "_pct.png")

ggsave(count_plot_file, plot = count_plot, width = 8, height = 6, dpi = 300)
ggsave(percentage_plot_file, plot = percentage_plot, width = 8, height = 6, dpi = 300)

print(paste("Count plot saved to", count_plot_file))
print(paste("Percentage plot saved to", percentage_plot_file))