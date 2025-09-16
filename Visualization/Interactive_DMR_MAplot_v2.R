#!/usr/bin/env Rscript

libs <- .libPaths("/home/ttm3567/63_tylert/DEseq_libpath_dir/")
options(repos = c(CRAN = "https://cloud.r-project.org"))


# ------------------------------------------------------------------------
# Interactive_DMR_MAplot.R
#
# Usage:
#   Rscript Interactive_DMR_MAplot.R <DMR_CSV> <DEG_CSV> <OUT_DIR>
#
#  - <DMR_CSV>   : path to DMR table (must contain meanMethy1, meanMethy2, diff.Methy, gene_symbol, etc.)
#  - <DEG_CSV>   : path to DEG results (must have Gene_Symbol and Regulation columns)
#  - <OUT_DIR>   : directory to write the resulting HTML file
# ------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Three arguments required:\n  1) DMR CSV\n  2) DEG CSV\n  3) Output directory\n",
       call. = FALSE)
}
dmr_file <- args[1]
deg_file <- args[2]
out_dir  <- args[3]

# Create output directory if needed
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Load libraries
suppressPackageStartupMessages({
  library(htmlwidgets)
  library(plotly)
  library(dplyr)
})

# 1) Read DMR table and compute meanAll
dmr_df <- read.csv(dmr_file, stringsAsFactors = FALSE) %>%
  mutate(
    meanAll = (meanMethy1 + meanMethy2) / 2
  )

# 2) Read DEG results
deg_df <- read.csv(deg_file, stringsAsFactors = FALSE) %>%
  select(Gene_Symbol, Regulation)

# 3) Join Regulation onto DMRs by gene_symbol
dmr_df <- dmr_df %>%
  left_join(deg_df, by = c("gene_symbol" = "Gene_Symbol")) %>%
  mutate(
    Regulation = ifelse(is.na(Regulation), "Other", Regulation),
    Regulation = factor(Regulation, levels = c("Downregulated", "Other", "Upregulated"))
  )

# 4) Build hover text
dmr_df <- dmr_df %>%
  mutate(
    hover_text = paste0(
      "chr: ", chr, "<br>",
      "start: ", start, "<br>",
      "end: ", end, "<br>",
      "length: ", length, "<br>",
      "nCG: ", nCG, "<br>",
      "meanMethy1: ", sprintf("%.4f", meanMethy1), "<br>",
      "meanMethy2: ", sprintf("%.4f", meanMethy2), "<br>",
      "diff.Methy: ", sprintf("%.4f", diff.Methy), "<br>",
      "areaStat: ", sprintf("%.2f", areaStat), "<br>",
      "feature_class: ", feature_class, "<br>",
      "nearest_gene: ", nearest_gene, "<br>",
      "distance_to_TSS: ", distance_to_TSS, "<br>",
      "in_CpG_island: ", in_CpG_island, "<br>",
      "Gene Symbol: ", gene_symbol, "<br>",
      "mean1: ", sprintf("%.3f", meanMethy1), "<br>",
      "mean2: ", sprintf("%.3f", meanMethy2), "<br>",
      "Δmeth: ", sprintf("%.3f", diff.Methy), "<br>",
      "Regulation: ", Regulation
    )
  )

# 5) Subsample if too many points
#set.seed(123)
#max_pts <- 10000
#dmr_sub <- if (nrow(dmr_df) > max_pts) {
#  sample_n(dmr_df, max_pts)
#} else {
#  dmr_df
#}

# 6) Create interactive MA‐plot
fig <- plot_ly(
  data      = dmr_df,
  x         = ~meanAll,
  y         = ~diff.Methy,
  type      = "scatter",
  mode      = "markers",
  hoverinfo = "text",
  text      = ~hover_text,
  color     = ~Regulation,
  colors    = c("blue", "lightgrey", "red"),
  marker    = list(size = 6, opacity = 0.6)
) %>%
  layout(
    title  = "Interactive MA‐Plot of DMRs",
    xaxis  = list(title = "Mean methylation ((group1+group2)/2)"),
    yaxis  = list(title = "Δ methylation (group2 – group1)"),
    legend = list(title = list(text = "<b>Regulation</b>"))
  )

# 7) Save HTML
out_file <- file.path(out_dir, "DMR_MA_plot.html")
saveWidget(fig, file = out_file, selfcontained = TRUE)

cat("✅ Done! Interactive MA‐plot saved to:\n  ", out_file, "\n", sep = "")
