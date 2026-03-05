#!/usr/bin/env Rscript

## Optional: set custom lib path first, then append existing ones
.libPaths(c("/home/ttm3567/63_tylert/DEseq_libpath_dir", .libPaths()))
#options(repos = c(CRAN = "https://cloud.r-project.org"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  ## DESeq2 not strictly required for this script, so you can drop it if you want
  ## library(DESeq2)
})

# ------------------------------------------------------------------------
# DMR_MAplot.R
#
# Usage:
#   Rscript DMR_MAplot.R <DMR_CSV> <DEG_CSV> <OUT_DIR>
#
#  - <DMR_CSV>   : path to DMR table (must contain meanMethy1, meanMethy2,
#                  diff.Methy, gene_symbol, etc.)
#  - <DEG_CSV>   : path to DEG results (must have Gene_Symbol, padj,
#                  and a 'log2FoldChange*' column)
#  - <OUT_DIR>   : directory to write the resulting figure + data
# ------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Three arguments required:\n",
    "  1) DMR CSV\n",
    "  2) DEG CSV\n",
    "  3) Output directory\n",
    "  4) Sample Name\n",
    call. = FALSE
  )
}

dmr_file <- args[1]
deg_file <- args[2]
out_dir  <- args[3]
sample_name <- args[4]

cat("DMR file :", dmr_file, "\n")
cat("DEG file :", deg_file, "\n")
cat("Out dir  :", out_dir, "\n")
cat("Sample Comp  :", sample_name, "\n")

# Create output directory if needed
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 1) Read DMR table and compute meanAll
dmr_df <- readr::read_csv(dmr_file, show_col_types = FALSE) %>%
  mutate(
    meanAll = (meanMethy1 + meanMethy2) / 2,
    gene_symbol = trimws(gene_symbol)
  )

# 2) Read DEG results
deg_df <- readr::read_csv(deg_file, show_col_types = FALSE) %>%
  mutate(Gene_Symbol = trimws(Gene_Symbol))

## Find the log2FC column, e.g. "log2FoldChange(MFplus_PMNs_stim_with_TNF/MF)"
lfc_col <- grep("^log2FoldChange", colnames(deg_df), value = TRUE)[1]

if (is.na(lfc_col)) {
  stop("Could not find a log2FoldChange* column in the DEG file.")
}

## Classify Up/Down/Other within the DEG table
padj_thr <- 0.10
lfc_thr  <- 0.10

deg_df <- deg_df %>%
  mutate(
    log2FC = as.numeric(.data[[lfc_col]]),
    padj   = as.numeric(padj),
    signif = !is.na(padj) & !is.na(log2FC) & padj <= padj_thr,
    Regulation = dplyr::case_when(
      signif & log2FC >=  lfc_thr ~ "Upregulated",
      signif & log2FC <= -lfc_thr ~ "Downregulated",
      TRUE                       ~ "Other"
    ),
    Regulation = factor(
      Regulation,
      levels = c("Downregulated", "Other", "Upregulated")
    )
  )

## Sanity check
cat("DEG Regulation counts:\n")
print(table(deg_df$Regulation, useNA = "ifany"))

## Overlap sanity check
overlap_genes <- intersect(unique(dmr_df$gene_symbol),
                           unique(deg_df$Gene_Symbol))
cat("Number of overlapping genes between DMR and DEG tables:",
    length(overlap_genes), "\n")
cat("Example overlap genes:\n")
print(head(overlap_genes, 20))

# 3) Join DEG info onto DMRs by gene_symbol
deg_for_join <- deg_df %>%
  select(Gene_Symbol, log2FC, padj, Regulation) %>%
  rename(deg_symbol = Gene_Symbol)

dmr_df <- dmr_df %>%
  left_join(deg_for_join, by = c("gene_symbol" = "deg_symbol")) %>%
  mutate(
    Regulation = dplyr::case_when(
      Regulation %in% c("Upregulated", "Downregulated") ~ as.character(Regulation),
      TRUE                                              ~ "Other"
    ),
    Regulation = factor(
      Regulation,
      levels = c("Downregulated", "Other", "Upregulated")
    )
  )

## 4) Sanity check on merged object
cat("DMR+DEG merged Regulation counts:\n")
print(table(dmr_df$Regulation, useNA = "ifany"))

dmr_plot_df <- dmr_df

# 5) Choose top 10 Up + top 10 Down genes to label
label_candidates <- dmr_plot_df %>%
  filter(Regulation %in% c("Upregulated", "Downregulated"),
         !is.na(gene_symbol),
         gene_symbol != "") %>%
  group_by(gene_symbol, Regulation) %>%
  slice_max(order_by = abs(diff.Methy), n = 1, with_ties = FALSE) %>%
  ungroup()

top_up <- label_candidates %>%
  filter(Regulation == "Upregulated") %>%
  slice_max(order_by = diff.Methy, n = 10, with_ties = FALSE)

top_down <- label_candidates %>%
  filter(Regulation == "Downregulated") %>%
  slice_min(order_by = diff.Methy, n = 10, with_ties = FALSE)

labels_df <- bind_rows(top_up, top_down)

# 6) Build ggplot MA-plot
p <- ggplot(
  dmr_plot_df,
  aes(x = meanAll, y = diff.Methy, color = Regulation)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_color_manual(
    name   = "Regulation",
    breaks = c("Downregulated", "Other", "Upregulated"),
    labels = c("Downregulated", "Other", "Upregulated"),
    values = c(
      "Downregulated" = "blue",
      "Other"         = "grey80",
      "Upregulated"   = "red"
    )
  ) +
  labs(
    title = "DMR MA-plot",
    x = "Mean methylation ((group1 + group2) / 2)",
    y = expression(Delta * " methylation (group2 - group1)")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(colour = "black")
  )

# Add labeled genes (if any)
if (nrow(labels_df) > 0) {
  p <- p +
    geom_point(
      data = labels_df,
      aes(x = meanAll, y = diff.Methy),
      size = 1.5,
      shape = 21,
      stroke = 0.3,
      color = "black",
      fill  = NA
    ) +
    ggrepel::geom_text_repel(
      data = labels_df,
      aes(x = meanAll, y = diff.Methy, label = gene_symbol),
      size          = 3,
      box.padding   = 0.3,
      point.padding = 0.2,
      max.overlaps  = Inf,
      segment.size  = 0.25
    )
}

# 7) Save merged DMR+DEG table used for plotting
dmr_data_out <- file.path(out_dir, paste0(sample_name,"DMR_MA_plot_data.csv"))
readr::write_csv(dmr_plot_df, dmr_data_out)

# 8) Save figure (PNG + PDF)
out_png <- file.path(out_dir, paste0(sample_name,"DMR_MA_plot_ggplot.png"))
out_pdf <- file.path(out_dir, paste0(sample_name,"DMR_MA_plot_ggplot.pdf"))

ggsave(out_png, plot = p, width = 6, height = 4.5, dpi = 300)
ggsave(out_pdf, plot = p, width = 6, height = 4.5)

cat("✅ Done! MA-plot saved to:\n  ",
    out_png, "\n  ",
    out_pdf, "\n",
    "Data table saved to:\n  ",
    dmr_data_out, "\n",
    sep = "")
