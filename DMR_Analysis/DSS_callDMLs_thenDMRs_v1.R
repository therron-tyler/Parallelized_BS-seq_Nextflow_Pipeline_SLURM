#!/usr/bin/env Rscript
libs <- .libPaths("/home/ttm3567/63_tylert/DEseq_libpath_dir/")

library(purrr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(readr)
library(BiocManager)
library(optparse)
library(DSS)

option_list <- list(
  make_option(c("-g","--groups"), type="character", help="CSV of Sample,Group"),
  make_option(     "--grp1"  , type="character", help="First group name (or comma-separated list)"),
  make_option(     "--grp2"  , type="character", help="Second group name (or comma-separated list)"),
  make_option(c("-b","--bsobj"), type="character", help="RDS file containing BSseq object"),
  make_option(c("-o","--out"),   type="character", default="DML", help="Output prefix")
)

opt <- parse_args(OptionParser(option_list=option_list))

## 1) read the mapping file
df.grp <- read.csv(opt$groups, stringsAsFactors=FALSE)

## 2) parse grp1/grp2 as one or more comparisons
grp1.vec <- strsplit(opt$grp1, ",")[[1]] |> trimws()
grp2.vec <- strsplit(opt$grp2, ",")[[1]] |> trimws()

if (length(grp1.vec) != length(grp2.vec)) {
  stop("Number of grp1 entries (", length(grp1.vec),
       ") does not match number of grp2 entries (", length(grp2.vec), ").")
}

n.comp <- length(grp1.vec)
message("Number of comparisons to run: ", n.comp)

## 3) load your BSseq object once
BSobj <- readRDS(opt$bsobj)
message("BSobj samples (first 10): ", paste(head(colnames(BSobj), 10), collapse=", "))

## 4) loop over each comparison
for (i in seq_len(n.comp)) {
  grp1.name <- grp1.vec[i]
  grp2.name <- grp2.vec[i]

  message("\n==========================")
  message("Running comparison ", i, " of ", n.comp, ": ",
          grp1.name, " vs ", grp2.name)
  message("==========================")

  # sample lists for this comparison
  grp1.samples <- df.grp$Sample[df.grp$Group == grp1.name]
  grp2.samples <- df.grp$Sample[df.grp$Group == grp2.name]

  if (length(grp1.samples) == 0 || length(grp2.samples) == 0) {
    stop("One of your groups has zero samples for comparison ",
         grp1.name, " vs ", grp2.name,
         ". Check --grp1/--grp2 against the CSV.")
  }

  message("Group 1 (", grp1.name, "): ", paste(grp1.samples, collapse=", "))
  message("Group 2 (", grp2.name, "): ", paste(grp2.samples, collapse=", "))

  # sanity checks
  stopifnot(all(grp1.samples %in% colnames(BSobj)))
  stopifnot(all(grp2.samples %in% colnames(BSobj)))
  stopifnot(length(unique(grp1.samples)) == length(grp1.samples))
  stopifnot(length(unique(grp2.samples)) == length(grp2.samples))

  # output prefix:
  #   - single comparison -> keep original behavior
  #   - multiple comparisons -> append comparison label
  if (n.comp == 1) {
    out_prefix_i <- opt$out
  } else {
    comp_label <- paste0(grp1.name, "_vs_", grp2.name)
    # sanitize a bit for filenames
    comp_label <- gsub("[^A-Za-z0-9._-]+", "_", comp_label)
    out_prefix_i <- paste0(opt$out, "_", comp_label)
  }

  ## 5) run the DML test
  dml.res <- DMLtest(
    BSobj,
    group1        = grp2.samples,
    group2        = grp1.samples,
    smoothing     = TRUE,
    smoothing.span= 500,
    equal.disp    = FALSE,
    ncores        = 12
  )

  ## 6) call DMLs
  dmls <- callDML(dml.res, delta = 0.05, p.threshold = 1e-3)

  ## 7) write out DML tables
  write.csv(
    dml.res,
    sprintf("%s_fullDMLtest.csv", out_prefix_i),
    row.names = FALSE
  )

  write.csv(
    dmls,
    sprintf("%s_significant_DMLs.csv", out_prefix_i),
    row.names = FALSE
  )

  message("Wrote:\n - ", out_prefix_i, "_fullDMLtest.csv\n - ",
          out_prefix_i, "_significant_DMLs.csv")

  ## 8) call DMRs
  dmrs <- callDMR(
    dml.res,
    delta      = 0.05,
    p.threshold= 1e-3,
    minlen     = 50,
    minCG      = 3,
    dis.merge  = 100,
    pct.sig    = 0.5
  )

  write_csv(dmrs, sprintf("%s_significant_DMRs.csv", out_prefix_i))
  message("Wrote:\n - ", out_prefix_i, "_significant_DMRs.csv")

  # optional peek for this comparison
  message("Top DMRs for ", grp1.name, " vs ", grp2.name, ":")
  print(head(dmrs))
  message("Number of DMRs: ", nrow(dmrs))
}
