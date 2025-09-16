libs <- .libPaths("/home/ttm3567/63_tylert/DEseq_libpath_dir/")

library(purrr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(readr)
library(BiocManager)

#install.packages("optparse")
library(optparse)
#BiocManager::install("DSS")
library(DSS)

option_list <- list(
  make_option(c("-g","--groups"), type="character", help="CSV of Sample,Group"),
  make_option(     "--grp1"  , type="character", help="First group name (exact match)"),
  make_option(     "--grp2"  , type="character", help="Second group name"),
  make_option(c("-b","--bsobj"), type="character", help="RDS file containing BSseq object"),
  make_option(c("-o","--out"),   type="character", default="DML", help="Output prefix")
)

opt <- parse_args(OptionParser(option_list=option_list))

# 1) read the mapping file
df.grp <- read.csv(opt$groups, stringsAsFactors=FALSE)

# 2) build two character vectors of sample names
grp1.samples <- df.grp$Sample[df.grp$Group == opt$grp1]
grp2.samples <- df.grp$Sample[df.grp$Group == opt$grp2]

if(length(grp1.samples)==0 || length(grp2.samples)==0){
  stop("One of your groups has zero samples. Check --grp1/--grp2 against the CSV.")
}

message("Group 1 (", opt$grp1, "): ", paste(grp1.samples, collapse=", "))
message("Group 2 (", opt$grp2, "): ", paste(grp2.samples, collapse=", "))

# 3) load your BSseq object
BSobj <- readRDS(opt$bsobj)

# 4) run the DML test
dml.res <- DMLtest(
  BSobj,
  group1        = grp1.samples,
  group2        = grp2.samples,
  smoothing     = TRUE,
  smoothing.span= 500,
  equal.disp    = FALSE,
  ncores        = 12
)

# 5) call DMLs at p < 1e-5
dmls <- callDML(dml.res, delta = 0.2, p.threshold=1e-3)

# 6) write out
write.csv(dml.res,
          sprintf("%s_fullDMLtest.csv", opt$out),
          row.names=FALSE)
write.csv(dmls,
          sprintf("%s_significant_DMLs.csv", opt$out),
          row.names=FALSE)

message("Wrote:\n - ", opt$out, "_fullDMLtest.csv\n - ", opt$out, "_significant_DMLs.csv")

# 5) take a look
head(dmls)
nrow(dmls)  # how many DMLs did you get?

# (Optionally) you can then call DMRs:
dmrs <- callDMR(dml.res,
                delta=0.2,
                p.threshold = 1e-3,
                minlen      = 100,
                minCG       = 5,
                dis.merge   = 300,
                pct.sig     = 0.5)
head(dmrs)

write_csv(dmrs, sprintf("%s_significant_DMRs.csv", opt$out))
message("Wrote:\n - ", opt$out, "_significant_DMRs.csv")
