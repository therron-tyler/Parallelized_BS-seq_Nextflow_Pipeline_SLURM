libs <- .libPaths("/home/ttm3567/63_tylert/DEseq_libpath_dir/")

# Load libraries
suppressMessages(library(optparse))
suppressMessages(library(DSS))
suppressMessages(library(bsseq))

# Define CLI options
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = NULL,
              help = "Directory with .cov_for_DSS.tsv files", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output RDS file path for BSseq object", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$input_dir) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("❌ Both --input_dir and --output must be specified.", call. = FALSE)
}

# Find input files
files <- list.files(opt$input_dir, pattern = "\\.cov_for_DSS\\.tsv$", full.names = TRUE)
if (length(files) == 0) {
  stop("❌ No .cov_for_DSS.tsv files found in specified input directory.")
}

cat("📂 Found", length(files), "files\n")

# Read all files
dat <- lapply(files, read.table, header = TRUE)
names(dat) <- sub("\\.cov_for_DSS\\.tsv$", "", basename(files))

# Create BSseq object
cat("🧬 Constructing BSseq object...\n")
BSobj <- makeBSseqData(dat, sampleNames = names(dat))

# Save to RDS
cat("💾 Saving BSseq object to:", opt$output, "\n")
saveRDS(BSobj, opt$output)
cat("✅ Done.\n")
