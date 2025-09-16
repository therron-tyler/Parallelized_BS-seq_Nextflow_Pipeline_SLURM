libs <- .libPaths("/home/ttm3567/63_tylert/DEseq_libpath_dir/")

#BiocManager::install("genomation")
#BiocManager::install("GenomicRanges")
#BiocManager::install("rtracklayer")

library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(optparse)

option_list <- list(
		    make_option(c("-d","--dmr_csv"),  type="character", help="CSV of significant DMRs"),
		    make_option(c("-r","--ref_seq"),  type="character", help="RefSeq reference for annotating Genes to DMRs"),
		    make_option(c("-o","--out"),      type="character", help="Output prefix"))

# Parse the command-line arguments
opt <- parse_args(OptionParser(option_list=option_list))
dmr_file <- opt$dmr_csv
RefSeq <- opt$ref_seq
output_prefix <- opt$out

# 7.2) Read the DMR table

dmr_df   <- read.csv(dmr_file, stringsAsFactors=FALSE)

dmr_gr <- GRanges(
  seqnames = dmr_df$chr,
  ranges   = IRanges(start=dmr_df$start, end=dmr_df$end),
  strand   = "*",
  mcols    = dmr_df[, c("length","nCG","meanMethy1","meanMethy2","diff.Methy","areaStat")]
)

# 7.3) Read a genome‐wide RefSeq BED12 (for hg19 here, adjust path to your file):
#tx_bed     <- "mm39.ncbiRefSeq.bed"
tx_bed <- RefSeq

gene_parts <- readTranscriptFeatures(tx_bed,
                                     remove.unusual=TRUE,
                                     up.flank=2000, 
                                     down.flank=0,
                                     unique.prom=TRUE)

# 7.4) Annotate DMRs to promoter/exon/intron/intergenic
annot_by_parts <- annotateWithGeneParts(target=dmr_gr,
                                        feature=gene_parts,
                                        strand=FALSE,
                                        intersect.chr=FALSE)

# 7.5) Extract the “members” matrix and assign feature_class
mem_mat <- getMembers(annot_by_parts)
feature_col <- apply(mem_mat, 1, function(r) {
  idx <- which(r == 1)
  if (length(idx)==1) colnames(mem_mat)[idx] else NA_character_
})
dmr_df$feature_class <- feature_col

# 7.6) Pull out nearest‐TSS/gene and distance
#dist_df <- annot_by_parts@dist.to.TSS
#dmr_df$nearest_gene    <- dist_df$TSS_name
#dmr_df$distance_to_TSS <- dist_df$distance

# 2) use getAssociationWithTSS() to get a data.frame with one row per DMR
tss_info <- getAssociationWithTSS(annot_by_parts)

print(head(tss_info))
str(tss_info)

nearest_gene_vector    <- rep(NA_character_, nrow(dmr_df))
distance_to_TSS_vector <- rep(NA_real_, nrow(dmr_df))

# (d) Fill in only those indices that actually appear in tss_info$target.row:
nearest_gene_vector[tss_info$target.row]    <- tss_info$feature.name
distance_to_TSS_vector[tss_info$target.row] <- tss_info$dist.to.feature

# (e) Now assign back into dmr_df:
dmr_df$nearest_gene    <- nearest_gene_vector
dmr_df$distance_to_TSS <- distance_to_TSS_vector

# 7.7) (Optionally) annotate CpG islands if you have a BED
cgi_gr <- rtracklayer::import(tx_bed)
anno_cgi <- annotateWithFeature(target=dmr_gr, feature=cgi_gr,
                                 strand=FALSE, extend=0,
                                 feature.name="CpGisland", intersect.chr=FALSE)
dmr_df$in_CpG_island <- as.logical(getMembers(anno_cgi))    

# take the RefSeq names -> change them to EntrzIDs -> change to Gene Symbols

all_refseq_ids <- unique( na.omit(dmr_df$nearest_gene) )
length(all_refseq_ids)

# 2) Strip off the trailing “.number” so that “NM_001371004.1” → “NM_001371004”:
all_refseq_nover <- sub("\\.\\d+$", "", all_refseq_ids)

refseq2entrez <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys   = all_refseq_nover,
  keytype = "REFSEQ",
  columns = c("ENTREZID")
)

head(refseq2entrez)

unique_entrez <- unique(refseq2entrez$ENTREZID)
entrez2symbol <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = unique_entrez,
  keytype = "ENTREZID",
  columns = c("SYMBOL")
)

head(entrez2symbol)

refseq2symbol <- merge(
  refseq2entrez,
  entrez2symbol,
  by="ENTREZID",
  all.x = TRUE,     # keep all REFSEQ entries, even if no SYMBOL
  stringsAsFactors = FALSE
)
# Now refseq2symbol has columns: ENTREZID | REFSEQ | SYMBOL
head(refseq2symbol)

# (b) Turn that into a named vector: REFSEQ → SYMBOL
refseq_to_symbol_vec <- setNames(
  refseq2symbol$SYMBOL,
  nm = refseq2symbol$REFSEQ
)

gene_symbol_vector <- rep(NA_character_, nrow(dmr_df))

# Wherever you have a non‐NA nearest_gene (RefSeq), look up its symbol:
nearest_vec <- dmr_df$nearest_gene
# Because nearest_vec is still the versioned IDs (e.g. “NM_001371004.1”), strip off “.1”
nearest_nover <- sub("\\.\\d+$", "", nearest_vec)

which_have_refseq <- which(!is.na(nearest_nover))

gene_symbol_vector[ which_have_refseq ] <-
  refseq_to_symbol_vec[ nearest_nover[which_have_refseq] ]

# Now add it:
dmr_df$gene_symbol <- gene_symbol_vector

write.csv(dmr_df, file=paste0(output_prefix,".csv"), row.names=FALSE)
