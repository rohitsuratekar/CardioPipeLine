# DESeq2 analysis
#
# Usage : Rscript deseq2.R <arguments>
#
# Arguments:
# -c, --counts : [String] Full path of count data file. (mandatory)
# -m, --meta : [String] Full path of the metadata file. (mandatory)
# -d, --design : [String] Design Formula (mandatory)
# -t, --test : [String] Stat Test used in the analysis (optional, default=Wald)
# -r, --reduced : [String] Rdecued design formula for LRT test (mandatory
# when test is LRT)
# -n, --contrast : [String] Contrast for the LFC Shrinkage. This should be
# exactly 3 strings seperated by '\t'. If provided LFC statistics is
# calculated else it will not. (optional)
# -a, --alpha : [Float] P-value threshold for the FDR calculation.
# (optional, default=0.1)
# -o, --out : [String] Full path of the output file. (optional,
# default=deseq2)
#
# Example:
# Rscript deseq2.R -c /path/to/count -m /path/to/meta -d "~ condition"
#

library("optparse") # Generate arguments
library("data.table") # To properly format output file
suppressPackageStartupMessages(library("DESeq2")) # For DESeq2 analysis

opt_counts <- optparse::make_option(c("-c", "--counts"), help =
  "Full path of count data file")

opt_meta <- optparse::make_option(c("-m", "--meta"), help = "Full path of
the metadata file")

opt_design <- optparse::make_option(c("-d", "--design"), help = "Design
formula ")

opt_test <- optparse::make_option(c("-t", "--test"), help = "Stat Test used in
the analysis")

opt_reduced <- optparse::make_option(c("-r", "--reduced"), help = "Rdecued
design formula for LRT test")

opt_out <- optparse::make_option(c("-o", "--out"), help = "Full path of the
output file.", default = "deseq2")

opt_contrast <- optparse::make_option(c("-n", "--contrast"), help = "Contrast
 for the LFC Shrinkage")

opt_alpha <- optparse::make_option(c("-a", "--alpha"), help = "P-value
threshold for the FDR calculation.")

opt_parser <- optparse::OptionParser(option_list = list(opt_counts, opt_out,
                                                        opt_meta, opt_reduced,
                                                        opt_design, opt_test,
                                                        opt_contrast, opt_alpha))

# Parse the options
opt <- optparse::parse_args(opt_parser)

# Convert the files into dataframes or tables for the further analysis
counts <- read.csv(opt$counts, header = TRUE)
rownames(counts) <- counts$gene_id
counts$gene_id <- NULL

meta <- read.csv(opt$meta, header = TRUE)
rownames(meta) <- meta$samples
meta$samples <- NULL

# Sanity check
rn <- rownames(meta)
cn <- colnames(counts)
if (!all(rn == cn) || length(rn) != length(cn)) {
  stop("Row names in metadata file and column names in count files should be
 in same order.")
}

# Create the DESeq data-set
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                      colData = meta,
                                      design = formula(opt$design))

# Perform analysis based on type of test
if (opt$test == "LRT") {
  dds <- DESeq2::DESeq(dds,
                       test = opt$test,
                       reduced = formula(opt$reduced))
}else {
  dds <- DESeq2::DESeq(dds)
}


# Normalized counts
normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
# Save Normalized counts
normalized_counts <- as.data.frame(normalized_counts)
data.table::setDT(normalized_counts, keep.rownames = "gene_id")
write.csv(normalized_counts,
          file = paste0(opt$out, ".normalized.csv"),
          row.names = FALSE)


# Perform LFC Shrinking if contrast is provided
if (!is.null(opt$contrast)) {

  contrast <- unlist(strsplit(opt$contrast, "\t"))

  # Get results
  res <- DESeq2::results(dds, contrast = contrast, alpha = opt$alpha)


  lfc_res <- DESeq2::lfcShrink(dds, contrast = contrast, res = res)
  lfc_df <- as.data.frame(lfc_res)
  data.table::setDT(lfc_df, keep.rownames = "gene_id")
  write.csv(lfc_df,
            file = paste0(opt$out, ".result_lfc.csv"),
            row.names = FALSE)

} else {
  # Get results without contrast and alpha
  res <- DESeq2::results(dds)

}

# Save results
res_df <- as.data.frame(res)
data.table::setDT(res_df, keep.rownames = "gene_id")
write.csv(res_df,
          file = paste0(opt$out, ".result.csv"),
          row.names = FALSE)

