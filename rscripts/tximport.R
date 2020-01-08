# Generates count matrix using tximport
#
# Usage : Rscript tximport.R <arguments>
#
# Arguments:
# -m, --method : [String] Method used for generating transcript alignment.
# Currently supports salmon, kallisto and stringtie (mandatory)
# -f, --file : [String] Full path of input file. (mandatory)
# -g, --gft : [String] Full path of the GTF file. (mandatory)
# -s, --srr : [String] SRR ID which will be used to make column name in the
# count matrix. If not provided, input filename will be used. (optional)
# -o, --out : [String] Full path of the output file. (optional,
# default=countmatrix.csv)
#
# Example:
# Rscript tximport.R -m salmon -f /path/to/quant.sf -g /path/to/gtf
#

library("optparse") # Generate arguments
library("tximport") # For generating count matrix
library("data.table") # To properly format output file
library("readr") # Reading csv properly
# Following packages needed to generate tx2gene file
suppressPackageStartupMessages(library('AnnotationDbi'))
suppressPackageStartupMessages(library('GenomicFeatures'))

# Make arguments
opt_method <- optparse::make_option(c("-m", "--method"), help = "Method used
 for generating transcript alignment. Currently supports Salmon and Kallisto")

opt_file <- optparse::make_option(c("-f", "--file"), help =
  "Full path of Quantification file from upstram analyis (e.g. quant.sf for
  salmon)")

opt_gtf <- optparse::make_option(c("-g", "--gtf", help = "Full
path of the GTF file."))

opt_srr <- optparse::make_option(c("-s", "--srr"), help = "SRR ID which will be
 used to make column name in the count matrix. If not provided, input
 filename will be used", default = NULL)

opt_out <- optparse::make_option(c("-o", "--out"), help = "Full path of the
output file.", default = "countmatrix.csv")

opt_parser <- optparse::OptionParser(option_list = list(opt_method,
                                                        opt_file, opt_srr,
                                                        opt_gtf, opt_out))

# Parse the options
opt <- optparse::parse_args(opt_parser)

# Generate tx2gene object
if (opt$method == "stringtie") {
  tmp <- readr::read_tsv(opt$file)
  tx2gene <- tmp[, c("t_name", "gene_id")]
} else {
  tx <- GenomicFeatures::makeTxDbFromGFF(file = opt$gtf, format = "gtf",
                                         organism = "Danio rerio")

  k <- AnnotationDbi::keys(tx, keytype = "GENEID")
  tx_df <- AnnotationDbi::select(tx, keys = k, keytype = "GENEID", columns =
    "TXNAME")
  tx2gene <- tx_df[, 2:1]
}

# # Generate transcript level summerization with tximport
txi <- tximport::tximport(opt$file,
                          type = opt$method,
                          tx2gene = tx2gene,
                          ignoreTxVersion = TRUE)

df <- as.data.frame(txi$counts)

# If SRR ID is provided, change the column name.
if (!is.null(opt$srr)) {
  colnames(df) <- opt$srr
} else {
  colnames(df) <- "sample"
}

# Set index column name before saving to the file
data.table::setDT(df, keep.rownames = "gene_id")

# Save to the file
write.csv(df, file = opt$out, row.names = FALSE)