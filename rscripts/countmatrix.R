# Generates count matrix from sorted BAM files and given GTF file
#
# Usage : Rscript countmatrix.R <arguments>
#
# Arguments:
# -f, --file : [String] Full path of input BAM file. (mandatory)
# -g, --gft : [String] Full path of the GTF file. (mandatory)
# -s, --srr : [String] SRR ID which will be used to make column name in the
# count matrix. If not provided, input filename will be used. (optional)
# -o, --out : [String] Full path of the output file. (optional,
# default=countmatrix.csv)
# -p, --paired : [Integer] Is current input file is paired-end? Use 1 for yes
# and 0 for no (optional, default=1)
# -t, --threads : [Integer] Number of threads (optional, default=1)
#
# Example:
# Rscript countmatrix.R -f /path/to/bam -g /path/to/gtf
#

library("optparse") # Need for getting proper arguments
library("Rsubread") # Need for featureCount function
library("data.table") # To properly format output file

# Generate options
opt_bam <- optparse::make_option(c("-f", "--file"),
                                 help = "Full path of input BAM file.")

opt_gtf <- optparse::make_option(c("-g", "--gtf", help = "Full
path of the GTF file."))

opt_srr <- optparse::make_option(c("-s", "--srr"), help = "SRR ID which will be
 used to make column name in the count matrix. If not provided, input
 filename will be used", default = NULL)

opt_out <- optparse::make_option(c("-o", "--out"), help = "Full path of the
output file.", default = "countmatrix.csv")

opt_paired <- optparse::make_option(c("-p", "--paired"), type = "integer",
                                    help = "Is current input file is paired-end? Use 1 for yes [default] and 0 for no
.", default = 1)

opt_threads <- optparse::make_option(c("-t", "--threads"), type = "integer",
                                     help = "Number of threads [Default = 1]", default = 1)

opt_parser <- optparse::OptionParser(option_list = list(opt_bam, opt_gtf,
                                                        opt_paired, opt_srr,
                                                        opt_out, opt_threads))

# Parse the options
opt <- optparse::parse_args(opt_parser)

paried_end <- TRUE
if (opt$paired == 0) {
  paired_end <- FALSE
}

# Create the featureCounts object
fc <- Rsubread::featureCounts(
  opt$file,
  annot.ext = opt$gtf,
  isGTFAnnotationFile = TRUE,
  isPairedEnd = paried_end,
  nthreads = opt$threads)

# Convert to DataFrame
df <- as.data.frame(fc$counts)

# If SRR ID is provided, change the column name.
# Else, input filename will be used as a column name
if (!is.null(opt$srr)) {
  colnames(df) <- opt$srr
}

# Set index column name before saving to the file
data.table::setDT(df, keep.rownames = "gene_id")

# Save to the file
write.csv(df, file = opt$out, row.names = FALSE)

# Create statistics file
stat <- as.data.frame(fc$stat)
log_file <- paste0(opt$out, ".stat")
# If SRR ID is provided, change the column name.
# Else, input filename will be used as a column name
if (!is.null(opt$srr)) {
  colnames(stat)[2] <- opt$srr
}
# Save log to the file
write.csv(stat, file = log_file, row.names = FALSE)
