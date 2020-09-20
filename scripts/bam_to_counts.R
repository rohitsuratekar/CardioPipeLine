# CardioPipeLine (c) 2020
# Author: Rohit Suratekar, ZDG Lab, IIMCB
#
# Generates count files for every sample

suppressPackageStartupMessages(library("data.table"))
library("Rsubread")


fc <- Rsubread::featureCounts(
  files = snakemake@input$bam,
  annot.ext = snakemake@input$gtf,
  isGTFAnnotationFile = TRUE,
  isPairedEnd = as.logical(snakemake@params$is_paired),
  nthreads = snakemake@threads)

colnames(fc$counts) <- snakemake@wildcards$SRR_ID
df <- data.frame(fc$counts)
data.table::setDT(df, keep.rownames = "gene_id")

fname <- paste(snakemake@wildcards$BASE, "deseq2", "counts", snakemake@wildcards$SRR_ID, sep = "/")
fname <- paste0(fname, "/", snakemake@wildcards$SRR_ID, ".star.counts")
write.csv(df, file = fname, row.names = FALSE)
write.csv(fc$stat, file = paste0(fname, ".stat"), row.names = FALSE)



