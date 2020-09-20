# CardioPipeLine (c) 2020
# Author: Rohit Suratekar, ZDG Lab, IIMCB
#
# Generates count files for every sample

library("tximport")
suppressPackageStartupMessages(library("AnnotationDbi"))

tx <- GenomicFeatures::makeTxDbFromGFF(file = snakemake@input$gtf,
                                       format = "gtf",
                                       organism = snakemake@params$animal)

k <- AnnotationDbi::keys(tx, keytype = "GENEID")
tx_df <- AnnotationDbi::select(tx,
                               keys = k,
                               keytype = "GENEID",
                               columns = "TXNAME")
tx2gene <- tx_df[, 2:1]

txi <- tximport::tximport(snakemake@input$file,
                          type = snakemake@params$method,
                          tx2gene = tx2gene,
                          ignoreTxVersion = TRUE)

colnames(txi$counts) <- snakemake@wildcards$SRR_ID
df <- data.frame(txi$counts)
data.table::setDT(df, keep.rownames = "gene_id")

fname <- paste(snakemake@wildcards$BASE, "deseq2", "counts",
               snakemake@wildcards$SRR_ID, sep = "/")
fname <- paste0(fname, "/", snakemake@wildcards$SRR_ID, ".",
                snakemake@params$method, ".counts")
write.csv(df, file = fname, row.names = FALSE)
