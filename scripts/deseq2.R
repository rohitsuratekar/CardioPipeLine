# CardioPipeLine (c) 2020
# Author: Rohit Suratekar, ZDG Lab, IIMCB
#
# Functions related to DESeq2 analysis

start_time <- Sys.time()

suppressPackageStartupMessages(library("DESeq2"))

# We can relevel our conditions initially so that we do not have to worry
# about them later on
conds <- factor(snakemake@params$conditions)
conds <- relevel(conds, snakemake@params$reference)

colData <- data.frame(sample = snakemake@params$samples,
                      condition = conds)

df <- read.csv(snakemake@input$combined)
rownames(df) <- df$gene_id
df$gene_id <- NULL

dds <- DESeq2::DESeqDataSetFromMatrix(countData = df,
                                      colData = colData,
                                      design = formula
                                      (snakemake@params$design))


# Start the DE analysis
dds <- DESeq2::DESeq(dds)

outputs <- NULL

# Get results from each combination of conditions
for (g in unique(snakemake@params$conditions)) {
  if (g != snakemake@params$reference) {
    contrast <- c(snakemake@params$column, g, snakemake@params$reference)
    res <- results(dds, contrast = contrast)
    fn <- paste(snakemake@wildcards$BASE, "deseq2", "analysis",
                snakemake@wildcards$METHOD, snakemake@wildcards$METHOD,
                sep = "/")
    name <- paste(fn, g, "vs", snakemake@params$reference, sep = "_")
    name <- paste0(name, ".csv")
    df <- as.data.frame(res)
    data.table::setDT(df, keep.rownames = "gene_id")
    write.csv(df, file = name, row.names = FALSE)
    outputs <- c(outputs, name)
  }
}

# Generate file for snakamake by renaming the input file
summary <- paste(snakemake@wildcards$BASE, "deseq2", "analysis",
                 snakemake@wildcards$METHOD, "analysis.log",
                 sep = "/")


fileConn <- file(summary)
fc <- "Analysis \t: DESeq2"
fc <- c(fc, paste0("Start Time \t: ", start_time))
fc <- c(fc, paste0("Finish Time \t: ", Sys.time()))
fc <- c(fc, paste0("Samples \t: ", paste0(snakemake@params$samples, collapse
  = ",")))
fc <- c(fc, paste0("Conditions \t: ", paste0(snakemake@params$conditions,
                                             collapse = ",")))
fc <- c(fc, paste0("Design Formula \t: ", snakemake@params$design))
fc <- c(fc, paste0("Reference \t: ", snakemake@params$reference))
fc <- c(fc, paste0("Outputs \t: ", paste(outputs, sep = "")))
writeLines(fc, fileConn)
close(fileConn)

