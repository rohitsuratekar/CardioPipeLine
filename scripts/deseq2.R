# CardioPipeLine (c) 2020
# Author: Rohit Suratekar, ZDG Lab, IIMCB
#
# Main DESeq2 analysis script

#library("yaml")
#suppressPackageStartupMessages(library("DESeq2"))
#suppressPackageStartupMessages(library("data.table"))

# Read the sample files and other input
data <- read.csv(snakemake@input$samples)
gtf <- snakemake@input$gtf
mode <- snakemake@params$mode
base_folder <- snakemake@wildcards$BASE
files_names <- snakemake@input$files
# Simple helper functions

get_files <- function(srr) {
  if (mode == "star") {
    file_name <- paste0(srr,)
    return(paste(base_folder, "bams", ""))
  }
}


sample <- NULL
conditions <- NULL
files <- NULL
for (row in 1:nrow(data)) {
  sample_name <- data[row, 'run']
  sample <- c(sample, sample_name)
  conditions <- c(conditions, data[row, 'condition'])
  for (r in files_names) {
    if (grepl(sample_name, r)) {
      files <- c(files, r)
    }
  }
}

names(files) <- sample  # Give column names are respective sample names
colData <- data.frame(sample, conditions) # Create the dataframee

print(files)
print(colData)

