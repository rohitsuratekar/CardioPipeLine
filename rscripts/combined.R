# Combined analysis of DESeq2

library("optparse") # For parsing arguments
library("jsonlite") # To read the json file (configuration file)


# Get all the functions from the function file
# Be careful about the path of the following file.
# It will change based on how you are runnig this script. Currently, it
# assumes that it is run from the python pipeline
source("rscripts/functions.R")


opt_meta <- optparse::make_option(c("-c", "--config"), help = "Full path of the
configuration file.")

opt_parser <- optparse::OptionParser(option_list = list(opt_meta))

# Parse the options
opt <- optparse::parse_args(opt_parser)

data <- jsonlite::fromJSON(opt$config, flatten = TRUE)

all_runs <- NULL
all_files <- NULL
for (s in data$samples) {
  for (run in names(data$input[[s]])) {
    in_file <- data$input[[s]][[run]]
    all_runs <- c(all_runs, run)
    all_files <- c(all_files, in_file)
  }
}

names(all_files) <- all_runs # Tximport will assign these as column names


get_txi <- function() {
  # Get the tx2gene object
  tx2gene <- get_tx2gene(method = data$method,
                         all_files = all_files,
                         gtf = data$gtf)

  # Generate the object
  txi <- generate_txi(method = data$method,
                      tx2gene = tx2gene,
                      all_files = all_files)

  # Save the count matrix for further downstream processing
  file_name <- paste0(data$output, "/", data$method, ".counts.csv")
  save_counts(txi$counts, file_name)

  return(txi)

}

get_meta <- function() {
  meta <- read.csv(data$meta, header = TRUE)
  rownames(meta) <- meta$samples
  meta$samples <- NULL

  return(meta)
}

perform_regular <- function() {

  txi <- get_txi()

  meta <- get_meta()

  # Sanity check
  rn <- rownames(meta)
  cn <- colnames(txi$counts)
  if (!all(rn == cn) || length(rn) != length(cn)) {
    stop("Row names in metadata file and column names in count files should be
 in same order.")
  }

  dds <- perform_deseq2_with_txi(txi, data, meta)
  res <- collect_deseq2_results(dds, data)

  save_transformed_data(dds, data)

  return(res)
}


if (data$method == "star") {
  counts <- generate_count_matrix(all_files, all_runs, data)
  dds <- perform_deseq2_with_counts(counts, data, get_meta())
  res <- collect_deseq2_results(dds, data)
  save_transformed_data(dds, data)
}else {
  res <- perform_regular()
}
