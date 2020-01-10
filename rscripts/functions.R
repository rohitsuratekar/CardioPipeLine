# All functions will go here

library("data.table") # To properly format output file
library("readr") # Reading csv properly
library("Rsubread") # Need for featureCount function
# For importing and generating count matrixes
library("tximport")
suppressPackageStartupMessages(library('AnnotationDbi'))
suppressPackageStartupMessages(library('GenomicFeatures'))
# For DESeq2
suppressPackageStartupMessages(library("DESeq2"))


get_tx2gene <- function(method, all_files, gtf) {
  # Generate tx2gene object
  if (method == "stringtie") {
    # In StringTie this is handled differently.
    # Pick any file to generate tx2gene file
    # We will select first file for sanity check
    tmp <- readr::read_tsv(all_files[1])
    return(tmp[, c("t_name", "gene_id")])
  } else {
    tx <- GenomicFeatures::makeTxDbFromGFF(file = gtf,
                                           format = "gtf",
                                           organism = "Danio rerio")

    k <- AnnotationDbi::keys(tx, keytype = "GENEID")
    tx_df <- AnnotationDbi::select(tx,
                                   keys = k,
                                   keytype = "GENEID",
                                   columns = "TXNAME")
    return(tx_df[, 2:1])
  }
}

save_counts <- function(d, file_name) {
  temp_df <- as.data.frame(d)
  data.table::setDT(temp_df, keep.rownames = "gene_id")
  write.csv(temp_df, file = file_name, row.names = FALSE)
}


generate_txi <- function(method, tx2gene, all_files) {
  # Generate the tximport object which can be used in further analysis
  txi <- tximport::tximport(all_files,
                            type = method,
                            tx2gene = tx2gene,
                            ignoreTxVersion = TRUE)

  return(txi)
}

perform_deseq2_with_txi <- function(txi, config, meta) {
  dds <- DESeq2::DESeqDataSetFromTximport(txi = txi,
                                          colData = meta,
                                          design = formula(config$design))
  return(dds)
}

perform_deseq2_with_counts <- function(counts, config, meta) {
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                        colData = meta,
                                        design = formula(config$design))

  return(dds)
}

collect_deseq2_results <- function(dds, data) {

  dds <- DESeq2::DESeq(dds)

  # Perform LFC Shrinking if contrast is provided
  if (!is.null(data$contrast)) {

    # Get results
    res <- DESeq2::results(dds, contrast = data$contrast, alpha = data$alpha)


    lfc_res <- DESeq2::lfcShrink(dds, contrast = data$contrast, res = res)

    # Save the result
    file_name <- paste0(data$output, "/", data$method, ".result_lfc.csv")
    save_counts(lfc_res, file_name)

    final_res <- lfc_res

  } else {
    # Get results without contrast and alpha
    res <- DESeq2::results(dds)
    final_res <- res
  }

  # Save the result
  file_name <- paste0(data$output, "/", data$method, ".result.csv")
  save_counts(lfc_res, file_name)
  # Save the result
  file_name <- paste0(data$output, "/", data$method, ".result.csv")
  save_counts(lfc_res, file_name)

  return(final_res)

}


generate_count_matrix <- function(all_files, all_runs, data) {

  if (sum(data$paired) == length(all_files)) {
    paried_end <- TRUE
  } else if (sum(data$paired) == 0) {
    paried_end <- FALSE
  } else {
    stop("Mixed data types found in current samples. Current pipeline
    accepts either all Paired end or all Single end")
  }

  fc <- Rsubread::featureCounts(
    all_files,
    annot.ext = data$gtf,
    isGTFAnnotationFile = TRUE,
    isPairedEnd = paried_end,
    nthreads = data$threads)

  df <- as.data.frame(fc$counts)
  colnames(df) <- all_runs
  name <- paste0(data$output, "/", data$method, ".counts.csv")
  save_counts(df, name)
  stat <- paste0(data$output, "/", data$method, ".counts.csv.stat")
  write.csv(fc$stat, file = stat, row.names = FALSE)

  return(df)
}