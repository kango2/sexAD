make_fastq_sample_sheet <- function(filepaths,
                                    output_csv = "samples.csv",
                                    strandedness = "auto") {
  
  # Ensure character vector
  filepaths <- as.character(filepaths)
  
  # Keep only FASTQ files
  filepaths <- filepaths[grepl("\\.fastq\\.gz$", filepaths)]
  
  # Extract filename only
  filenames <- basename(filepaths)
  
  # Identify read direction
  read_type <- ifelse(grepl("_R1_", filenames), "R1",
                      ifelse(grepl("_R2_", filenames), "R2", NA))
  
  if (any(is.na(read_type))) {
    stop("Some files do not contain _R1_ or _R2_ in the name.")
  }
  
  # Extract sample + lane (everything before _R1_001.fastq.gz etc.)
  sample_name <- sub("(_R[12]_001\\.fastq\\.gz)$", "", filenames)
  
  # Create data frame
  df <- data.frame(
    sample = sample_name,
    filepath = filepaths,
    read = read_type,
    stringsAsFactors = FALSE
  )
  
  # Reshape to wide format (R1 and R2 columns)
  df_wide <- reshape(
    df,
    idvar = "sample",
    timevar = "read",
    direction = "wide"
  )
  
  # Clean column names
  colnames(df_wide) <- gsub("filepath\\.", "fastq_", colnames(df_wide))
  
  # Add strandedness column
  df_wide$strandedness <- strandedness
  
  # Reorder columns
  df_wide <- df_wide[, c("sample", "fastq_R1", "fastq_R2", "strandedness")]
  colnames(df_wide)[2:3] <- c("fastq_1", "fastq_2")
  
  # Write CSV
  write.csv(df_wide, output_csv, row.names = FALSE, quote = FALSE)
  
  return(df_wide)
}
