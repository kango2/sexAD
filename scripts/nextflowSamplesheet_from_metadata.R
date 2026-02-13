nextflowSamplesheet_from_metadata <- function(
    fastq_paths,
    metadata_csv,
    output_prefix = "samplesheet",
    strandedness = "auto",
    split_by_sex = TRUE
) {
  
  ## -----------------------------
  ## 1. Read + clean metadata
  ## -----------------------------
  meta <- read.csv(metadata_csv, stringsAsFactors = FALSE, check.names = FALSE)
  
  if (!all(c("Case", "Gender") %in% colnames(meta))) {
    stop("Metadata must contain 'Case' and 'Gender' columns.")
  }
  
  meta$Case_clean <- meta$Case
  meta$Case_clean <- sub(" .*", "", meta$Case_clean)   # remove trailing text
  meta$Case_clean <- gsub("/", "_", meta$Case_clean)  # replace slash
  meta$Case_clean <- trimws(meta$Case_clean)
  meta$Gender <- trimws(meta$Gender)
  
  ## -----------------------------
  ## 2. Process FASTQs
  ## -----------------------------
  fastq_paths <- fastq_paths[grepl("\\.fastq\\.gz$", fastq_paths)]
  filenames <- basename(fastq_paths)
  
  case_id <- sub("^([^_]+_[^_]+)_.*", "\\1", filenames)
  lane <- sub(".*_(L[0-9]{3})_R[12]_001\\.fastq\\.gz$", "\\1", filenames)
  
  read_type <- ifelse(grepl("_R1_", filenames), "R1",
                      ifelse(grepl("_R2_", filenames), "R2", NA))
  
  if (any(is.na(read_type))) {
    stop("Some FASTQ files missing _R1_/_R2_ designation.")
  }
  
  sample_name <- paste(case_id, lane, sep = "_")
  
  df <- data.frame(
    sample = sample_name,
    case = case_id,
    read = read_type,
    filepath = fastq_paths,
    stringsAsFactors = FALSE
  )
  
  ## -----------------------------
  ## 3. Match metadata
  ## -----------------------------
  df$Gender <- meta$Gender[match(df$case, meta$Case_clean)]
  
  missing_cases <- unique(df$case[is.na(df$Gender)])
  if (length(missing_cases) > 0) {
    warning("Cases not found in metadata:\n",
            paste(missing_cases, collapse = ", "))
  }
  
  ## -----------------------------
  ## 4. Reshape to wide
  ## -----------------------------
  df_wide <- reshape(
    df,
    idvar = c("sample", "case", "Gender"),
    timevar = "read",
    direction = "wide"
  )
  
  colnames(df_wide) <- gsub("filepath\\.", "fastq_", colnames(df_wide))
  df_wide$strandedness <- strandedness
  
  df_wide <- df_wide[, c("sample", "fastq_R1", "fastq_R2", "strandedness", "Gender")]
  colnames(df_wide)[2:3] <- c("fastq_1", "fastq_2")
  
  ## -----------------------------
  ## 5. Write outputs
  ## -----------------------------
  
  # Always write combined (with Gender column)
  write.csv(df_wide,
            paste0(output_prefix, "_combined.csv"),
            row.names = FALSE, quote = FALSE)
  
  if (split_by_sex) {
    
    females <- subset(df_wide, Gender == "Female")
    males   <- subset(df_wide, Gender == "Male")
    
    # Remove Gender column for split outputs
    females <- females[, c("sample", "fastq_1", "fastq_2", "strandedness")]
    males   <- males[, c("sample", "fastq_1", "fastq_2", "strandedness")]
    
    write.csv(females,
              paste0(output_prefix, "_females.csv"),
              row.names = FALSE, quote = FALSE)
    
    write.csv(males,
              paste0(output_prefix, "_males.csv"),
              row.names = FALSE, quote = FALSE)
  }
  
  return(invisible(df_wide))
}

