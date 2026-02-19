nextflowContrast_from_metadata <- function(
    fastq_paths,
    metadata_csv,
    output_file = "contrast_sheet.csv",
    condition_label = "control"
) {
  
  ## -----------------------------
  ## 1. Read + clean metadata
  ## -----------------------------
  meta <- read.csv(metadata_csv, stringsAsFactors = FALSE, check.names = FALSE)
  
  if (!all(c("Case", "Gender", "Group") %in% colnames(meta))) {
    stop("Metadata must contain 'Case', 'Gender', and 'Group' columns.")
  }
  
  meta$Case_clean <- meta$Case
  meta$Case_clean <- sub(" .*", "", meta$Case_clean)   # remove trailing text (e.g. " F")
  meta$Case_clean <- gsub("/", "_", meta$Case_clean)  # replace slash
  meta$Case_clean <- trimws(meta$Case_clean)
  
  meta$Gender <- trimws(meta$Gender)
  meta$Group  <- trimws(meta$Group)
  
  ## -----------------------------
  ## 2. Process FASTQs
  ## -----------------------------
  fastq_paths <- fastq_paths[grepl("\\.fastq\\.gz$", fastq_paths)]
  filenames <- basename(fastq_paths)
  
  case_id <- sub("^([^_]+_[^_]+)_.*", "\\1", filenames)
  lane    <- sub(".*_(L[0-9]{3})_R[12]_001\\.fastq\\.gz$", "\\1", filenames)
  
  sample_name <- paste(case_id, lane, sep = "_")
  
  df <- data.frame(
    sample = sample_name,
    case   = case_id,
    stringsAsFactors = FALSE
  )
  
  ## Remove duplicate R1/R2 rows (keep one per sample)
  df <- df[!duplicated(df$sample), ]
  
  ## -----------------------------
  ## 3. Match metadata
  ## -----------------------------
  df$gender <- meta$Gender[match(df$case, meta$Case_clean)]
  df$group  <- meta$Group[match(df$case, meta$Case_clean)]
  
  missing_cases <- unique(df$case[is.na(df$gender)])
  if (length(missing_cases) > 0) {
    warning("Cases not found in metadata:\n",
            paste(missing_cases, collapse = ", "))
  }
  
  ## -----------------------------
  ## 4. Build contrast sheet
  ## -----------------------------
  contrast_df <- data.frame(
    file      = df$sample,
    condition = condition_label,
    group     = df$group,
    gender    = df$gender,
    stringsAsFactors = FALSE
  )
  
  ## -----------------------------
  ## 5. Write output
  ## -----------------------------
  write.csv(
    contrast_df,
    output_file,
    row.names = FALSE,
    quote = FALSE
  )
  
  return(invisible(contrast_df))
}
