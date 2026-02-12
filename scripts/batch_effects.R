batch_1<- contrast_2merged_2
batch_1$file 
for(i in 1:nrow(batch_1)){
  batch_1[i,1] <- paste0(batch_1[i,1], "_batch2") 
}



for(things in contrast_2merged_2$file){
  things <- paste0(things, "_batch2")
}
contrast_2merged_2
contrastBatchMerged <- rbind(contrastmerged,batch_1)
contrastBatchMerged$file

# male, female batch (1)
# male, female batch (2)

quantsfmerged
quantsf_2merged_2
quantsfmerged_batch <- c(quantsfmerged,quantsf_2merged_2)

file.list <- c("/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out_run2/workdirectory/star_salmon/tx2gene.tsv",
               "/g/data/xl04/eh8642/RNAseq_run/Male_RNAseq_out_run2/workdirectory/star_salmon/tx2gene.tsv", 
               "/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out/workdirectory/star_salmon/tx2gene.tsv",
               "/g/data/xl04/eh8642/RNAseq_run/Male_RNAseq_out/workdirectory/star_salmon/tx2gene.tsv")

tx2gene_batch <- map_dfr(file.list, read_tsv) %>%
  unique()
colnames(tx2gene_batch) <- c("transcript_id","gene_name","gene_id")


txi_batch <- tximport(
  quantsfmerged_batch,
  type = "salmon",
  txOut = FALSE,
  tx2gene = tx2gene_batch,
  countsFromAbundance = "no"  # <-- critical
)


dds_batch <- DESeqDataSetFromTximport(
  txi_batch,
  colData = contrastBatchMerged,
  design = ~ batch
)

gene_names_batch <- tx2gene_batch |>
  dplyr::select(gene_id, gene_name) |>
  dplyr::distinct()

rowData(dds_batch)$gene_name <- gene_names_batch$gene_name[
  match(rownames(dds_batch), gene_names_batch$gene_id)
]

rowData(dds_batch)$gene_name <- ifelse(
  is.na(rowData(dds_batch)$gene_name),
  rownames(dds_batch),
  rowData(dds_batch)$gene_name
)

rowData(dds_batch)$gene_name <- make.unique(rowData(dds_batch)$gene_name)


smallestGroupSize <- 20
keep <- rowSums(counts(dds_batch) >= 10) >= smallestGroupSize
dds_batch <- dds_batch[keep,]


dds_batch$batch <- factor(dds_batch$batch, levels = c("batch_1","batch_2"))

dds_batch <- DESeq(dds_batch)
res_batch <- results(dds_batch)

#res_2 <- res_2ults(dds_batch, contrast_2=c("condition","control","treated"))
resultsNames(dds_batch)
dds_batch$batch
resLFC_batch <- lfcShrink(dds_batch, coef="batch_batch_2_vs_batch_1", type="apeglm")

res_ordered_batch <- res_batch[order(res_batch$pvalue),]

sum(res_batch$padj < 0.1, na.rm=TRUE)

res_batch.05 <- results(dds_batch, alpha=0.05)
summary(res_batch.05)
sum(res_205$padj < 0.05, na.rm=TRUE)

plotMA(res_2, ylim=c(-2,2))
plot_MA <- plotMA(res_2LFC, ylim=c(-2,2))

legend("topright", # Position keyword
       legend = levels(res_2),
       col = c("blue", "grey"),
       pch = 19, # Must match the symbol used in plot()
       title = "Group")


#idx <- identify(res_2$baseMean, res_2$log2FoldChange)
#rownames(res_2)[idx]

plotCounts(dds_batch, gene=which.min(res_batch$padj), intgroup="batch")
plotDispEsts(dds_batch)

res_ordered_batch <- as.data.frame(res_ordered_batch) |>
  unique() |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::mutate(
    gene_name = rowData(dds_batch)$gene_name[
      match(gene_id, rownames(dds_batch))
    ]
  ) |>
  dplyr::relocate(gene_name, .before = gene_id)

res_2Sig <- subset(res_2Ordered, padj < 0.05)
res_2Sig <- as.data.frame(res_2Sig) |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::mutate(
    gene_name = rowData(dds_batch)$gene_name[
      match(gene_id, rownames(dds_batch))
    ]
  ) |>
  dplyr::relocate(gene_name, .before = gene_id)


# PCA plot 
vsd <- vst(dds_batch, blind = FALSE)
plotPCA(vsd, intgroup="group")

write.csv(as.data.frame(res_ordered_batch), 
          file="/g/data/xl04/eh8642/RNAseq_run/deSeq_results/Batch_effects_p_vals.csv")





################################################################
################################################################
######################CHATGPT CODE##############################
################################################################
################################################################


run_deseq_batch <- function(
    coldata_batch1,
    coldata_batch2,
    quants_batch1,
    quants_batch2,
    tx2gene_files,
    batch_labels = c("batch_1", "batch_2"),
    min_count = 10,
    min_samples = 20,
    alpha = 0.05,
    lfc_coef = "batch_batch_2_vs_batch_1",
    output_csv = NULL
) {
  
  suppressPackageStartupMessages({
    library(DESeq2)
    library(tximport)
    library(dplyr)
    library(readr)
    library(purrr)
    library(tibble)
  })
  
  # -----------------------------
  # 1. Prepare colData
  # -----------------------------
  
  coldata_batch1$batch <- batch_labels[1]
  coldata_batch2$batch <- batch_labels[2]
  
  coldata_merged <- dplyr::bind_rows(coldata_batch1, coldata_batch2)
  coldata_merged$batch <- factor(
    coldata_merged$batch,
    levels = batch_labels
  )
  
  # -----------------------------
  # 2. Merge quant files
  # -----------------------------
  
  quants_merged <- c(quants_batch1, quants_batch2)
  
  # -----------------------------
  # 3. Build tx2gene
  # -----------------------------
  
  tx2gene <- purrr::map_dfr(tx2gene_files, readr::read_tsv) |>
    dplyr::distinct()
  
  colnames(tx2gene) <- c("transcript_id", "gene_name", "gene_id")
  
  # -----------------------------
  # 4. tximport
  # -----------------------------
  
  txi <- tximport(
    quants_merged,
    type = "salmon",
    txOut = FALSE,
    tx2gene = tx2gene,
    countsFromAbundance = "no"
  )
  
  # -----------------------------
  # 5. Build DESeq object
  # -----------------------------
  
  dds <- DESeqDataSetFromTximport(
    txi,
    colData = coldata_merged,
    design = ~ batch
  )
  
  # Attach gene names
  gene_names <- tx2gene |>
    dplyr::select(gene_id, gene_name) |>
    dplyr::distinct()
  
  rowData(dds)$gene_name <- gene_names$gene_name[
    match(rownames(dds), gene_names$gene_id)
  ]
  
  rowData(dds)$gene_name[is.na(rowData(dds)$gene_name)] <- rownames(dds)
  rowData(dds)$gene_name <- make.unique(rowData(dds)$gene_name)
  
  # -----------------------------
  # 6. Filter low counts
  # -----------------------------
  
  keep <- rowSums(counts(dds) >= min_count) >= min_samples
  dds <- dds[keep, ]
  
  # -----------------------------
  # 7. Run DESeq
  # -----------------------------
  
  dds <- DESeq(dds)
  
  res <- results(dds, alpha = alpha)
  res_lfc <- lfcShrink(dds, coef = lfc_coef, type = "apeglm")
  
  # -----------------------------
  # 8. Order + tidy results
  # -----------------------------
  
  res_ordered <- as.data.frame(res[order(res$pvalue), ]) |>
    tibble::rownames_to_column("gene_id") |>
    dplyr::mutate(
      gene_name = rowData(dds)$gene_name[
        match(gene_id, rownames(dds))
      ]
    ) |>
    dplyr::relocate(gene_name, .before = gene_id)
  
  res_sig <- subset(res_ordered, padj < alpha)
  
  # -----------------------------
  # 9. Optional write to disk
  # -----------------------------
  
  if (!is.null(output_csv)) {
    write.csv(res_ordered, file = output_csv, row.names = FALSE)
  }
  
  # -----------------------------
  # 10. Return structured output
  # -----------------------------
  
  return(list(
    dds = dds,
    results = res,
    results_lfc = res_lfc,
    ordered_results = res_ordered,
    significant = res_sig
  ))
}


















