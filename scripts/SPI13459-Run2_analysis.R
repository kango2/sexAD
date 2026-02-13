library(BiocManager)
library(tidyverse)
library(DESeq2)
library(tximport)
library(tidyverse)
library("vsn")
library("pheatmap")
library("apeglm")
library(org.Hs.eg.db)
library(AnnotationDbi)


#### Intersex analysis ##### 
Femaleinputnames_2 <- list.files("/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out_run2/workdirectory/star_salmon", pattern="*.bai") %>%
  {. <- gsub(".markdup.sorted.bam.bai","",.);.}
maleinputnames_2 <- list.files("/g/data/xl04/eh8642/RNAseq_run/Male_RNAseq_out_run2/workdirectory/star_salmon", pattern="*.bai") %>%
  {. <- gsub(".markdup.sorted.bam.bai","",.);.}

Femalequantsf_2 <- list.files("/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out_run2/workdirectory/star_salmon", pattern="*.bai",full.names = T) %>%
  {. <- gsub(".markdup.sorted.bam.bai","/quant.sf",.);.}
malequantsf_2 <- list.files("/g/data/xl04/eh8642/RNAseq_run/Male_RNAseq_out_run2/workdirectory/star_salmon", pattern="*.bai",full.names = T) %>%
  {. <- gsub(".markdup.sorted.bam.bai","/quant.sf",.);.}

quantsf_2merged_2 <- NULL

quantsf_2merged_2 <- c(malequantsf_2,Femalequantsf_2)
quantsf_2merged_2

Femalecontrast_2 <- as_tibble(Femaleinputnames_2)
Femalecontrast_2$group <- as.character(rep(0,nrow(Femalecontrast_2)))
Femalecontrast_2$sex <- rep("F",nrow(Femalecontrast_2))

for(i in 1:nrow(Femalecontrast_2)){
  if(grepl("Control",Femalecontrast_2[i,1])==TRUE){
    Femalecontrast_2[i,2] <- "control"
  }else{
    Femalecontrast_2[i,2] <- "affected"
  }
}

malecontrast_2 <- as_tibble(maleinputnames_2)
malecontrast_2$group <- as.character(rep(0,nrow(malecontrast_2)))
malecontrast_2$sex <- rep("M",nrow(malecontrast_2))

for(i in 1:nrow(malecontrast_2)){
  if(grepl("Control",malecontrast_2[i,1])==TRUE){
    malecontrast_2[i,2] <- "control"
  }else{
    malecontrast_2[i,2] <- "affected"
  }
}

contrast_2merged_2 <- rbind(malecontrast_2,Femalecontrast_2)
colnames(contrast_2merged_2) <- c("file","condition","group")
contrast_2merged_2$condition <- factor(contrast_2merged_2$condition,levels=c("control","affected"))
contrast_2merged_2$group <- factor(contrast_2merged_2$group,levels=c("M","F"))
contrast_2merged_2$batch <- rep("batch_2",nrow(contrast_2merged_2))

Femaletx2gene_2 <- read_tsv(
  "/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out_run2/workdirectory/star_salmon/tx2gene.tsv",
  col_names = TRUE
)
Maletx2gene_2 <- read_tsv(
  "/g/data/xl04/eh8642/RNAseq_run/Male_RNAseq_out_run2/workdirectory/star_salmon/tx2gene.tsv",
  col_names = TRUE
)
file.list <- c("/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out_run2/workdirectory/star_salmon/tx2gene.tsv",
               "/g/data/xl04/eh8642/RNAseq_run/Male_RNAseq_out_run2/workdirectory/star_salmon/tx2gene.tsv")

merged_2tx2gene <- map_dfr(file.list, read_tsv) %>%
  unique()
colnames(merged_2tx2gene) <- c("transcript_id","gene_name","gene_id")

merged_2txi <- tximport(
  quantsf_2merged_2,
  type = "salmon",
  txOut = FALSE,
  tx2gene = merged_2tx2gene,
  countsFromAbundance = "no"  # <-- critical
)


dds_2 <- DESeqDataSetFromTximport(
  merged_2txi,
  colData = contrast_2merged_2,
  design = ~ group
)

gene_names_2 <- merged_2tx2gene |>
  dplyr::select(gene_id, gene_name) |>
  dplyr::distinct()

rowData(dds_2)$gene_name <- gene_names_2$gene_name[
  match(rownames(dds_2), gene_names_2$gene_id)
]

rowData(dds_2)$gene_name <- ifelse(
  is.na(rowData(dds_2)$gene_name),
  rownames(dds_2),
  rowData(dds_2)$gene_name
)

rowData(dds_2)$gene_name <- make.unique(rowData(dds_2)$gene_name)


smallestGroupSize <- 20
keep <- rowSums(counts(dds_2) >= 10) >= smallestGroupSize
keep
dds_3 <- dds_2[keep,]
dds_2$
dds_2$condition <- factor(dds_2$condition, levels = c("control","affected"))

dds_2 <- DESeq(dds_2)
res_2 <- results(dds_2)
#res_2 <- res_2ults(dds_2, contrast_2=c("condition","control","treated"))
resultsNames(dds_2)

res_2LFC <- lfcShrink(dds_2, coef="group_F_vs_M", type="apeglm")

res_2Ordered <- res_2[order(res_2$pvalue),]
res_2Ordered

sum(res_2$padj < 0.1, na.rm=TRUE)

res_205 <- results(dds_2, alpha=0.05)
summary(res_205)
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

plotCounts(dds_2, gene=which.min(res_2$padj), intgroup="group")
plotDispEsts(dds_2)

res_2Ordered <- as.data.frame(res_2Ordered) |>
  unique() |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::mutate(
    gene_name = rowData(dds_2)$gene_name[
      match(gene_id, rownames(dds_2))
    ]
  ) |>
  dplyr::relocate(gene_name, .before = gene_id)

res_2Sig <- subset(res_2Ordered, padj < 0.05)
res_2Sig <- as.data.frame(res_2Sig) |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::mutate(
    gene_name = rowData(dds_2)$gene_name[
      match(gene_id, rownames(dds_2))
    ]
  ) |>
  dplyr::relocate(gene_name, .before = gene_id)



write.csv(as.data.frame(res_2Ordered), 
          file="/g/data/xl04/eh8642/RNAseq_run/deSeq_results/F_M_rnaseq_results_batch2.csv")
write.csv(as.data.frame(res_2Sig), 
          file="/g/data/xl04/eh8642/RNAseq_run/deSeq_results/F_M_rnaseq_results_p0.05_batch2.csv")
vsd <- vst(dds_2, blind = FALSE)
plotPCA(vsd)
##################### Female comparison #####################

Femaletxi <- tximport(
  Femalequantsf_2,
  type = "salmon",
  txOut = FALSE,
  tx2gene = Femaletx2gene,
  countsFromAbundance = "no"  # <-- critical
)

Femalecontrast_2$group <- factor(Femalecontrast_2$group,levels=c("control","affected"))

Female_dds_2 <- DESeqDataSetFromTximport(
  Femaletxi,
  colData = Femalecontrast_2,
  design = ~ group
)

gene_names_2_Female <- Femaletx2gene |>
  dplyr::select(gene_id, gene_name) |>
  dplyr::distinct()

rowData(Female_dds_2)$gene_name <- gene_names_2_Female$gene_name[
  match(rownames(Female_dds_2), gene_names_2_Female$gene_id)
]

rowData(Female_dds_2)$gene_name <- ifelse(
  is.na(rowData(Female_dds_2)$gene_name),
  rownames(Female_dds_2),
  rowData(Female_dds_2)$gene_name
)

rowData(Female_dds_2)$gene_name <- make.unique(rowData(Female_dds_2)$gene_name)

smallestGroupSize <- 10
keep <- rowSums(counts(Female_dds_2) >= 10) >= smallestGroupSize
Female_dds_2 <- Female_dds_2[keep,]


Female_dds_2$group <- factor(Female_dds_2$group, levels = c("control","affected"))
Female_dds_2 <- DESeq(Female_dds_2)
res_2_Female <- res_2ults(Female_dds_2)
res_2_Female  <- res_2ults(Female_dds_2, contrast_2=c("group","control","affected"))
res_2ultsNames(Female_dds_2)

res_2LFC_Female <- lfcShrink(Female_dds_2, coef="group_affected_vs_control", type="apeglm")
res_2Ordered_Female <- res_2_Female[order(res_2_Female$pvalue),]
res_2Ordered_Female

sum(res_2_Female$padj < 0.1, na.rm=TRUE)
Femaleres_205 <- res_2ults(Female_dds_2, alpha=0.05)

sum(Femaleres_205$padj < 0.05, na.rm=TRUE)

plotMA(res_2_Female, ylim=c(-2,2))
plot_MA <- plotMA(res_2LFC_Female,ylim=c(-2,2))

legend("topright", # Position keyword
       legend = levels(res_2),
       col = c("blue", "grey"),
       pch = 19, # Must match the symbol used in plot()
       title = "Group")


#Female_idx <- identify(res_2_Female$baseMean, res_2_Female$log2FoldChange)
#rownames(res_2_Female)[Female_idx]

plotCounts(Female_dds_2, gene=which.min(res_2_Female$padj), intgroup="group") 

res_2Ordered_Female <- as.data.frame(res_2Ordered_Female) |>
  unique() |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::mutate(
    gene_name = rowData(Female_dds_2)$gene_name[
      match(gene_id, rownames(Female_dds_2))
    ]
  ) |>
  dplyr::relocate(gene_name, .before = gene_id)

res_2Sig_Female <- subset(res_2Ordered_Female, padj < 0.05)
res_2Sig_Female <- as.data.frame(res_2Sig_Female) |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::mutate(
    gene_name = rowData(Female_dds_2)$gene_name[
      match(gene_id, rownames(Female_dds_2))
    ]
  ) |>
  dplyr::relocate(gene_name, .before = gene_id)

plotDispEsts(Female_dds_2)

write.csv(as.data.frame(res_2Ordered_Female), 
          file="/g/data/xl04/eh8642/RNAseq_run/deSeq_res_2ults/F_control_affected_rnaseq_res_2ults.csv")
write.csv(as.data.frame(res_2Sig_Female), 
          file="/g/data/xl04/eh8642/RNAseq_run/deSeq_res_2ults/F_control_affected_rnaseq_res_2ults_p0.05.csv")


##################### Male comparison #####################

maletxi <- tximport(
  malequantsf_2,
  type = "salmon",
  txOut = FALSE,
  tx2gene = Maletx2gene,
  countsFromAbundance = "no"  # <-- critical
)

malecontrast_2$group <- factor(malecontrast_2$group,levels=c("control","affected"))

male_dds_2 <- DESeqDataSetFromTximport(
  maletxi,
  colData = malecontrast_2,
  design = ~ group
)

gene_names_2_male <- maletx2gene |>
  dplyr::select(gene_id, gene_name) |>
  dplyr::distinct()

rowData(male_dds_2)$gene_name <- gene_names_2_male$gene_name[
  match(rownames(male_dds_2), gene_names_2_male$gene_id)
]

rowData(male_dds_2)$gene_name <- ifelse(
  is.na(rowData(male_dds_2)$gene_name),
  rownames(male_dds_2),
  rowData(male_dds_2)$gene_name
)

rowData(male_dds_2)$gene_name <- make.unique(rowData(male_dds_2)$gene_name)

smallestGroupSize <- 10
keep <- rowSums(counts(male_dds_2) >= 10) >= smallestGroupSize
male_dds_2 <- male_dds_2[keep,]


male_dds_2$group <- factor(male_dds_2$group, levels = c("control","affected"))
male_dds_2 <- DESeq(male_dds_2)
res_2_male <- res_2ults(male_dds_2)
res_2_male  <- res_2ults(male_dds_2, contrast_2=c("group","control","affected"))
res_2ultsNames(male_dds_2)

res_2LFC_male <- lfcShrink(male_dds_2, coef="group_affected_vs_control", type="apeglm")
res_2Ordered_male <- res_2_male[order(res_2_male$pvalue),]
res_2Ordered_male

sum(res_2_male$padj < 0.1, na.rm=TRUE)
maleres_205 <- res_2ults(male_dds_2, alpha=0.05)

sum(maleres_205$padj < 0.05, na.rm=TRUE)

plotMA(res_2_male, ylim=c(-2,2))
plot_MA <- plotMA(res_2LFC_male,ylim=c(-2,2))

legend("topright", # Position keyword
       legend = levels(res_2),
       col = c("blue", "grey"),
       pch = 19, # Must match the symbol used in plot()
       title = "Group")


#male_idx <- identify(res_2_male$baseMean, res_2_male$log2FoldChange)
#rownames(res_2_male)[male_idx]

plotCounts(male_dds_2, gene=which.min(res_2_male$padj), intgroup="group") 

res_2Ordered_male <- as.data.frame(res_2Ordered_male) |>
  unique() |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::mutate(
    gene_name = rowData(male_dds_2)$gene_name[
      match(gene_id, rownames(male_dds_2))
    ]
  ) |>
  dplyr::relocate(gene_name, .before = gene_id)

res_2Sig_male <- subset(res_2Ordered_male, padj < 0.05)
res_2Sig_male <- as.data.frame(res_2Sig_male) |>
  unique() |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::mutate(
    gene_name = rowData(male_dds_2)$gene_name[
      match(gene_id, rownames(male_dds_2))
    ]
  ) |>
  dplyr::relocate(gene_name, .before = gene_id)

plotDispEsts(male_dds_2)

write.csv(as.data.frame(res_2Ordered_male), 
          file="/g/data/xl04/eh8642/RNAseq_run/deSeq_res_2ults/M_control_affected_rnaseq_res_2ults.csv")
write.csv(as.data.frame(res_2Sig_male), 
          file="/g/data/xl04/eh8642/RNAseq_run/deSeq_res_2ults/M_control_affected_rnaseq_res_2ults_p0.05.csv")




##################### Multi-factor #####################

dds_2MF <- dds_2
design(dds_2MF) <- formula(~ condition + group)
dds_2MF <- DESeq(dds_2MF)
res_2MF <- res_2ults(dds_2MF)


# this gives log2(n + 1)
ntd <- normTransform(dds_2MF)
#colnames(ntd) <- contrast_2merged_2$file
meanSdPlot(assay(ntd))
select <- order(rowMeans(counts(dds_2MF,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_2)[,c("condition","group")])

mat <- assay(vst(dds_2MF))[select, ]
rownames(mat) <- rowData(dds_2)$gene_name[select]
# ensure uniqueness
rownames(mat) <- make.unique(rownames(mat))

pheatmap(mat, cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df, show_colnames = TRUE)



ntd <- normTransform(dds_2)
#colnames(ntd) <- contrast_2merged_2$file
meanSdPlot(assay(ntd))
select <- order(rowMeans(counts(dds_2,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_2)[,"group"])
colnames(df) <- c("group")
mat <- assay(vst(dds_2))[select, ]
rownames(mat) <- rowData(dds_2)$gene_name[select]
# ensure uniqueness
rownames(mat) <- make.unique(rownames(mat))

pheatmap(mat, cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df, show_colnames = TRUE)





























