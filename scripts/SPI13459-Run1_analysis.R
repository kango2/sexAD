library(BiocManager)
library(tidyverse)
library(tximport)
library(tidyverse)
library("vsn")
library("pheatmap")
library("apeglm")
library(org.Hs.eg.db)
library(AnnotationDbi)


female_rawcounts <- read_tsv("/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out/workdirectory/star_salmon/salmon.merged.gene_counts.tsv")
female_rawcounts
female_salmon <- tximport("/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out/workdirectory/star_salmon/salmon.merged.gene_counts.tsv",
                          type="salmon", 
                          tx2gene="/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out/workdirectory/star_salmon/tx2gene.tsv")


female_controls <- list.files("/g/data/xl04/hrp561/adrna/analyses/alnbam", pattern="*Control_F_", full.names = F) %>%
  {. <- gsub(".bam.bai", "",.);.} %>%
  {. <- gsub(".bam", "",.);.} %>%
  unique()
female_treatment <- list.files("/g/data/xl04/hrp561/adrna/analyses/alnbam", pattern="*Familial_F_", full.names = F) %>%
  {. <- gsub(".bam.bai", "",.);.} %>%
  {. <- gsub(".bam", "",.);.} %>%
  unique()

female_contrasts<- read_csv("/g/data/xl04/eh8642/RNAseq_run/script_template/Female_contrasts.csv")
female_contrasts$condition <- as.factor(female_contrasts$condition)
female_contrasts$group <- as.factor(female_contrasts$group)

female_rawcountsSelect <- female_rawcounts[, female_contrasts$file]

rownames(female_rawcountsSelect) <- female_rawcounts$gene_id
all(female_contrasts$file %in% colnames(female_rawcounts))

femaleReadin <- DESeqDataSetFromTximport("/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out/workdirectory/star_salmon/salmon.merged.gene_counts.tsv", 
                                         female_contrasts, 
                                         ~condition)


dirs <- list.dirs("/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out/workdirectory/star_salmon")
dirs
controls <- grep("Control",dirs)
dirs_controls <- dirs[controls]
samples <- c("Control_F_1_L001","Control_F_1_L002","Familial_F_10_L001","Familial_F_10_L002")
tests <- c("/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out/workdirectory/star_salmon/Control_F_1_L001",
           "/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out/workdirectory/star_salmon/Control_F_1_L002",  
           "/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out/workdirectory/star_salmon/Familial_F_10_L001",
           "/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out/workdirectory/star_salmon/Familial_F_10_L002")
quantsf <- c()
for(dirs in tests){
  quantsf <- c(quantsf,paste0(dirs,"/quant.sf"))
}
names(quantsf) <- samples

tx2gene <- read_tsv(
  "/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out/workdirectory/star_salmon/tx2gene.tsv",
  col_names = TRUE
)

txi <- tximport(
  quantsf,
  type = "salmon",
  txOut = FALSE,
  tx2gene = tx2gene,
  countsFromAbundance = "no"  # <-- critical
)
female_contrasts <- female_contrasts[c(1,2,11,12),]

dds <- DESeqDataSetFromTximport(
  txi,
  colData = female_contrasts,
  design = ~ condition
)
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$condition <- factor(dds$condition, levels = c("control","treated"))

dds <- DESeq(dds)
res <- results(dds)
res <- results(dds, contrast=c("condition","control","treated"))
resultsNames(dds)

resLFC <- lfcShrink(dds, coef="condition_treated_vs_control", type="apeglm")
resLFC


#### Intersex analysis ##### 
femaleinputnames <- list.files("/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out/workdirectory/star_salmon", pattern="*.bai") %>%
  {. <- gsub(".markdup.sorted.bam.bai","",.);.}

maleinputnames <- list.files("/g/data/xl04/eh8642/RNAseq_run/Male_RNAseq_out/workdirectory/star_salmon", pattern="*.bai") %>%
  {. <- gsub(".markdup.sorted.bam.bai","",.);.}

femalequantsf <- list.files("/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out/workdirectory/star_salmon", pattern="*.bai",full.names = T) %>%
  {. <- gsub(".markdup.sorted.bam.bai","/quant.sf",.);.}
malequantsf <- list.files("/g/data/xl04/eh8642/RNAseq_run/Male_RNAseq_out/workdirectory/star_salmon", pattern="*.bai",full.names = T) %>%
  {. <- gsub(".markdup.sorted.bam.bai","/quant.sf",.);.}

quantsfmerged <- merge(malequantsf,femalequantsf)
quantsfmerged <- c(malequantsf,femalequantsf)

femalecontrast <- as_tibble(femaleinputnames)
femalecontrast$group <- as.character(rep(0,nrow(femalecontrast)))
femalecontrast$sex <- rep("F",nrow(femalecontrast))

for(i in 1:nrow(femalecontrast)){
  if(grepl("Control",femalecontrast[i,1])==TRUE){
    femalecontrast[i,2] <- "control"
  }else{
    femalecontrast[i,2] <- "affected"
  }
}

malecontrast <- as_tibble(maleinputnames)
malecontrast$group <- as.character(rep(0,nrow(malecontrast)))
malecontrast$sex <- rep("M",nrow(malecontrast))

for(i in 1:nrow(malecontrast)){
  if(grepl("Control",malecontrast[i,1])==TRUE){
    malecontrast[i,2] <- "control"
  }else{
    malecontrast[i,2] <- "affected"
  }
}

contrastmerged <- rbind(malecontrast,femalecontrast)
colnames(contrastmerged) <- c("file","condition","group")
contrastmerged$condition <- factor(contrastmerged$condition,levels=c("control","affected"))
contrastmerged$group <- factor(contrastmerged$group,levels=c("M","F"))
contrastmerged$batch <- rep("batch_1", nrow(contrastmerged))

Femaletx2gene <- read_tsv(
  "/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out/workdirectory/star_salmon/tx2gene.tsv",
  col_names = TRUE
)
Maletx2gene <- read_tsv(
  "/g/data/xl04/eh8642/RNAseq_run/Male_RNAseq_out/workdirectory/star_salmon/tx2gene.tsv",
  col_names = TRUE
)
file.list <- c("/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out/workdirectory/star_salmon/tx2gene.tsv",
               "/g/data/xl04/eh8642/RNAseq_run/Male_RNAseq_out/workdirectory/star_salmon/tx2gene.tsv")

mergedtx2gene <- map_dfr(file.list, read_tsv) %>%
  unique()
colnames(mergedtx2gene) <- c("transcript_id","gene_name","gene_id")

Mergedtxi <- tximport(
  quantsfmerged,
  type = "salmon",
  txOut = FALSE,
  tx2gene = mergedtx2gene,
  countsFromAbundance = "no"  # <-- critical
)

dds <- DESeqDataSetFromTximport(
  Mergedtxi,
  colData = contrastmerged,
  design = ~ group
)

gene_names <- mergedtx2gene |>
  dplyr::select(gene_id, gene_name) |>
  dplyr::distinct()

rowData(dds)$gene_name <- gene_names$gene_name[
  match(rownames(dds), gene_names$gene_id)
]

rowData(dds)$gene_name <- ifelse(
  is.na(rowData(dds)$gene_name),
  rownames(dds),
  rowData(dds)$gene_name
)

rowData(dds)$gene_name <- make.unique(rowData(dds)$gene_name)


smallestGroupSize <- 20
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$condition <- factor(dds$condition, levels = c("control","affected"))

dds <- DESeq(dds)
res <- results(dds)
#res <- results(dds, contrast=c("condition","control","treated"))
resultsNames(dds)

resLFC <- lfcShrink(dds, coef="group_F_vs_M", type="apeglm")

resOrdered <- res[order(res$pvalue),]
resOrdered

sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

plotMA(res, ylim=c(-2,2))
plot_MA <- plotMA(resLFC, ylim=c(-2,2))

legend("topright", # Position keyword
       legend = levels(res),
       col = c("blue", "grey"),
       pch = 19, # Must match the symbol used in plot()
       title = "Group")


#idx <- identify(res$baseMean, res$log2FoldChange)
#rownames(res)[idx]

plotCounts(dds, gene=which.min(res$padj), intgroup="group")


resOrdered <- as.data.frame(resOrdered) |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::mutate(
    gene_name = rowData(dds)$gene_name[
      match(gene_id, rownames(dds))
    ]
  ) |>
  dplyr::relocate(gene_name, .before = gene_id)

resSig <- subset(resOrdered, padj < 0.05)
resSig <- as.data.frame(resSig) |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::mutate(
    gene_name = rowData(dds)$gene_name[
      match(gene_id, rownames(dds))
    ]
  ) |>
  dplyr::relocate(gene_name, .before = gene_id)



write.csv(as.data.frame(resOrdered), 
          file="/g/data/xl04/eh8642/RNAseq_run/deSeq_results/F_M_rnaseq_results.csv")
write.csv(as.data.frame(resSig), 
          file="/g/data/xl04/eh8642/RNAseq_run/deSeq_results/F_M_rnaseq_results_p0.05.csv")


##################### Female comparison #####################

femaletxi <- tximport(
  femalequantsf,
  type = "salmon",
  txOut = FALSE,
  tx2gene = Femaletx2gene,
  countsFromAbundance = "no"  # <-- critical
)

femalecontrast$group <- factor(femalecontrast$group,levels=c("control","affected"))

Female_dds <- DESeqDataSetFromTximport(
  femaletxi,
  colData = femalecontrast,
  design = ~ group
)

gene_names_female <- Femaletx2gene |>
  dplyr::select(gene_id, gene_name) |>
  dplyr::distinct()

rowData(Female_dds)$gene_name <- gene_names_female$gene_name[
  match(rownames(Female_dds), gene_names_female$gene_id)
]

rowData(Female_dds)$gene_name <- ifelse(
  is.na(rowData(Female_dds)$gene_name),
  rownames(Female_dds),
  rowData(Female_dds)$gene_name
)

rowData(Female_dds)$gene_name <- make.unique(rowData(Female_dds)$gene_name)

smallestGroupSize <- 10
keep <- rowSums(counts(Female_dds) >= 10) >= smallestGroupSize
Female_dds <- Female_dds[keep,]


Female_dds$group <- factor(Female_dds$group, levels = c("control","affected"))
Female_dds <- DESeq(Female_dds)
res_female <- results(Female_dds)
res_female  <- results(Female_dds, contrast=c("group","control","affected"))
resultsNames(Female_dds)

resLFC_female <- lfcShrink(Female_dds, coef="group_affected_vs_control", type="apeglm")
resOrdered_female <- res_female[order(res_female$pvalue),]
resOrdered_female

sum(res_female$padj < 0.1, na.rm=TRUE)
femaleres05 <- results(Female_dds, alpha=0.05)

sum(femaleres05$padj < 0.05, na.rm=TRUE)

plotMA(res_female, ylim=c(-2,2))
plot_MA <- plotMA(resLFC_female,ylim=c(-2,2))

legend("topright", # Position keyword
       legend = levels(res),
       col = c("blue", "grey"),
       pch = 19, # Must match the symbol used in plot()
       title = "Group")


#female_idx <- identify(res_female$baseMean, res_female$log2FoldChange)
#rownames(res_female)[female_idx]

plotCounts(Female_dds, gene=which.min(res_female$padj), intgroup="group") 

resOrdered_female <- as.data.frame(resOrdered_female) |>
  unique() |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::mutate(
    gene_name = rowData(Female_dds)$gene_name[
      match(gene_id, rownames(Female_dds))
    ]
  ) |>
  dplyr::relocate(gene_name, .before = gene_id)

resSig_female <- subset(resOrdered_female, padj < 0.05)
resSig_female <- as.data.frame(resSig_female) |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::mutate(
    gene_name = rowData(Female_dds)$gene_name[
      match(gene_id, rownames(Female_dds))
    ]
  ) |>
  dplyr::relocate(gene_name, .before = gene_id)

plotDispEsts(Female_dds)

write.csv(as.data.frame(resOrdered_female), 
          file="/g/data/xl04/eh8642/RNAseq_run/deSeq_results/F_control_affected_rnaseq_results.csv")
write.csv(as.data.frame(resSig_female), 
          file="/g/data/xl04/eh8642/RNAseq_run/deSeq_results/F_control_affected_rnaseq_results_p0.05.csv")


##################### Male comparison #####################

maletxi <- tximport(
  malequantsf,
  type = "salmon",
  txOut = FALSE,
  tx2gene = Maletx2gene,
  countsFromAbundance = "no"  # <-- critical
)

malecontrast$group <- factor(malecontrast$group,levels=c("control","affected"))

male_dds <- DESeqDataSetFromTximport(
  maletxi,
  colData = malecontrast,
  design = ~ group
)

gene_names_male <- maletx2gene |>
  dplyr::select(gene_id, gene_name) |>
  dplyr::distinct()

rowData(male_dds)$gene_name <- gene_names_male$gene_name[
  match(rownames(male_dds), gene_names_male$gene_id)
]

rowData(male_dds)$gene_name <- ifelse(
  is.na(rowData(male_dds)$gene_name),
  rownames(male_dds),
  rowData(male_dds)$gene_name
)

rowData(male_dds)$gene_name <- make.unique(rowData(male_dds)$gene_name)

smallestGroupSize <- 10
keep <- rowSums(counts(male_dds) >= 10) >= smallestGroupSize
male_dds <- male_dds[keep,]


male_dds$group <- factor(male_dds$group, levels = c("control","affected"))
male_dds <- DESeq(male_dds)
res_male <- results(male_dds)
res_male  <- results(male_dds, contrast=c("group","control","affected"))
resultsNames(male_dds)

resLFC_male <- lfcShrink(male_dds, coef="group_affected_vs_control", type="apeglm")
resOrdered_male <- res_male[order(res_male$pvalue),]
resOrdered_male

sum(res_male$padj < 0.1, na.rm=TRUE)
maleres05 <- results(male_dds, alpha=0.05)

sum(maleres05$padj < 0.05, na.rm=TRUE)

plotMA(res_male, ylim=c(-2,2))
plot_MA <- plotMA(resLFC_male,ylim=c(-2,2))

legend("topright", # Position keyword
       legend = levels(res),
       col = c("blue", "grey"),
       pch = 19, # Must match the symbol used in plot()
       title = "Group")


#male_idx <- identify(res_male$baseMean, res_male$log2FoldChange)
#rownames(res_male)[male_idx]

plotCounts(male_dds, gene=which.min(res_male$padj), intgroup="group") 

resOrdered_male <- as.data.frame(resOrdered_male) |>
  unique() |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::mutate(
    gene_name = rowData(male_dds)$gene_name[
      match(gene_id, rownames(male_dds))
    ]
  ) |>
  dplyr::relocate(gene_name, .before = gene_id)

resSig_male <- subset(resOrdered_male, padj < 0.05)
resSig_male <- as.data.frame(resSig_male) |>
  unique() |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::mutate(
    gene_name = rowData(male_dds)$gene_name[
      match(gene_id, rownames(male_dds))
    ]
  ) |>
  dplyr::relocate(gene_name, .before = gene_id)

plotDispEsts(male_dds)

write.csv(as.data.frame(resOrdered_male), 
          file="/g/data/xl04/eh8642/RNAseq_run/deSeq_results/M_control_affected_rnaseq_results.csv")
write.csv(as.data.frame(resSig_male), 
          file="/g/data/xl04/eh8642/RNAseq_run/deSeq_results/M_control_affected_rnaseq_results_p0.05.csv")




##################### Multi-factor #####################

ddsMF <- dds
design(ddsMF) <- formula(~ condition + group)
ddsMF <- DESeq(ddsMF)
resMF <- results(ddsMF)


# this gives log2(n + 1)
ntd <- normTransform(ddsMF)
#colnames(ntd) <- contrastmerged$file
meanSdPlot(assay(ntd))
select <- order(rowMeans(counts(ddsMF,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","group")])

mat <- assay(vst(ddsMF))[select, ]
rownames(mat) <- rowData(dds)$gene_name[select]
# ensure uniqueness
rownames(mat) <- make.unique(rownames(mat))

pheatmap(mat, cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df, show_colnames = TRUE)



ntd <- normTransform(dds)
#colnames(ntd) <- contrastmerged$file
meanSdPlot(assay(ntd))
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,"group"])
colnames(df) <- c("group")
mat <- assay(vst(dds))[select, ]
rownames(mat) <- rowData(dds)$gene_name[select]
# ensure uniqueness
rownames(mat) <- make.unique(rownames(mat))

pheatmap(mat, cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df, show_colnames = TRUE)



















































