# libraries
library(dplyr)
library(ggfortify)
#library(umapr)
library(edgeR)
library(data.table)
library(stringr)

# read count data
count_data <- read.csv("../../data/RNA_Seq/genes/Counts/primary_tumor/tcga_coadread_htseq_count_aliquot_tumor.csv", header = T, row.names = 1)

# read genome metrics
metrics <- read.csv("../../results/sample_genome_metrics.csv", header = T, row.names = 1)

sample_id <- function(x){
  return(paste(strsplit(as.character(x), "-")[[1]][1:4], collapse="-"))
}

# match samples genome metrics and RNA seq count data
count_data <- count_data %>% mutate(sample_ids = sapply(aliquot_ids, FUN = sample_id))

overlaps_samples <- intersect(count_data$sample_ids, metrics$sample_ids)

metrics <- metrics[match(overlaps_samples, metrics$sample_ids),]

count_data <- count_data[match(overlaps_samples, count_data$sample_ids),]

# write out metrics data with samples with all data
write.csv(metrics, "../../results/sample_genome_metrics_matched.csv")

count_data <- count_data[!duplicated(count_data$sample_ids),]

# change rownames to aliquot ids
rownames(count_data) <- count_data$sample_ids

count_meta_data <- subset(count_data, select = c(aliquot_ids, sample_ids))

count_data <- subset(count_data, select = -c(aliquot_ids, sample_ids, 
                                             X__no_feature, X__ambiguous, 
                                             X__too_low_aQual, X__not_aligned, 
                                             X__alignment_not_unique))

count_data <- count_data[ , which(apply(count_data, 2, var) != 0)]

#pc_tumor <- prcomp(count_data, center = T, scale. = T)
#autoplot(pc_tumor)

# embedding <- umap(count_data)
# embedding %>% 
#   ggplot(aes(UMAP1, UMAP2)) + geom_point()

# reorient data for downstream analysis
countData <- transpose(count_data)
colnames(countData) <- rownames(count_data)
rownames(countData) <- colnames(count_data)

# strip version nr from ensemble id
rownames(countData) <- str_replace(rownames(countData), pattern = ".[0-9]+$",
                                   replacement = "")

# normalize count data
dgeObj <- DGEList(countData)
dgeObj <- calcNormFactors(dgeObj, method = "TMM")

# transform count data
tmm <- cpm(countData)
tmm <- as.data.frame(tmm)

pc_data <- transpose(tmm)
rownames(pc_data) <- colnames(tmm)
colnames(pc_data) <- rownames(tmm)

#pca_tumor <- prcomp(pc_data, center = T, scale. = T)
#autoplot(pca_tumor)

# embedding <- umap(pc_data)
# embedding %>% 
#   ggplot(aes(UMAP1, UMAP2)) + geom_point()

count_normalized_data <- transpose(tmm)
rownames(count_normalized_data) <- colnames(tmm)
colnames(count_normalized_data) <- rownames(tmm)

dir.create("../../results/RNA_Seq/genes/Counts", recursive = T)

write.csv(count_normalized_data, file = "../../results/RNA_Seq/genes/Counts/normalized_RNA_Seq.csv")
