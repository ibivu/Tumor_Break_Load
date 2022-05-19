library(dplyr)
library(tibble)
library(plyr)
library(tidyr)

# Read mutation data
mut_data_coad <- read.table(file = "../../data/Mutation_data/MuTect2/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf", header = T, sep = "\t", quote = "") 
mut_data_read <- read.table(file = "../../data/Mutation_data/MuTect2/TCGA.READ.mutect.faa5f62a-2731-4867-a264-0e85b7074e87.DR-10.0.somatic.maf", header = T, sep = "\t", quote = "") 

mutation_data <- rbind(mut_data_coad, mut_data_read)

genes_of_interest <- c("POLE", "POLD1")

mutation_goi <- mutation_data[mutation_data$Hugo_Symbol %in% genes_of_interest,]

non_silent_mutation_goi <- mutation_goi[mutation_goi$Variant_Classification != "Silent",]

mutation_profile <- mutation_goi %>% group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% 
  dplyr::summarise(n = n()) %>% spread(Hugo_Symbol, n) %>% 
  replace(is.na(.), 0) %>% column_to_rownames(var = "Tumor_Sample_Barcode")

# fill mutation profile with samples with no mutations for genes of interest
new_samples <- setdiff(mutation_data$Tumor_Sample_Barcode, rownames(mutation_profile))

zero_df <- matrix(0, nrow = length(new_samples), ncol = length(colnames(mutation_profile)))
zero_df <- data.frame(zero_df)
colnames(zero_df) <- colnames(mutation_profile)
rownames(zero_df) <- new_samples

mutation_profile <- rbind(mutation_profile, zero_df)

colnames(mutation_profile) <- paste(colnames(mutation_profile), "mut", sep="_")

mutation_profile <- data.frame(mutation_profile > 0)

# read genomic metrics
genome_metrics <- read.csv("../../results/sample_genome_metrics.csv", header = T, row.names = 1)

sample_id <- function(x){
  return(paste(strsplit(as.character(x), "-")[[1]][1:4], collapse="-"))
}

mutation_profile <- mutation_profile %>% mutate(sample_ids = sapply(rownames(mutation_profile), FUN = sample_id))

genome_metrics <- cbind(genome_metrics, mutation_profile[match(genome_metrics$sample_ids, mutation_profile$sample_ids),])

write.csv(genome_metrics, "../../results/sample_genome_metrics.csv")
