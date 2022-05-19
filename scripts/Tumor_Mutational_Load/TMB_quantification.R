library(maftools)
library(dplyr)
library(ggplot2)

# read masked copy number variation
outputdirTumor <- "../../results/SNP6/masked_copy_number_segment/primary_tumor/"

mut_data_coad <- read.table(file = "../../data/Mutation_data/MuTect2/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf", header = T, sep = "\t", quote = "") 
mut_data_read <- read.table(file = "../../data/Mutation_data/MuTect2/TCGA.READ.mutect.faa5f62a-2731-4867-a264-0e85b7074e87.DR-10.0.somatic.maf", header = T, sep = "\t", quote = "") 

mutation_data <- rbind(mut_data_coad, mut_data_read)

coadread <- read.maf(mutation_data)

# calculate TMB
mut_tmb <- tmb(coadread, logScale = T, captureSize = 100)

sample_id <- function(x){
  return(paste(strsplit(as.character(x), "-")[[1]][1:4], collapse="-"))
}

mut_tmb <- mut_tmb %>% mutate(sample_ids = sapply(Tumor_Sample_Barcode, FUN = sample_id))

# READ TBL data
Tumor_Break_Load <- read.csv(paste(outputdirTumor, "tumor_break_load.csv", sep=""), header = T, row.names = 1)

Tumor_Break_Load <- Tumor_Break_Load %>% mutate(sample_ids = sapply(sample, FUN = sample_id))

# add TMB
Tumor_Break_Load$TMB <- mut_tmb[match(Tumor_Break_Load$sample_ids, mut_tmb$sample_ids), ]$total

Tumor_Break_Load <- Tumor_Break_Load[!duplicated(Tumor_Break_Load$sample_ids),]

Tumor_Break_Load <- Tumor_Break_Load[!is.na(Tumor_Break_Load$TMB),]

#ggplot(Tumor_Break_Load, aes(x = count, y = TMB, color = msi_status)) + geom_point()
#ggplot(Tumor_Break_Load, aes(x = aneuploidy_score, y = TMB, color = msi_status)) + geom_point()
#ggplot(Tumor_Break_Load, aes(x = fraction_abr, y = TMB, color = msi_status)) + geom_point()

write.csv(Tumor_Break_Load, "../../results/sample_genome_metrics.csv")
