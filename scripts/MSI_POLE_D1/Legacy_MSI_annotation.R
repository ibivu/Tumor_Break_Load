library(TCGAbiolinks)
source("../Break_calling/cnd_functions.R")
library(dplyr)

outputdir <- "../../results/SNP6/masked_copy_number_segment/primary_tumor/"
Tumor_Break_Load <- read.csv(paste(outputdir, "tumor_break_load.csv", sep=""), header = T, row.names = 1)

Tumor_Break_Load <- Tumor_Break_Load %>% mutate(barcode = sapply(sample, FUN = barcode))

# Query Legacy GDC portal to obtain MSI annotation
query <- GDCquery(project = c("TCGA-COAD", "TCGA-READ"),
                  data.category = "Other",
                  legacy = T,
                  access = "open",
                  data.type = "Auxiliary test",
                  barcode = Tumor_Break_Load$barcode)

GDCdownload(query)

msi_results <- GDCprepare_clinic(query, "msi")

Tumor_Break_Load$msi_status <- msi_results[match(Tumor_Break_Load$barcode, 
                                         msi_results$bcr_patient_barcode), 
                                   "mononucleotide_and_dinucleotide_marker_panel_analysis_status"]

write.csv(Tumor_Break_Load, paste(outputdir, "tumor_break_load.csv", sep=""))
