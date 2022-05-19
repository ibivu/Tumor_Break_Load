# compare the DFS for CMS 2 and 4

# Libraries
library(survival)
library(survminer)
library(plyr)
library(TCGAbiolinks)
library(stringr)
library(data.table)
library(tidyverse)
library(ggplot2)
library(lmvar)
library(caret)
library(ggpmisc)
library(MASS)
library(readxl)
library(pROC)

# result output directories
outputDir = "../../results/survival_analysis/CMS/"
dir.create(outputDir)

# read data, clinical data from liu et.al 2018 cell and genomic metrics
clinical_data = read_excel("../../data/mmc1.xlsx", sheet = "TCGA-CDR")

genome_metrics = read.csv("../../results/sample_genome_metrics_tbl_class.csv", 
                          header = T, row.names = 1)

# read CMS classifciaton data
cms <- read.table("../../data/CMS/cms_labels_public_all.txt", header = T)
cms <- cms[cms$dataset == "tcga",] # select only TCGA samples

# process survival data
survival_data = subset(clinical_data[which(clinical_data$bcr_patient_barcode %in% unique(genome_metrics$barcode)),], 
                       select = c(bcr_patient_barcode, OS.time, OS, DSS, DSS.time, 
                                  DFI, DFI.time, PFI, PFI.time))

survival_data = survival_data[!duplicated(survival_data$bcr_patient_barcode),]
survival_data$OS.time <- as.numeric(as.character(survival_data$OS.time))
survival_data$DSS.time <- as.numeric(as.character(survival_data$DSS.time))
survival_data$DFI.time <- as.numeric(as.character(survival_data$DFI.time))
survival_data$PFI.time <- as.numeric(as.character(survival_data$PFI.time))
survival_data[survival_data == "[Not Available]"] <- NA
survival_data$OS <- as.numeric(survival_data$OS)
survival_data$DSS <- as.numeric(survival_data$DSS)
survival_data$DFI <- as.numeric(survival_data$DFI)
survival_data$PFI <- as.numeric(survival_data$PFI)

# add cms classification to TBL data frame
survival_data$cms <- cms[match(as.character(survival_data$bcr_patient_barcode), 
                     as.character(cms$sample)), "CMS_final_network_plus_RFclassifier_in_nonconsensus_samples"]

# annotate with staging
query = GDCquery(project = c("TCGA-COAD", "TCGA-READ"),
                 data.category = "Clinical",
                 file.type = "xml",
                 barcode = survival_data$bcr_patient_barcode)
GDCdownload(query)
clinical = GDCprepare_clinic(query, clinical.info = "stage_event")

survival_data$pathological_stage = clinical[match(survival_data$bcr_patient_barcode, 
                                                 clinical$bcr_patient_barcode), 
                                           "pathologic_stage"]

survival_data$pathological_stage = factor(str_replace(survival_data$pathological_stage, 
                                                     pattern = "[ABC]", replacement = ""))

survival_filtered = survival_data[survival_data$pathological_stage %in% 
                               c("Stage I", "Stage II", "Stage III", "Stage IV"),]

survival_cms = survival_filtered[survival_filtered$pathological_stage %in% c("Stage II", "Stage III") & 
                                   survival_filtered$cms %in% c("CMS2", "CMS4"),]

sfit = survfit(Surv(DFI.time, DFI) ~ cms, data = survival_cms)
cox_TBL <- coxph(Surv(DFI.time, DFI) ~ factor(cms, levels = c("CMS4", "CMS2")), 
                 data = survival_cms)
s_tbl <- summary(cox_TBL)
p_tbl = surv_pvalue(sfit)
coldata <- survival_cms[!is.na(survival_cms$DFI),]$cms
my_xlab <- paste(levels(factor(coldata)), " (N=", table(coldata), ")", sep="")
ggsurv <- ggsurvplot(sfit, conf.int = F, pval = F, risk.table = F, 
                     palette = c("dodgerblue2", "orchid2"),
                     title = "TCGA MSS CRC Stage 2 & 3",
                     risk.table.height = .15)
ggsurv$plot <- ggsurv$plot  + 
  scale_color_discrete(name = "", labels = my_xlab) +
  ggplot2::annotate(
    "text",
    x = 4000, y = .3,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl$pval, 3), "\nHR = ", round(exp(cox_TBL$coefficients), 3)," (",  
                  round(s_tbl$conf.int[3], 3), " - ", round(s_tbl$conf.int[4], 3), ")", sep=""),
    size = 9
  ) + xlab("Time (Days)") + ylab("Disease Free Survival (DFS)") + theme_classic(base_size = 21) + 
  theme(axis.text = element_text(size=18), plot.title = element_text(hjust = 0.5), legend.position = "top")
pdf(paste(outputDir, "survival_analysis_DFI_MSS_23_CMS_2_4.pdf", sep=""))
ggsurv
dev.off()
ggsurv

sfit = survfit(Surv(DFI.time, DFI) ~ cms, data = survival_filtered[survival_filtered$pathological_stage %in% c("Stage II", "Stage III"),])
cox_TBL <- coxph(Surv(DFI.time, DFI) ~ factor(cms, levels = c("CMS1", "CMS2", "CMS3", "CMS4")), 
                 data = survival_filtered[survival_filtered$pathological_stage %in% c("Stage II", "Stage III"),])
s_tbl <- summary(cox_TBL)
p_tbl = surv_pvalue(sfit)
coldata <- survival_filtered[!is.na(survival_filtered$DFI) & survival_filtered$pathological_stage %in% c("Stage II", "Stage III"),]$cms
my_xlab <- paste(levels(factor(coldata)), " (N=", table(coldata), ")", sep="")
ggsurv <- ggsurvplot(sfit, conf.int = F, pval = F, risk.table = F, 
                     palette = c("dodgerblue2", "orchid2"),
                     title = "TCGA MSS CRC Stage 2 & 3",
                     risk.table.height = .15)
ggsurv$plot <- ggsurv$plot + theme_classic(base_size = 21) + theme(axis.text = element_text(size=18),
                                   plot.title = element_text(hjust = 0.5), legend.position = "top", legend.text = element_text(size=10)) + 
  scale_color_discrete(name = "", labels = my_xlab) +
  ggplot2::annotate(
    "text",
    x = 4000, y = .3,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl$pval, 3), sep=""),
    size = 9
  ) + xlab("Time (Days)") + ylab("Disease Free Survival (DFS)")
pdf(paste(outputDir, "survival_analysis_DFI_MSS_23_CMS_all.pdf", sep=""))
ggsurv
dev.off()
ggsurv