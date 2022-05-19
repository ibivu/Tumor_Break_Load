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
outputDir = "../../results/survival_analysis/pathway/"
dir.create(outputDir)

barcode <- function(x){
  # Convert TCGA identifier to TCGA barcode (TCGA case identifier)
  return(paste(strsplit(as.character(x), "-")[[1]][1:3], collapse="-"))
}

# read data, clinical data from liu et.al 2018 cell and genomic metrics
clinical_data = read_excel("../../data/mmc1.xlsx", sheet = "TCGA-CDR")

genome_metrics = read.csv("../../results/sample_genome_metrics_tbl_class.csv", 
                          header = T, row.names = 1)

# read curated pathway alteration info (Sanchez-Vega,F. et al. (2018) Oncogenic Signaling Pathways in The Cancer Genome Atlas. Cell, 173, 321-337.e10.)
pathway_info <- read_excel("../../data/1-s2.0-S0092867418303593-mmc4.xlsx", sheet = "Pathway level", na = c("NA"))

pathway_info <- pathway_info %>% mutate(barcode = sapply(SAMPLE_BARCODE, FUN = barcode))
pathway_info$barcode <- as.character(pathway_info$barcode)

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
survival_data$TP53_pathway <- pathway_info[match(as.character(survival_data$bcr_patient_barcode), 
                               as.character(pathway_info$barcode)),]$TP53
survival_data$RAS_RTK_pathway <- pathway_info[match(as.character(survival_data$bcr_patient_barcode), 
                                                    as.character(pathway_info$barcode)),]$`RTK RAS`

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

survival_local = survival_data[survival_data$pathological_stage %in% 
                                 c("Stage II", "Stage III"),]

survival_local$TP53_pathway = revalue(as.character(survival_local$TP53_pathway), c("0" = "WT", "1" = "Altered"))

sfit = survfit(Surv(DFI.time, DFI) ~ TP53_pathway, data = survival_local[!is.na(survival_local$TP53_pathway),])
cox_TBL <- coxph(Surv(DFI.time, DFI) ~ factor(TP53_pathway, levels = c("WT", "Altered")), 
                 data = survival_local[!is.na(survival_local$TP53_pathway),])
s_tbl <- summary(cox_TBL)
p_tbl = surv_pvalue(sfit)
coldata <- survival_local[!is.na(survival_local$DFI) & !is.na(survival_local$TP53_pathway),]$TP53_pathway
my_xlab <- paste(levels(factor(coldata)), " (N=", table(coldata), ")", sep="")
ggsurv <- ggsurvplot(sfit, conf.int = F, pval = F, risk.table = F, 
                     palette = c("dodgerblue2", "orchid2"),
                     title = "TCGA MSS CRC Stage 2 & 3",
                     risk.table.height = .15)
ggsurv$plot <- ggsurv$plot + theme_classic(base_size = 21) + theme(axis.text = element_text(size=18),
                                   plot.title = element_text(hjust = 0.5), legend.position = "top",
                                   legend.text = element_text(size=16)) + 
  scale_color_discrete(name = "TP53 pathway", labels = my_xlab) +
  ggplot2::annotate(
    "text",
    x = 4000, y = .3,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl$pval, 3), "\nHR = ", round(exp(cox_TBL$coefficients), 3)," (",  
                  round(s_tbl$conf.int[3], 3), " - ", round(s_tbl$conf.int[4], 3), ")", sep=""),
    size = 8
  ) + xlab("Time (Days)") + ylab("Disease Free Survival (DFS)")
pdf(paste(outputDir, "survival_analysis_DFI_MSS_23_TP53_pathway.pdf", sep=""))
ggsurv
dev.off()
ggsurv

survival_local$RAS_RTK_pathway = revalue(as.character(survival_local$RAS_RTK_pathway), c("0" = "WT", "1" = "Altered"))

sfit = survfit(Surv(DFI.time, DFI) ~ RAS_RTK_pathway, data = survival_local[!is.na(survival_local$RAS_RTK_pathway),])
cox_TBL <- coxph(Surv(DFI.time, DFI) ~ factor(RAS_RTK_pathway, levels = c("WT", "Altered")), 
                 data = survival_local[!is.na(survival_local$RAS_RTK_pathway),])
s_tbl <- summary(cox_TBL)
p_tbl = surv_pvalue(sfit)
coldata <- survival_local[!is.na(survival_local$DFI) & !is.na(survival_local$RAS_RTK_pathway),]$RAS_RTK_pathway
my_xlab <- paste(levels(factor(coldata)), " (N=", table(coldata), ")", sep="")
ggsurv <- ggsurvplot(sfit, conf.int = F, pval = F, risk.table = F, 
                     palette = c("dodgerblue2", "orchid2"),
                     title = "TCGA MSS CRC Stage 2 & 3",
                     risk.table.height = .15)
ggsurv$plot <- ggsurv$plot + theme_classic(base_size = 21) + theme(axis.text = element_text(size=18), 
                                   plot.title = element_text(hjust = 0.5), legend.position = "top", 
                                   legend.text = element_text(size=16)) + 
  scale_color_discrete(name = "RAS RTK pathway", labels = my_xlab) +
  ggplot2::annotate(
    "text",
    x = 4000, y = .3,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl$pval, 3), "\nHR = ", round(exp(cox_TBL$coefficients), 3)," (",  
                  round(s_tbl$conf.int[3], 3), " - ", round(s_tbl$conf.int[4], 3), ")", sep=""),
    size = 8
  ) + xlab("Time (Days)") + ylab("Disease Free Survival (DFS)")
pdf(paste(outputDir, "survival_analysis_DFI_MSS_23_RAS_RTK_pathway.pdf", sep=""))
ggsurv
dev.off()
ggsurv