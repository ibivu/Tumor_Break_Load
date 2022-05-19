# libraries
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

outputDir = "../../results/survival_analysis/"
dir.create(outputDir)

# read data, clinical data from liu et.al 2018 cell
clinical_data = read_excel("../../data/mmc1.xlsx", sheet = "TCGA-CDR")
genome_metrics = read.csv("../../results/sample_genome_metrics_tbl_class.csv", header = T, row.names = 1)

sample_id <- function(x){
  return(paste(strsplit(as.character(x), "-")[[1]][1:4], collapse="-"))
}

length(intersect(genome_metrics$barcode, clinical_data$bcr_patient_barcode))

# process survival data
survival_data = subset(clinical_data[which(clinical_data$bcr_patient_barcode %in% unique(genome_metrics$barcode)),], 
                       select = c(bcr_patient_barcode, OS.time, OS, DSS, DSS.time, DFI, DFI.time, PFI, PFI.time))

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

# plot overall survival of cohort
sfit = survfit(Surv(OS.time, OS)~1, data = survival_data)
plt = ggsurvplot(sfit)
plot.new()
print(plt, newpage = F)

# annotate TBL variables to the survival data
survival_data$TBL = genome_metrics[match(survival_data$bcr_patient_barcode, genome_metrics$barcode), "count"]
survival_data$TBL_status = genome_metrics[match(survival_data$bcr_patient_barcode, genome_metrics$barcode), "tbl_predicted"] %>% revalue(c("abr"="High", "wt"="Low"))
survival_data$TBL_status = factor(survival_data$TBL_status, levels = c("Low", "High"))
survival_data$TBL_predefined = genome_metrics[match(survival_data$bcr_patient_barcode, genome_metrics$barcode), "TBL_state"] %>% revalue(c("H" = "High", "I" = "Medium", "L" = "Low")) 
survival_data$TBL_predefined = factor(survival_data$TBL_predefined, levels = c("Low", "Medium", "High"))

# plot DFS of the cohort
fit <- survfit(Surv(DFI.time, DFI) ~ 1, data = survival_data)
ggsurvplot(
  fit,
  risk.table = T,
  pval = T,
  conf.int = F
)

# predicted TBL state (DFI)
sfit_TBL_status = survfit(Surv(DFI.time, DFI) ~ TBL_status, data = survival_data)
plt = ggsurvplot(sfit_TBL_status)
plot.new()
print(plt, newpage = F)

# predefined TBL states (DFI)
surf_model = survival_data[survival_data$TBL_predefined != "Medium",]
surf_model$TBL_predefined = factor(surf_model$TBL_predefined)
sfit_TBL_predefined = survfit(Surv(DFI.time, DFI) ~ TBL_predefined, data = surf_model)
plt = ggsurvplot(sfit_TBL_predefined)
plot.new()
print(plt, newpage = F)

cox_TBL <- coxph(Surv(DFI.time, DFI) ~ TBL_status, data = survival_data)
s_tbl <- summary(cox_TBL)
p_tbl = surv_pvalue(sfit_TBL_status)

cox_TBL_pre <- coxph(Surv(DFI.time, DFI) ~ TBL_predefined, data = surf_model)
s_tb_pre <- summary(cox_TBL_pre)
p_tbl_pre = surv_pvalue(sfit_TBL_predefined)

ggsurv <- ggsurvplot(sfit_TBL_status, conf.int = F, pval = F, risk.table = T,
                     legend.labs = c("LOW", "HIGH"), legend.title = "TBL predicted",
                     palette = c("dodgerblue2", "orchid2"),
                     title = "Kaplan-Meier Curve for TBL MSS colorectal cancer Stage 1-4",
                     risk.table.height = .15)
ggsurv$plot <- ggsurv$plot +
  ggplot2::annotate(
    "text",
    x = Inf, y = Inf,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl$pval, 3), "\nHR = ", round(exp(cox_TBL$coefficients), 3)," (",  round(s_tbl$conf.int[3], 3), " - ", round(s_tbl$conf.int[4], 3), ")", sep=""),
    size = 5
  ) + xlab("Time (Months)") + ylab("Disease Free Survival (DFS)")
pdf(paste(outputDir, "survival_analysis_TBL_MSS_all_DFS.pdf", sep=""))
ggsurv
dev.off()

ggsurv <- ggsurvplot(sfit_TBL_predefined, conf.int = F, pval = F, risk.table = T,
                     legend.labs = c("LOW", "HIGH"), legend.title = "TBL predefined",
                     palette = c("dodgerblue2", "orchid2"),
                     title = "Kaplan-Meier Curve for TBL MSS colorectal cancer Stage 1-4",
                     risk.table.height = .15)
ggsurv$plot <- ggsurv$plot +
  ggplot2::annotate(
    "text",
    x = Inf, y = Inf,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl_pre$pval, 3), "\nHR = ", round(exp(cox_TBL_pre$coefficients), 3)," (",  round(s_tb_pre$conf.int[3], 3), " - ", round(s_tb_pre$conf.int[4], 3), ")", sep=""),
    size = 5
  ) + xlab("Time (Months)") + ylab("Disease Free Survival (DFS)")
pdf(paste(outputDir, "survival_analysis_TBL_predefined_MSS_all_DFS.pdf", sep=""))
ggsurv
dev.off()

pdf(paste(outputDir, "survival_analysis_TBL_combined_MSS_all_DFS.pdf", sep=""))
ggsurvplot_combine(list(predicted = sfit_TBL_status, predefined = sfit_TBL_predefined),
                   pval = T, 
                   data=list(predicted = survival_data, predefined = surf_model),
                   legend.labs = c("Low predicted", "High predicted", "Low predefined", "High predefined"),
                   palette = c(alpha("dodgerblue2", 0.8), alpha("orchid2", 0.8), alpha("green", 0.1), alpha("red", 0.1)),
                   legend = "bottom")
plot.new()
print(plt, newpage = F)
dev.off()
plt

sfit_TBL_status = survfit(Surv(DFI.time, DFI) ~ TBL_status, data = 
                              survival_data)

plt = ggsurvplot(sfit_TBL_status, conf.int = F, pval = T, risk.table = T)
plot.new()
print(plt, newpage = F)

cox_TBL <- coxph(Surv(DFI.time, DFI) ~ TBL_status, data = survival_data)
s_tbl <- summary(cox_TBL)
p_tbl = surv_pvalue(sfit_TBL_status)

coldata <- survival_data[!is.na(survival_data$DFI),]$TBL_status
my_xlab <- paste(levels(factor(coldata)), " (N=", table(coldata), ")", sep="")
ggsurv <- ggsurvplot(sfit_TBL_status, conf.int = F, pval = F, risk.table = F, 
                     legend = c(0.75, 0.4),
                     palette = c("dodgerblue2", "orchid2"),
                     title = "TCGA MSS CRC Stage 1-4")
ggsurv$plot <- ggsurv$plot + theme(legend.text = element_text(size=12), legend.title = element_text(size=14),
                                   plot.title = element_text(hjust = 0.5)) + 
  scale_color_discrete(name = "TBL expressional phenotype", labels = my_xlab, direction = -1) + 
  ggplot2::annotate(
    "text",
    x = Inf, y = 0.25,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl$pval, 3), "\nHR = ", 
                  round(exp(cox_TBL$coefficients), 3)," (",  
                  round(s_tbl$conf.int[3], 3), " - ", 
                  round(s_tbl$conf.int[4], 3), ")", sep=""),
    size = 5
  ) + xlab("Time (Days)") + ylab("Disease Free Survival (DFS)")

pdf(paste(outputDir, "survival_analysis_DFI_MSS_all.pdf", sep=""), width = 8, height = 8)
ggsurv
dev.off()
ggsurv


# Annotate with stage information
query = GDCquery(project = c("TCGA-COAD", "TCGA-READ"),
                 data.category = "Clinical",
                 file.type = "xml",
                 barcode = genome_metrics$barcode)
GDCdownload(query)
clinical = GDCprepare_clinic(query, clinical.info = "stage_event")

genome_metrics$pathological_stage = clinical[match(genome_metrics$barcode, clinical$bcr_patient_barcode), "pathologic_stage"]

genome_metrics$pathological_stage = factor(str_replace(genome_metrics$pathological_stage, pattern = "[ABC]", replacement = ""))

genomic_filtered = genome_metrics[genome_metrics$pathological_stage %in% c("Stage I", "Stage II", "Stage III", "Stage IV"),]

genomic_filtered = genomic_filtered[!duplicated(genomic_filtered$sample_ids),]

rownames(genomic_filtered) = genomic_filtered$sample_ids

survival_filtered = survival_data[survival_data$bcr_patient_barcode %in% genomic_filtered$barcode,]
survival_filtered$stage <- genomic_filtered[match(survival_filtered$bcr_patient_barcode, genomic_filtered$barcode), "pathological_stage"]

sfit_TBL_status = survfit(Surv(DFI.time, DFI) ~ TBL_status, data = 
                            survival_filtered[survival_filtered$stage %in% c("Stage I", "Stage II", "Stage III", "Stage IV"),])
summary(sfit_TBL_status)


plt = ggsurvplot(sfit_TBL_status, conf.int = F, pval = T)
plot.new()
print(plt, newpage = F)

cox_TBL <- coxph(Surv(DFI.time, DFI) ~ TBL_status, data = 
                   survival_filtered[survival_filtered$stage %in% c("Stage I", "Stage II", "Stage III", "Stage IV"),])
s_tbl <- summary(cox_TBL)
p_tbl = surv_pvalue(sfit_TBL_status)

ggsurv <- ggsurvplot(sfit_TBL_status, conf.int = F, pval = F, risk.table = T,
                     legend.labs = c("LOW", "HIGH"), legend.title = "TBL predicted",
                     palette = c("dodgerblue2", "orchid2"),
                     title = "Kaplan-Meier Curve for TBL MSS colorectal cancer Stage 1-4",
                     risk.table.height = .15)
ggsurv$plot <- ggsurv$plot +
  ggplot2::annotate(
    "text",
    x = Inf, y = 0.25,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl$pval, 3), "\nHR = ", round(exp(cox_TBL$coefficients), 3)," (",  round(s_tbl$conf.int[3], 3), " - ", round(s_tbl$conf.int[4], 3), ")", sep=""),
    size = 5
  ) + xlab("Time (Days)") + ylab("Disease Free Interval (DFI)")
pdf(paste(outputDir, "survival_analysis_DFS_MSS_all.pdf", sep=""), width = 15, height = 15)
ggsurv
dev.off()
ggsurv

sfit_TBL_status = survfit(Surv(OS.time, OS) ~ TBL_status, data = 
                            survival_filtered[survival_filtered$stage %in% c("Stage I", "Stage II", "Stage III", "Stage IV"),])
summary(sfit_TBL_status)

plt = ggsurvplot(sfit_TBL_status, conf.int = F, pval = T)
plot.new()
print(plt, newpage = F)

cox_TBL <- coxph(Surv(OS.time, OS) ~ TBL_status, data = 
                   survival_filtered[survival_filtered$stage %in% c("Stage I", "Stage II", "Stage III", "Stage IV"),])
s_tbl <- summary(cox_TBL)
p_tbl = surv_pvalue(sfit_TBL_status)

ggsurv <- ggsurvplot(sfit_TBL_status, conf.int = F, pval = F, risk.table = T,
                     legend.labs = c("LOW", "HIGH"), legend.title = "TBL predicted",
                     palette = c("dodgerblue2", "orchid2"),
                     title = "Kaplan-Meier Curve for TBL MSS colorectal cancer Stage 1-4",
                     risk.table.height = .15)
ggsurv$plot <- ggsurv$plot +
  ggplot2::annotate(
    "text",
    x = Inf, y = 0.25,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl$pval, 3), "\nHR = ", round(exp(cox_TBL$coefficients), 3)," (",  round(s_tbl$conf.int[3], 3), " - ", round(s_tbl$conf.int[4], 3), ")", sep=""),
    size = 5
  ) + xlab("Time (Days)") + ylab("Overall Survival (OS)")
pdf(paste(outputDir, "survival_analysis_OS_MSS_all.pdf", sep=""), width = 15, height = 15)
ggsurv
dev.off()
ggsurv

sfit_TBL_status23 = survfit(Surv(DFI.time, DFI) ~ TBL_status, data = 
                             survival_filtered[survival_filtered$stage %in% c("Stage II", "Stage III"),])

plt = ggsurvplot(sfit_TBL_status23, conf.int = F, pval = T, risk.table = T)
plot.new()
print(plt, newpage = F)

cox_TBL <- coxph(Surv(DFI.time, DFI) ~ TBL_status, data = 
                   survival_filtered[survival_filtered$stage %in% c("Stage II", "Stage III"),])
s_tbl <- summary(cox_TBL)
p_tbl = surv_pvalue(sfit_TBL_status23)

coldata <- survival_filtered[survival_filtered$stage %in% c("Stage II", "Stage III") & !is.na(survival_filtered$DFI),]$TBL_status
my_xlab <- paste(levels(factor(coldata)), " (N=", table(coldata), ")", sep="")
ggsurv <- ggsurvplot(sfit_TBL_status23, conf.int = F, pval = F, risk.table = F, legend = c(0.8, 0.4),
                     palette = c("dodgerblue2", "orchid2"),
                     title = "TCGA MSS CRC Stage 2 and 3")
ggsurv$plot <- ggsurv$plot + theme_classic(base_size = 21) + theme(axis.text = element_text(size=18),
                                   plot.title = element_text(hjust = 0.5), legend.position = "top",
                                   legend.text = element_text(size = 14)) + 
  scale_color_discrete(name = "TBL expressional phenotype", labels = my_xlab, direction = -1) + 
  ggplot2::annotate(
    "text",
    x = Inf, y = 0.25,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl$pval, 3), "\nHR = ", round(exp(cox_TBL$coefficients), 3)," (",  round(s_tbl$conf.int[3], 3), " - ", round(s_tbl$conf.int[4], 3), ")", sep=""),
    size = 8
  ) + xlab("Time (Days)") + ylab("Disease Free Interval (DFI)")

pdf(paste(outputDir, "survival_analysis_DFS_MSS_23.pdf", sep=""), width = 8, height = 8)
ggsurv
dev.off()
ggsurv

sfit_TBL_status2 = survfit(Surv(DFI.time, DFI) ~ TBL_status, data = 
                            survival_filtered[survival_filtered$stage %in% c("Stage II"),])

plt = ggsurvplot(sfit_TBL_status2, conf.int = F, pval = T)
plot.new()
print(plt, newpage = F)

cox_TBL <- coxph(Surv(DFI.time, DFI) ~ TBL_status, data = 
                   survival_filtered[survival_filtered$stage %in% c("Stage II"),])
s_tbl <- summary(cox_TBL)
p_tbl = surv_pvalue(sfit_TBL_status2)

ggsurv <- ggsurvplot(sfit_TBL_status2, conf.int = F, pval = F, risk.table = T,
                     legend.labs = c("LOW", "HIGH"), legend.title = "TBL predicted",
                     palette = c("dodgerblue2", "orchid2"),
                     title = "Kaplan-Meier Curve for TBL MSS colorectal cancer Stage 2",
                     risk.table.height = .15)
ggsurv$plot <- ggsurv$plot +
  ggplot2::annotate(
    "text",
    x = Inf, y = 0.25,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl$pval, 3), "\nHR = ", round(exp(cox_TBL$coefficients), 3)," (",  round(s_tbl$conf.int[3], 3), " - ", round(s_tbl$conf.int[4], 3), ")", sep=""),
    size = 5
  ) + xlab("Time (Days)") + ylab("Disease Free Interval (DFI)")

pdf(paste(outputDir, "survival_analysis_DFS_MSS_2.pdf", sep=""), width = 15, height = 15)
ggsurv
dev.off()

sfit_TBL_status3 = survfit(Surv(DFI.time, DFI) ~ TBL_status, data = 
                             survival_filtered[survival_filtered$stage %in% c("Stage III"),])

plt = ggsurvplot(sfit_TBL_status3, conf.int = F, pval = T)
plot.new()
print(plt, newpage = F)

cox_TBL <- coxph(Surv(DFI.time, DFI) ~ TBL_status, data = 
                   survival_filtered[survival_filtered$stage %in% c("Stage III"),])
s_tbl <- summary(cox_TBL)
p_tbl = surv_pvalue(sfit_TBL_status3)

ggsurv <- ggsurvplot(sfit_TBL_status3, conf.int = F, pval = F, risk.table = T,
                     legend.labs = c("LOW", "HIGH"), legend.title = "TBL predicted",
                     palette = c("dodgerblue2", "orchid2"),
                     title = "Kaplan-Meier Curve for TBL MSS colorectal cancer Stage 3",
                     risk.table.height = .15)
ggsurv$plot <- ggsurv$plot +
  ggplot2::annotate(
    "text",
    x = Inf, y = 0.25,
    vjust = 1, hjust = 1,
    label = paste("Log-rank p = ", round(p_tbl$pval, 3), "\nHR = ", round(exp(cox_TBL$coefficients), 3)," (",  round(s_tbl$conf.int[3], 3), " - ", round(s_tbl$conf.int[4], 3), ")", sep=""),
    size = 5
  ) + xlab("Time (Days)") + ylab("Disease Free Interval (DFI)")

pdf(paste(outputDir, "survival_analysis_DFS_MSS_3.pdf", sep=""), width = 15, height = 15)
ggsurv
dev.off()
