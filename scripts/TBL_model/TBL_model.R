# load libraries
library(caret)
library(doParallel)
library(stringr)
library(PRROC)
library(ROCit)      
library(tidyverse)
library(plyr)
library(biomaRt)
library(org.Hs.eg.db)
library(TCGAbiolinks)
library(AnnotationDbi)
library(DESeq2)
library(ggridges)
library(pheatmap)
library(ggpubr)
library(gridExtra)
source("model_helper_functions.R")
set.seed(0202)

outputDir <- "../../results/TBL_models/"
dir.create(outputDir, recursive = T)

#####################################
##### Preprocessing of the data #####
#####################################

# read data
rna_data <- read.csv("../../results/RNA_Seq/genes/Counts/normalized_RNA_Seq.csv", header = T, row.names = 1)
genomic_metrics <- read.csv("../../results/sample_genome_metrics_matched.csv", header = T, row.names = 1)

# RNA gene annotation
gene_ids <- data.frame(ensembl = colnames(rna_data))

gene_ids$entrez <- mapIds(org.Hs.eg.db,
               keys=as.character(gene_ids$ensembl), 
               column="ENTREZID",
               keytype="ENSEMBL",
               multiVals="first")

gene_ids$symbol <- mapIds(org.Hs.eg.db,
                          keys=as.character(gene_ids$ensembl), 
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")

gene_ids$gene_description <- mapIds(org.Hs.eg.db,
                                    keys=as.character(gene_ids$ensembl), 
                                    column="GENENAME",
                                    keytype="ENSEMBL",
                                    multiVals="first")

# filter out genes with no Entrez ID
gene_ids <- gene_ids[!is.na(gene_ids$entrez),]

rna_data <- rna_data[,gene_ids$ensembl]

# annotate MSI and POLE/D1 as one variable
genomic_metrics$msi_pole <- as.character(genomic_metrics$msi_status)
genomic_metrics[which(genomic_metrics$msi_pole %in% c("MSI-L", "Indeterminate")), "msi_pole"] <- "MSS"
genomic_metrics[is.na(genomic_metrics$msi_pole), "msi_pole"] <- "-"
genomic_metrics[which((genomic_metrics$POLE_mut == T | genomic_metrics$POLD1_mut == T) & genomic_metrics$msi_pole != "MSI-H"), "msi_pole"] <- "POLE/D1"
genomic_metrics$msi_pole <- as.factor(genomic_metrics$msi_pole)

# remove MSI and POLE/D1 mutated samples
genomic_filtered <- genomic_metrics[genomic_metrics$msi_pole == "MSS",]

genomic_filtered <- genomic_filtered[!duplicated(genomic_filtered$sample_ids),]

rownames(genomic_filtered) <- genomic_filtered$sample_ids

rna_data <- rna_data[rownames(genomic_filtered),]

############################################
####### Group selection for training #######
############################################

# determine thresholds to indicate high bp count and low bp count
quantile_marks <- quantile(genomic_filtered$count)
high_mark <- quantile_marks[["75%"]]
low_mark <- quantile_marks[["25%"]]

genomic_filtered$TBL_state <- "I"
genomic_filtered[genomic_filtered$count < low_mark, "TBL_state"] <- "L"
genomic_filtered[genomic_filtered$count > high_mark, "TBL_state"] <- "H"
genomic_filtered$TBL_state <- factor(genomic_filtered$TBL_state)

TBL_filter <- genomic_filtered[genomic_filtered$TBL_state != "I",]
rna_filter <- rna_data[as.character(TBL_filter$sample_ids),]

###########################
###### Model training #####
###########################

# create labels for model
gene_bp_status <- factor(TBL_filter[match(rownames(rna_filter), TBL_filter$sample_ids),]$TBL_state)
names(gene_bp_status) <- rownames(rna_filter)
gene_bp_status <- gene_bp_status %>% revalue(c("H"="abr", "L"="wt"))

set.seed(203)

# run model training only if there are enough bps
splits <- createDataPartition(gene_bp_status, times = 1, p = 0.65, list = TRUE)
train_labels <- gene_bp_status[splits$Resample1]
test_labels <- gene_bp_status[-splits$Resample1]

# subsample rna-seq data
RNA_train <- rna_filter[names(train_labels),]
RNA_test <- rna_filter[names(test_labels),]

results <- run_model(RNA_train, train_labels, method = "rf", folds = 5,
                     tune_length = 10, features = 300, 
                     seq_length = 25, iteration_length = 10,
                     cores = 30, min_features = 10)

ggplot(results$best_features) + xlim(10, 300)
ggplot(results$best_model)

prediction <- predict(results$best_model, RNA_test)
conf_results <- confusionMatrix(data = prediction, reference = test_labels)

prediction <- predict(results$best_model, RNA_test, type = "prob")

roc_bin <- rocit(score = prediction$abr, class = test_labels, negref = "wt")

pdf(paste(outputDir, "TBL_MSS", "_ROC.pdf", sep=""), width = 12, height = 12)
plot(roc_bin, YIndex = F, 
     values = F, legend = F)
legend("bottomright", col = c(1,alpha("grey50")),
       c(paste("Empirical ROC curve (AUC: ", format(roc_bin$AUC, digits = 3), ")", sep=""), "Chance line"), 
       lty= c("solid", "dashed"), lwd = 2)
dev.off()

library(pROC)
myRoc <- roc(predictor = prediction$abr, response = factor(test_labels, 
                                                           levels = c("wt", "abr")), 
             positive = "abr")
plot(myRoc, print.auc = T, print.auc.x = 0.2, print.auc.y = 0.2, las=1, bty="l",
     cex.lab = 1.5, xlim = c(1,0), ylim = c(0,1))

roc_data <- data.frame(obs = prediction, label = as.character(test_labels))
roc_data$label <- as.numeric(revalue(roc_data$label, c(wt = 0, abr = 1)))

library(plotROC)
pdf(paste(outputDir, "TBL_MSS_ROC.pdf", sep=""))
g <- ggplot(roc_data, aes(m=obs.abr, d=label)) +
  geom_roc(n.cuts=0) +
  coord_equal() + 
  style_roc(xlab = "1 - Specificity (FPR)", ylab = "Sensitivity (TPR)", 
            theme = theme_classic, major.breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  geom_abline(linetype = "dashed") + theme_classic(base_size = 21) + 
  theme(axis.text = element_text(size = 18))

g + annotate("text", x=0.75, y=0.25, 
             label=paste("AUC =", round((calc_auc(g))$AUC, 4)),
             size = 7)
dev.off()
  
control_models <- results$best_model
control_performance <- results$max_performance 
control_performance_list <- results$performance_list
final_control_performance <- roc_bin$AUC
final_control_conf_table <- conf_results

saveRDS(results, paste(outputDir, "TBL_model.rds", sep=""))

#########################################
##### Model prediction on all data ######
#########################################
  
# reclassify all samples using the predictive model
full_prediction <- predict(results$best_model, rna_data)
names(full_prediction) <- rownames(rna_data)

genomic_filtered$tbl_predicted <- full_prediction[match(genomic_filtered$sample_ids, names(full_prediction))]

# plot reclassified TBL distribution
pdf(paste(outputDir, "TBL_prediction_distribution.pdf", sep=""))
p <- ggplot(genomic_filtered, aes(x=count, fill = tbl_predicted, colour=tbl_predicted)) + 
  geom_histogram(aes(y=..density..), binwidth = 5, alpha=.5, position = "identity") +
  geom_density(alpha=.3) + xlab("Tumor Break Load") + 
  scale_fill_discrete(name = "TBL expressional profile", 
                      labels = c("Low", "High"), direction = -1) + 
  scale_colour_discrete(name = "TBL expressional profile", 
                        labels = c("Low", "High"), direction = -1) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme_classic(base_size = 21) +  theme(legend.position = c(0.7, 0.3), 
                                       legend.title = element_text(), 
                                       legend.title.align = 0.5, 
                                       axis.text = element_text(size = 18))
p
dev.off()

column <- "tbl_predicted"
coldata <- genomic_filtered[,column]
my_xlab <- paste(levels(factor(coldata, levels = c("wt", "abr"), labels = c("Low", "High"))), "\n(N=", table(coldata), ")", sep="")
combinations <- combn(length(levels(factor(coldata))), 2, simplify = F)
plt_tbl <- ggplot(genomic_filtered, aes(x=tbl_predicted, y=log(count), color = tbl_predicted)) + 
  geom_violin(lwd = 1) + ylab("ln(TBL)") + xlab("") +
  scale_color_discrete(direction = -1) +
  scale_x_discrete(labels = my_xlab) +
  geom_boxplot(varwidth = T, alpha = 0.2, width = 0.1, lwd = 1) +
  stat_summary(fun = median, geom = "text", 
               aes(label=sprintf("%1.1f", ..y..)),
               position = position_nudge(x=0.12), size = 8) +
  stat_compare_means(comparisons = list(c(1,2)), 
                     label.y = c(6.2), size = 8, label = "p.signif") + 
  theme_classic(base_size = 21) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 18)) 
plt_tbl

pdf(paste(outputDir, "full_TBL_prediction_distribution.pdf", sep=""))
p + 
  annotation_custom(ggplotGrob(plt_tbl),
                    xmin = 180, xmax = 340, ymin = 0.018, ymax = 0.035)
dev.off()

write.csv(genomic_filtered, "../../results/sample_genome_metrics_tbl_class.csv")

################################################################
####### Model prediction of MSI or hypermutators samples #######
################################################################

# prepare rna data
# read data
rna_all <- read.csv("../../results/RNA_Seq/genes/Counts/normalized_RNA_Seq.csv", 
                    header = T, row.names = 1)

# RNA gene annotation
gene_ids <- data.frame(ensembl = colnames(rna_all))

gene_ids$entrez <- mapIds(org.Hs.eg.db,
                          keys=as.character(gene_ids$ensembl), 
                          column="ENTREZID",
                          keytype="ENSEMBL",
                          multiVals="first")

gene_ids$symbol <- mapIds(org.Hs.eg.db,
                          keys=as.character(gene_ids$ensembl), 
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")

gene_ids$gene_description <- mapIds(org.Hs.eg.db,
                                    keys=as.character(gene_ids$ensembl), 
                                    column="GENENAME",
                                    keytype="ENSEMBL",
                                    multiVals="first")

# filter out genes with no Entrez ID
gene_ids <- gene_ids[!is.na(gene_ids$entrez),]

rna_all <- rna_all[,gene_ids$ensembl]

# select MSI or hyper mutatoters
genomic_hyper <- genomic_metrics[genomic_metrics$msi_pole != "MSS",]

genomic_hyper <- genomic_hyper[!duplicated(genomic_hyper$sample_ids),]

rownames(genomic_hyper) <- genomic_hyper$sample_ids

rna_hyper <- rna_all[rownames(genomic_hyper),]
rna_msi <- rna_all[rownames(genomic_hyper[genomic_hyper$msi_pole == "MSI-H",]),]

hyper_prediction <- predict(results$best_model, rna_hyper)
names(hyper_prediction) <- rownames(rna_hyper)
msi_prediction <- predict(results$best_model, rna_msi)
names(msi_prediction) <- rownames(rna_msi)

genomic_hyper$tbl_predicted_hyper <- hyper_prediction[match(genomic_hyper$sample_ids, names(hyper_prediction))]
genomic_msi <- genomic_hyper[genomic_hyper$msi_pole == "MSI-H",]
genomic_msi$tbl_predicted_msi <- msi_prediction[match(genomic_msi$sample_ids, names(msi_prediction))]

# plot hyper mutated TBL distribution
pdf(paste(outputDir, "TBL_hyper_prediction_distribution.pdf", sep=""))
p <- ggplot(genomic_hyper, aes(x=count, fill = tbl_predicted_hyper, colour=tbl_predicted_hyper)) + 
  geom_histogram(aes(y=..density..), binwidth = 5, alpha=.5, position = "identity") +
  geom_density(alpha=.3) + xlab("Tumor Break Load") + 
  scale_fill_discrete(name = "TBL predictions", labels = c("High", "Low")) + 
  scale_colour_discrete(name = "TBL predictions", labels = c("High", "Low")) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme_classic(base_size = 21) + 
  ggtitle("TBL predictions Hyper Mutators stage 1-4") +
  theme(axis.text = element_text(size = 18))
p
dev.off()

column <- "tbl_predicted_hyper"
coldata <- genomic_hyper[,column]
my_xlab <- paste(levels(factor(coldata, levels = c("abr", "wt"), labels = c("High", "Low"))), "\n(N=", table(coldata), ")", sep="")
combinations <- combn(length(levels(factor(coldata))), 2, simplify = F)
plt_tbl <- ggplot(genomic_hyper, aes(x=tbl_predicted_hyper, y=log(count), color = tbl_predicted_hyper)) + 
  geom_violin(lwd = 1) + ylab("ln(TBL)") + xlab("") +
  scale_x_discrete(labels = my_xlab) +
  geom_boxplot(varwidth = T, alpha = 0.2, width = 0.1, lwd = 1) +
  stat_summary(fun = median, geom = "text", 
               aes(label=sprintf("%1.1f", ..y..)),
               position = position_nudge(x=0.12), size = 8) +
  stat_compare_means(comparisons = list(c(1,2)), 
                     label.y = c(6.2), size = 8, label = "p.signif") + 
  theme_classic(base_size = 21) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size = 18)) 
plt_tbl

pdf(paste(outputDir, "TBL_prediction_distribution_hyper.pdf", sep=""))
p + 
  annotation_custom(ggplotGrob(plt_tbl),
                    xmin = 140, xmax = 240, ymin = 0.018, ymax = 0.035)
dev.off()

pdf(paste(outputDir, "TBL_only_hyper_prediction_distribution.pdf", sep=""))
p <- ggplot(genomic_hyper[genomic_hyper$msi_pole != "MSI-H",], aes(x=count, fill = tbl_predicted_hyper, colour=tbl_predicted_hyper)) + 
  geom_histogram(aes(y=..density..), binwidth = 5, alpha=.5, position = "identity") +
  geom_density(alpha=.3) + xlab("Tumor Break Load") + scale_fill_discrete(name = "TBL predictions", labels = c("High", "Low")) + scale_colour_discrete(name = "TBL predictions", labels = c("High", "Low")) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme_classic(base_size = 21) + 
  ggtitle("TBL predictions only Hyper Mutators stage 1-4") + 
  theme(axis.text = element_text(size = 18))
p
dev.off()

column <- "tbl_predicted_hyper"
coldata <- genomic_hyper[genomic_hyper$msi_pole != "MSI-H",column]
my_xlab <- paste(levels(factor(coldata, levels = c("abr", "wt"), labels = c("High", "Low"))), "\n(N=", table(coldata), ")", sep="")
combinations <- combn(length(levels(factor(coldata))), 2, simplify = F)
plt_tbl <- ggplot(genomic_hyper[genomic_hyper$msi_pole != "MSI-H",], aes(x=tbl_predicted_hyper, y=log(count), color = tbl_predicted_hyper)) + 
  geom_violin(lwd = 1) + ylab("ln(TBL)") + xlab("") +
  scale_x_discrete(labels = my_xlab) +
  geom_boxplot(varwidth = T, alpha = 0.2, width = 0.1, lwd = 1) +
  stat_summary(fun = median, geom = "text", 
               aes(label=sprintf("%1.1f", ..y..)),
               position = position_nudge(x=0.12), size = 8) +
  stat_compare_means(comparisons = list(c(1,2)), 
                     label.y = c(6.2), size = 8, label = "p.signif") + 
  theme_classic(base_size = 21) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 18)) 
plt_tbl

pdf(paste(outputDir, "TBL_prediction_distribution_only_hyper.pdf", sep=""))
p + 
  annotation_custom(ggplotGrob(plt_tbl),
                    xmin = 140, xmax = 240, ymin = 0.025, ymax = 0.05)
dev.off()

pdf(paste(outputDir, "TBL_pol_prediction_distribution.pdf", sep=""))
p <- ggplot(genomic_hyper[genomic_hyper$msi_pole == "POLE/D1",], aes(x=count, fill = tbl_predicted_hyper, colour=tbl_predicted_hyper)) + 
  geom_histogram(aes(y=..density..), binwidth = 5, alpha=.5, position = "identity") +
  geom_density(alpha=.3) + xlab("Tumor Break Load") + scale_fill_discrete(name = "TBL predictions", labels = c("High", "Low")) + scale_colour_discrete(name = "TBL predictions", labels = c("High", "Low")) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme_classic(base_size = 21) + ggtitle("TBL predictions POLE/D1 stage 1-4") +
  theme(axis.text = element_text(size = 18))
p
dev.off()

column <- "tbl_predicted_hyper"
coldata <- genomic_hyper[genomic_hyper$msi_pole == "POLE/D1",column]
my_xlab <- paste(levels(factor(coldata, levels = c("abr", "wt"), labels = c("High", "Low"))), "\n(N=", table(coldata), ")", sep="")
combinations <- combn(length(levels(factor(coldata))), 2, simplify = F)
plt_tbl <- ggplot(genomic_hyper[genomic_hyper$msi_pole != "MSI-H",], aes(x=tbl_predicted_hyper, y=log(count), color = tbl_predicted_hyper)) + 
  geom_violin(lwd = 1) + ylab("ln(TBL)") + xlab("") +
  scale_x_discrete(labels = my_xlab) +
  geom_boxplot(varwidth = T, alpha = 0.2, width = 0.1, lwd = 1) +
  stat_summary(fun = median, geom = "text", 
               aes(label=sprintf("%1.1f", ..y..)),
               position = position_nudge(x=0.12), size = 8) +
  stat_compare_means(comparisons = list(c(1,2)), 
                     label.y = c(6.2), size = 8, label = "p.signif") + theme_classic(base_size=21) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size = 18)) 
plt_tbl

pdf(paste(outputDir, "TBL_prediction_distribution_pol.pdf", sep=""))
p + 
  annotation_custom(ggplotGrob(plt_tbl),
                    xmin = 80, xmax = 140, ymin = 0.035, ymax = 0.06)
dev.off()

# plot MSI mutated TBL distribution
pdf(paste(outputDir, "TBL_msi_prediction_distribution.pdf", sep=""))
p <- ggplot(genomic_msi, aes(x=count, fill = tbl_predicted_msi, colour=tbl_predicted_msi)) + 
  geom_histogram(aes(y=..density..), binwidth = 5, alpha=.5, position = "identity") +
  geom_density(alpha=.3) + xlab("Tumor Break Load") + scale_fill_discrete(name = "TBL predictions", labels = c("High", "Low")) + scale_colour_discrete(name = "TBL predictions", labels = c("High", "Low")) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme_classic(base_size=21) + ggtitle("TBL predictions MSI-H stage 1-4") +
  theme(axis.text = element_text(size = 18), legend.position = "top")
p
dev.off()

column <- "tbl_predicted_msi"
coldata <- genomic_msi[,column]
my_xlab <- paste(levels(factor(coldata, levels = c("abr", "wt"), labels = c("High", "Low"))), "\n(N=", table(coldata), ")", sep="")
combinations <- combn(length(levels(factor(coldata))), 2, simplify = F)
plt_tbl <- ggplot(genomic_msi, aes(x=tbl_predicted_msi, y=log(count), color = tbl_predicted_msi)) + 
  geom_violin(lwd = 1) + ylab("ln(TBL)") + xlab("") +
  scale_x_discrete(labels = my_xlab) +
  geom_boxplot(varwidth = T, alpha = 0.2, width = 0.1, lwd = 1) +
  stat_summary(fun = median, geom = "text", 
               aes(label=sprintf("%1.1f", ..y..)),
               position = position_nudge(x=0.2), size = 3, color="black") +
  stat_compare_means(comparisons = list(c(1,2)), 
                     label.y = c(6.2), size = 2.5, label = "p.signif") + theme_classic(base_size=14) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size = 11)) 
plt_tbl

pdf(paste(outputDir, "TBL_prediction_distribution_msi.pdf", sep=""))
p + 
  annotation_custom(ggplotGrob(plt_tbl),
                    xmin = 70, xmax = 140, ymin = 0.020, ymax = 0.06)
dev.off()

rownames(genomic_metrics) <- as.character(genomic_metrics$sample_ids)
genomic_metrics$TBL_prediction <- "-"
genomic_metrics[as.character(genomic_filtered$sample_ids), "TBL_prediction"] <- as.character(genomic_filtered$tbl_predicted)
genomic_metrics[as.character(genomic_hyper$sample_ids), "TBL_prediction"] <- as.character(genomic_hyper$tbl_predicted_hyper)

p <- ggplot(genomic_metrics[genomic_metrics$TBL_prediction != "-",], 
            aes(x = count, y = log10(TMB), color = TBL_prediction, shape = msi_pole)) + 
  geom_point() + xlab("Tumor Break Load (TBL)") + ylab("Tumor Mutational Burden (log10(TMB))") + 
  theme_classic(base_size=21) + 
  scale_color_discrete(name = "TBL expr. phen.", labels = c("High", "Low")) +
  labs(shape = "MSI status") +
  theme(axis.text = element_text(size = 18))

pdf("../../results/TBL_models/TBL_TMB_classes.pdf", width = 8, height = 6)
p
dev.off()

write.csv(genomic_metrics, "../../results/sample_genome_metrics_tbl_class_all.csv")

#####################################################################
##### Association of tumor staging and the TBL classification #######
#####################################################################

# create table of stage predicted TBL state counts
query <- GDCquery(project = c("TCGA-COAD", "TCGA-READ"),
                  data.category = "Clinical",
                  file.type = "xml",
                  barcode = genomic_filtered$barcode)
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "stage_event")

genomic_filtered$pathological_stage <- clinical[match(genomic_filtered$barcode, clinical$bcr_patient_barcode), "pathologic_stage"]

genomic_filtered$pathological_stage <- factor(str_replace(genomic_filtered$pathological_stage, pattern = "[ABC]", replacement = ""))

genomic_filtered <- genomic_filtered[genomic_filtered$pathological_stage %in% c("Stage I", "Stage II", "Stage III", "Stage IV"),]

stage_overview <- genomic_filtered %>% group_by(pathological_stage, tbl_predicted) %>% dplyr::summarise(Freq = n())

# test for assessment of independence between the pathological stage and the TBL classification
# there is a significant association between stage and TBL (0.00847)
fisher.test(genomic_filtered$tbl_predicted, genomic_filtered$pathological_stage)

pdf(paste(outputDir, "predicted_TBL_state_in_stages_perc.pdf", sep=""))
ggplot(stage_overview, aes(x=pathological_stage, y=Freq, fill=tbl_predicted)) + 
  geom_bar(position = "fill", stat="identity") + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_discrete(name = "TBL predictions", labels = c("High", "Low")) + 
  ylab("Percentage") + xlab("Pathological stage") +
  scale_x_discrete(labels=c("Stage I" = "I", "Stage II" = "II", "Stage III" = "III", "Stage IV" = "IV")) +
  theme_minimal(base_size = 21) + theme(axis.text = element_text(size = 18))
dev.off()

pdf(paste(outputDir, "predicted_TBL_state_in_stages.pdf", sep=""))
ggplot(stage_overview, aes(x=pathological_stage, y=Freq, fill=tbl_predicted)) + 
  geom_bar(stat="identity") +
  scale_fill_discrete(name = "TBL predictions", labels = c("High", "Low")) + ylab("Frequency") +
  theme_minimal(base_size = 21) + theme(axis.text = element_text(size = 18))
dev.off()

save.image(paste(outputDir, "TBL_rfe_model_all_mss.RData", sep=""))

# randomization performance

TBL_random <- genomic_filtered

random_list <- 1:100

random_results <- NULL
random_results[random_list] <- list(NULL)

random_performance <- NULL
random_performance[random_list] <- list(NULL)

random_performance_list <- NULL
random_performance_list[random_list] <- list(NULL)

final_random_performance <- NULL
final_random_performance[random_list] <- list(NULL)

final_random_roc <- NULL
final_random_roc[random_list] <- list(NULL)

final_random_conf_table <- NULL
final_random_conf_table[random_list] <- list(NULL)

for (i in 1:100){

  set.seed(i)

  # randomize labels
  TBL_random$TBL_state <- sample(TBL_random$TBL_state)

  TBL_filter <- TBL_random[TBL_random$TBL_state != "I",]
  rna_filter <- rna_data[as.character(TBL_filter$sample_ids),]

  # create labels for model
  gene_bp_status <- factor(TBL_filter[match(rownames(rna_filter), TBL_filter$sample_ids),]$TBL_state)
  names(gene_bp_status) <- rownames(rna_filter)

  set.seed(203)

  # run model training only if there are enough bps
  splits <- createDataPartition(gene_bp_status, times = 1, p = 0.70, list = TRUE)
  train_labels <- gene_bp_status[splits$Resample1]
  test_labels <- gene_bp_status[-splits$Resample1]

  # subsample rna-seq data
  RNA_train <- rna_data[names(train_labels),]
  RNA_test <- rna_data[names(test_labels),]

  results.random <- run_model(RNA_train, train_labels, method = "rf", folds = 5,
                              tune_length = 10, features = 300,
                              seq_length = 25, iteration_length = 1,
                              cores = 50, min_features = 10)

  prediction <- predict(results.random$best_model, RNA_test)
  conf_results.random <- confusionMatrix(data = prediction, reference = test_labels)

  prediction <- predict(results.random$best_model, RNA_test, type = "prob")

  roc_results.random <- rocit(score = prediction$H, class = test_labels,
                             negref = "L", method = "bin")

  random_results[[i]] <- results.random
  random_performance[[i]] <- results.random$max_performance
  random_performance_list[[i]] <- results.random$performance_list
  final_random_performance[[i]] <- roc_results.random$auc
  final_random_roc[[i]] <- roc_results.random
  final_random_conf_table[[i]] <- conf_results.random

}

roc_bin <- rocit(score = prediction$H, class = test_labels, negref = "L")

plot(roc_bin, YIndex = F,
     values = F, legend = T)
for (line in final_random_roc) {
  lines(line$TPR~line$FPR,
        col = alpha("grey50", 0.3), lty = "dashed")
}

test <- 0
test_list <- c()
for (line in final_random_roc) {
  test <- test + line$AUC
  test_list <- c(test_list, line$AUC)
}
test <- test/100
sd(test_list)

sens <- NULL
sens[random_list] <- list(NULL)
spec <- NULL
spec[random_list] <- list(NULL)

for (i in 1:100) {
  sens[[i]] <- final_random_roc[[i]]$sensitivities
  spec[[i]] <- final_random_roc[[i]]$specificities
}

roc_data <- data.frame(sens = sens, spec = spec)

ggplot(roc_data, aes(x = sens, y = spec)) + geom_point() + geom_smooth() + 
  scale_x_continuous(trans = "reverse") + geom_abline(slope=1, intercept = 1, 
                                                      linetype = "dashed", 
                                                      alpha=0.7, color = "grey") + 
  coord_equal()



