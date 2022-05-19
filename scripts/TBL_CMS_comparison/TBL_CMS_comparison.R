# read CMS classifciaton data
cms <- read.table("../../data/CMS/cms_labels_public_all.txt", header = T)
cms <- cms[cms$dataset == "tcga",] # select only TCGA samples

# read TBL classification
tbl <- read.csv("../../results/sample_genome_metrics_tbl_class.csv", header = T)

# add cms classification to TBL data frame
tbl$cms <- cms[match(as.character(tbl$barcode), as.character(cms$sample)), "CMS_final_network_plus_RFclassifier_in_nonconsensus_samples"]

# label NAs NOLBL
tbl[is.na(tbl$cms),"cms"] <- "NOLBL"

# map TBL levels to high and low
tbl$tbl_predicted <- factor(tbl$tbl_predicted, levels = c("abr", "wt"), labels = c("TBL high", "TBL low"))

barplot(table(tbl[,c("tbl_predicted", "cms")]), legend = T)

# TBL and the CMS classification is not independent and are related
# significance of 1.154e-05
fisher.test(tbl$tbl_predicted, tbl$cms)

write.csv(table(tbl[,c("tbl_predicted", "cms")]), "tbl_cms.csv")
