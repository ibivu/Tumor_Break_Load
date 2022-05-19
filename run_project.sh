Rscript -e "rmarkdown::render('scripts/Break_calling/break_point_detection.Rmd')"
Rscript scripts/Break_calling/segment_profile.R
Rscript scripts/Tumor_Break_Load/Tumor_Break_Load_quantification.R
Rscript scripts/MSI_POLE_D1/Legacy_MSI_annotation.R 
Rscript scripts/CNV_score/copy_number_variation_score.R
Rscript scripts/Tumor_Mutational_Load/TMB_quantification.R
Rscript scripts/MSI_POLE_D1/POL_annotation.R
Rscript -e "rmarkdown::render('scripts/Tumor_Break_Load/TBL_analysis.Rmd')"
Rscript scripts/RNA_normalization/RNA_Seq_normalization.R
Rscript scripts/MSI_model/MSI_model.R
Rscript scripts/TBL_model/TBL_model.R
Rscript -e "rmarkdown::render('scripts/FGA_model/FGA_model.Rmd')"
Rscript -e "rmarkdown::render('scripts/Pathway_genetic_alteration_landscape/curated_pathway_tbl_association.Rmd')"
Rscript scripts/survival_analysis/TBL_survival_analysis_curated_data.R
Rscript scripts/survival_analysis/TBL_proc_threshold_survival_analysis.R
Rscript scripts/survival_analysis/curated_pathway_survival_analysis.R
Rscript scripts/survival_analysis/CMS_survival_analysis.R
Rscript scripts/TBL_CMS_comparison/TBL_CMS_comparison.R
Rscript validation/Orsetti/scripts/tbl_clinical_validation_proc.R
Rscript -e "rmarkdown::render('validation/SaraLahoz/scripts/breakpoint_detection_dataSara.Rmd')"