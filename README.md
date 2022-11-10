# Tumor Break Load
Supporting data and code for the paper titled "Tumor Break Load is a biologically relevant feature of genomic instability with prognostic value in colorectal cancer"

# Abstract
**Background**

Over 80% of colorectal cancers (CRC) exhibit chromosomal instability (CIN). Chromosome segregation errors and double-strand break repair defects lead to somatic copy number aberrations (SCNAs) and chromosomal rearrangement-associated structural variants (SVs), respectively. We hypothesize that the number of SVs is a distinct feature of genomic instability and propose a new measure to quantify SVs: the tumor break load (TBL).

**Aim** 

To characterize the biological impact and clinical relevance of TBL in CRC.

**Data and Methods** 

Disease-free survival (DFS) and SCNA data were obtained from The Cancer Genome Atlas (TCGA) and two independent CRC studies. TBL was defined as the sum of SCNA-associated SVs. TCGA RNA gene expression data of microsatellite stable (MSS) CRC samples were used to train an RNA-based TBL classifier. Dichotomized DNA-based TBL data were used for survival analysis.

**Results**

There is large variation in TBL across CRC samples (median: 47; range: 0-337) with poor correlation to the number of point mutations (tumor mutational burden, R2: 0.057) and abundancy of SCNAs (fraction of genome altered, R2: 0.072). TBL-status could be classified with high accuracy (AUC: 0.88; p<0.01), illustrating its impact on tumor biology. High TBL was associated with higher risk of disease recurrence in 85 stage II-III MSS CRCs from TCGA (HR: 6.1; p = 0.007), which was further validated in independent series of 57 untreated stage II-III (HR: 4.1; p = 0.012) and 74 untreated stage II MSS CRCs (HR: 2.4; p = 0.01).

**Conclusion**

TBL is a prognostic biomarker in patients with non-metastatic MSS CRC.

# Project
The code provided in this repository describes three aspects introduced in the paper:
  - The calculation of the tumor break load and its characteristics in colorectal cancer.
  - Identifying if TBL has impact on CRC biology. 
  - Clinical impact of the tumor break load in CRC.
  
run_project.sh will run the whole analysis workflow to reproduce the results presented in the paper.

Data and results files can be retrieved from codeocean (https://codeocean.com/capsule/1605813/tree) including the model files.

# Requirements
With TBL_Dockerfile.dockerfile the required environment to run the current workflow can be recreated.

# Contact
dr. Sanne Abeln, s.abeln@vu.nl

MSc Soufyan Lakbir, s.lakbir@vu.nl
