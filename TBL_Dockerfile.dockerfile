# hash:sha256:9f269e4ca98d3a6ad9bb6da9b6fecf773162ab36db2f5b47badd9df4315e9308
FROM registry.codeocean.com/codeocean/r-studio:1.4.1106-r4.0.5-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        cmake=3.10.2-1ubuntu2.18.04.2 \
        libnlopt-dev=2.4.2+dfsg-4 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("PRROC", "1.3.1")' \
    && Rscript -e 'remotes::install_version("ROCit", "2.1.1")' \
    && Rscript -e 'remotes::install_version("RcppEigen", "0.3.3.9.2")' \
    && Rscript -e 'remotes::install_version("abind", "1.4-5")' \
    && Rscript -e 'remotes::install_version("carData", "3.0-5")' \
    && Rscript -e 'remotes::install_version("caret", "6.0-92")' \
    && Rscript -e 'remotes::install_version("corrplot", "0.92")' \
    && Rscript -e 'remotes::install_version("cowplot", "1.1.1")' \
    && Rscript -e 'remotes::install_version("devtools", "2.4.3")' \
    && Rscript -e 'remotes::install_version("fmsb", "0.7.3")' \
    && Rscript -e 'remotes::install_version("ggfortify", "0.4.14")' \
    && Rscript -e 'remotes::install_version("ggpmisc", "0.4.6")' \
    && Rscript -e 'remotes::install_version("ggprism", "1.0.3")' \
    && Rscript -e 'remotes::install_version("ggrepel", "0.9.1")' \
    && Rscript -e 'remotes::install_version("ggridges", "0.5.3")' \
    && Rscript -e 'remotes::install_version("ggsci", "2.9")' \
    && Rscript -e 'remotes::install_version("ggsignif", "0.6.3")' \
    && Rscript -e 'remotes::install_version("ggtext", "0.1.1")' \
    && Rscript -e 'remotes::install_version("ggthemes", "4.2.4")' \
    && Rscript -e 'remotes::install_version("lmvar", "1.5.2")' \
    && Rscript -e 'remotes::install_version("maptools", "1.1-4")' \
    && Rscript -e 'remotes::install_version("maxstat", "0.7-25")' \
    && Rscript -e 'remotes::install_version("minqa", "1.2.4")' \
    && Rscript -e 'remotes::install_version("pheatmap", "1.0.12")' \
    && Rscript -e 'remotes::install_version("plotROC", "2.2.1")' \
    && Rscript -e 'remotes::install_version("randomForest", "4.7-1")' \
    && Rscript -e 'remotes::install_version("sampling", "2.9")' \
    && Rscript -e 'remotes::install_version("survMisc", "0.5.6")' \
    && Rscript -e 'remotes::install_version("survminer", "0.4.9")'

RUN Rscript -e 'install.packages("BiocManager")'
RUN Rscript -e 'BiocManager::install(c( \
        "AnnotationDbi", \
        "BiocVersion", \
        "ComplexHeatmap", \
        "CopyNumberPlots", \
        "DESeq2", \
        "GenomicRanges", \
        "MASS", \
        "PRROC", \
        "RColorBrewer", \
        "ROCit", \
        "TCGAbiolinks", \
        "biomaRt", \
        "caret", \
        "circlize", \
        "data.table", \
        "doParallel", \
        "dplyr", \
        "edgeR", \
        "fmsb", \
        "ggbio", \
        "ggfortify", \
        "ggplot2", \
        "ggpmisc", \
        "ggprism", \
        "ggpubr", \
        "ggridges", \
        "ggthemes", \
        "gridExtra", \
        "karyoploteR", \
        "knitr", \
        "lmvar", \
        "maftools", \
        "org.Hs.eg.db", \
        "pROC", \
        "pheatmap", \
        "plotROC", \
        "plyr", \
        "readxl", \
        "reshape2", \
        "scales", \
        "stringr", \
        "survival", \
        "survminer", \
        "tibble", \
        "tidyr", \
        "tidyverse", \
        "umapr" \
    ))' # Original versions: 1.56.2 3.14.0 2.10.0 1.10.0 1.34.0 1.46.1 7.3-57 1.3.1 1.1-3 2.1.1 2.22.4 2.50.3 6.0-92 0.4.14 1.14.2 1.0.17 1.0.9 3.36.0 0.7.3 1.42.0 0.4.14 3.3.6 0.4.6 1.0.3 0.4.0 0.5.3 4.2.4 2.3 1.20.3 1.39 1.5.2 2.10.05 3.14.0 1.18.0 1.0.12 2.2.1 1.8.7 1.4.0 1.4.4 1.2.0 1.4.0 3.3-1 0.4.9 3.1.7 1.2.0 1.3.1 latest
