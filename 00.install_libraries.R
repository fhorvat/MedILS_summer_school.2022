### CRAN libraries
libraries_cran <- c("dplyr", "readr", "stringr", "tibble", "ggplot2", "magrittr", 
                    "RColorBrewer", "plotly", "htmlwidgets", "gprofiler2")

# install if not installed
libraries_cran_to_install <- libraries_cran[!(libraries_cran %in% installed.packages()[, "Package"])]
if(length(libraries_cran_to_install)) install.packages(libraries_cran_to_install)


### Bioconductor libraries
libraries_bio <- c("DESeq2", "pheatmap", "EnhancedVolcano")

# first install BiocManager
if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# install if not installed
libraries_bio_to_install <- libraries_bio[!(libraries_bio %in% installed.packages()[, "Package"])]
if(length(libraries_bio_to_install)) BiocManager::install(libraries_bio_to_install)

