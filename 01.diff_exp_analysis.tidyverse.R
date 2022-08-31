### INFO: differential expression analysis script for MedILS summer school 2022
### DATE: Thu Aug 25 19:35:32 2022
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/PhD/MedILS_summer_school/practical_lecture")

######################################################## LIBRARIES
# data wrangling libraries
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(ggplot2)
library(magrittr)

# biology stuff libraries
library(DESeq2)
library(gprofiler2)

# make stuff pretty libraries
library(RColorBrewer)
library(pheatmap)
library(EnhancedVolcano)

# make stuff clicky libraries
library(plotly)
library(htmlwidgets)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set input path 
inpath <- getwd()

# create and set outpath
outpath <- file.path(getwd(), "results")
dir.create(outpath)

# counts table path
counts_path <- file.path(inpath, "Stein_2015_PLoSGenet_GSE57514.ensembl.93.GRCm38.p6.20180919.UCSCseqnames.counts.txt")

## some other ways to get path to file
# counts_path <- list.files(inpath, ".*\\.counts\\.txt")
# counts_path <- "C:/Users/fhorvat/Dropbox/Bioinfo/PhD/MedILS_summer_school/practical_lecture/Stein_2015_PLoSGenet_GSE57514.ensembl.93.GRCm38.p6.20180919.UCSCseqnames.counts.txt"
# counts_path <- "https://s3-elixir.cloud.e-infra.cz/fhorvat/public/other_data/Stein_2015_PLoSGenet_GSE57514.ensembl.93.GRCm38.p6.20180919.UCSCseqnames.counts.txt"

# sample table path
sample_table_path <- file.path(inpath, "Stein_2015_PLoSGenet_GSE57514.sampleTable.csv")

######################################################## READ DATA
# read counts from featureCounts
counts_tb <-
  read_delim(counts_path, delim = "\t", comment = "#") %>%
  set_colnames(., basename(colnames(.)))

# read sample table
sample_table <- read_csv(sample_table_path)

######################################################## MAIN CODE
### prepare tables
# prepare sample table for SummarizedExperiment colData()
sample_table_dds <-
  sample_table %>%
  mutate(genotype = factor(genotype, levels = c("Ago2_WT", "Ago2_KO", "Ago2_CI",
                                                "Dicer_WT", "Dicer_KO"))) %>% 
  as.data.frame(.) %>%
  set_rownames(., .$sample_id) 

# counts table
se <-
  counts_tb %>%
  dplyr::select(-c(Chr:Length)) %>%
  dplyr::rename(gene_id = Geneid) %>%
  mutate_if(is.numeric, round, digits = 0) %>%
  as.data.frame(.) %>%
  set_rownames(., .$gene_id) %>%
  dplyr::select(-gene_id) %>%
  as.matrix(.)


### create SummarizedExperiment object
# filter data to include only chosen stage, set colData()
se <- SummarizedExperiment(list(counts = se))
colnames(se) <- str_remove_all(colnames(se), "\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam$")
se <- se[, colnames(se)[match(rownames(sample_table_dds), colnames(se))]]

# check if colnames of assay match rownames of sample table DataFrame
if(all(colnames(se) == rownames(sample_table_dds))){
  
  # set sample table as colData
  colData(se) <- DataFrame(sample_table_dds)
  
}else{
  
  # stop script with warning
  stop("Columns in assay are not matching row of sample table. Please check your data annotation")
  
}

### transform to DESeqDataSet => class similar to SummarizedExperiment
dds <- DESeqDataSet(se, design = ~genotype)


####### EXPLORATORY ANALYSIS
if(TRUE){
  
  ### data for exploration = rlog transformed counts
  rlog_df <- rlog(dds)
  
  
  ### PCA plot
  # plot PCA
  DESeq2::plotPCA(rlog_df, intgroup = "genotype")
  
  
  ### heatmap of the sample-to-sample distances
  # calculate distance
  dist_df <-
    rlog_df %>%
    assay(.) %>% 
    t(.) %>%
    dist(.)
  
  # make matrix
  dist_matrix <- as.matrix(dist_df)
  colnames(dist_matrix) <- NULL
  
  # annotation data.frame
  annotation_df <-
    sample_table_dds %>%
    dplyr::select(genotype)
  
  # rownames annotation
  annotation_rownames <-
    rownames(dist_matrix) %>%
    str_remove_all(., "^s_|\\.PE$|\\.SE$") %>%
    str_replace_all(., "_", " ")
  
  # plot
  pheatmap(dist_matrix,
           clustering_distance_rows = dist_df,
           clustering_distance_cols = dist_df,
           col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
           annotation_row = annotation_df,
           labels_row = annotation_rownames,
           # file = file.path(outpath, "plot.distance.Stein_2015.png"),
           height = 10,
           width = 14)
  
}


### run main DESeq2 function
# DESeq
dds_deseq <- DESeq(dds)

# get results, shrink logFC
dds_shrink <- lfcShrink(dds_deseq, coef = "genotype_Ago2_KO_vs_Ago2_WT")

# get results table
results_tb <-
  dds_shrink %>%
  as_tibble(., rownames = "gene_id") %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(comparison = str_c("Ago2_KO", "_vs_", "Ago2_WT"))

# set significance cutoffs - p-adjusted and log2 fold change values
padj_cut <- 0.01
lfc_cut <- 2

# get only significant results
results_tb_sign <-
  results_tb  %>%
  dplyr::filter(padj <= padj_cut, 
                abs(log2FoldChange) >= lfc_cut)


### write results
if(TRUE){
  
  ### write results as simple .csv
  readr::write_csv(x = results_tb, file = file.path(outpath, "results.Stein_2015.Ago2_KO_vs_Ago2_WT.all.csv"))
  readr::write_csv(x = results_tb_sign, file = file.path(outpath, "results.Stein_2015.Ago2_KO_vs_Ago2_WT.significant.csv"))
  
  
  ### write results as fancy .xlsx
  # open workbook
  results_wb <- openxlsx::createWorkbook()
  
  # add sheets
  openxlsx::addWorksheet(wb = results_wb, sheetName = "Ago2_KO_vs_Ago2_WT.all")
  openxlsx::addWorksheet(wb = results_wb, sheetName = "Ago2_KO_vs_Ago2_WT.significant")
  
  # write data to sheets
  openxlsx::writeData(wb = results_wb, sheet = "Ago2_KO_vs_Ago2_WT.all", x = results_tb)
  openxlsx::writeData(wb = results_wb, sheet = "Ago2_KO_vs_Ago2_WT.significant", x = results_tb_sign)
  
  # save to file
  openxlsx::saveWorkbook(wb = results_wb,
                         file = file.path(outpath, "results.Stein_2015.Ago2_KO_vs_Ago2_WT.xlsx"),
                         overwrite = TRUE)
}


### Visualization of results
## plot counts for individual gene
# get the gene with the smallest p-adjusted value
plotCounts(dds, gene = which.min(results_tb$padj), intgroup = "genotype")

# plot the most up-regulated gene
plotCounts(dds, gene = which.max(results_tb$log2FoldChange), intgroup = "genotype")

# plot the most down-regulated gene
plotCounts(dds, gene = which.max(-results_tb$log2FoldChange), intgroup = "genotype")


### plot heatmap of the counts matrix
# get the top 20 genes with the largest mean expression
top20_genes_mean_expression <- order(rowMeans(counts(dds_deseq, normalized = TRUE)), decreasing = TRUE)[1:20]

# get annotation 
annot_df <- as.data.frame(colData(dds)[, c("sample_id", "genotype")])

# get the heatmap
pheatmap(assay(rlog_df)[top20_genes_mean_expression, ], 
         cluster_rows = FALSE, 
         show_rownames = TRUE,
         cluster_cols = FALSE, 
         annotation_col = annot_df)


### plots
# plot MA
DESeq2::plotMA(dds_shrink, ylim = c(-4, 4), alpha = padj_cut)

# Volcano plots
EnhancedVolcano(dds_shrink, 
                lab = rownames(dds_shrink), 
                x = "log2FoldChange", y = "pvalue", 
                pCutoff = padj_cut,
                FCcutoff = lfc_cut,
                pointSize = 3.0)


### functional enrichment
# prepare gene list
geneList <- results_tb_sign$log2FoldChange
names(geneList) <- results_tb_sign$gene_id
geneList <- sort(geneList, decreasing = TRUE)

# do the enrichment
gostres <- gost(query = names(geneList), 
                organism = "mmusculus", ordered_query = TRUE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)

# get results
goterms_results_tb <- gostres$result

# save
readr::write_csv(x = goterms_results_tb, file = file.path(outpath, "results.GO_terms.Stein_2015.Ago2_KO_vs_Ago2_WT.csv"))


## visualize interactively
# create object
interactive_go_plot <- gostplot(gostres, capped = TRUE, interactive = TRUE)

# save as html widget
htmlwidgets::saveWidget(plotly::as_widget(interactive_go_plot),
                        file = file.path(outpath, "plot.GO_plot.Stein_2015.Ago2_KO_vs_Ago2_WT.shrink.interactive.html"),
                        selfcontained = T)

# visualize as plot with table
p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
publish_gostplot(p, highlight_terms = c("GO:0032502", "HP:0100626"), 
                 width = NA, height = NA, filename = NULL)

# visualize as table
publish_gosttable(gostres, highlight_terms = gostres$result[c(1:2, 10, 33), ],
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)


### Interactive analysis using Ideal 
## user guide: http://bioconductor.org/packages/release/bioc/vignettes/ideal/inst/doc/ideal-usersguide.html
## to install:
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ideal")
library("ideal")
ideal(dds_obj = dds_deseq, res_obj = dds_shrink)

