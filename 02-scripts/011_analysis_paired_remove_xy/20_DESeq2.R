#!/usr/bin/env Rscript
'02_DESeq2.R

Usage:
  02_DESeq2.R <sample_sheet> <count_table> --result_dir=<res_dir> --c1=<c1> --c2=<c2> [options]
  02_DESeq2.R --help

Arguments:
  <sample_sheet>                CSV file with the sample annotations.
  <count_table>                 TSV file with the read counts

Mandatory options:
  --result_dir=<res_dir>        Output directory
  --c1=<c1>                     Contrast level 1 (perturbation). Needs to be contained in condition_col.
  --c2=<c2>                     Contrast level 2 (baseline). Needs to be contained in condition_col.

Optional options:
  --condition_col=<cond_col>    Column in sample annotation that contains the condition [default: group]
  --covariate_formula=<formula> Formula to model additional covariates (need to be columns in the samplesheet)
                              that will be appended to the formula built from `condition_col`.
                              E.g. `+ age + sex`. Per default, no covariates are modelled.
  --sample_col=<sample_col>     Column in sample annotation that contains the sample names
                                (needs to match the colnames of the count table). [default: sample]
  --batch_col=<batch_col>       Optional: column in sample annotation that contains the batch [default: batch]
  --plot_title=<title>          Title shown above plots. Is built from contrast per default.
  --prefix=<prefix>             Results file prefix. Is built from contrasts per default.
  --fdr_cutoff=<fdr>            False discovery rate for GO analysis and volcano plots [default: 0.1]
  --fc_cutoff=<log2 fc cutoff>  Fold change (log2) cutoff for volcano plots [default: 1]
  --gtf_file=<gtf>              Path to the GTF file used for featurecounts. If specified, a Biotype QC
                                will be performed.
  --gene_id_type=<id_type>      Type of the identifier in the `gene_id` column compatible with AnnotationDbi [default: ENSEMBL]
  --n_cpus=<n_cpus>             Number of cores to use for DESeq2 [default: 1]
  --skip_gsea                   Skip Gene-Set-Enrichment-Analysis step
  --genes_of_interest=<genes>   File (tsv) containing a list of genes to highlight in the volcano plot (column must be named
                                "gene_name").
                                If an optional column named "group" is present, each gene will be associated with the corresponding
                                gene group and a separate volcanon plot for each gene group will be generated (e.g. cytokines).
  --organism=<human|mouse>      Ensebml annotation db [default: human]
  --save_workspace              Save R workspace for this analysis [default: FALSE]
' -> doc

library("conflicted")
library("docopt")
#arguments <- docopt(doc, version = "0.1")

# print(arguments)

library("BiocParallel")
library("DESeq2")
library("IHW")
library("ggplot2")
library("pcaExplorer")
library("topGO")
library("clusterProfiler")
library("ReactomePA")
library("writexl")
library("readr")
library("dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("count", "dplyr")
library("EnhancedVolcano")
library("ggpubr")
library("tibble")
library("stringr")
library("ggrepel")
conflict_prefer("paste", "base")
conflict_prefer("rename", "dplyr")
library("biomaRt")
remove_ensg_version = function(x) gsub("\\.[0-9]*$", "", x)


#sampleAnnotationCSV <- arguments$sample_sheet
#readCountFile <- arguments$count_table
#result_dir = arguments$result_dir
#c1 = arguments$c1
#c2 = arguments$c2
#dir.create(result_dir, recursive=TRUE, showWarnings=FALSE)
#paired_grp <- arguments$paired_grp
#
##prefix and plot title
#prefix <- arguments$prefix
#plot_title <- arguments$plot_title
#
## Sample information and contrasts
#nfcore = arguments$nfcore
#cond_col = arguments$condition_col
#sample_col = arguments$sample_col
#contrast = c(cond_col,c1, c2)
#gene_id_type = arguments$gene_id_type
#covariate_formula = arguments$covariate_formula
#remove_batch_effect = arguments$remove_batch_effect
#batch_col = arguments$batch_col
#
##Cutoff
#fdr_cutoff = as.numeric(arguments$fdr_cutoff)
#fc_cutoff = as.numeric(arguments$fc_cutoff)
#
##GTF for Biotype QC
#gtf_file = arguments$gtf_file
#
## Other
#n_cpus = as.numeric(arguments$n_cpus)
#skip_gsea = arguments$skip_gsea
#genes_of_interest = arguments$genes_of_interest
#
##set organism (human or mouse)
#organism = arguments$organism
#
##save R workspace
#save_ws = arguments$save_workspace
#
#
# DEBUG parameters 
sampleAnnotationCSV="/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/012_LUAD/deseq2_out/pb_cell_type/tumor_vs_normal/B_cell/samplesheet_B_cell.csv"
readCountFile="/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/012_LUAD/deseq2_out/pb_cell_type/tumor_vs_normal/B_cell/counts_B_cell.csv"
resDir="/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/012_LUAD/deseq2_out/out/" 
resDir_plot="/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/012_LUAD/deseq2_out/out/figures/" 

c1="tumor_primary"
c2="normal_adjacent"
cond_col ="origin"
contrast = c(cond_col, c1, c2)
organism="human"
n_cpus = 8
plot_title="DESEQ2"
prefix = "nsclc_gender_B_cell"
gene_id_type="ENSEMBL"
sample_col="sample"
covariate_formula=""
fdr_cutoff=0.1
fc_cutoff=1
paired_grp="donor_id"


# Reading the Annotation sample csv file
sampleAnno <- read.csv(sampleAnnotationCSV,row.names=1)
rownames(sampleAnno) <- gsub("-","_",rownames(sampleAnno))
rownames(sampleAnno) <- gsub(" ","_",rownames(sampleAnno))
sampleAnno$sample <- rownames(sampleAnno)
sampleAnno <- sampleAnno[,-1]

if(is.null(paired_grp)) {
  design_formula <- as.formula(paste0("~", cond_col, covariate_formula))
} else {
  design_formula <- as.formula(paste0("~", paired_grp , " +", cond_col, covariate_formula))
}


# Reading the Count matrix tsv file
count_mat <- read_csv(readCountFile)
colnames(count_mat)[1] <- "gene_id"

# colnames(count_mat)[1] <- "gene_id"

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- c(count_mat$gene_id) 

gene_symbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = ensembl_ids,
                      mart = ensembl)
colnames(gene_symbols)[1] <- "gene_id"
colnames(gene_symbols)[2] <- "gene_name"

# Merge the gene symbols with the original dataframe 'counts'
count_mat <- merge(count_mat, gene_symbols, by.x = "gene_id", by.y = "gene_id", all.x = TRUE)

ensg_to_genesymbol = count_mat %>% select(gene_id, gene_name)
count_mat <- count_mat[, -ncol(count_mat)]

count_mat <- count_mat |>
  column_to_rownames("gene_id") |>
  round()


colnames(count_mat) <- gsub("-", "_", colnames(count_mat))


# Check if names are the same
if (!all(rownames(sampleAnno) %in% colnames(count_mat))) {
  print('Row names in sample annotation and column names in count matrix are not the same')
  break
}

# Start processing
#design_formula <- as.formula(paste0("~", cond_col, covariate_formula))

dds <- DESeqDataSetFromMatrix(countData = round(count_mat),
                              colData = sampleAnno,
                              design = design_formula)

## Keep only genes where  >= 10 reads per sample condition in total
keep <- rowSums(counts(collapseReplicates(dds, dds[[cond_col]]))) >= 10
dds <- dds[keep, ]

# run DESeq
dds.res <- DESeq(dds, parallel = (n_cpus > 1))

# Define contrast names
contrasts <- list(c(cond_col, c1, c2))
names(contrasts) <- sprintf("%s_vs_%s", c1, c2)

# ## IHW

# Use of IHW for p value adjustment of DESeq2 results
resSHRINK  <- lfcShrink(dds.res, contrast= contrast, type ="normal") #specifying "normal" because "apeglm" need coef instead of contrast. 
resSHRINK <- as.data.frame(resSHRINK) |>
  rownames_to_column(var = "gene_id") |>
  as_tibble() |>
  arrange(padj) |>
  rename(lfcSE_shrink = lfcSE) |>
  rename(log2FoldChange_shrink = log2FoldChange)

resIHW <- lapply(names(contrasts), function(name) {
  contrast = contrasts[[name]]
  results(dds.res, filterFun = ihw, contrast = contrast, altHypothesis="greaterAbs",lfcThreshold = 1) |>
    as_tibble(rownames = "gene_id") |>
    mutate(comparison = name) |>
    arrange(pvalue)
}) |> bind_rows()


resIHW <- left_join(resIHW, select(resSHRINK, c(gene_id,log2FoldChange_shrink,lfcSE_shrink)), by="gene_id")


# Add gene names
resIHW_gene_name = resIHW %>%  inner_join(ensg_to_genesymbol, by=c("gene_id"))

# Filter for adjusted p-value < fdr_cutoff
resIHWsig <- resIHW_gene_name %>% filter(padj < fdr_cutoff)

# significant genes as DE gene FDR < fdr_cutoff & abs(logfoldchange) > fc_cutoff , all genes as background
resIHWsig_fc <- resIHWsig %>% filter(abs(log2FoldChange) > fc_cutoff)


# write results to TSV files
write_csv(resIHW_gene_name,  paste0(resDir,"/",prefix, "_all_genes_DESeq2_result.csv"))
write_csv(resIHWsig,  paste0(resDir,"/",prefix, "_sig_genes_DESeq2_result.csv"))
write_csv(resIHWsig_fc,  paste0(resDir,"/",prefix, "_sig_fc_genes_DESeq2_result.csv"))
resIHW_gene_name_tumor = read_csv("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/deseq2_out/corrected_ds/nsclc_gender_tumor_all_genes_DESeq2_result.csv")

resIHW_gene_name_normal = read_csv("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/deseq2_out/corrected_ds/nsclc_gender_normal_all_genes_DESeq2_result.csv")

resIHW_gene_name_all = read_csv("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/deseq2_out/corrected_ds/nsclc_gender_all_all_genes_DESeq2_result.csv")


my_list = c("RPS4Y1",  "XIST"  ,  "ANGPTL4", "CST6"  ,  "SFTPA2",  "CTSE" ,   "MAL2",    "SFTA3" )
v <- EnhancedVolcano(resIHW_gene_name,
                     lab = resIHW_gene_name$gene_id,
                     x = 'log2FoldChange',
                     y = 'padj',
                     pointSize = 1.0,
                     #labSize = 3.0,
                     pCutoff = fdr_cutoff,
                     #drawConnectors = TRUE,
                     FCcutoff = fc_cutoff,
                     caption = paste0("fold change cutoff: ", round(2**fc_cutoff, 1), ", adj.p-value cutoff: ", fdr_cutoff))

v

ggsave(paste0(resDir,prefix,"_volcano_plot.png"), plot = v, width = 6, height = 12, units = "in")


