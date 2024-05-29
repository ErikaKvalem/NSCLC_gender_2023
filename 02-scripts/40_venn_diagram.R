#install.packages("VennDiagram")
require(data.table)
# Load library
library(VennDiagram)

# load ggvenn package 
library("ggvenn") 


# load gplots package 
library("gplots") 

input_path = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/deseq2_out/"
chrom = read_csv("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/input/adata_var_nsclc_chrom.csv")
colnames(chrom)[2] <- "gene_id"

normal_deg = read_csv(paste0(input_path, "nsclc_gender_normal_all_genes_DESeq2_result.csv"))
tumor_deg = read_csv(paste0(input_path, "nsclc_gender_tumor_all_genes_DESeq2_result.csv"))
tumor_and_normal_deg = read_csv(paste0(input_path, "nsclc_gender_tumor_and_normal_all_genes_DESeq2_result.csv"))

# Assuming ch_df and normal are your dataframes
normal_deg_merged <- merge(normal_deg, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)
tumor_deg_merged <- merge(tumor_deg, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)
tumor_and_normal_deg_merged <- merge(tumor_and_normal_deg, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)

# Assuming normal is your dataframe
normal_genes_X <- normal_deg_merged$gene_id[normal_deg_merged$chromosome_name == "X"]
normal_genes_X <-na.omit(normal_genes_X)
normal_genes_Y <- normal_deg_merged$gene_id[normal_deg_merged$chromosome_name == "Y"]
normal_genes_Y <-na.omit(normal_genes_Y)
normal_genes_autosomal <- normal_deg_merged$gene_id[!(normal_deg_merged$chromosome_name == "X" | normal_deg_merged$chromosome_name == "Y")]
normal_genes_autosomal <- na.omit(normal_genes_autosomal)

tumor_genes_X <- tumor_deg_merged$gene_id[tumor_deg_merged$chromosome_name == "X"]
tumor_genes_X <-na.omit(tumor_genes_X)
tumor_genes_Y <- tumor_deg_merged$gene_id[tumor_deg_merged$chromosome_name == "Y"]
tumor_genes_Y <-na.omit(tumor_genes_Y)
tumor_genes_autosomal <- tumor_deg_merged$gene_id[!(tumor_deg_merged$chromosome_name == "X" | tumor_deg_merged$chromosome_name == "Y")]
tumor_genes_autosomal <- na.omit(tumor_genes_autosomal)

tumor_and_normal_genes_X <- tumor_and_normal_deg_merged$gene_id[tumor_and_normal_deg_merged$chromosome_name == "X"]
tumor_and_normal_genes_X <-na.omit(tumor_and_normal_genes_X)
tumor_and_normal_genes_Y <- tumor_and_normal_deg_merged$gene_id[tumor_and_normal_deg_merged$chromosome_name == "Y"]
tumor_and_normal_genes_Y <-na.omit(tumor_and_normal_genes_Y)
tumor_and_normal_genes_autosomal <- tumor_and_normal_deg_merged$gene_id[!(tumor_and_normal_deg_merged$chromosome_name == "X" | tumor_and_normal_deg_merged$chromosome_name == "Y")]
tumor_and_normal_genes_autosomal <- na.omit(tumor_and_normal_genes_autosomal)



########### 2 circles 
# use list as input 
X <-list(normal_x = c(normal_genes_X),tumor_x= c(tumor_genes_X)) 
Y <-list(normal_y = c(normal_genes_Y),tumor_y= c(tumor_genes_Y)) 
autosomal <- list(normal_autosomal = c(normal_genes_autosomal),tumor_autosomal= c(tumor_genes_autosomal)) 
# create venn diagram and display all the sets 
ggvenn(X)
ggvenn(Y)
ggvenn(autosomal)

############ 3 circles 

# creating venn diagram by passing  
# list as direct parameter 
venn(list(normal_x=c(normal_genes_X),tumor_x=c(tumor_genes_X),tumor_and_normal_x=c(tumor_and_normal_genes_X)))
venn(list(normal_y=c(normal_genes_Y),tumor_y=c(tumor_genes_Y),tumor_and_normal_y=c(tumor_and_normal_genes_Y)))
venn(list(normal_autosomal=c(normal_genes_autosomal),tumor_autosomal=c(tumor_genes_autosomal),tumor_and_normal_autosomal=c(tumor_and_normal_genes_autosomal)))
