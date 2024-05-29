# library
library(ggplot2)

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
normal_genes_Y <- normal_deg_merged$gene_id[normal_deg_merged$chromosome_name == "Y"]
normal_genes_autosomal <- normal_deg_merged$gene_id[!(normal_deg_merged$chromosome_name == "X" | normal_deg_merged$chromosome_name == "Y")]

tumor_genes_X <- tumor_deg_merged$gene_id[tumor_deg_merged$chromosome_name == "X"]
tumor_genes_Y <- tumor_deg_merged$gene_id[tumor_deg_merged$chromosome_name == "Y"]
tumor_genes_autosomal <- tumor_deg_merged$gene_id[!(tumor_deg_merged$chromosome_name == "X" | tumor_deg_merged$chromosome_name == "Y")]

tumor_and_normal_genes_X <- tumor_and_normal_deg_merged$gene_id[tumor_and_normal_deg_merged$chromosome_name == "X"]
tumor_and_normal_genes_Y <- tumor_and_normal_deg_merged$gene_id[tumor_and_normal_deg_merged$chromosome_name == "Y"]
tumor_and_normal_genes_autosomal <- tumor_and_normal_deg_merged$gene_id[!(tumor_and_normal_deg_merged$chromosome_name == "X" | tumor_and_normal_deg_merged$chromosome_name == "Y")]


# create a dataset
origin <- c(rep("normal" , 3) , rep("tumor" , 3) , rep("tumor_and_normal" , 3))
condition <- rep(c("DEG X chrom" , "DEG Y chrom","DEG autosome") , 3)
n_deg <- c(length(normal_genes_X),length(normal_genes_Y), length(normal_genes_autosomal),length(tumor_genes_X),length(tumor_genes_Y), length(tumor_genes_autosomal),length(tumor_and_normal_genes_X),length(tumor_and_normal_genes_Y), length(tumor_and_normal_genes_autosomal))
data <- data.frame(origin,condition,n_deg)

# Stacked
ggplot(data, aes(fill=condition, y=n_deg, x=origin)) + 
  geom_bar(position="dodge", stat="identity")
