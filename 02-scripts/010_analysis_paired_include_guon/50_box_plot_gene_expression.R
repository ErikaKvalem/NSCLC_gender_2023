library(readr)
library(ggplot2)
library(ggsignif)
library(forcats)
library(tidyr)

# #####################

log1p_norm_counts = read.table("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/input/log1p_norm_counts.csv",sep= ",", header=TRUE, row.names=1)
samplesheet = read_csv("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/input/samplesheet.csv")

resDir = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/figures/violin_plot/"
input_path = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/deseq2_out/corrected_ds/"
chrom = read_csv("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/input/adata_var_nsclc_chrom.csv")
colnames(chrom)[2] <- "gene_id"
normal_deg = read_csv(paste0(input_path, "nsclc_gender_normal_sig_fc_genes_DESeq2_result.csv"))
tumor_deg = read_csv(paste0(input_path, "nsclc_gender_tumor_sig_fc_genes_DESeq2_result.csv"))
tumor_and_normal_deg = read_csv(paste0(input_path, "nsclc_gender_all_sig_fc_genes_DESeq2_result.csv"))


# Convert log1p_norm_counts to dataframe
df_nc <- as.data.frame(log1p_norm_counts)
#################################################### NORMAL
# List of gene names
index_list <- head(normal_deg$gene_id,15)
index_list_name <- head(normal_deg$gene_name,15)

# Initialize an empty dataframe to store the result
result <- data.frame()
# Initialize an empty dataframe to store gene names and p-values
p_values_df <- data.frame(gene_name = character(), p_val = numeric())

# Iterate over each gene name in index_list
for (i in seq_along(index_list)) {
  # Filter dataframe based on gene name
  df_nc_filtered <- df_nc[index_list[i], ]
  
  # Transpose the filtered dataframe
  df_t <- t(df_nc_filtered)
  
  # Convert transposed dataframe to data frame
  df_t <- as.data.frame(df_t)
  
  # Rename the column to 'log1p_norm'
  colnames(df_t)[colnames(df_t) == index_list[i]] <- 'log1p_norm'
  
  # Add 'gene_name' column
  df_t$gene_id <- as.character(index_list[i])
  df_t$gene_name <- as.character(index_list_name[i])
  df_t$sample <- row.names(df_t)
  df_t <- merge(df_t, samplesheet[, c("sample", "sex", "origin")], by = "sample")
  # Filter for male
  male_df <- subset(df_t, sex == "male")
  
  # Filter for female
  female_df <- subset(df_t, sex == "female")
  
  # Calculate annotation
  anno = t.test(male_df$log1p_norm, female_df$log1p_norm)
  pval = anno$p.value
  df_t$pval = as.numeric(pval)
  
  # Append the result
  result <- rbind(result, df_t)

}

######################################################## TUMOR
# List of gene names
index_list_tumor <- head(tumor_deg$gene_id,15)
index_list_name_tumor <- head(tumor_deg$gene_name,30)

# Initialize an empty dataframe to store the result
result2 <- data.frame()

# Iterate over each gene name in index_list
for (i in seq_along(index_list_tumor)) {
  # Filter dataframe based on gene name
  df_nc_filtered <- df_nc[index_list_tumor[i], ]
  
  # Transpose the filtered dataframe
  df_t <- t(df_nc_filtered)
  
  # Convert transposed dataframe to data frame
  df_t <- as.data.frame(df_t)
  
  # Rename the column to 'log1p_norm'
  colnames(df_t)[colnames(df_t) == index_list_tumor[i]] <- 'log1p_norm'
  
  # Add 'gene_name' column
  df_t$gene_id <- as.character(index_list_tumor[i])
  df_t$gene_name <- as.character(index_list_name_tumor[i])
  df_t$sample <- row.names(df_t)
  df_t <- merge(df_t, samplesheet[, c("sample", "sex", "origin")], by = "sample")
  # Filter for male
  male_df <- subset(df_t, sex == "male")
  
  # Filter for female
  female_df <- subset(df_t, sex == "female")
  
  # Calculate annotation
  anno = t.test(male_df$log1p_norm, female_df$log1p_norm)
  pval = anno$p.value
  df_t$pval = as.numeric(pval)
  
  # Append the result
  result2 <- rbind(result2, df_t)

}


####################
result <- result %>%
  unite(gn_sex, gene_name, sex, remove=FALSE) 

result2 <- result2 %>%
  unite(gn_sex, gene_name, sex, remove=FALSE) 
result2$dataset <- substr(result2$sample, 1, 5)

result$condition ="normal_deg"
result2$condition ="tumor_deg"

results_merged <- rbind(result,result2)
####################


r1 <- result %>%
  unite(gn_sex, gene_name, sex) %>%
  ggplot(aes(x = fct_inorder(gn_sex), y = log1p_norm, fill=gn_sex)) +
  geom_violin() +
  geom_violin(trim=FALSE)+
  scale_fill_manual(values=c("red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue"))+
  guides(fill = FALSE)+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab(NULL)

r1

# Save the plot as an image file (e.g., PNG)
ggsave(paste0(resDir,"head_normal_deg_violin_plot.jpg"), plot = r1, width = 8, height = 10)  # Adjust width and height as needed


r2 <- result2 %>%
  unite(gn_sex, gene_name, sex) %>%
  ggplot(aes(x = fct_inorder(gn_sex), y = log1p_norm, fill=gn_sex)) +
  geom_violin() +
  geom_point(position = position_jitterdodge(seed = 1, dodge.width = 0.9))+
 # scale_color_manual(values = c("Govei" = "pink", "Guo_Z" = "green","He_Fa"="purple", "Kim_L"="brown","Lambr"="red","Leade":"yellow" )) +  # Add a manual color scale for datasets
  
  #Govei
  geom_violin(trim=FALSE)+
  scale_fill_manual(values=c("red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue"))+
ggsignif::geom_signif(
    comparisons = list(
      unique(result2$gn_sex)[1:2],
      unique(result2$gn_sex)[3:4],
      unique(result2$gn_sex)[5:6],
      unique(result2$gn_sex)[7:8],
      unique(result2$gn_sex)[9:10],
      unique(result2$gn_sex)[11:12],
      unique(result2$gn_sex)[13:14],
      unique(result2$gn_sex)[15:16]),y_position = c(11,11,11,11,11,11,11,11)
  )+
  guides(fill = FALSE)+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab(NULL)
r2
ggsave(paste0(resDir,"head_tumor_deg_violin_plot.jpg"), plot = r2, width = 8, height = 10)  # Adjust width and height as needed

# ###################################



p <- results_merged %>%
  unite(gn_sex, gene_name, sex) %>%
  ggplot(aes(x = fct_inorder(gn_sex), y = log1p_norm, fill=gn_sex)) +
  geom_violin() +
  geom_violin(trim=FALSE)+
  scale_fill_manual(values=c("red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue","red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue"))+
  #ggsignif::geom_signif(
  #  comparisons = list(
  #    unique(result$gn_sex)[1:2],
  #    unique(result$gn_sex)[3:4],
  #    unique(result$gn_sex)[5:6],
  #    unique(result$gn_sex)[7:8],
  #    unique(result$gn_sex)[9:10],
  #    unique(result$gn_sex)[11:12],c("AKR1C3_male","AKR1C3_female"),
  #    c("KRT19_male","KRT19_female"),
  #    c("KRT8_male","KRT8_female"),
  #    c("KRT7_male","KRT7_female"),
  #    c("KRT17_male","KRT17_female")),
  #  y_position = c(11,11,11,11,11,11,11,11,11,11,11)
  #)+
  guides(fill = FALSE)+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab(NULL)
# Use vars() to supply faceting variables:
p <- p + facet_wrap(vars(condition), scales = "free_y")

p
#ggsave(paste0(resDir,"head_tumor_and_normal_deg_violin_plot.jpg"), plot = p, width = 20, height = 10)  # Adjust width and height as needed
