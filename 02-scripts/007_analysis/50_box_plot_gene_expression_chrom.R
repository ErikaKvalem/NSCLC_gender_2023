library(readr)
library(ggplot2)
library(ggsignif)
library(forcats)

######################

resDir = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/figures/"
input_path = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/deseq2_out/"
chrom = read_csv("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/input/adata_var_nsclc_chrom.csv")
colnames(chrom)[2] <- "gene_id"


log1p_norm_counts = read.table("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/input/log1p_norm_counts.csv",sep= ",", header=TRUE, row.names=1)
samplesheet = read_csv("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/input/samplesheet.csv")

normal_deg = read_csv(paste0(input_path, "nsclc_gender_normal_sig_fc_genes_DESeq2_result.csv"))
tumor_deg = read_csv(paste0(input_path, "nsclc_gender_tumor_sig_fc_genes_DESeq2_result.csv"))
tumor_and_normal_deg = read_csv(paste0(input_path, "nsclc_gender_tumor_and_normal_sig_fc_genes_DESeq2_result.csv"))

normal_deg_all = read_csv(paste0(input_path, "nsclc_gender_normal_all_genes_DESeq2_result.csv"))
tumor_deg_all = read_csv(paste0(input_path, "nsclc_gender_tumor_all_genes_DESeq2_result.csv"))
tumor_and_normal_deg_all = read_csv(paste0(input_path, "nsclc_gender_tumor_and_normal_all_genes_DESeq2_result.csv"))





################################################################# CHROM
normal_deg_merged <- merge(normal_deg_all, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)
tumor_deg_merged <- merge(tumor_deg_all, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)
tumor_and_normal_deg_merged <- merge(tumor_and_normal_deg_all, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)

# Remove rows with NA values in a specific column
normal_deg_merged <- normal_deg_merged[complete.cases(normal_deg_merged$gene_name), ]
tumor_deg_merged <- tumor_deg_merged[complete.cases(tumor_deg_merged$gene_name), ]
tumor_and_normal_deg_merged <- tumor_and_normal_deg_merged[complete.cases(tumor_and_normal_deg_merged$gene_name), ]

normal_genes_X <- normal_deg_merged$gene_name[normal_deg_merged$chromosome_name == "X"]
#normal_genes_X <-na.omit(normal_genes_X)
normal_genes_Y <- normal_deg_merged$gene_name[normal_deg_merged$chromosome_name == "Y"]
#normal_genes_Y <-na.omit(normal_genes_Y)
normal_genes_autosomal <- normal_deg_merged$gene_name[!(normal_deg_merged$chromosome_name == "X" | normal_deg_merged$chromosome_name == "Y")]
#normal_genes_autosomal <- na.omit(normal_genes_autosomal)

tumor_genes_X <- tumor_deg_merged$gene_name[tumor_deg_merged$chromosome_name == "X"]
#tumor_genes_X <-na.omit(tumor_genes_X)
tumor_genes_Y <- tumor_deg_merged$gene_name[tumor_deg_merged$chromosome_name == "Y"]
#tumor_genes_Y <-na.omit(tumor_genes_Y)
tumor_genes_autosomal <- tumor_deg_merged$gene_name[!(tumor_deg_merged$chromosome_name == "X" | tumor_deg_merged$chromosome_name == "Y")]
#tumor_genes_autosomal <- na.omit(tumor_genes_autosomal)

tumor_and_normal_genes_X <- tumor_and_normal_deg_merged$gene_name[tumor_and_normal_deg_merged$chromosome_name == "X"]
#tumor_and_normal_genes_X <-na.omit(tumor_and_normal_genes_X)
tumor_and_normal_genes_Y <- tumor_and_normal_deg_merged$gene_name[tumor_and_normal_deg_merged$chromosome_name == "Y"]
#tumor_and_normal_genes_Y <-na.omit(tumor_and_normal_genes_Y)
tumor_and_normal_genes_autosomal <- tumor_and_normal_deg_merged$gene_name[!(tumor_and_normal_deg_merged$chromosome_name == "X" | tumor_and_normal_deg_merged$chromosome_name == "Y")]
#tumor_and_normal_genes_autosomal <- na.omit(tumor_and_normal_genes_autosomal)

################################################
# Convert log1p_norm_counts to dataframe
df_nc <- as.data.frame(log1p_norm_counts)

###################################################### X CHROM normal
# List of gene names

normal_genes_X_id <- normal_deg_merged$gene_id[normal_deg_merged$chromosome_name == "X"]
normal_genes_X_name <- normal_deg_merged$gene_name[normal_deg_merged$chromosome_name == "X"]

index_list <- normal_genes_X_id[1:6]
index_list_name <- normal_genes_X_name[1:6]

# Initialize an empty dataframe to store the result
resultx <- data.frame()
# Initialize an empty dataframe to store gene names and p-values
#p_values_df <- data.frame(gene_name = character(), p_val = numeric())

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
  #male_df <- subset(df_t, sex == "male")
  
  # Filter for female
  #female_df <- subset(df_t, sex == "female")
  
  # Calculate annotation
  #anno = t.test(male_df$log1p_norm, female_df$log1p_norm)
  #pval = anno$p.value
  #df_t$pval = as.numeric(pval)
  
  # Append the result
  resultx <- rbind(resultx, df_t)
  
}
resultx <- resultx %>%
  unite(gn_sex, gene_name, sex, remove=FALSE) 

resultx$condition ="normal_deg_x"
rx <- resultx %>%
  unite(gn_sex, gene_name, sex) %>%
  ggplot(aes(x = fct_inorder(gn_sex), y = log1p_norm, fill=gn_sex)) +
  geom_violin() +
  geom_violin(trim=FALSE)+
  scale_fill_manual(values=c("red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue"))+
  ggsignif::geom_signif(
    comparisons = list(
      unique(resultx$gn_sex)[1:2],
      unique(resultx$gn_sex)[3:4],
      unique(resultx$gn_sex)[5:6],
      unique(resultx$gn_sex)[7:8],
      unique(resultx$gn_sex)[9:10],
      unique(resultx$gn_sex)[11:12]),y_position = c(11,11,11,11,11,11)
  )+
  guides(fill = FALSE)+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab(NULL)

rx

ggsave(paste0(resDir,"head_normal_X_deg_violin_plot.jpg"), plot = rx, width = 8, height = 10)  # Adjust width and height as needed
###################################################### NORMAL Y CHROM

# List of gene names

normal_genes_y_id <- normal_deg_merged$gene_id[normal_deg_merged$chromosome_name == "Y"]
normal_genes_y_name <- normal_deg_merged$gene_name[normal_deg_merged$chromosome_name == "Y"]

index_list <- normal_genes_y_id[1:6]
index_list_name <- normal_genes_y_name[1:6]

# Initialize an empty dataframe to store the result
resulty <- data.frame()
# Initialize an empty dataframe to store gene names and p-values
#p_values_df <- data.frame(gene_name = character(), p_val = numeric())

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
  #male_df <- subset(df_t, sex == "male")
  
  # Filter for female
  #female_df <- subset(df_t, sex == "female")
  
  # Calculate annotation
  #anno = t.test(male_df$log1p_norm, female_df$log1p_norm)
  #pval = anno$p.value
  #df_t$pval = as.numeric(pval)
  
  # Append the result
  resulty<- rbind(resulty, df_t)
  
}
resulty <- resulty %>%
  unite(gn_sex, gene_name, sex, remove=FALSE) 

resulty$condition ="normal_deg_y"
ry<- resulty %>%
  unite(gn_sex, gene_name, sex) %>%
  ggplot(aes(x = fct_inorder(gn_sex), y = log1p_norm, fill=gn_sex)) +
  geom_violin() +
  geom_violin(trim=FALSE)+
  scale_fill_manual(values=c("red","blue")) +
  ggsignif::geom_signif(
    comparisons = list(
      unique(resulty$gn_sex)[1:2]
    ),y_position = c(11)
  )+
  guides(fill = FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab(NULL)

ry

ggsave(paste0(resDir,"head_normal_Y_deg_violin_plot.jpg"), plot = ry, width = 8, height = 10)  # Adjust width and height as needed
############################
############################
############################
###################################################### X CHROM tumor
# List of gene names

tumor_genes_X_id <- tumor_deg_merged$gene_id[tumor_deg_merged$chromosome_name == "X"]
tumor_genes_X_name <- tumor_deg_merged$gene_name[tumor_deg_merged$chromosome_name == "X"]

index_list <- tumor_genes_X_id[1:6]
index_list_name <- tumor_genes_X_name[1:6]

# Initialize an empty dataframe to store the result
resultxt <- data.frame()
# Initialize an empty dataframe to store gene names and p-values
#p_values_df <- data.frame(gene_name = character(), p_val = numeric())

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
  #male_df <- subset(df_t, sex == "male")
  
  # Filter for female
  #female_df <- subset(df_t, sex == "female")
  
  # Calculate annotation
  #anno = t.test(male_df$log1p_norm, female_df$log1p_norm)
  #pval = anno$p.value
  #df_t$pval = as.numeric(pval)
  
  # Append the result
  resultxt <- rbind(resultxt, df_t)
  
}
resultxt <- resultxt %>%
  unite(gn_sex, gene_name, sex, remove=FALSE) 

resultxt$condition ="tumor_deg_x"
rxt <- resultxt %>%
  unite(gn_sex, gene_name, sex) %>%
  ggplot(aes(x = fct_inorder(gn_sex), y = log1p_norm, fill=gn_sex)) +
  geom_violin() +
  geom_violin(trim=FALSE)+
  scale_fill_manual(values=c("red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue"))+
  ggsignif::geom_signif(
    comparisons = list(
      unique(resultxt$gn_sex)[1:2],
      unique(resultxt$gn_sex)[3:4],
      unique(resultxt$gn_sex)[5:6],
      unique(resultxt$gn_sex)[7:8],
      unique(resultxt$gn_sex)[9:10],
      unique(resultxt$gn_sex)[11:12]),y_position = c(11,11,11,11,11,11)
  )+
  guides(fill = FALSE)+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab(NULL)

rxt

ggsave(paste0(resDir,"head_tumor_X_deg_violin_plot.jpg"), plot = rxt, width = 8, height = 10)  # Adjust width and height as needed
###################################################### NORMAL Y tumor

# List of gene names

tumor_genes_y_id <- tumor_deg_merged$gene_id[tumor_deg_merged$chromosome_name == "Y"]
tumor_genes_y_name <- tumor_deg_merged$gene_name[tumor_deg_merged$chromosome_name == "Y"]

index_list <- tumor_genes_y_id
index_list_name <- tumor_genes_y_name

# Initialize an empty dataframe to store the result
resultyt <- data.frame()
# Initialize an empty dataframe to store gene names and p-values
#p_values_df <- data.frame(gene_name = character(), p_val = numeric())

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
  #df_t <- subset(df_t, origin == "tumor_primary")

  # Append the result
  resultyt<- rbind(resultyt, df_t)
  
}
resultyt <- resultyt %>%
  unite(gn_sex, gene_name, sex, remove=FALSE) 

resulty$condition ="normal_deg_y"
ryt<- resultyt %>%
  unite(gn_sex, gene_name, sex) %>%
  ggplot(aes(x = fct_inorder(gn_sex), y = log1p_norm, fill=gn_sex)) +
  geom_violin() +
  geom_violin(trim=FALSE)+
  scale_fill_manual(values=c("red","blue")) +
  ggsignif::geom_signif(
    comparisons = list(
      unique(resultyt$gn_sex)[1:2]
    ),y_position = c(11)
  )+
  guides(fill = FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab(NULL)

ryt

ggsave(paste0(resDir,"head_tumor_Y_deg_violin_plot.jpg"), plot = ryt, width = 8, height = 10)  # Adjust width and height as needed

########################################################################################## facet wrap normal and tumor x 

results_mergedx <- rbind(resultx,resultxt)

#NORMAL AND TUMOR 
pX <- results_mergedx %>%
  unite(gn_sex, gene_name, sex) %>%
  ggplot(aes(x = fct_inorder(gn_sex), y = log1p_norm, fill=gn_sex)) +
  geom_violin() +
  geom_violin(trim=FALSE)+
  scale_fill_manual(values=c("red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue","red","blue", "red","blue", "red","blue", "red","blue", "red","blue"))+
  ggsignif::geom_signif(
    comparisons = list(
      unique(results_mergedx$gn_sex)[1:2],
      unique(results_mergedx$gn_sex)[3:4],
      unique(results_mergedx$gn_sex)[5:6],
      unique(results_mergedx$gn_sex)[7:8],
      unique(results_mergedx$gn_sex)[9:10],
      unique(results_mergedx$gn_sex)[11:12]),
    y_position = c(11,11,11,11,11,11,11,11,11,11,11)
  )+
  guides(fill = FALSE)+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab(NULL)
# Use vars() to supply faceting variables:
pX <- pX + facet_wrap(vars(condition), scales = "free")
pX
ggsave(paste0(resDir,"head_normal_and_tumor_x_deg_violin_plot.jpg"), plot = pX, width = 8, height = 10)  # Adjust width and height as needed
####################################################################################### facet wrap normal and tumor y 
results_mergedy <- rbind(resulty,resultyt)

pY<- results_mergedy %>%
  unite(gn_sex, gene_name, sex) %>%
  ggplot(aes(x = fct_inorder(gn_sex), y = log1p_norm, fill=gn_sex)) +
  geom_violin() +
  geom_violin(trim=FALSE)+
  scale_fill_manual(values=c("red","blue")) +
  ggsignif::geom_signif(
    comparisons = list(
      unique(results_mergedy$gn_sex)[1:2]
    ),y_position = c(11)
  )+
  guides(fill = FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab(NULL)

# Use vars() to supply faceting variables:
pY <- pY + facet_wrap(vars(condition), scales = "free")
pY
ggsave(paste0(resDir,"head_normal_and_tumor_y_deg_violin_plot.jpg"), plot = pY, width = 8, height = 10)  # Adjust width and height as needed