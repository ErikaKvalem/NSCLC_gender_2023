library(readr)
library(ggplot2)
library(ggsignif)
library(forcats)
library(tidyr)
library(ggpubr)
library(rstatix)

library(dplyr)
# #####################
log1p_norm_counts = read.table("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/input/log1p_norm_counts.csv",sep= ",", header=TRUE, row.names=1)
samplesheet = read_csv("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/input/samplesheet.csv")
counts = read.table("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/input/counts.csv",sep= ",", header=TRUE, row.names=1)

resDir = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/figures/violin_plot/"
input_path = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/deseq2_out/corrected_ds/"
chrom = read_csv("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/input/adata_var_nsclc_chrom.csv")
colnames(chrom)[2] <- "gene_id"
normal_deg = read_csv(paste0(input_path, "nsclc_gender_normal_all_genes_DESeq2_result.csv"))
tumor_deg = read_csv(paste0(input_path, "nsclc_gender_tumor_all_genes_DESeq2_result.csv"))
tumor_and_normal_deg = read_csv(paste0(input_path, "nsclc_gender_all_all_genes_DESeq2_result.csv"))


# Convert log1p_norm_counts to dataframe
df_nc <- as.data.frame(counts)
#################################################### NORMAL
# List of gene names
index_list <- head(normal_deg$gene_id,8)
index_list_name <- head(normal_deg$gene_name,8)

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
  
  df_t$log1p_norm <- scale(df_t$log1p_norm, center = TRUE, scale = TRUE)
  
  # Apply the log(x + 1) transformation
  df_t$log1p_norm <- log1p(df_t$log1p_norm)
  
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
  #anno = t.test(male_df$log1p_norm, female_df$log1p_norm)
  anno = wilcox.test(male_df$log1p_norm, female_df$log1p_norm, alternative = "two.sided",  correct = TRUE)
  padjval = anno$p.value
  df_t$padjval = as.numeric(padjval)
  
  # Append the result
  result <- rbind(result, df_t)

}

######################################################## TUMOR
# List of gene names
index_list_tumor <- head(tumor_deg$gene_id,8)
index_list_name_tumor <- head(tumor_deg$gene_name,8)

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
  df_t$log1p_norm <- scale(df_t$log1p_norm, center = TRUE, scale = TRUE)
  
  # Apply the log(x + 1) transformation
  df_t$log1p_norm <- log1p(df_t$log1p_norm)
  
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
  #anno = t.test(male_df$log1p_norm, female_df$log1p_norm)
  anno = wilcox.test(male_df$log1p_norm, female_df$log1p_norm, alternative = "two.sided",  correct = TRUE)
  padjval = anno$p.value
  df_t$padjval = as.numeric(padjval)
  
  # Append the result
  result2 <- rbind(result2, df_t)

}



####################
result <- result %>%
  unite(gn_sex, gene_name, sex, remove=FALSE) 

result2 <- result2 %>%
  unite(gn_sex, gene_name, sex, remove=FALSE) 
result2$dataset <- substr(result2$sample, 1, 10)
result$dataset <- substr(result$sample, 1, 10)

result$condition ="normal_deg"
result2$condition ="tumor_deg"

results_merged <- rbind(result,result2)


############################################# GGPUBR

#conflicts_prefer(rstatix::filter)
## Perform Wilcoxon tests and collect p-values
test_results <- result %>%
  group_by(gene_name) %>%
  wilcox_test(log1p_norm ~ sex) %>%
  ungroup()

# Calculate the covariate (mean expression level) for each gene
covariate <- result %>%
  group_by(gene_name) %>%
  summarize(mean_expression = mean(log1p_norm)) %>%
  pull(mean_expression)

# Perform IHW adjustment
ihw_result <- ihw(test_results$p, covariate, alpha = 0.1)

# Extract adjusted p-values
test_results$p.adj <- adj_pvalues(ihw_result)
stat.test <- test_results |> add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.1, 1),
                 symbols = c("****", "***", "**", "*", "ns"))

stat.test$p.scient <- format(stat.test$p.adj, scientific = TRUE, digits = 3)

stat.test <- stat.test %>% add_xy_position(x = "sex")
stat.test.normal <- stat.test

#stat.test <- result |>
#  group_by(gene_name)  |>
#  wilcox_test(log1p_norm ~ sex) |>
#  adjust_pvalue(method = "fdr") |>
#  add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.1, 1),
#                   symbols = c("****", "***", "**", "*", "ns"))
#stat.test$p.scient <- format(stat.test$p.adj, scientific = TRUE, digits = 3)
#
#stat.test <- stat.test %>% add_xy_position(x = "sex")
#stat.test.normal <- stat.test

p  <-  ggviolin(result, x = "sex", y = "log1p_norm", trim=FALSE, title="Top 20 Normal DE genes corrected") +  
  # scale_color_manual(values = c("male" = "blue", "female" = "red")) +  # Customize colors for "sex"
  facet_wrap(~ gene_name) +stat_pvalue_manual(
    stat.test, bracket.nudge.y = -1, hide.ns = FALSE,
    label = "{p.signif}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  geom_jitter(height = 0, width = 0.1, aes(colour = dataset)) 
p



ggsave(paste0("/home/kvalem/myScratch/projects/2023/NSCLC_gender_2023/05-results/010_analysis_paired_include_guon/","normal_deg_violin_plot_corrected_all_ihw_top8.jpg"), plot = p, width = 10, height = 10)  # Adjust width and height as needed

test_results <- result2 %>%
  group_by(gene_name) %>%
  wilcox_test(log1p_norm ~ sex) %>%
  ungroup()

# Calculate the covariate (mean expression level) for each gene
covariate <- result2 %>%
  group_by(gene_name) %>%
  summarize(mean_expression = mean(log1p_norm)) %>%
  pull(mean_expression)

# Perform IHW adjustment
ihw_result <- ihw(test_results$p, covariate, alpha = 0.1)

# Extract adjusted p-values
test_results$p.adj <- adj_pvalues(ihw_result)
stat.test <- test_results |> add_significance("p", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.1, 1),
                                              symbols = c("****", "***", "**", "*", "ns"))

stat.test$p.scient <- format(stat.test$p.adj, scientific = TRUE, digits = 3)

stat.test <- stat.test %>% add_xy_position(x = "sex")
stat.test.tumor <- stat.test


#stat.test <- result2 |>
#  group_by(gene_name)  |>
#  wilcox_test(log1p_norm ~ sex) |>
#  adjust_pvalue(method ="fdr") |>
#  add_significance("p.adj", cutpoints = c(0, 1e-04, 0.001, 0.01, 0.1, 1),
#                   symbols = c("****", "***", "**", "*", "ns"))
#stat.test$p.scient <- format(stat.test$p.adj, scientific = TRUE, digits = 3)
#
#stat.test <- stat.test %>% add_xy_position(x = "sex")
#stat.test.tumor <- stat.test

p2<-  ggviolin(result2, x = "sex", y = "log1p_norm", trim=FALSE,  title="Top 20 Tumor DE genes corrected") +  
  # scale_color_manual(values = c("male" = "blue", "female" = "red")) +  # Customize colors for "sex"
  facet_wrap(~ gene_name) +stat_pvalue_manual(
    stat.test, bracket.nudge.y = -1, hide.ns = FALSE,
    label =  "{p.signif}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  geom_jitter(height = 0, width = 0.1, aes(colour = dataset)) 
p2

stat.test.df <- as.data.frame(stat.test)
stat.test.df  %>% filter(row(gene_name) == "RPS4Y1")
p2<-  ggviolin(result2, x = "sex", y = "log1p_norm", trim=FALSE,  title="Top 20 Tumor DE genes corrected") +  
  # scale_color_manual(values = c("male" = "blue", "female" = "red")) +  # Customize colors for "sex"
  #facet_wrap(~ gene_name) +
  stat_pvalue_manual(
    stat.test, bracket.nudge.y = -1, hide.ns = FALSE,
    label =  "{p.signif}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  geom_jitter(height = 0, width = 0.1, aes(colour = dataset)) 
p2
#ggsave(paste0("/home/kvalem/myScratch/projects/2023/NSCLC_gender_2023/05-results/010_analysis_paired_include_guon/","tumor_deg_violin_plot_corrected_all_ihw_top8.jpg"), plot = p2, width = 10, height = 10)  # Adjust width and height as needed

#####################################

my_comparisons <- list(c("male", "female"))
result2 <- result2 %>%
  filter(grepl("MAL2", gene_name))

ggplot(df_t, aes(x=sex,y=gene_name, fill=sex))+
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.1,fill="white")+
    theme(legend.position="none")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 25))+
    theme(axis.text.y = element_text(size = 25),
        axis.title=element_text(size=25))+
    theme(axis.title.x = element_blank())+
   # ylab('PD-L1')+
    scale_fill_manual(values=c("#D11141", "#414487FF"))+
    stat_compare_means(comparisons = my_comparisons, label.y = 6,size=8,label ="p.signif") +
    stat_compare_means(label = "p.format",size=7,label.y = 6.35,label.x = 1.05 )#+
    #ylim(0,7)




