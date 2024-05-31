
devtools::install_github("hrbrmstr/hrbrthemes")

resDir = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/figures/"
input_path = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/deseq2_out/"
chrom = read_csv("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/input/adata_var_nsclc_chrom.csv")
colnames(chrom)[2] <- "gene_id"
counts = read.table("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/input/counts.csv",sep= ",", header=TRUE, row.names=1)
samplesheet = read_csv("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/input/samplesheet.csv")
normal_deg = read_csv(paste0(input_path, "nsclc_gender_normal_sig_fc_genes_DESeq2_result.csv"))
tumor_deg = read_csv(paste0(input_path, "nsclc_gender_tumor_sig_fc_genes_DESeq2_result.csv"))
tumor_and_normal_deg = read_csv(paste0(input_path, "nsclc_gender_tumor_and_normal_sig_fc_genes_DESeq2_result.csv"))

# Assuming ch_df and normal are your dataframes
normal_deg_merged <- merge(normal_deg, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)
tumor_deg_merged <- merge(tumor_deg, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)
tumor_and_normal_deg_merged <- merge(tumor_and_normal_deg, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)

first_normal_deg_merged <- head(normal_deg_merged)
first_tumor_deg_merged <- head(tumor_deg_merged)
first_tumor_and_normal_deg_merged <- head(tumor_and_normal_deg_merged)



genes_of_interest = first_normal_deg_merged$gene_id
filtered_gene_counts <- counts[counts$gene_id %in% genes_of_interest, ]
rownames(filtered_gene_counts) <- filtered_gene_counts$gene_id # Set the values in the "gene_id" column as row names


transposed_gene_counts <- t(filtered_gene_counts) # Transpose the gene counts data
row_index <- which(rownames(transposed_gene_counts) == "gene_id") # Find the index of the row with the specified name
transposed_gene_counts <- transposed_gene_counts[-row_index, ] # Remove the row with the specified name

transposed_gene_counts<- as.data.frame(transposed_gene_counts)
sample_list <- rownames(transposed_gene_counts)
transposed_gene_counts$sample <- sample_list


# Reset row names
rownames(transposed_gene_counts) <- NULL
merged_data <- merge(transposed_gene_counts, samplesheet, by.x = "sample", by.y = "sample")

# Log normalization
log_norm <- function(x) {
  log(x + 1)  # Adding 1 to avoid log(0)
}

normalized_log <- as.data.frame(lapply(transposed_gene_counts, log_norm))
print("Log normalized dataframe:")

library(tidyverse)
library(hrbrthemes)
# Plot
# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)

# create a dataset
data <- data.frame(
  name=c( rep("A",500), rep("B",500), rep("B",500), rep("C",20), rep('D', 100)  ),
  value=c( rnorm(500, 10, 5), rnorm(500, 13, 1), rnorm(500, 18, 1), rnorm(20, 25, 4), rnorm(100, 12, 1) )
)

# Plot
merged_data %>%
  ggplot( aes(x=as.factor(origin), y=ENSG00000001497, fill=sex)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")


# A really basic boxplot.
ggplot(mtcars, aes(x=as.factor(cyl), y=mpg)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  xlab("cyl")


