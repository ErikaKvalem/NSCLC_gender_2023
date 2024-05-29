#install.packages("VennDiagram")
require(data.table)

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



# Load the required library
library(VennDiagram)

# Define the sets for the Venn diagram
  normal <- normal_deg_merged$gene_id
  tumor <- tumor_deg_merged$gene_id
  tumor_and_normal <- tumor_and_normal_deg_merged$gene_id

# Create the Venn diagram
#venn.plot(list(set1 = normal, set2 = tumor, set3 = tumor_and_normal), category.names = c("normal", "tumor", "tumor_and_normal"))

# Create the Venn diagram
  venn_object <- venn.diagram(
  x = list(
    normal = normal,
    tumor = tumor,
    tumor_and_normal = tumor_and_normal
  ),
  category.names = c("normal", "tumor", "tumor_and_normal"),
  filename = NULL,  # To prevent saving as a file
  output = TRUE,    # To output the plot
  imagetype = "png" # Specify image type if saving as a file
)
# Display the plot
grid.newpage()
grid.draw(venn_object)
