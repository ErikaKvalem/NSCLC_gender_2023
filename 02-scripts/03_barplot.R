# library
library(ggplot2)
resDir = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/figures/"
input_path = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/deseq2_out/"
chrom = read_csv("/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/input/adata_var_nsclc_chrom.csv")
colnames(chrom)[2] <- "gene_id"
  

normal_deg = read_csv(paste0(input_path, "nsclc_gender_normal_all_genes_DESeq2_result.csv"))
tumor_deg = read_csv(paste0(input_path, "nsclc_gender_tumor_all_genes_DESeq2_result.csv"))
tumor_and_normal_deg= read_csv(paste0(input_path, "nsclc_gender_tumor_and_normal_all_genes_DESeq2_result.csv"))

normal_deg_sig = read_csv(paste0(input_path, "nsclc_gender_normal_sig_genes_DESeq2_result.csv"))
tumor_deg_sig = read_csv(paste0(input_path, "nsclc_gender_tumor_sig_genes_DESeq2_result.csv"))
tumor_and_normal_deg_sig = read_csv(paste0(input_path, "nsclc_gender_tumor_and_normal_sig_genes_DESeq2_result.csv"))

normal_deg_sig_fc = read_csv(paste0(input_path, "nsclc_gender_normal_sig_fc_genes_DESeq2_result.csv"))
tumor_deg_sig_fc = read_csv(paste0(input_path, "nsclc_gender_tumor_sig_fc_genes_DESeq2_result.csv"))
tumor_and_normal_deg_sig_fc = read_csv(paste0(input_path, "nsclc_gender_tumor_and_normal_sig_fc_genes_DESeq2_result.csv"))

#################################### ALL
normal_deg_merged <- merge(normal_deg, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)
tumor_deg_merged <- merge(tumor_deg, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)
tumor_and_normal_deg_merged <- merge(tumor_and_normal_deg, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)


normal_genes_X <- normal_deg_merged$gene_id[normal_deg_merged$chromosome_name == "X"]
normal_genes_Y <- normal_deg_merged$gene_id[normal_deg_merged$chromosome_name == "Y"]
normal_genes_autosomal <- normal_deg_merged$gene_id[!(normal_deg_merged$chromosome_name == "X" | normal_deg_merged$chromosome_name == "Y")]

tumor_genes_X <- tumor_deg_merged$gene_id[tumor_deg_merged$chromosome_name == "X"]
tumor_genes_Y <- tumor_deg_merged$gene_id[tumor_deg_merged$chromosome_name == "Y"]
tumor_genes_autosomal <- tumor_deg_merged$gene_id[!(tumor_deg_merged$chromosome_name == "X" | tumor_deg_merged$chromosome_name == "Y")]

tumor_and_normal_genes_X <- tumor_and_normal_deg_merged$gene_id[tumor_and_normal_deg_merged$chromosome_name == "X"]
tumor_and_normal_genes_Y <- tumor_and_normal_deg_merged$gene_id[tumor_and_normal_deg_merged$chromosome_name == "Y"]
tumor_and_normal_genes_autosomal <- tumor_and_normal_deg_merged$gene_id[!(tumor_and_normal_deg_merged$chromosome_name == "X" | tumor_and_normal_deg_merged$chromosome_name == "Y")]


origin <- c(rep("normal" , 3) , rep("tumor" , 3) , rep("tumor_and_normal" , 3))
condition <- rep(c("DEG X chrom" , "DEG Y chrom","DEG autosome") , 3)
n_deg <- c(length(normal_genes_X),length(normal_genes_Y), length(normal_genes_autosomal),length(tumor_genes_X),length(tumor_genes_Y), length(tumor_genes_autosomal),length(tumor_and_normal_genes_X),length(tumor_and_normal_genes_Y), length(tumor_and_normal_genes_autosomal))
data <- data.frame(origin,condition,n_deg)

# Stacked
p <- ggplot(data, aes(fill=condition, y=n_deg, x=origin)) + 
  geom_bar(position="dodge", stat="identity") 
ggsave(filename = paste0(resDir,"deg_stacked_plot.png"), plot = p)
# Adding text labels on top of each bar
p <- p + geom_text(aes(label=n_deg), position=position_dodge(width=0.9), vjust=-0.5, size=3)




#ggsave( paste0(resDir,"deg_stacked_plot.png"), plot = p, width = 20, height = 10) 
########################################SIG

normal_deg_merged <- merge(normal_deg_sig, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)
tumor_deg_merged <- merge(tumor_deg_sig, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)
tumor_and_normal_deg_merged <- merge(tumor_and_normal_deg_sig, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)


normal_genes_X <- normal_deg_merged$gene_id[normal_deg_merged$chromosome_name == "X"]
normal_genes_Y <- normal_deg_merged$gene_id[normal_deg_merged$chromosome_name == "Y"]
normal_genes_autosomal <- normal_deg_merged$gene_id[!(normal_deg_merged$chromosome_name == "X" | normal_deg_merged$chromosome_name == "Y")]

tumor_genes_X <- tumor_deg_merged$gene_id[tumor_deg_merged$chromosome_name == "X"]
tumor_genes_Y <- tumor_deg_merged$gene_id[tumor_deg_merged$chromosome_name == "Y"]
tumor_genes_autosomal <- tumor_deg_merged$gene_id[!(tumor_deg_merged$chromosome_name == "X" | tumor_deg_merged$chromosome_name == "Y")]

tumor_and_normal_genes_X <- tumor_and_normal_deg_merged$gene_id[tumor_and_normal_deg_merged$chromosome_name == "X"]
tumor_and_normal_genes_Y <- tumor_and_normal_deg_merged$gene_id[tumor_and_normal_deg_merged$chromosome_name == "Y"]
tumor_and_normal_genes_autosomal <- tumor_and_normal_deg_merged$gene_id[!(tumor_and_normal_deg_merged$chromosome_name == "X" | tumor_and_normal_deg_merged$chromosome_name == "Y")]

origin <- c(rep("normal" , 3) , rep("tumor" , 3) , rep("tumor_and_normal" , 3))
condition <- rep(c("DEG X chrom" , "DEG Y chrom","DEG autosome") , 3)
n_deg <- c(length(normal_genes_X),length(normal_genes_Y), length(normal_genes_autosomal),length(tumor_genes_X),length(tumor_genes_Y), length(tumor_genes_autosomal),length(tumor_and_normal_genes_X),length(tumor_and_normal_genes_Y), length(tumor_and_normal_genes_autosomal))
data <- data.frame(origin,condition,n_deg)

# Stacked
p2 <- ggplot(data, aes(fill=condition, y=n_deg, x=origin)) + 
  geom_bar(position="dodge", stat="identity") 
ggsave(filename = paste0(resDir,"deg_stacked_plot.png"), plot = p)
# Adding text labels on top of each bar
p2 <- p2 + geom_text(aes(label=n_deg), position=position_dodge(width=0.9), vjust=-0.5, size=3)




ggsave( paste0(resDir,"deg_sig_stacked_plot.png"), plot = p2, width = 20, height = 10) 
####################################### 

normal_deg_merged <- merge(normal_deg_sig_fc, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)
tumor_deg_merged <- merge(tumor_deg_sig_fc, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)
tumor_and_normal_deg_merged <- merge(tumor_and_normal_deg_sig_fc, chrom[, c("gene_id", "chromosome_name")], by = "gene_id", all.x = TRUE)


normal_genes_X <- normal_deg_merged$gene_id[normal_deg_merged$chromosome_name == "X"]
normal_genes_Y <- normal_deg_merged$gene_id[normal_deg_merged$chromosome_name == "Y"]
normal_genes_autosomal <- normal_deg_merged$gene_id[!(normal_deg_merged$chromosome_name == "X" | normal_deg_merged$chromosome_name == "Y")]

tumor_genes_X <- tumor_deg_merged$gene_id[tumor_deg_merged$chromosome_name == "X"]
tumor_genes_Y <- tumor_deg_merged$gene_id[tumor_deg_merged$chromosome_name == "Y"]
tumor_genes_autosomal <- tumor_deg_merged$gene_id[!(tumor_deg_merged$chromosome_name == "X" | tumor_deg_merged$chromosome_name == "Y")]

tumor_and_normal_genes_X <- tumor_and_normal_deg_merged$gene_id[tumor_and_normal_deg_merged$chromosome_name == "X"]
tumor_and_normal_genes_Y <- tumor_and_normal_deg_merged$gene_id[tumor_and_normal_deg_merged$chromosome_name == "Y"]
tumor_and_normal_genes_autosomal <- tumor_and_normal_deg_merged$gene_id[!(tumor_and_normal_deg_merged$chromosome_name == "X" | tumor_and_normal_deg_merged$chromosome_name == "Y")]

origin <- c(rep("normal" , 3) , rep("tumor" , 3) , rep("tumor_and_normal" , 3))
condition <- rep(c("DEG X chrom" , "DEG Y chrom","DEG autosome") , 3)
n_deg <- c(length(normal_genes_X),length(normal_genes_Y), length(normal_genes_autosomal),length(tumor_genes_X),length(tumor_genes_Y), length(tumor_genes_autosomal),length(tumor_and_normal_genes_X),length(tumor_and_normal_genes_Y), length(tumor_and_normal_genes_autosomal))
data <- data.frame(origin,condition,n_deg)

# Stacked
p3<- ggplot(data, aes(fill=condition, y=n_deg, x=origin)) + 
  geom_bar(position="dodge", stat="identity") 
ggsave(filename = paste0(resDir,"deg_stacked_plot.png"), plot = p)
# Adding text labels on top of each bar
p3 <- p3 + geom_text(aes(label=n_deg), position=position_dodge(width=0.9), vjust=-0.5, size=3)




ggsave( paste0(resDir,"deg_sig_fc_stacked_plot.png"), plot = p3, width = 20, height = 10) 

