library("conflicted")
library("docopt")
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
remove_ensg_version = function(x) gsub("\\.[0-9]*$", "", x)
library("enrichplot")
organism = "human"
gene_id_type = "ENSEMBL"
if (organism == "human") {
  anno_db = "org.Hs.eg.db"
  org_kegg = "hsa"
  org_reactome = "human"
  org_wp = "Homo sapiens"
} else if (organism == "mouse") {
  anno_db = "org.Mm.eg.db"
  org_kegg = "mmu"
  org_reactome = "mouse"
  org_wp = "Mus musculus"
} else {
  msg <- paste0("Organism not implemented: ", organism)
  stop(msg)
}
library(anno_db, character.only = TRUE)


require(data.table)
################### INTESTINE 
############################################################################### 
resDir = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/deseq2_out/corrected_ds/GSEA/ORA/" 


inputDir = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/deseq2_out/corrected_ds/GSEA"
data_normal<-as.data.frame(fread(paste0(inputDir,"/nsclc_gender_normal_GSEA_KEGG.tsv")))
data_tumor <- as.data.frame(fread(paste0(inputDir,"/nsclc_gender_tumor_GSEA_KEGG.tsv")))
data_tumor_and_normal <- as.data.frame(fread(paste0(inputDir,"/nsclc_gender_tumor_and_normal_GSEA_KEGG.tsv")))
data_normal$condition = "normal"
data_tumor$condition = "tumor"
data_tumor_and_normal$condition = "tumor_normal"



#data_normal <- data_normal %>%
#  filter(p.adjust < 0.05)
#data_tumor <- data_tumor %>%
#  filter(p.adjust < 0.05)
#data_tumor_and_normal <- data_tumor_and_normal %>%
#  filter(p.adjust < 0.05)

total <- rbind(head(data_normal,10),head(data_tumor,10),head(data_tumor_and_normal,10))
total$Description <- gsub("- M.*", "", total$Description)

p <- ggplot(total, aes(NES,Description,fill = NES)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free", space = "free") +  
  labs(y = "KEGG Pathways", color = "black") # Change the y-axis label to "Pathways"
p2 <- ggplot(total, aes(NES,Description,fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free", space = "free") +  
  labs(y = "KEGG Pathways", color = "black") # Change the y-axis label to "Pathways"


# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
p2 <- p2 + scale_fill_gradient(low = "blue", high = "red")
# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired
p2 <- p2 + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired
# Print the plot
print(p)
print(p2)

# Save the plot as an image file (e.g., PNG)
ggsave(paste0(resDir, "KEGG_pathways_facet_grid.png"), plot = p, width = 8, height = 10)  # Adjust width and height as needed

ggsave(paste0(resDir,"KEGG_pathways_padjust_facet_grid.png"), plot = p2, width = 8, height = 10)  # Adjust width and height as needed


###################################################################
data_normal<-as.data.frame(fread(paste0(inputDir,"/nsclc_gender_normal_ORA_GO_BP.tsv")))
data_tumor <- as.data.frame(fread(paste0(inputDir,"/nsclc_gender_tumor_ORA_GO_BP.tsv")))
data_tumor_and_normal <- as.data.frame(fread(paste0(inputDir,"/nsclc_gender_tumor_and_normal_ORA_GO_BP.tsv")))
data_normal$condition = "normal"
data_tumor$condition = "tumor"
data_tumor_and_normal$condition = "tumor_normal"

#data_normal <- data_normal %>%
#  filter(p.adjust < 0.05)
#data_tumor <- data_tumor %>%
#  filter(p.adjust < 0.05)
#data_tumor_and_normal <- data_tumor_and_normal %>%
#  filter(p.adjust < 0.05)



total <- rbind(head(data_normal,10),head(data_tumor,10),head(data_tumor_and_normal,10))
total$Description <- gsub("- M.*", "", total$Description)

p <- ggplot(total, aes(Count,Description,fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free", space = "free") +  
  labs(y = "ORA GO BP", color = "black") # Change the y-axis label to "Pathways"


# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired

# Print the plot
print(p)

# Save the plot as an image file (e.g., PNG)
ggsave(paste0(resDir,"ORA_GO_BP_facet_grid.png"), plot = p, width = 8, height = 10)  # Adjust width and height as needed

########################################################
data_normal<-as.data.frame(fread(paste0(inputDir,"/nsclc_gender_normal_ORA_GO_CC.tsv")))
data_tumor <- as.data.frame(fread(paste0(inputDir,"/nsclc_gender_tumor_ORA_GO_CC.tsv")))
data_tumor_and_normal <- as.data.frame(fread(paste0(inputDir,"/nsclc_gender_tumor_and_normal_ORA_GO_CC.tsv")))
data_normal$condition = "normal"
data_tumor$condition = "tumor"
data_tumor_and_normal$condition = "tumor_normal"



total <- rbind(head(data_normal,10),head(data_tumor,10),head(data_tumor_and_normal,10))
total$Description <- gsub("- M.*", "", total$Description)

p <- ggplot(total, aes(Count,Description,fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free", space = "free") +  
  labs(y = "ORA GO CC", color = "black") # Change the y-axis label to "Pathways"


# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired

# Print the plot
print(p)

# Save the plot as an image file (e.g., PNG)
ggsave(paste0(resDir,"ORA_GO_CC_facet_grid.png"), plot = p, width = 8, height = 10)  # Adjust width and height as needed

################################################### CANCER 

resDir = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/deseq2_out/corrected_ds/GSEA/ORA/" 


inputDir = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/deseq2_out/corrected_ds/GSEA"
data_normal<-as.data.frame(fread(paste0(inputDir,"/nsclc_gender_normal_GSEA_KEGG.tsv")))
data_tumor <- as.data.frame(fread(paste0(inputDir,"/nsclc_gender_tumor_GSEA_KEGG.tsv")))
data_tumor_and_normal <- as.data.frame(fread(paste0(inputDir,"/nsclc_gender_tumor_and_normal_GSEA_KEGG.tsv")))
data_normal$condition = "normal"
data_tumor$condition = "tumor"
data_tumor_and_normal$condition = "tumor_normal"



total <- rbind(data_normal,data_tumor,data_tumor_and_normal)


subset_df <- total[grepl("cancer", total$Description, ignore.case = TRUE), ]

# p plots counts and padjust 
p <- ggplot(subset_df, aes(NES,Description,fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free", space = "free") +  
  labs(y = "KEGG pathways", color = "black") # Change the y-axis label to "Pathways"

p2 <- ggplot(subset_df, aes(NES,Description,fill = NES)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free", space = "free") +  
  labs(y = "KEGG pathways", color = "black") # Change the y-axis label to "Pathways"



# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
p2 <- p2 + scale_fill_gradient(low = "blue", high = "red")

# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired
p2 <- p2 + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired
# Print the plot
print(p)
print(p2)
# Save the plot as an image file (e.g., PNG)
ggsave(paste0(resDir, "KEGG_pathways_facet_grid_cancer_nes_padj.png"), plot = p, width = 8, height = 10)  # Adjust width and height as needed
ggsave(paste0(resDir, "KEGG_pathways_facet_grid_cancer_nes.png"), plot = p2, width = 8, height = 10)  # Adjust width and height as needed

############################################# GSEA GO  BP 
resDir = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/deseq2_out/corrected_ds/GSEA/ORA/" 


inputDir = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/deseq2_out/corrected_ds/GSEA"
data_normal<-as.data.frame(fread(paste0(inputDir,"/nsclc_gender_normal_GSEA_GO_BP.tsv")))
data_tumor <- as.data.frame(fread(paste0(inputDir,"/nsclc_gender_tumor_GSEA_GO_BP.tsv")))
data_tumor_and_normal <- as.data.frame(fread(paste0(inputDir,"/nsclc_gender_tumor_and_normal_GSEA_GO_BP.tsv")))
data_normal$condition = "normal"
data_tumor$condition = "tumor"
data_tumor_and_normal$condition = "tumor_normal"



#data_normal <- data_normal %>%
#  filter(p.adjust < 0.05)
#data_tumor <- data_tumor %>%
#  filter(p.adjust < 0.05)
#data_tumor_and_normal <- data_tumor_and_normal %>%
#  filter(p.adjust < 0.05)

total <- rbind(head(data_normal,10),head(data_tumor,10),head(data_tumor_and_normal,10))
total$Description <- gsub("- M.*", "", total$Description)



p <- ggplot(total, aes(NES,Description,fill = NES)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free", space = "free") +  
  labs(y = "GSEA GO BP", color = "black") # Change the y-axis label to "Pathways"
p2 <- ggplot(total, aes(NES,Description,fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free", space = "free") +  
  labs(y = "GSEA GO BP", color = "black") # Change the y-axis label to "Pathways"


# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
p2 <- p2 + scale_fill_gradient(low = "blue", high = "red")
# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired
p2 <- p2 + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired
# Print the plot
print(p)
print(p2)

# Save the plot as an image file (e.g., PNG)
ggsave(paste0(resDir, "GSEA__GO_BP_pathways_facet_grid.png"), plot = p, width = 8, height = 10)  # Adjust width and height as needed

ggsave(paste0(resDir,"GSEA_GO_BP_pathways_padjust_facet_grid.png"), plot = p2, width = 8, height = 10)  # Adjust width and height as needed


############################################# GSEA GO  BP  immune system 

resDir = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/deseq2_out/corrected_ds/GSEA/ORA/" 


inputDir = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/010_analysis_paired_include_guon/tables/deseq2_out/corrected_ds/GSEA"
data_normal<-as.data.frame(fread(paste0(inputDir,"/nsclc_gender_normal_GSEA_GO_BP.tsv")))
data_tumor <- as.data.frame(fread(paste0(inputDir,"/nsclc_gender_tumor_GSEA_GO_BP.tsv")))
data_tumor_and_normal <- as.data.frame(fread(paste0(inputDir,"/nsclc_gender_tumor_and_normal_GSEA_GO_BP.tsv")))
data_normal$condition = "normal"
data_tumor$condition = "tumor"
data_tumor_and_normal$condition = "tumor_normal"

data_normal <- data_normal[grepl("immune", data_normal$Description, ignore.case = TRUE), ]
data_tumor <- data_tumor[grepl("immune", data_tumor$Description, ignore.case = TRUE), ]
data_tumor_and_normal <- data_tumor_and_normal[grepl("immune", data_tumor_and_normal$Description, ignore.case = TRUE), ]




#data_normal <- data_normal %>%
#  filter(p.adjust < 0.05)
#data_tumor <- data_tumor %>%
#  filter(p.adjust < 0.05)
#data_tumor_and_normal <- data_tumor_and_normal %>%
#  filter(p.adjust < 0.05)

total <- rbind(head(data_normal,10),head(data_tumor,10),head(data_tumor_and_normal,10))
total$Description <- gsub("- M.*", "", total$Description)

# Wrap the Description text to limit its length
total$Description <- str_wrap(total$Description, width = 40)  # Adjust width as needed

# Truncate long labels
total$Description <- ifelse(nchar(total$Description) > 50, paste(str_sub(total$Description, end = 30), "..."), total$Description)


p <- ggplot(total, aes(NES,Description,fill = NES)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free", space = "free") +  
  labs(y = "GSEA GO BP", color = "black") # Change the y-axis label to "Pathways"
p2 <- ggplot(total, aes(NES,Description,fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free", space = "free") +  
  labs(y = "GSEA GO BP", color = "black") # Change the y-axis label to "Pathways"


# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
p2 <- p2 + scale_fill_gradient(low = "blue", high = "red")
# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired
p2 <- p2 + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired
# Print the plot
(p)
print(p2)

# Save the plot as an image file (e.g., PNG)
ggsave(paste0(resDir, "GSEA__GO_BP_pathways_facet_grid_immune.png"), plot = p, width = 8, height = 10)  # Adjust width and height as needed

ggsave(paste0(resDir,"GSEA_GO_BP_pathways_padjust_facet_grid_immune.png"), plot = p2, width = 8, height = 10)  # Adjust width and height as needed

