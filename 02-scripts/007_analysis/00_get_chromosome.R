

library(biomaRt)

path ="/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/"
file="adata_var_nsclc.csv"

adata_var_nsclc <- read.csv(paste0(path,file))

#Get gene names annotation

biolist <- as.data.frame(listMarts())
ensembl=useMart("ensembl")
esemblist <- as.data.frame(listDatasets(ensembl))
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

t2g<-getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version",'chromosome_name','start_position','end_position'), mart = ensembl)

my_ids <- data.frame(adata_var_nsclc$gene_id)

colnames(my_ids)[1] <- "ensembl_gene_id"

adata_var_nsclc_chrom <- merge(my_ids, t2g, by= 'ensembl_gene_id')


write.csv(adata_var_nsclc_chrom, paste0(path,"adata_var_nsclc_chrom.csv"))
