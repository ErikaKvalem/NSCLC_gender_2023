#!/bin/bash
#SBATCH --job-name=gender_deseq2
#SBATCH --output=gender_deseq2.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=06:00:00

sbatch 02_DESeq2.R /data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/input/samplesheet.csv /data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/input/counts.csv --result_dir=/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/deseq2_out --c1=male --c2=female --sample_col=sample --prefix=all --condition_col=sex --covariate_formula="+origin + disease +  tumor_stage"  --fdr_cutoff=0.1 --fc_cutoff=1
 
sbatch 02_DESeq2.R /data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/input/samplesheet_normal.csv /data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/input/counts_normal.csv --result_dir=/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/deseq2_out --c1=male --c2=female --sample_col=sample --prefix=normal --condition_col=sex --covariate_formula=""  --fdr_cutoff=0.1 --fc_cutoff=1
 
sbatch 02_DESeq2.R /data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/input/samplesheet_tumor.csv /data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/input/counts_tumor.csv --result_dir=/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp/out/007_re_analysis/tables/deseq2_out --c1=male --c2=female --sample_col=sample --prefix=tumor --condition_col=sex --covariate_formula=""  --fdr_cutoff=0.1 --fc_cutoff=1
 
