args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
chr <- args[2]
cell <- args[3]

.libPaths(c(.libPaths(),"/share/data/R4.3_lib/library","/share/home/zhaiwl/R/x86_64-pc-linux-gnu-library/4.3"))
library("data.table")
library("coloc")
library("argparse")
library("stringr")
library("BEDMatrix")
library("dplyr")
library("readxl")
library(arrow)

study ='ROSMAP'
i = as.integer(str_split_fixed(chr,'chr',n=2)[,2])

# res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/results/'
# cells = c("Endo", "Exc", "Inh", "Micro", "Oligo", "OPC", "Astro")
# 
# cpg_info = fread('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cg_BED.bed.gz')
# 
# for (cell in cells) {
#   res = fread(paste0(res_dir,'CBMAP_',cell,'_imqtl_1e-5.txt'),header=T)
#   length(unique(res$phenotype_id))
#   cpg_df = data.frame('cpg'=unique(res$phenotype_id),'chr'=NA)
#   cpg_df$chr = cpg_info[match(cpg_df$cpg, cpg_info$pid),`#Chr`]
#   for (i in 1:22) {
#     cpg_by_chr = cpg_df$cpg[which(cpg_df$chr == i)]
#     write.table(cpg_by_chr,file=paste0('/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/coloc/tmp/','CBMAP_',cell,'_chr',i,'_cpg.txt'),quote=F,col.names = F,row.names = F,sep='\t')
#   }
# }

current_cpgs <- read.table(paste0('/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/coloc/tmp/',study,'_',cell,'_chr',i,'_cpg.txt'),header=F)$V1
if (trait == 'BMI') { # for quasi-continuous phenotypes
  if (study == 'CBMAP') {
    ### load GWAS data
    GWAS = fread("/data/shared_data/neuropsych_GWAS/EAS/processed/BMI/BMI_processed.txt")[,c('rsid', 'effect_allele', 'reference_allele', 'beta', 'se', 'Frq')]
    GWAS$varbeta = GWAS$se ^ 2; GWAS = GWAS[,-c('se')]
    colnames(GWAS) = c('snp', 'ALT', 'REF', 'beta', 'MAF', 'varbeta')
    GWAS = distinct(GWAS, snp, .keep_all= TRUE)
    GWAS$MAF = ifelse(GWAS$MAF == 0, 0.0001, GWAS$MAF)
    GWAS$MAF = ifelse(GWAS$MAF == 1, 0.9999, GWAS$MAF)
    N = 158284
    ### load mQTL data
    sd_cpg = read.table('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cpg_sdY.txt',header=T)
    coloc_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/coloc/CBMAP/'
    mqtl_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/results/res_chr/'
    qtl_bim = fread('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/sample_geno.bim')[,c('V2', 'V4', 'V5', 'V6')]
    cis_mqtl = as.data.table(read_parquet(paste0(mqtl_res_dir, cell, '_chr', i, '.cis_qtl_pairs.', i, '.parquet'),col_select = c('phenotype_id','b_g','b_g_se','variant_id')))
    cis_mqtl <- cis_mqtl[phenotype_id %in% current_cpgs]
    cis_mqtl = unique(cis_mqtl)
    colnames(cis_mqtl) = c('cpg', 'beta', 'varbeta', 'snp')
    cis_mqtl$varbeta = cis_mqtl$varbeta ^ 2
    colnames(qtl_bim) = c('snp', 'position', 'ALT', 'REF')
    cis_mqtl = merge(cis_mqtl, qtl_bim, by='snp')
  } else if (study == 'ROSMAP') {
    ### load GWAS data
    GWAS = fread("/data/shared_data/neuropsych_GWAS/EUR/processed/BMI/BMI_processed.txt")[,c('rsid', 'effect_allele', 'reference_allele', 'beta', 'se', 'Freq_Tested_Allele')]
    GWAS$varbeta = GWAS$se ^ 2; GWAS = GWAS[,-c('se')]
    colnames(GWAS) = c('snp', 'ALT', 'REF', 'beta', 'MAF', 'varbeta')
    GWAS = distinct(GWAS, snp, .keep_all= TRUE)
    GWAS$MAF = ifelse(GWAS$MAF == 0, 0.0001, GWAS$MAF)
    GWAS$MAF = ifelse(GWAS$MAF == 1, 0.9999, GWAS$MAF)
    N = 694649
    ### load mQTL data
    sd_cpg = read.table('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/ROSMAP_mQTL/mQTL_mapping/b37/QTLtools_input/cpg_sdY.txt',header=T)
    coloc_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/coloc/ROSMAP/'
    mqtl_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/results/res_chr/'
    qtl_bim = fread('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/ROSMAP_mQTL/mQTL_mapping/b37/QTLtools_input/pca_selection/geno_pca/sample_geno.bim')[,c('V2', 'V4', 'V5', 'V6')]
    cis_mqtl = as.data.table(read_parquet(paste0(mqtl_res_dir, 'ROSMAP_', cell, '_chr', i, '.cis_qtl_pairs.', i, '.parquet'),col_select = c('phenotype_id','b_g','b_g_se','variant_id')))
    cis_mqtl <- cis_mqtl[phenotype_id %in% current_cpgs]
    cis_mqtl = unique(cis_mqtl)
    colnames(cis_mqtl) = c('cpg', 'beta', 'varbeta', 'snp')
    cis_mqtl$varbeta = cis_mqtl$varbeta ^ 2
    colnames(qtl_bim) = c('snp', 'position', 'ALT', 'REF')
    cis_mqtl = merge(cis_mqtl, qtl_bim, by='snp')
  } 
} else if (trait %chin% c('AD','SCZ')) { # for binary phenotypes
  if (study == 'CBMAP') {
    ### load GWAS data
    if (trait == 'AD') {
      GWAS = fread("/data/shared_data/TPMI/result/meta_analysis/meta_results_1.tbl")[,c('MarkerName', 'Allele1', 'Allele2', 'Effect', 'StdErr')]
      colnames(GWAS) = c('rsid', 'effect_allele', 'reference_allele', 'beta', 'se')
    } else if (trait == 'SCZ') {
      GWAS = fread("/data/shared_data/neuropsych_GWAS/EAS/processed/SCZ/SCZ_processed2.txt")[,c('rsid', 'effect_allele', 'reference_allele', 'beta', 'se')]
    }
    GWAS$varbeta = GWAS$se ^ 2; GWAS = GWAS[,-c('se')]
    colnames(GWAS) = c('snp', 'ALT', 'REF', 'beta', 'varbeta')
    GWAS = distinct(GWAS, snp, .keep_all= TRUE)
    ### load mQTL data
    sd_cpg = read.table('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cpg_sdY.txt',header=T)
    coloc_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/coloc/CBMAP/'
    mqtl_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/results/res_chr/'
    qtl_bim = fread('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/sample_geno.bim')[,c('V2', 'V4', 'V5', 'V6')]
    cis_mqtl = as.data.table(read_parquet(paste0(mqtl_res_dir, cell, '_chr', i, '.cis_qtl_pairs.', i, '.parquet'),col_select = c('phenotype_id','b_g','b_g_se','variant_id')))
    cis_mqtl <- cis_mqtl[phenotype_id %in% current_cpgs]
    cis_mqtl = unique(cis_mqtl)
    colnames(cis_mqtl) = c('cpg', 'beta', 'varbeta', 'snp')
    cis_mqtl$varbeta = cis_mqtl$varbeta ^ 2
    colnames(qtl_bim) = c('snp', 'position', 'ALT', 'REF')
    cis_mqtl = merge(cis_mqtl, qtl_bim, by='snp')
  } else if (study == 'ROSMAP') {
    ### load GWAS data
    if (trait == 'AD') {
      GWAS = fread("/data/shared_data/neuropsych_GWAS/EUR/processed/AD/file1/AD2_processed.txt")[,c('rsid', 'effect_allele', 'reference_allele', 'beta', 'se')]
    } else if (trait == 'SCZ') {
      GWAS = fread("/data/shared_data/neuropsych_GWAS/EUR/processed/SCZ/file1/SCZfile1_processed.txt")[,c('rsid', 'effect_allele', 'reference_allele', 'beta', 'se')]
    }
    GWAS$varbeta = GWAS$se ^ 2; GWAS = GWAS[,-c('se')]
    colnames(GWAS) = c('snp', 'ALT', 'REF', 'beta', 'varbeta')
    GWAS = distinct(GWAS, snp, .keep_all= TRUE)
    ### load mQTL data
    sd_cpg = read.table('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/ROSMAP_mQTL/mQTL_mapping/b37/QTLtools_input/cpg_sdY.txt',header=T)
    coloc_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/coloc/ROSMAP/'
    mqtl_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/results/res_chr/'
    qtl_bim = fread('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/ROSMAP_mQTL/mQTL_mapping/b37/QTLtools_input/pca_selection/geno_pca/sample_geno.bim')[,c('V2', 'V4', 'V5', 'V6')]
    cis_mqtl = as.data.table(read_parquet(paste0(mqtl_res_dir, 'ROSMAP_', cell, '_chr', i, '.cis_qtl_pairs.', i, '.parquet'),col_select = c('phenotype_id','b_g','b_g_se','variant_id')))
    cis_mqtl <- cis_mqtl[phenotype_id %in% current_cpgs]
    cis_mqtl = unique(cis_mqtl)
    colnames(cis_mqtl) = c('cpg', 'beta', 'varbeta', 'snp')
    cis_mqtl$varbeta = cis_mqtl$varbeta ^ 2
    colnames(qtl_bim) = c('snp', 'position', 'ALT', 'REF')
    cis_mqtl = merge(cis_mqtl, qtl_bim, by='snp')
  } 
}

### merge data
dat_merge = merge(cis_mqtl, GWAS, by='snp', suffixes = c("_qtl","_gwas"));rm(cis_mqtl);rm(GWAS)
dat_merge$ALT_gwas = toupper(dat_merge$ALT_gwas); dat_merge$REF_gwas = toupper(dat_merge$REF_gwas)
dat_merge[dat_merge == ""] = NA
dat_merge = na.omit(dat_merge)
dat_merge <- dat_merge %>% filter((ALT_qtl == ALT_gwas & REF_qtl == REF_gwas) |
                                    (ALT_qtl == REF_gwas & REF_qtl == ALT_gwas))
dat_merge$beta_gwas <- ifelse(dat_merge$ALT_qtl == dat_merge$ALT_gwas, dat_merge$beta_gwas, -dat_merge$beta_gwas)

### run coloc
my.res.all = list()
cpg0 = unique(dat_merge$cpg)
for (j in 1:length(cpg0)) {
  cpg1 = cpg0[j]
  dat_merge1 = dat_merge[cpg == cpg1]
  my.res.all[[j]] = list()
  my.res.all[[j]]$cpg = cpg1
  # Minimum to perform coloc set to 5 GWAS-mQTL variants
  if (nrow(dat_merge1) < 5) {
    my.res.all[[j]]$summary = NA
    my.res.all[[j]]$result = NA
  } else {
    dataset1 = as.list(dat_merge1[,c('snp', 'beta_qtl', 'varbeta_qtl', 'position')])
    names(dataset1) = c('snp', 'beta', 'varbeta', 'position')
    dataset1$sdY = sd_cpg$sdY[which(sd_cpg$cg_id == cpg1)]
    if (trait %chin% c('AD','SCZ')) {
      dataset2 = as.list(dat_merge1[,c('snp', 'beta_gwas', 'varbeta_gwas')])
      dataset2$type = 'cc'
      names(dataset2) = c('snp', 'beta', 'varbeta', 'type')
    } else if (trait == 'BMI') {
      dataset2 = as.list(dat_merge1[,c('snp', 'beta_gwas', 'varbeta_gwas', 'MAF')])
      dataset2$type = 'quant'
      dataset2$N = N
      names(dataset2) = c('snp', 'beta', 'varbeta', 'MAF', 'type', 'N')
    }
    dataset1$type = 'quant'
    check_dataset(dataset1); check_dataset(dataset2)
  
    ### run coloc
    my.res <- coloc.abf(dataset1=dataset1, dataset2=dataset2)
    my.res.all[[j]]$summary = my.res$summary
    my.res.all[[j]]$result = my.res$results[,c('snp','position','SNP.PP.H4')]
  }
  my.res.all[[j]]$result = merge(my.res.all[[j]]$result, dat_merge1[,c('snp', 'ALT_qtl', 'REF_qtl')], by='snp')
}


if(!dir.exists(paste0(coloc_res_dir,trait,'/'))){
  dir.create(paste0(coloc_res_dir,trait,'/'))
}
save(my.res.all, file=paste0(coloc_res_dir, trait,'/', study, '_',cell,'_',trait,'_coloc_chr', i, '.RData'))

# library(readxl)
# library(data.table)
# traits <- c('AD','SCZ','BMI')
# cells <- c("Exc", "Inh", "Micro", "Oligo", "OPC", "Astro","Endo")
# chr_batches <- list(
#   chr1 = 1,
#   chr2 = 1,
#   chr3 = 1,
#   chr4 = 1,
#   chr5 = 1,
#   chr6 = 1,
#   chr7 = 1,
#   chr8 = 1,
#   chr9 = 1,
#   chr10 = 1,
#   chr11 = 1,
#   chr12 = 1,
#   chr13 = 1,
#   chr14 = 1,
#   chr15 = 1,
#   chr16 = 1,
#   chr17 = 1,
#   chr18 = 1,
#   chr19 = 1,
#   chr20 = 1,
#   chr21 = 1,
#   chr22 = 1
# )
# 
# lines <- c()
# for (trait in traits) {
#   for (chr in names(chr_batches)) {
#     for (cell in cells) {
#       lines <- c(lines, paste(trait, chr, cell, sep = ","))
#     }
#   }
# }
# 
# writeLines(lines, "/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/coloc/tmp/cbmap_task_list2.txt")
