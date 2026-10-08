args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
chr <- args[2]
cell <- args[3]

library("data.table")
library("coloc")
library("argparse")
library("stringr")
library("BEDMatrix")
library("dplyr")
library("readxl")
library(arrow)

study ='CBMAP'
i = as.integer(str_split_fixed(chr,'chr',n=2)[,2])

bbj_trait_all = read_excel('/data/shared_data/BBJ_GWAS/BBJ_trait.xlsx')
bbj_type = as.character(bbj_trait_all[match(trait, bbj_trait_all$Trait), "Status"])

current_cpgs <- read.table(paste0('/methylation/mQTL/celltype_mQTL/tensorqtl/coloc/tmp/',study,'_',cell,'_chr',i,'_cpg.txt'),header=F)$V1
if (bbj_type == 'conti') {
    ### load GWAS data
    GWAS = fread(paste0('/data/shared_data/BBJ_GWAS/hum0197.v3.BBJ.',trait,'.v1/GWASsummary_',trait,'_Japanese_SakaueKanai2020.auto.txt.gz'))
    N = as.integer(bbj_trait_all[match(trait, bbj_trait_all$Trait), 4])
    GWAS = GWAS[,c('SNP', 'ALLELE1', 'ALLELE0', 'BETA', 'SE', 'A1FREQ')]
    GWAS$varbeta = GWAS$SE ^ 2; GWAS = GWAS[,-c('SE')]
    colnames(GWAS) = c('snp', 'ALT', 'REF', 'beta', 'MAF', 'varbeta')
    GWAS = distinct(GWAS, snp, .keep_all= TRUE)
    ### load mQTL data
    sd_cpg = read.table('/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cpg_sdY.txt',header=T)
    coloc_res_dir = '/methylation/mQTL/celltype_mQTL/tensorqtl/coloc/CBMAP/'
    mqtl_res_dir = '/methylation/mQTL/celltype_mQTL/tensorqtl/results/res_chr/'
    qtl_bim = fread('/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/sample_geno.bim')[,c('V2', 'V4', 'V5', 'V6')]
    cis_mqtl = as.data.table(read_parquet(paste0(mqtl_res_dir, cell, '_chr', i, '.cis_qtl_pairs.', i, '.parquet'),col_select = c('phenotype_id','b_gi','b_gi_se','variant_id')))
    cis_mqtl <- cis_mqtl[phenotype_id %in% current_cpgs]
    cis_mqtl = unique(cis_mqtl)
    colnames(cis_mqtl) = c('cpg', 'beta', 'varbeta', 'snp')
    cis_mqtl$varbeta = cis_mqtl$varbeta ^ 2
    colnames(qtl_bim) = c('snp', 'position', 'ALT', 'REF')
    cis_mqtl = merge(cis_mqtl, qtl_bim, by='snp')
} else if (bbj_type == 'bin') {
    ### load GWAS data
    GWAS = fread(paste0('/data/shared_data/BBJ_GWAS/hum0197.v3.BBJ.',trait,'.v1/GWASsummary_',trait,'_Japanese_SakaueKanai2020.auto.txt.gz'))[,c('SNPID', 'Allele2', 'Allele1', 'BETA', 'SE')]
    GWAS$varbeta = GWAS$SE ^ 2; GWAS = GWAS[,-c('SE')]
    colnames(GWAS) = c('snp', 'ALT', 'REF', 'beta', 'varbeta')
    GWAS = distinct(GWAS, snp, .keep_all= TRUE)
    ### load mQTL data
    sd_cpg = read.table('/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cpg_sdY.txt',header=T)
    coloc_res_dir = '/methylation/mQTL/celltype_mQTL/tensorqtl/coloc/CBMAP/'
    mqtl_res_dir = '/methylation/mQTL/celltype_mQTL/tensorqtl/results/res_chr/'
    qtl_bim = fread('/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/sample_geno.bim')[,c('V2', 'V4', 'V5', 'V6')]
    cis_mqtl = as.data.table(read_parquet(paste0(mqtl_res_dir, cell, '_chr', i, '.cis_qtl_pairs.', i, '.parquet'),col_select = c('phenotype_id','b_gi','b_gi_se','variant_id')))
    cis_mqtl <- cis_mqtl[phenotype_id %in% current_cpgs]
    cis_mqtl = unique(cis_mqtl)
    colnames(cis_mqtl) = c('cpg', 'beta', 'varbeta', 'snp')
    cis_mqtl$varbeta = cis_mqtl$varbeta ^ 2
    colnames(qtl_bim) = c('snp', 'position', 'ALT', 'REF')
    cis_mqtl = merge(cis_mqtl, qtl_bim, by='snp')
 }


### merge data
dat_merge = merge(cis_mqtl, GWAS, by='snp', suffixes = c("_qtl","_gwas"));rm(cis_mqtl);rm(GWAS)
dat_merge[dat_merge == ""] = NA
dat_merge = na.omit(dat_merge)
dat_merge <- dat_merge %>% filter((ALT_qtl == ALT_gwas & REF_qtl == REF_gwas) |
                                    (ALT_qtl == REF_gwas & REF_qtl == ALT_gwas))
dat_merge$beta_gwas <- ifelse(dat_merge$ALT_qtl == dat_merge$ALT_gwas, dat_merge$beta_gwas, -dat_merge$beta_gwas)


### run coloc
my.res.all = list()
cpg0 = unique(dat_merge$cpg)
for (j in 1:length(unique(dat_merge$cpg))) {
  cpg1 = cpg0[j]
  dat_merge1 = dat_merge[cpg == cpg1]
  my.res.all[[j]] = list()
  # Minimum to perform coloc set to 5 GWAS-mQTL variants
  if (nrow(dat_merge1) < 5) {
    my.res.all[[j]]$cpg = cpg1
    my.res.all[[j]]$summary = NA
    my.res.all[[j]]$result = NA
  } else {
    if (study == 'PKUH6') {
      if (bbj_type == 'bin') {
        dataset1 = as.list(dat_merge1[,c('snp', 'beta_qtl', 'varbeta_qtl', 'position', 'MAF')])
        names(dataset1) = c('snp', 'beta', 'varbeta', 'position', 'MAF')
        dataset1$N = 750
        dataset2 = as.list(dat_merge1[,c('snp', 'beta_gwas', 'varbeta_gwas')])
        dataset2$type = 'cc'
        names(dataset2) = c('snp', 'beta', 'varbeta', 'type')
      } else if (bbj_type == 'quant') {
        dataset1 = as.list(dat_merge1[,c('snp', 'beta_qtl', 'varbeta_qtl', 'position', 'MAF_qtl')])
        names(dataset1) = c('snp', 'beta', 'varbeta', 'position', 'MAF')
        dataset1$N = 750
        dataset2 = as.list(dat_merge1[,c('snp', 'beta_gwas', 'varbeta_gwas', 'MAF_gwas')])
        dataset2$type = 'quant'
        dataset2$N = N
        names(dataset2) = c('snp', 'beta', 'varbeta', 'MAF', 'type', 'N')
      }
    } else if (study == 'ROSMAP' | study == 'CBMAP') {
      dataset1 = as.list(dat_merge1[,c('snp', 'beta_qtl', 'varbeta_qtl', 'position')])
      names(dataset1) = c('snp', 'beta', 'varbeta', 'position')
      dataset1$sdY = sd_cpg$sdY[which(sd_cpg$cg_id == cpg1)]
      if (bbj_type == 'bin') {
        dataset2 = as.list(dat_merge1[,c('snp', 'beta_gwas', 'varbeta_gwas')])
        dataset2$type = 'cc'
        names(dataset2) = c('snp', 'beta', 'varbeta', 'type')
      } else if (bbj_type == 'conti') {
        dataset2 = as.list(dat_merge1[,c('snp', 'beta_gwas', 'varbeta_gwas', 'MAF')])
        dataset2$type = 'quant'
        dataset2$N = N
        names(dataset2) = c('snp', 'beta', 'varbeta', 'MAF', 'type', 'N')
      }
      
    }
    dataset1$type = 'quant'
    check_dataset(dataset1); check_dataset(dataset2)
    
    ### run coloc
    my.res <- coloc.abf(dataset1=dataset1, dataset2=dataset2)
    my.res.all[[j]]$cpg = cpg1
    my.res.all[[j]]$summary = my.res$summary
    my.res.all[[j]]$result = my.res$results[,c('snp','position','SNP.PP.H4')]
    my.res.all[[j]]$result = merge(my.res.all[[j]]$result, dat_merge1[,c('snp', 'ALT_qtl', 'REF_qtl')], by='snp')
  }
}
if(!dir.exists(paste0(coloc_res_dir,trait,'/'))){
  dir.create(paste0(coloc_res_dir,trait,'/'))
}
save(my.res.all, file=paste0(coloc_res_dir, trait,'/', study, '_',cell,'_',trait,'_coloc_chr', i, '.RData'))

