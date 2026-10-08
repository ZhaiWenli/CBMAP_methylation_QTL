args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
chr <- args[2]
u <- as.integer(args[3])

cat(sprintf("Running coloc for chromosome: %s, batch: %d\n", chr, u))

library("data.table")
library("coloc")
library("argparse")
library("stringr")
library("BEDMatrix")
library("dplyr")

i = as.integer(str_split_fixed(chr,'chr',n=2)[,2])
print(u); print(i)
study = 'CBMAP'

if (study == 'CBMAP') {
  ### load GWAS data
  GWAS = fread("/TPMI/result/meta_analysis/meta_results_1.tbl")[,c('MarkerName', 'Allele1', 'Allele2', 'Effect', 'StdErr')]
  GWAS$varbeta = GWAS$StdErr ^ 2; GWAS = GWAS[,-c('StdErr')]
  colnames(GWAS) = c('snp', 'ALT', 'REF', 'beta', 'varbeta')
  GWAS = distinct(GWAS, snp, .keep_all= TRUE)
  ### load mQTL data
  sd_cpg = read.table('/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cpg_sdY.txt',header=T)
  coloc_res_dir = '/methylation/mQTL/coloc/result/cbmap/coloc/'
  mqtl_res_dir = '/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_result/chr_res/'
  qtl_bim = fread('/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/sample_geno.bim')[,c('V2', 'V4', 'V5', 'V6')]
  cis_mqtl = fread(paste0(mqtl_res_dir, 'qtltools_cis_nominal_chr', i, '.txt'))[,c('V1', 'V14', 'V15', 'V8')]
  cis_mqtl = unique(cis_mqtl)
  unique_cpg <- unique(cis_mqtl$V1)
  batch_size <- 4000
  num_batches <- ceiling(length(unique_cpg) / batch_size)
  start_idx <- (u - 1) * batch_size + 1
  end_idx <- min(u * batch_size, length(unique_cpg))
  current_cpgs <- unique_cpg[start_idx:end_idx]
  cis_mqtl <- cis_mqtl[V1 %in% current_cpgs]
  colnames(cis_mqtl) = c('cpg', 'beta', 'varbeta', 'snp')
  cis_mqtl$varbeta = cis_mqtl$varbeta ^ 2
  colnames(qtl_bim) = c('snp', 'position', 'ALT', 'REF')
  cis_mqtl = merge(cis_mqtl, qtl_bim, by='snp')
} else if (study == 'ROSMAP') {
  ### load GWAS data
  GWAS = fread("/neuropsych_GWAS/EUR/processed/AD/file1/AD2_processed.txt")[,c('rsid', 'effect_allele', 'reference_allele', 'beta', 'se')]
  GWAS$varbeta = GWAS$se ^ 2; GWAS = GWAS[,-c('se')]
  colnames(GWAS) = c('snp', 'ALT', 'REF', 'beta', 'varbeta')
  GWAS = distinct(GWAS, snp, .keep_all= TRUE)
  ### load mQTL data
  sd_cpg = read.table('/methylation/mQTL/mQTL_mapping/result/ROSMAP_mQTL/mQTL_mapping/b37/QTLtools_input/cpg_sdY.txt',header=T)
  coloc_res_dir = '/methylation/mQTL/coloc/result/rosmap/coloc/'
  mqtl_res_dir = '/methylation/mQTL/mQTL_mapping/result/ROSMAP_mQTL/mQTL_mapping/b37/QTLtools_result/chr_res/'
  qtl_bim = fread('/methylation/mQTL/mQTL_mapping/result/ROSMAP_mQTL/mQTL_mapping/b37/QTLtools_input/pca_selection/geno_pca/sample_geno.bim')[,c('V2', 'V4', 'V5', 'V6')]
  cis_mqtl = fread(paste0(mqtl_res_dir, 'qtltools_cis_nominal_chr', i, '.txt'))[,c('V1', 'V14', 'V15', 'V8')]
  cis_mqtl = unique(cis_mqtl)
  colnames(cis_mqtl) = c('cpg', 'beta', 'varbeta', 'snp')
  cis_mqtl$varbeta = cis_mqtl$varbeta ^ 2
  colnames(qtl_bim) = c('snp', 'position', 'ALT', 'REF')
  cis_mqtl = merge(cis_mqtl, qtl_bim, by='snp')
} else if (study == 'PKUH6') {
  ### load GWAS data
  GWAS = fread("/TPMI/result/meta_analysis/meta_results_1.tbl")[,c('MarkerName', 'Allele1', 'Allele2', 'Effect', 'StdErr')]
  GWAS$varbeta = GWAS$StdErr ^ 2; GWAS = GWAS[,-c('StdErr')]
  colnames(GWAS) = c('snp', 'ALT', 'REF', 'beta', 'varbeta')
  GWAS = distinct(GWAS, snp, .keep_all= TRUE)
  ### load mQTL data
  coloc_res_dir = '/methylation/mQTL/coloc/result/pkuh6/coloc/'
  mqtl_res_dir = '/methylationQTL/beijing_6_hospital/mqtl/'
  qtl_bim = fread('/methylation/mQTL/coloc/result/pkuh6/snp_anno.txt',header=T)[,c('rsid', 'pos', 'ALT', 'REF', 'MAF','ID')]
  colnames(qtl_bim) = c('snp', 'position', 'ALT', 'REF', 'MAF','ID')
  cis_mqtl = fread(paste0(mqtl_res_dir, 'cis.nominal.chr', i))[,c('V1', 'V14', 'V15', 'V8')]
  unique_cpg <- unique(cis_mqtl$V1)
  batch_size <- 4000
  num_batches <- ceiling(length(unique_cpg) / batch_size)
  start_idx <- (u - 1) * batch_size + 1
  end_idx <- min(u * batch_size, length(unique_cpg))
  current_cpgs <- unique_cpg[start_idx:end_idx]
  cis_mqtl <- cis_mqtl[V1 %in% current_cpgs]
  colnames(cis_mqtl) = c('cpg', 'beta', 'varbeta', 'id')
  cis_mqtl$varbeta = cis_mqtl$varbeta ^ 2
  cis_mqtl$snp = qtl_bim[match(cis_mqtl$id, qtl_bim$ID), 'snp']
  cis_mqtl = na.omit(cis_mqtl)
  cis_mqtl = distinct(cis_mqtl, id, cpg, .keep_all=T)
  cis_mqtl = merge(cis_mqtl, qtl_bim, by.x='id', by.y='ID')
  cis_mqtl = cis_mqtl[,-c('snp.x')]
  colnames(cis_mqtl)[5] = 'snp'
  cis_mqtl = distinct(cis_mqtl, snp, cpg, .keep_all=T)
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
  # Minimum to perform coloc set to 5 GWAS-mQTL variants
  if (nrow(dat_merge1) < 5) {
    my.res.all[[j]]$cpg = cpg1
    my.res.all[[j]]$summary = NA
    my.res.all[[j]]$result = NA
  } else {
    if (study == 'PKUH6') {
      dataset1 = as.list(dat_merge1[,c('snp', 'beta_qtl', 'varbeta_qtl', 'position', 'MAF')])
      names(dataset1) = c('snp', 'beta', 'varbeta', 'position', 'MAF')
      dataset1$N = 750
    } else if (study == 'ROSMAP') {
      dataset1 = as.list(dat_merge1[,c('snp', 'beta_qtl', 'varbeta_qtl', 'position')])
      names(dataset1) = c('snp', 'beta', 'varbeta', 'position')
      dataset1$sdY = sd_cpg$sdY[which(sd_cpg$cg_id == cpg1)]
    } else if (study == 'CBMAP') {
      dataset1 = as.list(dat_merge1[,c('snp', 'beta_qtl', 'varbeta_qtl', 'position')])
      names(dataset1) = c('snp', 'beta', 'varbeta', 'position')
      dataset1$sdY = sd_cpg$sdY[which(sd_cpg$cg_id == cpg1)]
    }
    dataset1$type = 'quant'
    dataset2 = as.list(dat_merge1[,c('snp', 'beta_gwas', 'varbeta_gwas')])
    dataset2$type = 'cc'
    names(dataset2) = c('snp', 'beta', 'varbeta', 'type')
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
if(!dir.exists(paste0(coloc_res_dir,trait,'/chr',i,'_',trait,'_batch/'))){
  dir.create(paste0(coloc_res_dir,trait,'/chr',i,'_',trait,'_batch/'))
}
save(my.res.all, file=paste0(coloc_res_dir,trait,'/chr',i,'_',trait,'_batch/', study, '_',trait,'_coloc_chr', i, '_batch', u, '.RData'))
