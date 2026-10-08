args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
chr <- args[2]
u <- as.integer(args[3])
study = args[4]

.libPaths(c(.libPaths(),"/share/data/R4.3_lib/library"))
library("data.table")
library("coloc")
library("argparse")
library("stringr")
library("BEDMatrix")
library("dplyr")
library("readxl")

i = as.integer(str_split_fixed(chr,'chr',n=2)[,2])
print(u); print(i)

unique_cpg = read.table(paste0('/data/projects/China_Brain_MultiOmics/methylation/mQTL/coloc/result/cbmap/coloc/tmp/',study,'_chr',i,'_cpg.txt'),header=F)$V1
batch_size <- 5000
start_idx <- (u - 1) * batch_size + 1
end_idx <- min(u * batch_size, length(unique_cpg))
current_cpgs <- unique_cpg[start_idx:end_idx]


if (study == 'CBMAP') {
  ### load GWAS data
  GWAS = fread("/data/shared_data/neuropsych_GWAS/EAS/processed/SCZ/SCZ_processed2.txt")[,c('rsid', 'effect_allele', 'reference_allele', 'beta', 'se')]
  GWAS$varbeta = GWAS$se ^ 2; GWAS = GWAS[,-c('se')]
  colnames(GWAS) = c('snp', 'ALT', 'REF', 'beta', 'varbeta')
  GWAS = distinct(GWAS, snp, .keep_all= TRUE)
  ### load mQTL data
  sd_cpg = read.table('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cpg_sdY.txt',header=T)
  coloc_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/coloc/result/cbmap/coloc/'
  mqtl_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_result/chr_res/'
  qtl_bim = fread('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/sample_geno.bim')[,c('V2', 'V4', 'V5', 'V6')]
  cis_mqtl = fread(paste0(mqtl_res_dir, 'qtltools_cis_nominal_chr', i, '.txt'))[,c('V1', 'V14', 'V15', 'V8')]
  cis_mqtl <- cis_mqtl[V1 %in% current_cpgs]
  cis_mqtl = unique(cis_mqtl)
  colnames(cis_mqtl) = c('cpg', 'beta', 'varbeta', 'snp')
  cis_mqtl$varbeta = cis_mqtl$varbeta ^ 2
  colnames(qtl_bim) = c('snp', 'position', 'ALT', 'REF')
  cis_mqtl = merge(cis_mqtl, qtl_bim, by='snp')
} else if (study == 'ROSMAP') {
  ### load GWAS data
  GWAS = fread("/data/shared_data/neuropsych_GWAS/EUR/processed/SCZ/file1/SCZfile1_processed.txt")[,c('rsid', 'effect_allele', 'reference_allele', 'beta', 'se')]
  GWAS$varbeta = GWAS$se ^ 2; GWAS = GWAS[,-c('se')]
  colnames(GWAS) = c('snp', 'ALT', 'REF', 'beta', 'varbeta')
  GWAS = distinct(GWAS, snp, .keep_all= TRUE)
  ### load mQTL data
  sd_cpg = read.table('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/ROSMAP_mQTL/mQTL_mapping/b37/QTLtools_input/cpg_sdY.txt',header=T)
  coloc_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/coloc/result/rosmap/coloc/'
  mqtl_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/ROSMAP_mQTL/mQTL_mapping/b37/QTLtools_result/chr_res/'
  qtl_bim = fread('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/ROSMAP_mQTL/mQTL_mapping/b37/QTLtools_input/pca_selection/geno_pca/sample_geno.bim')[,c('V2', 'V4', 'V5', 'V6')]
  cis_mqtl = fread(paste0(mqtl_res_dir, 'qtltools_cis_nominal_chr', i, '.txt'))[,c('V1', 'V14', 'V15', 'V8')]
  cis_mqtl <- cis_mqtl[V1 %in% current_cpgs]
  cis_mqtl = unique(cis_mqtl)
  colnames(cis_mqtl) = c('cpg', 'beta', 'varbeta', 'snp')
  cis_mqtl$varbeta = cis_mqtl$varbeta ^ 2
  colnames(qtl_bim) = c('snp', 'position', 'ALT', 'REF')
  cis_mqtl = merge(cis_mqtl, qtl_bim, by='snp')
} else if (study == 'PKUH6') {
  ### load GWAS data
  GWAS = fread("/data/shared_data/neuropsych_GWAS/EAS/processed/SCZ/SCZ_processed2.txt")[,c('rsid', 'effect_allele', 'reference_allele', 'beta', 'se')]
  GWAS$varbeta = GWAS$se ^ 2; GWAS = GWAS[,-c('se')]
  colnames(GWAS) = c('snp', 'ALT', 'REF', 'beta', 'varbeta')
  GWAS = distinct(GWAS, snp, .keep_all= TRUE)
  ### load mQTL data
  coloc_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/coloc/result/pkuh6/coloc/'
  mqtl_res_dir = '/data/shared_data/methylationQTL/beijing_6_hospital/mqtl/'
  qtl_bim = fread('/data/projects/China_Brain_MultiOmics/methylation/mQTL/coloc/result/pkuh6/snp_anno.txt',header=T)[,c('rsid', 'pos', 'ALT', 'REF', 'MAF','ID')]
  colnames(qtl_bim) = c('snp', 'position', 'ALT', 'REF', 'MAF','ID')
  cis_mqtl = fread(paste0(mqtl_res_dir, 'cis.nominal.chr', i))[,c('V1', 'V14', 'V15', 'V8')]
  cis_mqtl <- cis_mqtl[V1 %in% current_cpgs]
  cis_mqtl = unique(cis_mqtl)
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
      dataset2 = as.list(dat_merge1[,c('snp', 'beta_gwas', 'varbeta_gwas')])
      dataset2$type = 'cc'
      names(dataset2) = c('snp', 'beta', 'varbeta', 'type')
    } else if (study %chin% c('ROSMAP','CBMAP')) {
      dataset1 = as.list(dat_merge1[,c('snp', 'beta_qtl', 'varbeta_qtl', 'position')])
      names(dataset1) = c('snp', 'beta', 'varbeta', 'position')
      dataset1$sdY = sd_cpg$sdY[which(sd_cpg$cg_id == cpg1)]
      dataset2 = as.list(dat_merge1[,c('snp', 'beta_gwas', 'varbeta_gwas')])
      dataset2$type = 'cc'
      names(dataset2) = c('snp', 'beta', 'varbeta', 'type')
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


if(!dir.exists(paste0(coloc_res_dir,'SCZ/'))){
  dir.create(paste0(coloc_res_dir,'SCZ/'))
}
save(my.res.all, file=paste0(coloc_res_dir, 'SCZ/', study, '_SCZ_coloc_chr', i, '_batch', u, '.RData'))


# DNAm <- fread('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cg_BED.bed')
# geno <- BEDMatrix('/data/projects/China_Brain_MultiOmics/methylation/results/CBMAP_mQTL/QTLtools_input/geno_pca_5/sample_geno.bed')
# colnames(geno) <- str_split_fixed(colnames(geno), '_', n=2)[,1]
# sdY <- data.frame(cg_id=rep(NA,nrow(DNAm)),sdY=NA)
# sdY$cg_id <- DNAm$pid
# sdY$sdY <- apply(DNAm[,7:ncol(DNAm)], 1, sd)
# sdY=distinct(sdY)
# write.table(sdY, file='/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cpg_sdY.txt',col.names=T,row.names=F,quote=F)
# sdX <- data.frame(snp_id=rep(NA,ncol(geno)),sdX=NA)
# sdX$snp_id <- colnames(geno)
# sdX$sdX <- apply(geno, 2, sd, na.rm=T)
# sdX=distinct(sdX)
# trans_mqtl <- merge(trans_mqtl, sdY, by='cg_id')
# trans_mqtl <- merge(trans_mqtl, sdX, by.x='variants_id', by.y='snp_id')
# trans_mqtl$beta <- trans_mqtl$cor_coef * (trans_mqtl$sdY / trans_mqtl$sdX)
# trans_mqtl$varbeta <- ((trans_mqtl$sdY / (trans_mqtl$sdX * sqrt(args$qtl_sample_size))) * sqrt(1 / (1 - trans_mqtl$cor_coef ^2))) ^ 2
# trans_mqtl = trans_mqtl[which(trans_mqtl$variants_chr==22),]

# my.res.all=as.data.frame(t(c(NA,NA,NA,NA,NA,NA)))
# colnames(my.res.all)=c('nsnps','PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf')
# my.res.all.top=as.data.frame(t(c(NA,NA)))
# colnames(my.res.all.top)=c('topColoc_variant','SNP.PP.H4')
# PP4_susie_gwas_mqtl.all=as.data.frame(t(rep(NA, 10)))
# colnames(PP4_susie_gwas_mqtl.all)=c('nsnps','hit1','hit2','PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf','idx1','idx2')
# cpgs_tested=c()
# 
# t1 = Sys.time()
# for (cpg in as.character(cpgs)) {
#   print(cpg)
#   mQTL = subset(mQTLall, cg_id%in%cpg)
#   stats = merge(GWAS, mQTL, by.x='rsid', by.y='variants_id')
#   
#   # Minimum to perform coloc set to 50 GWAS-mQTL variants
#   if (nrow(stats) < 50) {next}
#   
#   if( args$type=="cc" ) {
#     # coloc.abf:
#     try(my.res <- coloc.abf(dataset1=list(beta=stats$beta.x, varbeta=stats$se ^ 2, snp=stats$rsid, type="cc"),
#                             dataset2=list(beta=stats$beta.y, varbeta=stats$varbeta, snp=stats$rsid, sdY=sdY[which(sdY$cg_id==cpg),'sdY'], type="quant")))
#     
#     # coloc.susie:
#     LD = cor(geno[,match(unlist(stats$rsid), colnames(geno))],use='pairwise.complete.obs')
#     s_gwas = runsusie(list(beta=stats$beta.x, varbeta=stats$se ^ 2, snp=stats$rsid, type="cc", LD=LD, N=args$cases, s=args$cases/args$gwas_sample_size))
#     s_qtl = runsusie(list(beta=stats$beta.y, varbeta=stats$varbeta, snp=stats$rsid, sdY=sdY[which(sdY$cg_id==cpg),'sdY'], type="quant", LD=LD, N=args$qtl_sample_size))
#     if(sum(dim(summary(s_gwas)$cs)[1]>0) & sum(dim(summary(s_qtl)$cs)[1]>0)){
#       susie_gwas_mqtl=coloc.susie(s_gwas,s_qtl)
#       PP4_susie_gwas_mqtl = susie_gwas_mqtl$summary
#       print(PP4_susie_gwas_mqtl)
#     }else{
#       PP4_susie_gwas_mqtl=NA
#     }
#   } else if( args$type=="quant" ) {
#     # coloc.abf:
#     try(my.res <- coloc.abf(dataset1=list(beta=stats$beta.x, varbeta=stats$se ^ 2, snp=stats$rsid, type="quant"),
#                             dataset2=list(beta=stats$beta.y, varbeta=stats$varbeta, snp=stats$rsid, sdY=sdY[which(sdY$cg_id==cpg),'sdY'], type="quant"),
#                             p1=args$p1,p2=args$p2,p12=args$p12))
#     
#     # coloc.susie:
#     LD = cor(geno[,match(unlist(stats$rsid), colnames(geno))],use='pairwise.complete.obs')
#     s_gwas = runsusie(list(beta=stats$beta.x, varbeta=stats$se ^ 2, snp=stats$rsid, type="quant", LD=LD, N=args$cases/args$gwas_sample_size))
#     s_qtl = runsusie(list(beta=stats$beta, varbeta=stats$varbeta, snp=stats$rsid, sdY=sdY[which(sdY$cg_id==cpg),'sdY'], type="quant", LD=LD, N=args$qtl_sample_size))
#     if(sum(dim(summary(s_gwas)$cs)[1]>0) & sum(dim(summary(s_qtl)$cs)[1]>0)){
#       susie_gwas_mqtl=coloc.susie(s_gwas,s_qtl)
#       PP4_susie_gwas_mqtl = susie_gwas_mqtl$summary
#       print(PP4_susie_gwas_mqtl)
#     }else{
#       PP4_susie_gwas_mqtl=NA
#     }
#   }
#   
#   if (!exists(deparse(substitute(my.res)))) {print(paste0(cpg," failed"));next}
#   cpgs_tested=c(cpgs_tested,cpg)
#   
#   
#   my.res$results=my.res$results[order(as.numeric(as.character(my.res$results$snp))),]
#   rownames(my.res$results)=my.res$results$snp
#   my.res$results=my.res$results[order(my.res$results$SNP.PP.H4,decreasing=T),]
#   variant=rownames(my.res$results)[1]
#   names(variant)='topColoc_variant'
#   SNP.PP.H4=my.res$results$SNP.PP.H4[1]
#   names(SNP.PP.H4)='SNP.PP.H4'
#   my.res.all.top=rbind(my.res.all.top,c(variant,SNP.PP.H4))
#   my.res.all=rbind(my.res.all,my.res$summary)	
#   rm(my.res)
#   
#   PP4_susie_gwas_mqtl.all = rbind(PP4_susie_gwas_mqtl.all,PP4_susie_gwas_mqtl)
# }
# 
# my.res.all=my.res.all[-1,]
# rownames(my.res.all)=as.character(cpgs_tested)
# my.res.all.top=my.res.all.top[-1,]
# rownames(my.res.all.top)=as.character(cpgs_tested)
# my.res.all=cbind(my.res.all,my.res.all.top)
# PP4_susie_gwas_mqtl.all = na.omit(PP4_susie_gwas_mqtl.all)
# write.table(my.res.all, file='/data/projects/China_Brain_MultiOmics/methylation/results/coloc/cbmap/cbmap_AD_coloc_chr22_res.txt', quote=F,sep='\t')
# write.table(PP4_susie_gwas_mqtl.all, file='/data/projects/China_Brain_MultiOmics/methylation/results/coloc/cbmap/cbmap_AD_susie_chr22_res.txt', quote=F,sep='\t')
# t2 = Sys.time()
# print(t2-t1)


### result summary
# library(purrr)
# 
# for (i in 1:22) {
#   if (study == 'CBMAP') {
#     coloc_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/coloc/result/cbmap/coloc/'
#     load(paste0(coloc_res_dir, study, '_SCZ_coloc_chr', i, '.RData'))
#   } else if (study == 'ROSMAP') {
#     coloc_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/coloc/result/rosmap/coloc/'
#     load(paste0(coloc_res_dir, study, '_SCZ_coloc_chr', i, '.RData'))
#   } else if (study == 'NSPT') {
#     coloc_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/coloc/result/nspt/coloc/'
#     load(paste0(coloc_res_dir, study, '_SCZ_coloc_chr', i, '.RData'))
#   }
#   filtered_list <- keep(my.res.all, ~ !all(is.na(.x$summary)))
#   my.res.all1 <- filtered_list
#   save(my.res.all1, file=paste0(coloc_res_dir, study, '_SCZ_coloc_chr', i, '.RData'))
# }
