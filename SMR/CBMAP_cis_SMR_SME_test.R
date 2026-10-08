######## mediation analysis SNP---Expression---Methylation ##########
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide the array index (n)")
}
i <- as.integer(args[1])

# .libPaths(c(.libPaths(),"/share/data/R4.3_lib/library","/share/home/zhaiwl/R/x86_64-pc-linux-gnu-library/4.3"))
library(data.table)
library(BEDMatrix)
library(stringr)
library(dplyr)
library(bruceR)
library(readxl)
library(rhdf5)


# loading data
# pairs to be analyzed
sem_res = fread('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/result/cis_SMR_SEM_res.txt',header=T)
chunk_size <- 1300
chunks <- split(sem_res, ceiling(seq_along(1:nrow(sem_res)) / chunk_size))
sem_res0 <- chunks[[i]];rm(sem_res);rm(chunks)

cpgs = unique(sem_res0$Outco_ID)
snps = unique(sem_res0$topSNP)
genes = unique(sem_res0$Expo_ID)
gene_ids = unlist(str_split_fixed(genes,"_",n=2)[,2])
gene_names = unlist(str_split_fixed(genes,"_",n=2)[,1])

# RNA-seq data
exp0 = fread('/data/projects/China_Brain_MultiOmics/humanBrain_RNAseq/CBMAP_RNAseq_protein_coding/RNAseq_process/RNAseq_data_qn_comb.txt',header=T)
exp = as.data.frame(t(exp0))
colnames(exp) = exp[1,]
exp = exp[-1,]
exp = apply(exp,2,as.numeric)
rownames(exp) = colnames(exp0)[2:ncol(exp0)]
rm(exp0)
# methylation data
methy <- h5read("/data/projects/China_Brain_MultiOmics/methylation/data/CBMAP/DNAm_processed/final_DNAm_invMdat.h5",'Mdat')
colnames(methy) = h5read("/data/projects/China_Brain_MultiOmics/methylation/data/CBMAP/DNAm_processed/final_DNAm_invMdat.h5",'colnames')
rownames(methy) = h5read("/data/projects/China_Brain_MultiOmics/methylation/data/CBMAP/DNAm_processed/final_DNAm_invMdat.h5",'rownames')
methy = methy[,cpgs]

# genotype data
geno = BEDMatrix('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/sample_geno.bed')
colnames(geno) = str_split_fixed(colnames(geno),"_",n=2)[,1]
rownames(geno) = str_split_fixed(rownames(geno),"_",n=2)[,1]
geno = geno[,snps]
idkey = as.data.frame(fread('/data/shared_data/China_Brain_MultiOmics/WGS/firstpass_20241126/WGS_firstpass_sample_id_info_20241127.txt'))
rownames(geno) = unlist(idkey[match(rownames(geno),idkey$WGS_sample_name),'sample_name'])
# covariants
cov = fread('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cov.txt',header=T)
cov = as.data.frame(t(cov))
colnames(cov) = cov[1,]
cov = cov[-1,]
cov[,c(1,4:ncol(cov))] = apply(cov[,c(1,4:ncol(cov))],2,as.numeric)
cov$sex = factor(cov$sex)
cov$bank = factor(cov$bank)
rownames(cov) = unlist(idkey[match(rownames(cov),idkey$WGS_sample_name),'sample_name'])
pmd_rin = read.csv('/data/shared_data/China_Brain_MultiOmics/sample_information/final/CBMAP_sample_info_1187_final_20250523.csv',header=T)
cov$PMD = unlist(pmd_rin[match(rownames(cov),pmd_rin$id),'PMD'])
cov$RIN = unlist(pmd_rin[match(rownames(cov),pmd_rin$id),'RIN'])

sample = intersect(rownames(exp),rownames(geno)) # 823 samples
exp = exp[sample,gene_ids] # 823 15621
geno = geno[sample,] # 823 283052
methy = methy[sample,] # 823 316834
cov = cov[sample,c(1,2,3,4,21,22)] # 823 22


##### mediation analysis #####
# library(ggplot2)
# library(data.table)

mediation_res = data.frame('cpg'=NA, 'gene'=NA, 'snp'=NA, 'Mediated.Prop'=NA, 'Indirect_beta'=NA, 'Indirect_P'=NA, 'Direct_beta'=NA, 'Direct_P'=NA, 'Total_beta'=NA, 'Total_P'=NA)
o = 0
for (j in 1:nrow(sem_res0)) {
  cpg = sem_res0$Outco_ID[j]
  gene = sem_res0$Expo_ID[j]
  gene_id = unlist(str_split_fixed(gene,"_",n=2)[,2])
  gene_name = unlist(str_split_fixed(gene,"_",n=2)[,1])
  snp = sem_res0$topSNP[j]
  dat = data.frame('methy'=methy[,cpg], 'expr'=exp[,gene_id], 'genotype'=geno[,snp],
                   'sex'=cov$sex, 'age'=cov$age, 'NeuN_pos'=cov$NeuN_pos, 'bank'=cov$bank, 'PMD'=cov$PMD, 'RIN'=cov$RIN)
  contcont <- PROCESS(data=dat, y='expr', x='genotype', meds='methy', covs=c('sex','age','NeuN_pos','bank','PMD','RIN'), nsim=1000, seed=1, digits=7)
  o <- o + 1
  mediation_res[o,1] = cpg
  mediation_res[o,2] = gene
  mediation_res[o,3] = snp
  mediation_res[o,4]= contcont$results[[1]]$mediation[1,1] / contcont$results[[1]]$mediation[3,1]
  mediation_res[o,5]= contcont$results[[1]]$mediation[1,1]
  mediation_res[o,6]= contcont$results[[1]]$mediation[1,'pval']
  mediation_res[o,7]= contcont$results[[1]]$mediation[2,1]
  mediation_res[o,8]= contcont$results[[1]]$mediation[2,'pval']
  mediation_res[o,9]= contcont$results[[1]]$mediation[3,1]
  mediation_res[o,10]= contcont$results[[1]]$mediation[3,'pval']
}
fwrite(mediation_res,file=paste0('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/result/SME_mediation_chunk_res/SME_mediation_chunk_',i,'.txt'),col.names=T,row.names=F,quote=F,sep='\t')


# mediation results summary
res = data.table()
files = list.files('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/result/SME_mediation_chunk_res',full.names = T)
for (file in files) {
  res0 = fread(file,header=T)
  res = rbind(res,res0)
}
res$fdr = p.adjust(res$Indirect_P,'BH')
sum(res$fdr<0.05) # 2008 P-threshold 0.0009669499
fwrite(res,file='/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/result/SME_mediation_res.txt',quote=F,sep='\t')

