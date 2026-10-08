args = as.numeric(commandArgs(TRUE))
print(args)
i = args


library(data.table)
library(stringr)
library(dplyr)

setwd('/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input/tmp')
smr_file_dir = '/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input/'
smr_res_dir = '/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/result/'

############ cis-eQTL input ###########
eqtl_res=fread('/humanBrain_RNAseq/CBMAP_RNAseq_protein_coding/QTLtools_res/qtltools_cis_nominal_main_p1.txt')
eqtl_res = eqtl_res[,c(1:12,14,16)]
eqtl_res=distinct(eqtl_res,pheno_id,variants_id,.keep_all=T)
geno_frq = fread('/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input/cis_eqtl_frq.frq')
snp = intersect(geno_frq$SNP, eqtl_res$variants_id)
eqtl_res = eqtl_res[which(eqtl_res$variants_id %chin% snp),]
geno_frq = geno_frq[which(geno_frq$SNP %chin% snp),]
fwrite(eqtl_res,file='cis_eqtl_qtltools.txt',col.names = F, row.names=F, quote=F,sep='\t')

smr_cmd <- paste0("smr ",
                  "--eqtl-summary ",smr_file_dir,"cis_eqtl_qtltools.txt ",
                  "--qtltools-nominal-format --make-besd ",
                  "--out ",smr_file_dir,"tmp/cis_eqtl_besd")
system(smr_cmd, wait=T)

# update epi file
epi=fread('cis_eqtl_besd.epi')
epi$V5 = str_split_fixed(epi$V2,'_',n=2)[,1]
fwrite(epi,file='cis_eqtl_besd.epi',col.names = F, row.names=F, quote=F,sep='\t')
# update esi file
esi = fread('cis_eqtl_besd.esi')
ind = esi$V2
esi = merge(esi, geno_frq[,2:5], by.x='V2', by.y='SNP', all.x=T)
esi = esi[,c(2,1,3,4,8:10)]
rownames(esi) = esi$V2
esi = esi[ind,]
fwrite(esi,file='cis_eqtl_besd.esi',col.names = F, row.names=F, quote=F,sep='\t')

smr_cmd <- paste0("smr ",
                  "--beqtl-summary ",smr_file_dir,"cis_eqtl_besd ",
                  "--update-esi ",smr_file_dir,"cis_eqtl_besd.esi")
system(smr_cmd, wait=T)
smr_cmd <- paste0("smr ",
                  "--beqtl-summary ",smr_file_dir,"cis_eqtl_besd ",
                  "--update-epi ",smr_file_dir,"cis_eqtl_besd.epi")
system(smr_cmd, wait=T)


############ cis-mQTL input ###########
# mqtl_res=fread(paste0('/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_result/chr_res/qtltools_cis_nominal_chr',i,'.txt'))
# mqtl_res = mqtl_res[,c(1:12,14,16)]
# mqtl_res = distinct(mqtl_res,V1,V8,.keep_all=T)
# geno_frq = fread('/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/snp.maf.frq')
# snp = intersect(geno_frq$SNP, mqtl_res$V8)
# mqtl_res = mqtl_res[which(mqtl_res$V8 %chin% snp),]
# geno_frq = geno_frq[which(geno_frq$SNP %chin% snp),]
# fwrite(mqtl_res,file=paste0(smr_file_dir,'cis_mqtl_chr/cis_mqtl_chr',i,'_qtltools.txt'),sep='\t',col.names = F, row.names=F, quote=F)
# 
# smr_cmd <- paste0("smr ",
#                   "--eqtl-summary ",smr_file_dir,"cis_mqtl_chr/cis_mqtl_chr", i, "_qtltools.txt ",
#                   "--qtltools-nominal-format --make-besd ",
#                   "--out ",smr_file_dir,"cis_mqtl_chr/cis_mqtl_chr",i,"_besd")
# system(smr_cmd, wait=T)
# 
# # update epi file
# epi=fread(paste0(smr_file_dir,"cis_mqtl_chr/cis_mqtl_chr",i,"_besd.epi"))
# epi$V5 = epi$V2
# fwrite(epi,file=paste0(smr_file_dir,"cis_mqtl_chr/cis_mqtl_chr",i,"_besd.epi"),col.names = F, row.names=F, quote=F,sep='\t')
# # update esi file
# esi = fread(paste0(smr_file_dir,"cis_mqtl_chr/cis_mqtl_chr",i,"_besd.esi"))
# esi = merge(esi, geno_frq[,2:5], by.x='V2', by.y='SNP', all.x=T)
# esi = esi[,c(2,1,3,4,8:10)]
# fwrite(esi,file=paste0(smr_file_dir,"cis_mqtl_chr/cis_mqtl_chr",i,"_besd.esi"),col.names = F, row.names=F, quote=F,sep='\t')
# 
# smr_cmd <- paste0("smr ",
#                   "--beqtl-summary ",smr_file_dir,"cis_mqtl_chr/cis_mqtl_chr",i,"_besd ",
#                   "--update-esi ",smr_file_dir,"cis_mqtl_chr/cis_mqtl_chr",i,"_besd.esi")
# system(smr_cmd, wait=T)
# smr_cmd <- paste0("smr ",
#                   "--beqtl-summary ",smr_file_dir,"cis_mqtl_chr/cis_mqtl_chr",i,"_besd ",
#                   "--update-epi ",smr_file_dir,"cis_mqtl_chr/cis_mqtl_chr",i,"_besd.epi")
# system(smr_cmd, wait=T)

########### Run cis-SMR ##########
smr_cmd <- paste0("/data/tools/SMR/smr-1.4.0-linux-x86_64/smr ",
                  "--bfile ","/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/sample_geno ",
                  "--beqtl-summary ",smr_file_dir,"cis_mqtl_chr/cis_mqtl_chr",i,"_besd ",
                  "--beqtl-summary ",smr_file_dir,"cis_eqtl_besd ",
                  "--thread-num 3 ",
                  "--out ",smr_res_dir,"cis_SMR/cis_SMR_chr",i)
system(smr_cmd, wait=T)
