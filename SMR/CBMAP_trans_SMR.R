args = as.numeric(commandArgs(TRUE))
print(args)
i = args


library(data.table)
library(dplyr)

setwd('/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input')
# cis_eqtl = fread('/humanBrain_RNAseq/CBMAP_RNAseq_protein_coding/QTLtools_res/qtltools_cis_nominal_main_p0.05.txt')   
# gene = data.frame('id'=cis_eqtl$pheno_id,
#                      'chr'=cis_eqtl$gene_chr,
#                      'start'=cis_eqtl$gene_start,
#                      'end'=cis_eqtl$gene_end)
# gene = distinct(gene,.keep_all = T)
# 
# trans_mqtl = fread("/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_result/transQTL_res.txt")   
# cpg = data.frame('id'=trans_mqtl$cg_id,
#                  'chr'=trans_mqtl$cg_chr,
#                  'start'=trans_mqtl$cg_start,
#                  'end'=trans_mqtl$cg_start)
# cpg = distinct(cpg,.keep_all = T)
# 
# trans_list <- lapply(1:nrow(gene), function(i) {
#   g <- gene[i, ]
#   
#   same_chr_trans <- cpg[
#     cpg$chr == g$chr &
#       (cpg$end < (g$start - 1e6) | cpg$start > (g$end + 1e6)), "id"
#   ]
#   
#   diff_chr_trans <- cpg[cpg$chr != g$chr, "id"]
#   
#   unique(c(same_chr_trans, diff_chr_trans))
# })
# names(trans_list) <- gene$id
# trans = unique(unlist(trans_list))
# write.table(trans,file='trans-cpg.txt',sep='\t',col.names=F,row.names=F,quote=F)

### run QTLtools for trans-CpGs
qtltools_file_dir = '/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/'
qtltools_file_res = '/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input/trans_tmp/'
trans_cpg = read.table('trans-cpg.txt')
cpg_id = trans_cpg[i,1]
write.table(cpg_id, file=paste0(qtltools_file_res,cpg_id,'.txt'),col.names=F,row.names=F,quote=F,sep='\t')
qtltools_cmd <- paste0("QTLtools trans ",
                       "--vcf ", qtltools_file_dir, "geno.vcf.gz ",
                       "--bed ", qtltools_file_dir, "cg_BED.bed.gz ", 
                       "--cov ", qtltools_file_dir, "cov.txt ",
                       "--include-phenotypes ", qtltools_file_res, cpg_id, ".txt ",
                       "--nominal --threshold 1 ",
                       "--window 1000000 ",
                       "--normal ",
                       "--out ", qtltools_file_res, cpg_id, ".txt")
system(qtltools_cmd, wait=T)

smr_file_dir = '/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input/'
smr_res_dir = '/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/result/trans_SMR/'
cpg_sd = fread('/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cpg_sdY.txt')
geno_frq = fread('/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/snp.maf.frq')  
geno_frq$sd <- sqrt(2 * geno_frq$MAF * (1 - geno_frq$MAF))
trans_res=fread(paste0(qtltools_file_res, cpg_id, ".txt.hits.txt.gz"),header=F)[,c(1,4,7,9)]
colnames(trans_res) = c('cg_id','SNP','p','cor_coef')

trans_res <- merge(trans_res, cpg_sd[, .(cg_id, sdY)], by = "cg_id", all.x = TRUE)
trans_res <- merge(trans_res, geno_frq[, .(SNP, A1, A2, MAF, sd)], by = "SNP", all.x = TRUE)
trans_res[, beta := cor_coef * (sdY / sd)]
trans_res[, z := sign(beta) * qnorm(1 - p/2)]
trans_res[, se := beta / z]
trans_res = trans_res[,c(1,6:8,10,12,3)]
colnames(trans_res) = c('SNP','A1','A2','freq','b','se','p')
trans_res$n = 1018
write.table(trans_res,file=paste0(qtltools_file_res, cpg_id, ".ma"),quote=F,col.names=T,row.names=F,sep='\t')

smr_cmd <- paste0(
  "smr ",
  "--bfile /methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/sample_geno ",
  "--gwas-summary ",qtltools_file_res, cpg_id, ".ma ",
  "--beqtl-summary ", smr_file_dir, "cis_eqtl_besd ",
  "--out ", smr_res_dir, cpg_id,"_trans_SMR ",
  "--thread-num 1")
system(smr_cmd, wait = TRUE)

setwd('/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input/trans_tmp')
files_to_delete <- list.files(pattern = cpg_id)
file.remove(files_to_delete)


