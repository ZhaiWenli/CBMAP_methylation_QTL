args = as.numeric(commandArgs(TRUE))
print(args)
i = args


library(data.table)
library(dplyr)

setwd('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input')
# cis_eqtl = fread('/data/projects/China_Brain_MultiOmics/humanBrain_RNAseq/CBMAP_RNAseq_protein_coding/QTLtools_res/qtltools_cis_nominal_main_p0.05.txt')   
# gene = data.frame('id'=cis_eqtl$pheno_id,
#                      'chr'=cis_eqtl$gene_chr,
#                      'start'=cis_eqtl$gene_start,
#                      'end'=cis_eqtl$gene_end)
# gene = distinct(gene,.keep_all = T)
# 
# trans_mqtl = fread("/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_result/transQTL_res.txt")   
# cpg = data.frame('id'=trans_mqtl$cg_id,
#                  'chr'=trans_mqtl$cg_chr,
#                  'start'=trans_mqtl$cg_start,
#                  'end'=trans_mqtl$cg_start)
# cpg = distinct(cpg,.keep_all = T)
# 
# trans_list <- lapply(1:nrow(gene), function(i) {
#   g <- gene[i, ]
#   
#   # 同染色体且距离 > 1Mb
#   same_chr_trans <- cpg[
#     cpg$chr == g$chr &
#       (cpg$end < (g$start - 1e6) | cpg$start > (g$end + 1e6)), "id"
#   ]
#   
#   # 不同染色体
#   diff_chr_trans <- cpg[cpg$chr != g$chr, "id"]
#   
#   unique(c(same_chr_trans, diff_chr_trans))
# })
# names(trans_list) <- gene$id
# trans = unique(unlist(trans_list))
# write.table(trans,file='trans-cpg.txt',sep='\t',col.names=F,row.names=F,quote=F)

#### Estimate the GRM from all the autosomal SNPs
# gcta_cmd <- paste0("/data/tools/gcta-1.94.1-linux-kernel-3-x86_64/gcta64 ",
#                    "--bfile /data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/sample_geno ",
#                    "--autosome --make-grm ",
#                    "--out /data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input/sample_geno")
# system(gcta_cmd, wait=T)

### make pheno files .phen
# trans_cpg = read.table('trans-cpg.txt')
# methy = fread('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cg_BED.bed.gz')  
# for (i in 1:nrow(trans_cpg)) {
#   cpg_id = trans_cpg[i,1]
#   methy_sub = data.frame('FID'=colnames(methy)[7:ncol(methy)], 'IID'=colnames(methy)[7:ncol(methy)], 
#                          'methy'=unlist(methy[which(methy$pid == cpg_id),7:ncol(methy)]))
#   write.table(methy_sub,file=paste0('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input/trans_tmp/',cpg_id,'.phen'),
#               col.names = F,row.names = F, quote = F, sep='\t')
# }

### make covar and qcovar file
# cov = fread('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cov.txt',header=T)  
# covar = data.table('FID'=colnames(cov)[2:ncol(cov)], 'IID'=colnames(cov)[2:ncol(cov)])
# covar$sex = unlist(cov[2,2:ncol(cov)])
# covar$bank = unlist(cov[3,2:ncol(cov)])
# write.table(covar, file='/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input/covar.covar',
#             quote = F, col.names = F,row.names = F, sep='\t')
# qcovar = data.table('FID'=colnames(cov)[2:ncol(cov)], 'IID'=colnames(cov)[2:ncol(cov)])
# cov_tmp = t(cov[c(1,4:20),])
# colnames(cov_tmp) = cov_tmp[1,]
# cov_tmp = cov_tmp[-1,]
# cov_tmp = as.data.frame(cov_tmp)
# cov_tmp$id = rownames(cov_tmp)
# qcovar = merge(qcovar, cov_tmp, by.x='IID', by.y='id')
# write.table(qcovar, file='/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input/qcovar.qcovar',
#             quote = F, col.names = F,row.names = F, sep='\t')

### run GCTA-MLMA for trans-CpGs
gcta_file_res = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input/trans_tmp/'
trans_cpg = read.table('trans-cpg.txt')
cpg_id = trans_cpg[i,1]
gcta_cmd <- paste0("/data/tools/gcta-1.94.1-linux-kernel-3-x86_64/gcta64 ",
                   "--mlma ",
                   "--bfile /data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/sample_geno ",
                   "--grm /data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input/sample_geno ",
                   "--pheno ", gcta_file_res, cpg_id, ".phen ",
                   "--covar /data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input/covar.covar ",
                   "--qcovar /data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input/qcovar.qcovar ",
                   "--mlma-no-preadj-covar ",
                   "--out ", gcta_file_res, cpg_id, " ",
                   "--thread-num 5"
                   )
system(gcta_cmd, wait=T)

# # 读取trans结果
smr_file_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input/'
smr_res_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/result/trans_SMR/'
trans_res=fread(paste0(gcta_file_res, cpg_id, ".mlma"),header=T)[,c(2,4:9)]
colnames(trans_res) = c('SNP','A1','A2','freq','b','se','p')
trans_res$n = 1018
write.table(trans_res,file=paste0(gcta_file_res, cpg_id, ".ma"),quote=F,col.names=T,row.names=F,sep='\t')

smr_cmd <- paste0(
  "/data/tools/SMR/smr-1.4.0-linux-x86_64/smr ",
  "--bfile /data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/sample_geno ",
  "--gwas-summary ",gcta_file_res, cpg_id, ".ma ",
  "--beqtl-summary ", smr_file_dir, "cis_eqtl_besd ",
  "--out ", smr_res_dir, cpg_id,"_trans_SMR ",
  "--thread-num 5")
system(smr_cmd, wait = TRUE)

# setwd('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/result/CBMAP_SMR/input/trans_tmp')
# files_to_delete <- list.files(pattern = cpg_id)
# file.remove(files_to_delete)

