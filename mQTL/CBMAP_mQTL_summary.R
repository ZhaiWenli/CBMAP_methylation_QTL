library(data.table)
library(dplyr)
library(stringr)

############################## Cis-mQTL ###################################
# calculate cis_Pvalue
cis_num <- 0
for (i in 1:22) {
  cis_snp_counts = fread(paste0('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/cissnp/cisSNP_chr',i,'.txt'),header=T)
  cis_num <- cis_num + sum(cis_snp_counts[,2])
}
cis_num
# 3757882127
cis_Pvalue <- 0.05/cis_num
cis_Pvalue
# [1] 1.330537e-11

# extract cis-mQTL
qtltool_file_res = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_result/chr_res/'
cis_res <- data.frame(cg_id=NA, cg_chr=NA, cg_start=NA, cg_end=NA, cg_strand=NA, num_variants=NA, distance=NA,
                      variants_id=NA, variants_chr=NA, variants_start=NA, variants_end=NA,
                      nominal_P=NA, r_squared=NA,beta=NA, se=NA,top_variant=NA)
for (i in 1:22) {
  qtltools_cis_nominal <- fread(paste0(qtltool_file_res,"qtltools_cis_nominal_chr",i,".txt")) %>% as.data.frame()
  colnames(qtltools_cis_nominal) <- c("cg_id", "cg_chr", "cg_start", "cg_end", "cg_strand", "num_variants", "distance",
                                      "variants_id", "variants_chr", "variants_start", "variants_end",
                                      "nominal_P", "r_squared","beta", "se","top_variant")
  qtltools_cis_nominal <- qtltools_cis_nominal[order(qtltools_cis_nominal$nominal_P),]
  qtltools_cis_nominal <- qtltools_cis_nominal[which(qtltools_cis_nominal$nominal_P < cis_Pvalue),]
  cis_res <- rbind(cis_res, qtltools_cis_nominal)
}  
cis_res <- cis_res[-1,]
dim(cis_res)
# [1] 20409142       16
length(unique(cis_res$variants_id))
# [1] 3733823
length(unique(cis_res$cg_id))
# [1] 171038
fwrite(cis_res, file=paste0(qtltool_file_res, "cisQTL_res.txt"), quote =F, sep = '\t', row.names = F, col.names = T)


############################## Trans-mQTL ###################################
# calculate trans_Pvalue
qtltool_file_res = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_result/chr_res/'
num_snp = dim(fread('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/sample_geno.bim',header=F))[1]
num_cg = 889881
trans_Pvalue <- 0.05 / (num_snp * num_cg - cis_num)
# [1] 9.960145e-15

# extract trans-mQTL
trans_res <- data.frame(cg_id=NA, cg_chr=NA, cg_start=NA, variants_id=NA, variants_chr=NA, variants_start=NA, nominal_P=NA, cor_coef=NA)
for (i in 1:22) {
  qtltools_trans_nominal <- fread(paste0(qtltool_file_res,"qtltools_trans_nominal_chr",i,".txt.hits.txt.gz"))[,c(1:7,9)] %>% as.data.frame()
  colnames(qtltools_trans_nominal) <- c("cg_id", "cg_chr", "cg_start", "variants_id", "variants_chr", "variants_start", "nominal_P", "cor_coef")
  qtltools_trans_nominal <- qtltools_trans_nominal[order(qtltools_trans_nominal$nominal_P),]
  qtltools_trans_nominal <- qtltools_trans_nominal[which(qtltools_trans_nominal$nominal_P < trans_Pvalue),]
  trans_res <- rbind(trans_res, qtltools_trans_nominal)
}  
trans_res <- trans_res[-1,]
dim(trans_res)
# [1] 850176      8
length(unique(trans_res$variants_id))
# [1] 444253
length(unique(trans_res$cg_id))
# [1] 7614
fwrite(trans_res, file=paste0(qtltool_file_res, "transQTL_res.txt"), quote =F, sep = '\t', row.names = F, col.names = T)





