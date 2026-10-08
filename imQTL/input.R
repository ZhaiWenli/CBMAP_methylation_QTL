library(data.table)
library(dplyr)


cell.pro=read.table('/methylation/results/ctp_analysis/CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt')[,1:7]
pdat <- read.csv('/methylation/data/CBMAP/DNAm_processed/phenotype.csv', header=T)
idkey = as.data.frame(fread('WGS_firstpass_sample_id_info_20241127.txt'))
rownames(cell.pro) <- pdat[match(rownames(cell.pro),pdat$barcode_id), 'Sample_Name']
sample = intersect(rownames(cell.pro),idkey$sample_name)
cell.pro = cell.pro[sample,]
rownames(cell.pro) <- idkey[match(rownames(cell.pro),idkey$sample_name), 'WGS_sample_name']


# methylation
for (i in 1:22) {
  methy = fread(paste0('/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cg_BED_chr/cg_BED_chr',i,'.bed.gz'))  
  sample = intersect(rownames(cell.pro),colnames(methy))
  methy0 = methy[,..sample]
  methy = cbind(methy[,1:4],methy0)
  methy = distinct(methy,pid,.keep_all = T)
  write.table(methy, paste0('/methylation/mQTL/celltype_mQTL/tensorqtl/input/methy_chr/methy_chr_',i,'.bed'), quote =F, sep = '\t', row.names = F, col.names = T)
  
  bgzip_cmd <- paste0('bgzip -f /methylation/mQTL/celltype_mQTL/tensorqtl/input/methy_chr/methy_chr_',i,'.bed && tabix -p bed /methylation/mQTL/celltype_mQTL/tensorqtl/input/methy_chr/methy_chr',i,'.bed.gz')  
  system(bgzip_cmd, wait=T)     #index
}
sample = data.frame('V1'=sample,'V2'=sample)
write.table(sample,file='/methylation/mQTL/celltype_mQTL/tensorqtl/input/sample.txt',sep='\t',col.names=F,row.names=F,quote=F)
# genotype pgen/pvar/psam format
plink_cmd <- paste0("plink2 --bfile /methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/sample_geno ",
                    "--make-pgen ",
                    "--keep /methylation/mQTL/celltype_mQTL/tensorqtl/input/sample.txt ",
                    "--out /methylation/mQTL/celltype_mQTL/tensorqtl/input/geno")                                              
system(plink_cmd, wait=T)                                                                               

# cov
sample = sample$V1
cov = fread('/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cov.txt',header=T)
cov.matrix = as.data.frame(t(cov[c(1:3,5:20),..sample]))
colnames(cov.matrix) = cov$id[-4]
cov.matrix[,c(1,2,4:ncol(cov.matrix))] = apply(cov.matrix[,c(1,2,4:ncol(cov.matrix))],2,as.numeric)
cov.matrix$bank = factor(cov.matrix$bank)
cov.matrix$sex = factor(cov.matrix$sex)
cov.matrix = t(cov.matrix)
write.table(cov.matrix,file='/methylation/mQTL/celltype_mQTL/tensorqtl/input/cov.txt',sep='\t',col.names = T,row.names = T, quote = F)  

# cell proportion
cell.pro = cell.pro[sample,1:7]
colnames(cell.pro) = c('Exc','Inh','Astro','Endo','Micro','Oligo','OPC')
for (i in 1:ncol(cell.pro)) {
  a = as.data.frame(cell.pro[,i])
  rownames(a) = rownames(cell.pro)
  write.table(a,file=paste0('/methylation/mQTL/celltype_mQTL/tensorqtl/input/',colnames(cell.pro)[i],'.txt'),col.names=F,row.names=T,quote=F,sep='\t')   
}
