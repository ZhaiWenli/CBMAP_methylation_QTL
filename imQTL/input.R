library(data.table)
library(dplyr)

# CBMAP
cell.pro=read.table('/data/projects/China_Brain_MultiOmics/methylation/results/ctp_analysis/CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt')[,1:7]
pdat <- read.csv('/data/projects/China_Brain_MultiOmics/methylation/data/CBMAP/DNAm_processed/phenotype.csv', header=T)
idkey = as.data.frame(fread('/data/shared_data/China_Brain_MultiOmics/WGS/firstpass_20241126/WGS_firstpass_sample_id_info_20241127.txt'))
rownames(cell.pro) <- pdat[match(rownames(cell.pro),pdat$barcode_id), 'Sample_Name']
sample = intersect(rownames(cell.pro),idkey$sample_name)
cell.pro = cell.pro[sample,]
rownames(cell.pro) <- idkey[match(rownames(cell.pro),idkey$sample_name), 'WGS_sample_name']


# methylation level
for (i in 1:22) {
  methy = fread(paste0('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cg_BED_chr/cg_BED_chr',i,'.bed.gz'))  
  sample = intersect(rownames(cell.pro),colnames(methy))
  methy0 = methy[,..sample]
  methy = cbind(methy[,1:4],methy0)
  methy = distinct(methy,pid,.keep_all = T)
  write.table(methy, paste0('/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/input/methy_chr/methy_chr_',i,'.bed'), quote =F, sep = '\t', row.names = F, col.names = T)
  
  bgzip_cmd <- paste0('bgzip -f /data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/input/methy_chr/methy_chr_',i,'.bed && tabix -p bed /data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/input/methy_chr/methy_chr',i,'.bed.gz')  
  system(bgzip_cmd, wait=T)     #index
}
sample = data.frame('V1'=sample,'V2'=sample)
write.table(sample,file='/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/input/sample.txt',sep='\t',col.names=F,row.names=F,quote=F)
# genotype pgen/pvar/psam format
plink_cmd <- paste0("/data/tools/plink_v2/plink2 --bfile /data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/sample_geno ",
                    "--make-pgen ",
                    "--keep /data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/input/sample.txt ",
                    "--out /data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/input/geno")                                              
system(plink_cmd, wait=T)                                                                               

# cov
sample = sample$V1
cov = fread('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cov.txt',header=T)
cov.matrix = as.data.frame(t(cov[c(1:3,5:20),..sample]))
colnames(cov.matrix) = cov$id[-4]
cov.matrix[,c(1,2,4:ncol(cov.matrix))] = apply(cov.matrix[,c(1,2,4:ncol(cov.matrix))],2,as.numeric)
cov.matrix$bank = factor(cov.matrix$bank)
cov.matrix$sex = factor(cov.matrix$sex)
cov.matrix = t(cov.matrix)
write.table(cov.matrix,file='/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/input/cov.txt',sep='\t',col.names = T,row.names = T, quote = F)  

# cell proportion
cell.pro = cell.pro[sample,1:7]
colnames(cell.pro) = c('Exc','Inh','Astro','Endo','Micro','Oligo','OPC')
for (i in 1:ncol(cell.pro)) {
  a = as.data.frame(cell.pro[,i])
  rownames(a) = rownames(cell.pro)
  write.table(a,file=paste0('/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/input/',colnames(cell.pro)[i],'.txt'),col.names=F,row.names=T,quote=F,sep='\t')   
}


# ROSMAP
cell.pro=read.table('/data/projects/China_Brain_MultiOmics/methylation/results/ctp_analysis/ROSMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt')[,1:7]
pdat=read.table('/data/projects/China_Brain_MultiOmics/methylation/data/ROSMAP_methylaition/phenotype.txt',header=T)
idkey = read.csv("/data/shared_data/ROSMAP/ROSMAP_IDkey.csv", header=T)
rownames(cell.pro) <- pdat[match(rownames(cell.pro),rownames(pdat)), 'projid']
cell.pro$gwas_id <- idkey[match(rownames(cell.pro),idkey$projid), 'gwas_id']


# methylation level
for (i in 1:22) {
  methy = fread(paste0('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/ROSMAP_mQTL/mQTL_mapping/b37/QTLtools_input/cg_BED_chr/cg_BED_chr',i,'.bed.gz'))  
  sample = intersect(colnames(methy),cell.pro$gwas_id) #582
  methy0 = methy[,..sample]
  methy = cbind(methy[,1:4],methy0)
  methy = distinct(methy,pid,.keep_all = T)
  write.table(methy, paste0('/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/input/methy_chr/ROSMAP_methy_chr_',i,'.bed'), quote =F, sep = '\t', row.names = F, col.names = T)
  
  bgzip_cmd <- paste0('bgzip -f /data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/input/methy_chr/ROSMAP_methy_chr_',i,'.bed && tabix -p bed /data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/input/methy_chr/ROSMAP_methy_chr',i,'.bed.gz')  
  system(bgzip_cmd, wait=T)     #index
}
sample = data.frame('V1'=sample,'V2'=sample)
write.table(sample,file='/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/input/ROSMAP_sample.txt',sep='\t',col.names=F,row.names=F,quote=F)
# genotype pgen/pvar/psam format
plink_cmd <- paste0("/data/tools/plink_v2/plink2 --bfile /data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/ROSMAP_mQTL/mQTL_mapping/b37/QTLtools_input/pca_selection/geno_pca/sample_geno ",
                    "--make-pgen ",
                    "--keep /data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/input/ROSMAP_sample.txt ",
                    "--out /data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/input/ROSMAP_geno")                                              
system(plink_cmd, wait=T)                                                                               

# cov
sample = sample$V1
cov = fread('/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/ROSMAP_mQTL/mQTL_mapping/b37/QTLtools_input/cov.txt',header=T)
cov.matrix = as.data.frame(t(cov[c(1:3,5:20),..sample]))
colnames(cov.matrix) = cov$id[-4]
rown = rownames(cov.matrix)
cov.matrix = apply(cov.matrix,2,as.numeric)
rownames(cov.matrix) = rown
cov.matrix = t(cov.matrix)
write.table(cov.matrix,file='/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/input/ROSMAP_cov.txt',sep='\t',col.names = T,row.names = T, quote = F)  

# cell proportion
cell.pro = cell.pro[which(cell.pro$gwas_id %chin% sample),]
rownames(cell.pro) = cell.pro$gwas_id
cell.pro = cell.pro[sample,1:7]
colnames(cell.pro) = c('Exc','Inh','Astro','Endo','Micro','Oligo','OPC')
for (i in 1:ncol(cell.pro)) {
  a = as.data.frame(cell.pro[,i])
  rownames(a) = rownames(cell.pro)
  write.table(a,file=paste0('/data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/input/ROSMAP_',colnames(cell.pro)[i],'.txt'),col.names=F,row.names=T,quote=F,sep='\t')   
}
