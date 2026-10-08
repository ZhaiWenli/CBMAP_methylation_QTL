##### QTLtools file prepare #####
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(optparse)
library(readxl)
library(rhdf5)
library(BEDMatrix)
library(stringr)

qtltool_input_dir = '/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/'


#### prepare cg_BED #Chr start end pid gid strand sampleid
# get cpg data
cg=h5read('/methylation/data/CBMAP/DNAm_processed/final_DNAm_invMdat.h5','Mdat')
colnames(cg) = h5read('/methylation/data/CBMAP/DNAm_processed/final_DNAm_invMdat.h5','colnames')
rownames(cg) = h5read('/methylation/data/CBMAP/DNAm_processed/final_DNAm_invMdat.h5','rownames')

# get genotype
geno = BEDMatrix("/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/genotype_processed/firstpass_20241126_03.bed")
rownames(geno) = str_split_fixed(rownames(geno), '_', n=2)[,2]

# match sample ID
idkey = as.data.frame(fread('WGS_firstpass_sample_id_info_20241127.txt'))
sampleid = intersect(rownames(cg), idkey$sample_name)
idkey = idkey[which(idkey$sample_name %chin% sampleid),]
geno = geno[idkey$WGS_sample_name,]
rownames(geno) = idkey[match(rownames(geno), idkey$WGS_sample_name), 'sample_name']
sampleid = intersect(rownames(cg), rownames(geno)) #1018 length

# cg quantification
cg_tmp = as.data.frame(t(cg[sampleid,]))
cg_tmp = cbind(rownames(cg_tmp),cg_tmp, stringsAsFactors=F)
colnames(cg_tmp)[1] = 'pid'
rownames(cg_tmp) = NULL
colnames(cg_tmp)[2:ncol(cg_tmp)] = idkey[match(colnames(cg_tmp)[2:ncol(cg_tmp)], idkey$sample_name),'WGS_sample_name']

cg_info <- read.csv('/methylation/935k_annotation/EPIC-8v2-0_A1.csv',header=F)
cg_info <- cg_info[-c(1:7),]
colnames(cg_info) <- cg_info[1,]
cg_info <- cg_info[-1,]
cg_info <- data.frame(pid=cg_info$Name, Chr=str_split_fixed(cg_info$CHR, 'chr', n=2)[,2], start=cg_info$MAPINFO, end=cg_info$MAPINFO, gid=cg_info$UCSC_RefGene_Name, strand=cg_info$Strand_FR)
cg_info$strand <- ifelse(cg_info$strand=='F', '+', cg_info$strand)
cg_info$strand <- ifelse(cg_info$strand=='R', '-', cg_info$strand)
cg_info[cg_info==""] <- NA
cg_info <- cg_info[which(cg_info$Chr %in% c(1:22)),]

cg_BED <- merge(cg_info, cg_tmp, by='pid', all.y=T)
cg_BED <- cg_BED[,c('Chr','start','end','pid','gid','strand', colnames(cg_tmp)[2:ncol(cg_tmp)])]
cg_BED$Chr=as.integer(cg_BED$Chr)
cg_BED$start=as.integer(cg_BED$start)
cg_BED <- cg_BED[order(cg_BED$Chr, cg_BED$start),]
cg_BED <- cg_BED[-which(is.na(cg_BED$start)==T),]
colnames(cg_BED)[1] <- '#Chr'
write.table(cg_BED, '/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/cg_BED.bed', quote =F, sep = '\t', row.names = F, col.names = T)

# group cg according to #Chr 
split_data <- split(cg_BED, cg_BED$'#Chr')
for (chr in names(split_data)) {
  chr_data <- split_data[[chr]]
  write.table(chr_data, paste0(qtltool_input_dir, 'cg_BED_chr/cg_BED_chr', chr, '.bed'), quote =F, sep = '\t', row.names = F, col.names = T)
  bgzip_cmd <- paste0('bgzip -f ',qtltool_input_dir,'cg_BED_chr/cg_BED_chr', chr, '.bed && tabix -p bed ',qtltool_input_dir,'cg_BED_chr/cg_BED_chr',chr,'.bed.gz')  
  system(bgzip_cmd, wait=T)   
}

#### prepare geno_CVF ####
#bfile to vcf
plink_cmd <- paste0("plink --bfile ",qtltool_input_dir,"geno_pca_20/sample_geno ",
                    "--recode vcf-fid ",
                    "--out ",qtltool_input_dir,"geno")                                        
system(plink_cmd, wait=T)                                                                               

#vcf to vcf.gz
bgzip_cmd <- paste0('bgzip -c ',qtltool_input_dir,'geno.vcf > ',
                    qtltool_input_dir,'geno.vcf.gz')                                     
system(bgzip_cmd, wait=T)                                                                                                      

#make index
tabix_cmd <- paste0('tabix -p vcf ',qtltool_input_dir,'geno.vcf.gz')                 
system(tabix_cmd, wait=T)


#### prepare COV file ####
#sex & bank & age
pdat = read.csv("methylation/data/CBMAP/DNAm_processed/phenotype.csv", header=T)
cov=data.frame(id=sampleid, age=NA, sex=NA, bank=NA, NeuN_pos=NA, gpca1=NA, gpca2=NA, gpca3=NA, gpca4=NA, gpca5=NA, gpca6=NA, gpca7=NA, gpca8=NA, gpca9=NA, gpca10=NA, mpca1=NA, mpca2=NA, mpca3=NA, mpca4=NA, mpca5=NA, mpca6=NA)
cov$sex = pdat[match(cov$id, pdat$Sample_Name), 'sex_male']
cov$bank = pdat[match(cov$id, pdat$Sample_Name), 'bank']
cov$age = pdat[match(cov$id, pdat$Sample_Name), 'age']

#neuron proportion
sh_celltype = read.table('/methylation/data/CBMAP/DNAm_processed/sh_DNAm_processed/celltype.txt', header=T)
zk_celltype = read.table('/methylation/data/CBMAP/DNAm_processed/zk_DNAm_processed/celltype.txt', header=T)
celltype = rbind(sh_celltype, zk_celltype)
pdat2 = read.csv('/methylation/data/CBMAP/DNAm_processed/phenotype.csv', header=T)
rownames(celltype) = pdat2[match(rownames(celltype), pdat2$barcode_id), 'Sample_Name']
celltype = celltype[sampleid,]
cov$NeuN_pos = celltype[match(cov$id, rownames(celltype)), 'NeuN_pos']

#top 5 geno_PCA
gpca = as.data.frame(fread('/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/geno_pca_20.eigenvec'))
gpca$V2 = idkey[match(gpca$V2, idkey$WGS_sample_name), 'sample_name']
cov$gpca1 = gpca[match(cov$id, gpca$V2), 'V3']
cov$gpca2 = gpca[match(cov$id, gpca$V2), 'V4']
cov$gpca3 = gpca[match(cov$id, gpca$V2), 'V5']
cov$gpca4 = gpca[match(cov$id, gpca$V2), 'V6']
cov$gpca5 = gpca[match(cov$id, gpca$V2), 'V7']
cov$gpca6 = gpca[match(cov$id, gpca$V2), 'V8']
cov$gpca7 = gpca[match(cov$id, gpca$V2), 'V9']
cov$gpca8 = gpca[match(cov$id, gpca$V2), 'V10']
cov$gpca9 = gpca[match(cov$id, gpca$V2), 'V11']
cov$gpca10 = gpca[match(cov$id, gpca$V2), 'V12']

#top 6 DNAm_PCA
mpca = as.data.frame(fread('/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/DNAm_pca_10/DNAm_pca_10.txt'))
cov$mpca1 = mpca[match(cov$id, mpca$V1), 'V2']
cov$mpca2 = mpca[match(cov$id, mpca$V1), 'V3']
cov$mpca3 = mpca[match(cov$id, mpca$V1), 'V4']
cov$mpca4 = mpca[match(cov$id, mpca$V1), 'V5']
cov$mpca5 = mpca[match(cov$id, mpca$V1), 'V6']
cov$mpca6 = mpca[match(cov$id, mpca$V1), 'V7']

cov=as.data.frame(t(cov))
colnames(cov)=cov[1,]
cov=cov[-1,]
cov = cbind(rownames(cov),cov)
colnames(cov)[1] = 'id'
colnames(cov)[2:ncol(cov)] = idkey[match(colnames(cov)[2:ncol(cov)], idkey$sample_name),'WGS_sample_name']
write.table(cov, paste0(qtltool_input_dir,'cov.txt'), quote =F, sep = '\t', row.names = F, col.names = T)
