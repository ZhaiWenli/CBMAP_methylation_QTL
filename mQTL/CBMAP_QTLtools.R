args = as.numeric(commandArgs(TRUE))
print(args)
i = args

##### QTLtools analysis #####
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(optparse)

# set work directory and file names
qtltool_file_dir = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/'
qtltool_file_res = '/data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_result/chr_res/'

# cis eQTL mapping
#nominal
qtl_cmd <- paste0("/data/tools/QTLtools/qtltools/bin/QTLtools cis ",
                  "--vcf ",qtltool_file_dir,"geno.vcf.gz ",
                  "--bed ",qtltool_file_dir,"cg_BED_chr/cg_BED_chr",i,".bed.gz ",
                  "--cov ",qtltool_file_dir,"cov.txt ",
                  "--nominal 1 ",
                  "--normal --std-err --seed 2025 ",
                  "--window 1000000 ",
                  "--out ",qtltool_file_res,"qtltools_cis_nominal_chr",i,".txt")
system(qtl_cmd, wait=T)

# trans eQTL mapping
#nominal
qtl_cmd <- paste0("/data/tools/QTLtools/qtltools/bin/QTLtools trans ",
                  "--vcf ",qtltool_file_dir,"geno.vcf.gz ",
                  "--bed ",qtltool_file_dir,"cg_BED_chr/cg_BED_chr",i,".bed.gz ",
                  "--cov ",qtltool_file_dir,"cov.txt ",
                  "--nominal --threshold 5e-8 ",
                  "--window 1000000 ",
                  "--normal ",
                  "--out ",qtltool_file_res,"qtltools_trans_nominal_chr",i,".txt")
system(qtl_cmd, wait=T)

