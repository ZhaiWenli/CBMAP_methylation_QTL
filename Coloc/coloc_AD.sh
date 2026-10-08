#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=100G
#SBATCH --job-name=coloc
#SBATCH --partition=mastercpu
#SBATCH --array=18-229%10
taskplugin=task/affinity

task=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" /methylation/mQTL/coloc/script/coloc_script/AD_task_cbmap_list.txt)
IFS=',' read trait chr batch <<< "$task"

echo "Running coloc for trait=$trait, chr=$chr, batch=$batch"
export PATH=/share/home/apps/R4.3.3/bin:$PATH
Rscript --vanilla -e '.libPaths(c("/usr/local/lib/R/site-library", .libPaths())); .libPaths(c("/usr/local/lib/R/site-library", .libPaths())); .libPaths(c("/usr/lib/R/library", .libPaths())); source("/methylation/mQTL/coloc/script/coloc_script/coloc_AD.R")' $trait $chr $batch
