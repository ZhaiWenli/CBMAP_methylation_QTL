#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=30G
#SBATCH --job-name=coloc
#SBATCH --partition=mastercpu
#SBATCH --array=1-461%40
taskplugin=task/affinity

task=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" /data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/coloc/tmp/cbmap_task_list2.txt)
IFS=',' read trait chr cell <<< "$task"

echo "Running coloc for trait=$trait, chr=$chr, cell=$cell"

Rscript /data/projects/China_Brain_MultiOmics/methylation/mQTL/celltype_mQTL/tensorqtl/script/coloc_AD_SCZ_BMI.R $trait $chr $cell
