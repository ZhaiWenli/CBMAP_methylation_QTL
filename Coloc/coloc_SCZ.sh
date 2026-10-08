#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --job-name=coloc
#SBATCH --partition=nodecpu
#SBATCH --array=0-26%13
taskplugin=task/affinity

task=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" /data/projects/China_Brain_MultiOmics/methylation/mQTL/coloc/script/coloc_script/T2D_task_rosmap_list.txt)
IFS=',' read trait chr batch <<< "$task"
study="ROSMAP"

echo "Running coloc for trait=$trait, chr=$chr, batch=$batch, study=$study"

Rscript /data/projects/China_Brain_MultiOmics/methylation/mQTL/coloc/script/coloc_script/coloc_SCZ.R $trait $chr $batch $study
