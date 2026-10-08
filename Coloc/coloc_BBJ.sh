#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=100G
#SBATCH --job-name=coloc
#SBATCH --partition=nodecpu
#SBATCH --array=0-22%6
taskplugin=task/affinity

task=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" /methylation/mQTL/coloc/script/coloc_script/bbj_task_list2.txt)
IFS=',' read trait chr batch <<< "$task"

echo "Running coloc for trait=$trait, chr=$chr, batch=$batch"

Rscript /methylation/mQTL/coloc/script/coloc_script/coloc_BBJ.R $trait $chr $batch
