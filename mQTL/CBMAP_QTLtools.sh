#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --job-name=mqtl
#SBATCH --partition=master
#SBATCH --array=1-22
taskplugin=task/affinity



Rscript /data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/script/CBMAP_mQTL_mapping/CBMAP_QTLtools.R ${SLURM_ARRAY_TASK_ID}
