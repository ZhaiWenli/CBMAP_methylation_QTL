#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=5
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=6G
#SBATCH --job-name=mqtl
#SBATCH --partition=all
#SBATCH --array=1-2
taskplugin=task/affinity

export LD_LIBRARY_PATH=/share/home/zhaiwl/gsl-2.7/lib:$LD_LIBRARY_PATH

Rscript /data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/script/SMR/CBMAP_trans_SMR.R ${SLURM_ARRAY_TASK_ID}
