#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH --job-name=medi
#SBATCH --partition=mastercpu
#SBATCH --array=1-80
taskplugin=task/affinity


Rscript /data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL-mCpG/script/SMR/archive/CBMAP_cis_SMR_SME_test.R ${SLURM_ARRAY_TASK_ID}
