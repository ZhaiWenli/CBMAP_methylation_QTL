#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=100G
#SBATCH --job-name=imqtl
#SBATCH --partition=mastercpu
#SBATCH --array=21
taskplugin=task/affinity


CHR_LIST=({1..22})
CELLTYPES=("Endo" "Exc" "Inh" "Micro" "Oligo" "OPC" "Astro")

CHR_INDEX=$((SLURM_ARRAY_TASK_ID % 22))
CELL_INDEX=$((SLURM_ARRAY_TASK_ID / 22))

CHR=${CHR_LIST[$CHR_INDEX]}
CELL=${CELLTYPES[$CELL_INDEX]}

GENO_PREFIX="/methylation/mQTL/celltype_mQTL/tensorqtl/input/geno"
PHENO_FILE="/methylation/mQTL/celltype_mQTL/tensorqtl/input/methy_chr/methy_chr_${CHR}.bed.gz"
COV_FILE="/methylation/mQTL/celltype_mQTL/tensorqtl/input/cov2.txt"
INTERACT_FILE="/methylation/mQTL/celltype_mQTL/tensorqtl/input/${CELL}.txt"

OUT_PREFIX="/methylation/mQTL/celltype_mQTL/tensorqtl/results/${CELL}_chr${CHR}"

echo "Running TensorQTL for cell type ${CELL}, chromosome ${CHR}"

source /conda/etc/profile.d/conda.sh
conda activate /conda/envs/tensorqtl
export R_HOME=/conda/envs/tensorqtl/lib/R

python -m tensorqtl ${GENO_PREFIX} ${PHENO_FILE} ${OUT_PREFIX} \
--covariates ${COV_FILE} \
--interaction ${INTERACT_FILE} \
--mode cis_nominal

conda deactivate
