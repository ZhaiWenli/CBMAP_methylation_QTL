#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --array=31
#SBATCH --job-name=predxc_s2
#SBATCH --partition=all
taskplugin=task/affinity


#SLURM_ARRAY_TASK_ID=1  #1 to N

main_dir="/data/projects/China_Brain_MultiOmics/methylation/mQTL/MA-FOCUS/db_files"

#step 2. generate .db and .cov file
#ml GCC OpenMPI R
Rscript /data/projects/China_Brain_MultiOmics/methylation/mQTL/MA-FOCUS/script/CBMAP_predixcan_step1.r \
        --generate_db_and_cov \
        --main_dir ${main_dir} \
        --batch_id ${SLURM_ARRAY_TASK_ID} \
        --plink_file_name /data/projects/China_Brain_MultiOmics/methylation/mQTL/mQTL_mapping/result/CBMAP_mQTL/QTLtools_input/geno_pca_20/sample_geno \
        --expression_file_name ${main_dir}/expression_file/CBMAP_cg_BED_residuals_batch_${SLURM_ARRAY_TASK_ID}.txt \
        --annotation_file_name ${main_dir}/annotation_file/CBMAP_cpg_anno_b38.txt \
        --output_file_name CBMAP_batch_${SLURM_ARRAY_TASK_ID}



