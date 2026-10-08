#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=42
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=22G
#SBATCH --job-name=predxc_s3
#SBATCH --partition=nodecpu
taskplugin=task/affinity


#SLURM_ARRAY_TASK_ID=1  #1 to N
trait="DrnkWk"

main_dir="/data/projects/China_Brain_MultiOmics/methylation/mQTL/MA-FOCUS"

#step 3. apply to gwas data (association test)
Rscript /data/projects/China_Brain_MultiOmics/methylation/mQTL/MA-FOCUS/script/CBMAP_predixcan_step1.r \
        --asso_test \
        --db_path ${main_dir}/db_files/CBMAP_output/CBMAP_Final_Combined.db \
        --cov_path ${main_dir}/db_files/CBMAP_output/CBMAP_Final_Combined.cov \
        --gwas_path ${main_dir}/GWAS_files/${trait}_EAS.txt \
        --gwas_variant_col SNP \
        --gwas_beta_col BETA \
        --gwas_se_col SE \
        --gwas_eff_allele_col A1 \
        --gwas_ref_allele_col A2 \
        --asso_out_path ${main_dir}/results/EAS_Predixcan/${trait}_MWAS.txt \
        --parallel
# If a '--parallel' flag is added, max(n-1,1) core(s) will be used for parallel association test, where n is the number of available cores.
# To specify the colname for the point estimate of GWAS effect size, use either "--gwas_or_col" or "gwas_beta_col"
# Use "--gwas_se_col" or "--gwas_p_col" to provide either the se(beta) or the p-value of beta.


