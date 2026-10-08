#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --array=151-197
#SBATCH --job-name=mafocus
#SBATCH --partition=all
taskplugin=task/affinity

module load R4.3.3
module load anaconda3
traits=('AD' 'SCZ' 'BMI' 'T2D' 'SmkInit' 'DrnkWk' 'CigDay' 'SmkCes' 'AgeSmk')
task_id=$SLURM_ARRAY_TASK_ID
trait_idx=$(( task_id / 22 ))
chr=$(( task_id % 22 + 1 ))
trait=${traits[$trait_idx]}
main_dir='/methylation/mQTL/MA-FOCUS'
source activate /anaconda3/envs/ma-focus
cd /methylation/mQTL/MA-FOCUS/results

focus finemap ${main_dir}/GWAS_files/${trait}_EUR.cleaned.sumstats.gz:${main_dir}/GWAS_files/${trait}_EAS.cleaned.sumstats.gz \
              ${main_dir}/LD_files/1000GP3_multiPop_allelesAligned/EUR/1000G.EUR.QC.allelesAligned.${chr}:${main_dir}/LD_files/1000GP3_multiPop_allelesAligned/EAS/1000G.EAS.QC.allelesAligned.${chr} \
              ${main_dir}/db_files/ROSMAP_output/ROSMAP_focus_aligned.db:${main_dir}/db_files/CBMAP_output/CBMAP_focus_aligned.db \
              --locations 38:EUR-EAS \
              --chr ${chr} \
              --prior-prob ${main_dir}/gencode_file/anno3.tsv \
              --out ${main_dir}/results/${trait}_mafocus.chr${chr}

echo "Finished task $task_id: $trait CHR $chr."
conda deactivate

