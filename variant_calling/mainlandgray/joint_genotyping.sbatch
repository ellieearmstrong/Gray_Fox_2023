#!/bin/bash
#SBATCH --job-name=joint_genotype
#SBATCH --partition=normal,hns,dpetrov
#SBATCH --time=36:00:00
#SBATCH -N 1
#SBATCH --array=1-38
#SBATCH --cpus-per-task=24
#SBATCH --mem=240000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ellieearmstrong@gmail.com

module load biology gatk/4.1.4.1

config=dbimport_config_wolf.txt
interval=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
scaffold=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
gatk --java-options '-Xmx250g' GenotypeGVCFs -R GCA_905319855.2_mCanLor1.2_genomic.fna -V gendb://graywolfdb_$scaffold -all-sites --output graywolf_mainland_allsites.${scaffold}.g.vcf -L $interval
