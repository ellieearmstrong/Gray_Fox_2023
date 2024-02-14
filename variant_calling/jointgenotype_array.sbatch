#!/bin/bash
#SBATCH --job-name=joint_genotype
#SBATCH --partition=normal,hns,dpetrov
#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --array=1-42
#SBATCH --cpus-per-task=12
#SBATCH --mem=120000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ellieearmstrong@gmail.com

module load biology gatk/4.1.4.1

config=dbimport_config.txt
interval=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
scaffold=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
gatk --java-options '-Xmx250g' GenotypeGVCFs -R greyfox-dovetail-final.fasta -V gendb://grayfoxdb_$scaffold -all-sites --output grayfox_allsites.${scaffold}.g.vcf -L $interval
