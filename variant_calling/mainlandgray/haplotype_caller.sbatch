#!/bin/bash
#SBATCH --job-name=gfwolf-hc
#SBATCH --partition=gpu,owners
#SBATCH --time=12:00:00
#SBATCH -C GPU_MEM:16GB
#SBATCH -N 1
#SBATCH -G 1
#SBATCH --array=3-42
#SBATCH --cpus-per-task=10
#SBATCH --mem=120000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ellieearmstrong@gmail.com

config=grayfox_wolfconfig.txt
output=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
vcf=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)
singularity exec --nv /oak/stanford/groups/dpetrov/armstrong/applications/parabricks-4.1.1-1.simg pbrun haplotypecaller --ref GCA_905319855.2_mCanLor1.2_genomic.fna --in-bam $output --out-variants $vcf --gvcf
