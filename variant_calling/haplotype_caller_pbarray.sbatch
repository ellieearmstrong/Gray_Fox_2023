#!/bin/bash
#SBATCH --job-name=af3-array
#SBATCH --partition=gpu,owners
#SBATCH --time=24:00:00
#SBATCH -C GPU_MEM:16GB
#SBATCH -N 1
#SBATCH -G 1
#SBATCH --array=1-42
#SBATCH --cpus-per-task=10
#SBATCH --mem=120000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ellieearmstrong@gmail.com

config=af_b3_config2.txt
input=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
output=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
singularity exec --nv /oak/stanford/groups/dpetrov/armstrong/applications/parabricks-4.1.1-1.simg pbrun haplotypecaller --ref ../GCA_018345385.1_ASM1834538v1_genomic.fna --in-bam $input --out-variants $output --gvcf
