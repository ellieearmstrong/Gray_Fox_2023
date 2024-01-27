#!/bin/bash
#SBATCH --job-name=genomicsdb
#SBATCH --partition=gpu,owners
#SBATCH --time=12:00:00
#SBATCH -C GPU_MEM:16GB
#SBATCH -N 1
#SBATCH -G 1
#SBATCH --array=2,38
#SBATCH --cpus-per-task=12
#SBATCH --mem=120000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ellieearmstrong@gmail.com

module load biology gatk/4.1.4.1
config=dbimport_config.txt
interval=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
scaffold=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
gatk GenomicsDBImport --java-options "-Xmx75g -Xms75g" --sample-name-map gf_samplemap.txt --genomicsdb-workspace-path grayfoxdb_$scaffold -L $interval --tmp-dir /scratch/users/elliea/
