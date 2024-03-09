#!/bin/bash
#SBATCH --job-name=dbarray
#SBATCH --time=16:00:00
#SBATCH --array=1-38
#SBATCH --partition=hns,normal,dpetrov
#SBATCH -n 12
#SBATCH --mem=120000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ellieearmstrong@gmail.com

module load biology gatk/4.1.4.1
config=dbimport_config_wolf.txt
interval=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
scaffold=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
gatk GenomicsDBImport --java-options "-Xmx75g -Xms75g" --sample-name-map gf_samplemap_gwolf.txt  --genomicsdb-workspace-path graywolfdb_$scaffold -L $interval --tmp-dir /scratch/users/elliea/
