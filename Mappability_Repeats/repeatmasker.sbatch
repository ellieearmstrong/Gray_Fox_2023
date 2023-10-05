#!/bin/bash
#SBATCH --job-name=repmask-gf
#SBATCH --time=48:00:00
#SBATCH --mem=240000
#SBATCH --partition=normal,hns
##SBATCH -n 24
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=ellieearmstrong@gmail.com


#Step 1: Build Database
#singularity run /oak/stanford/groups/dpetrov/armstrong/applications/singularity-containers/dfam-tetools-latest.sif BuildDatabase -name greyfox_repmask greyfox-dovetail-final.fasta

#Step 2: Repeat Modeler
#singularity run /oak/stanford/groups/dpetrov/armstrong/applications/singularity-containers/dfam-tetools-latest.sif RepeatModeler -database greyfox_repmask -threads 32 

#Step 3: Repeat Masker
singularity run /oak/stanford/groups/dpetrov/armstrong/applications/singularity-containers/dfam-tetools-latest.sif RepeatMasker greyfox-dovetail-final.fasta -lib RM_28818.WedAug160217552023/consensi.fa.classified -pa 24 --xsmall
