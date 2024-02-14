#! /bin/bash
#SBATCH --job-name=grayfox-sift
#SBATCH --time=6-00:00:00
#SBATCH --partition=dpetrov
#SBATCH --mem=240000
#SBATCH -n 24
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ellieearmstrong@gmail.com

source /oak/stanford/groups/dpetrov/armstrong/applications/miniconda3/etc/profile.d/conda.sh
#conda activate sift-test

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta.gz
gunzip uniref90.fasta.gz

#module load system gcc/12.1.0c
conda deactivate
module load perl/5.26.0
module load gcc/12.1.0
perl make-SIFT-db-all.pl -config /oak/stanford/groups/dpetrov/armstrong/gray-foxes/Gray_Fox_2023/SIFT_database/grayfox-sift4G.txt
