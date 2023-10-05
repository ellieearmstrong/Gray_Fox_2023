#!/bin/bash
#SBATCH --job-name=genmap_fox
#SBATCH --time=4:00:00
#SBATCH --mem=120000
#SBATCH --partition=dpetrov
#SBATCH -n 12
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ellieearmstrong@gmail.com

/oak/stanford/groups/dpetrov/armstrong/applications/genmap index -F GCF_000002285.3_CanFam3.1_genomic.fna -I canfam3_index

/oak/stanford/groups/dpetrov/armstrong/applications/genmap map -K 30 -E 2 -I canfam3_index/ -O canfam3_mappability -t -b

module load ruby
#process mappability scores
ruby /oak/stanford/groups/dpetrov/armstrong/applications/RatesTools/bin/filterGM.rb grayfox_mappability.bed 1.0 exclude > grayfox_mappability_genmap.1.0.bed
