#!/bin/bash
#SBATCH --job-name=genmap_cf3
#SBATCH --time=12:00:00
#SBATCH --mem=240000
#SBATCH --partition=normal,hns,dpetrov
#SBATCH -n 24
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ellieearmstrong@gmail.com

#/oak/stanford/groups/dpetrov/armstrong/applications/genmap index -F GCF_000002285.3_CanFam3.1_genomic.fna -I canfam3_index

#/oak/stanford/groups/dpetrov/armstrong/applications/genmap map -K 30 -E 2 -I canfam3_index/ -O canfam3_mappability -t -b

module load system ruby
#process mappability scores
#ruby /oak/stanford/groups/dpetrov/armstrong/applications/RatesTools/bin/filterGM.rb grayfox_mappability.bed 1.0 exclude > grayfox_mappability_genmap.1.0.bed

#ruby /oak/stanford/groups/dpetrov/armstrong/applications/RatesTools/bin/simplify_bed.rb canfam3_mappability_genmap.1.0.bed canfam3_repeats.bed > canfam3.repeats.mappability1.bed

ruby /oak/stanford/groups/dpetrov/armstrong/applications/RatesTools/bin/simplify_bed.rb grayfox_mappability_genmap.1.0.bed Grayfox_repeatmask.bed> grayfox.repeats.mappability1.bed
