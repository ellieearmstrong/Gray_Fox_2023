#!/bin/bash
#SBATCH --job-name=genmap_cf3
#SBATCH --time=12:00:00
#SBATCH --mem=240000
#SBATCH --partition=normal,hns,dpetrov
#SBATCH -n 24
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ellieearmstrong@gmail.com

<<<<<<< HEAD
/oak/stanford/groups/dpetrov/armstrong/applications/genmap index -F GCF_000002285.3_CanFam3.1_genomic.fna -I canfam3_index
/oak/stanford/groups/dpetrov/armstrong/applications/genmap index -F GCA_011100685.1_UU_Cfam_GSD_1.0_genomic.fna -I canfam4_index
/oak/stanford/groups/dpetrov/armstrong/applications/genmap index -F GCA_018345385.1_ASM1834538v1_genomic.fna -I arcticfox_index
/oak/stanford/groups/dpetrov/armstrong/applications/genmap index -F greyfox-dovetail-final.fasta -I grayfox_index

/oak/stanford/groups/dpetrov/armstrong/applications/genmap map -K 30 -E 2 -I canfam3_index/ -O canfam3_mappability -t -b
/oak/stanford/groups/dpetrov/armstrong/applications/genmap map -K 30 -E 2 -I canfam4_index/ -O canfam4_mappability -t -b
/oak/stanford/groups/dpetrov/armstrong/applications/genmap map -K 30 -E 2 -I arcticfox_index/ -O arcticfox_mappability -t -b
/oak/stanford/groups/dpetrov/armstrong/applications/genmap map -K 30 -E 2 -I grayfox_index/ -O grayfox_mappability -t -b

>>>>>>> 5cc5dd87536c0d3f9bdea7cd35f851666162dcec

module load system ruby
#process mappability scores
#ruby /oak/stanford/groups/dpetrov/armstrong/applications/RatesTools/bin/filterGM.rb grayfox_mappability.bed 1.0 exclude > grayfox_mappability_genmap.1.0.bed

#ruby /oak/stanford/groups/dpetrov/armstrong/applications/RatesTools/bin/simplify_bed.rb canfam3_mappability_genmap.1.0.bed canfam3_repeats.bed > canfam3.repeats.mappability1.bed

ruby /oak/stanford/groups/dpetrov/armstrong/applications/RatesTools/bin/simplify_bed.rb grayfox_mappability_genmap.1.0.bed Grayfox_repeatmask.bed> grayfox.repeats.mappability1.bed
