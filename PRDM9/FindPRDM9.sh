#!/bin/sh
#SBATCH --job-name=prdm9
#SBATCH --output=/home1/jazlynmo/sc1/logfiles/prdm9.out
#SBATCH --error=/home1/jazlynmo/sc1/logfiles/prdm9.err
#SBATCH --time=03:00:00
#SBATCH -p qcb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12000MB
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlynmo@usc.edu

#load modules
module purge
module load gcc/11.3.0 blast-plus/2.12.0 bedtools2/2.30.0

#create the blast database
#makeblastdb -in ../greyfox-dovetail-final.fasta -dbtype 'nucl' -hash_index -out greyfox-dovetail-final

#finding matching sequence
#blastn -query vulpes.fa -task blastn -db greyfox-dovetail-final -num_threads 4 -outfmt '6' > matching_greyfox_vulpes_int.fa
#blastn -query AW15.fa -task blastn -db greyfox-dovetail-final -num_threads 4 -outfmt '6' > matching_greyfox_wolf_int.fa

#pull the range from the reference
#bedtools getfasta -fi ../greyfox-dovetail-final.fasta -bed coords_PRDM9.bed >
greyfox-dovetail-final_PRMD9Scaffold_6__1_contigs__length_90720172_coords_65480835_65554467.fa
