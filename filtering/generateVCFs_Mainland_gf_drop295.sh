#!/bin/sh
#SBATCH --job-name=gf_filterGVCF
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/makevcfs_gf.out  #you will need to modify the path
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/makevcfs_gf.err  #you will need to modify the path
#SBATCH --time=12:00:00 #run time
#SBATCH -p hns,normal #the main partition
#SBATCH -n 24
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --mem 240000 #this is equivalent to 5G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlynmo@usc.edu

#load modules that you need
ml biology bcftools/1.16 htslib/1.16 vcftools/0.1.15

#Filter out repeat regions and low mappability regions and down get all sites file is biallelic sites that pass filters but may or may not be snps

#SECONDS=0
#bcftools view --max-alleles 2 --exclude-types indels /scratch/users/elliea/jazlyn-ellie/mainlandgf_gfref_merge_allsites_Jan24.g.vcf.gz | bcftools view -T ^/scratch/users/elliea/jazlyn-ellie/grayfox_mappability_genmap.1.0.bed | bcftools annotate --rename-chrs /scratch/users/elliea/jazlyn-ellie/grayfox_renameChroms_number.txt -Oz -o  /scratch/users/elliea/jazlyn-ellie/grayfox_filtered.renameChroms.Mainland.gvcf.gz
#echo $SECONDS

#SECONDS=0
#tabix -p vcf  /scratch/users/elliea/jazlyn-ellie/grayfox_filtered.renameChroms.Mainland.gvcf.gz
#echo $SECONDS

##Filtering based on depth and allele count
#SECONDS=0
#bcftools filter -i 'INFO/DP > 204 && INFO/DP < 501 && QUAL>=30'  /scratch/users/elliea/jazlyn-ellie/grayfox_filtered.renameChroms.Mainland.gvcf.gz -Oz -o  /scratch/users/elliea/jazlyn-ellie/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500.gvcf.gz
#echo $SECONDS

#SECONDS=0
#tabix -p vcf  /scratch/users/elliea/jazlyn-ellie/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500.gvcf.gz
#echo $SECONDS

#####Turn the files into snp only vcfs
#SECONDS=0
#bcftools view -m2 -M2 -v snps  /scratch/users/elliea/jazlyn-ellie/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500.gvcf.gz -Oz -o  /scratch/users/elliea/jazlyn-ellie/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500.vcf.gz
#echo $SECONDS

#SECONDS=0
#tabix -p vcf  /scratch/users/elliea/jazlyn-ellie/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500.vcf.gz
#echo $SECONDS

#Filter based on AN
#vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500.gvcf.gz --missing-site --out /scratch/users/elliea/jazlyn-ellie/grayfox_gvcf
#vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500.vcf.gz --missing-site --out /scratch/users/elliea/jazlyn-ellie/grayfox_vcf

#pull out sites with more than 22 missing GTs and remove them from files
#awk 'NR > 1 && $5 > 22 {print $1"\t"$2-1"\t"$2}' /scratch/users/elliea/jazlyn-ellie/grayfox_gvcf.lmiss > /scratch/users/elliea/jazlyn-ellie/grayfox_gvcf_ANlt60.bed
#awk 'NR > 1 && $5 > 22 {print $1"\t"$2-1"\t"$2}' /scratch/users/elliea/jazlyn-ellie/sum_stats/grayfox_vcf.lmiss > /scratch/users/elliea/jazlyn-ellie/grayfox_vcf_ANlt60.bed

#Filter VCF to exclude sites in the BED file
#bcftools view -T ^/scratch/users/elliea/jazlyn-ellie/grayfox_vcf_ANlt60.bed /scratch/users/elliea/jazlyn-ellie/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500.vcf.gz -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500_ANgr59.vcf.gz
#tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500_ANgr59.vcf.gz

#Filter gvcf to exclude sites in the BED file
bcftools view -T ^/scratch/users/elliea/jazlyn-ellie/grayfox_gvcf_ANlt60.bed /scratch/users/elliea/jazlyn-ellie/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500.gvcf.gz -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500_ANgr59.gvcf.gz
tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500_ANgr59.gvcf.gz
