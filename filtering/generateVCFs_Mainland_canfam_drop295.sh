#!/bin/sh
#SBATCH --job-name=bfilterGVCF
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/makevcfs_main.drop295.bcf.out  #you will need to modify the path 
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/makevcfs_main.drop295.bcf.err  #you will need to modify the path
#SBATCH --time=48:00:00 #two day run time
#SBATCH -p normal,hns #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --mem-per-cpu=5000MB #this is equivalent to 5G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlyn.mooney@gmail.com 

#load modules that you need

ml biology bcftools/1.16
ml biology samtools/1.16.1

####Each filtering step takes about 13 hours
#Filter out repeat regions and low mappability regions and down get all sites file is biallelic sites that pass filters but may or may not be snps
#Filter canfam3.1
#SECONDS=0
#bcftools view -M 2 /scratch/users/elliea/jazlyn-ellie/grayfox_2023/raw_data/batch3/sacks_cfam3_hc_gvcfs/batch3-mainland-canfam3.merged.allsites.vcf.gz| bcftools view -T ^/scratch/users/elliea/jazlyn-ellie/grayfox_2023/genmap_repmask/canfam3_mappability_genmap.1.0.bed -T ^/scratch/users/elliea/jazlyn-ellie/grayfox_2023/genmap_repmask/Canfam3.1_repeatmask.bed | bcftools annotate --rename-chrs /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/canfam3.1_renameChroms.txt -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.Mainland.gvcf.gz
#echo $SECONDS

#SECONDS=0
#tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.Mainland.gvcf.gz
#echo $SECONDS

#Filter canfam4
#SECONDS=0 
#bcftools view -M 2 /scratch/users/elliea/jazlyn-ellie/grayfox_2023/raw_data/canfam4_map/batch3/sacks_cfam4_hc/batch3-mainland-canfam4.merged.allsites.vcf.gz | bcftools view -T ^/scratch/users/elliea/jazlyn-ellie/grayfox_2023/genmap_repmask/canfam4_mappability_genmap.1.0.bed -T ^/scratch/users/elliea/jazlyn-ellie/grayfox_2023/genmap_repmask/Canfam4_repeatmask.bed | bcftools annotate --rename-chrs /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/canfam4_renameChroms_number.txt -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.Mainland.gvcf.gz
#echo $SECONDS

#SECONDS=0 
#tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.Mainland.gvcf.gz
#echo $SECONDS

#Remove the problematic individual and recode the info column
#delete SRR24465295 from individuals text then retry other filtering
#SECONDS=0 
#bcftools view /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.Mainland.gvcf.gz -s ^SRR24465295 | bcftools plugin fill-tags -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.Mainland.drop295.gvcf.gz
#echo $SECONDS

#SECONDS=0
#tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.Mainland.drop295.gvcf.gz
#echo $SECONDS

##Filtering based on depth and allele count
#SECONDS=0
#bcftools filter -e 'AC > 59 && INFO/DP > 204 && INFO/DP < 501' /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.Mainland.drop295.gvcf.gz -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.Mainland.drop295.ACgr25_DPgr165lt500.gvcf.gz
#echo $SECONDS

#SECONDS=0
#tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.Mainland.drop295.ACgr25_DPgr165lt500.gvcf.gz
#echo $SECONDS


##Turn the files into snp only vcfs
#SECONDS=0
#bcftools view -m2 -M2 -v snps /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.Mainland.drop295.ACgr25_DPgr165lt500.gvcf.gz -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.Mainland.drop295.ACgr25_DPgr165lt500.vcf.gz 
#echo $SECONDS

#SECONDS=0
#tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.Mainland.drop295.ACgr25_DPgr165lt500.vcf.gz
#echo $SECONDS

###Repeat with Canfam4
#Remove the problematic individual and recode the info column
#delete SRR24465295 from individuals text then retry other filtering
#SECONDS=0 
#bcftools view /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.Mainland.gvcf.gz -s ^SRR24465295 |  bcftools plugin fill-tags -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.Mainland.drop295.gvcf.gz
#echo $SECONDS

SECONDS=0
tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.Mainland.drop295.gvcf.gz
echo $SECONDS

##Filtering based on depth and allele count
SECONDS=0
bcftools filter -e 'AC < 61 || INFO/DP < 205 || INFO/DP > 500' /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.Mainland.drop295.gvcf.gz -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.Mainland.drop295.ACgr25_DPgr165lt500.gvcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.Mainland.drop295.ACgr25_DPgr165lt500.gvcf.gz
echo $SECONDS

#####Turn the files into snp only vcfs
SECONDS=0
bcftools view -m2 -M2 -v snps /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.Mainland.drop295.ACgr25_DPgr165lt500.gvcf.gz -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.Mainland.drop295.ACgr25_DPgr165lt500.vcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.Mainland.drop295.ACgr25_DPgr165lt500.vcf.gz
echo $SECONDS


