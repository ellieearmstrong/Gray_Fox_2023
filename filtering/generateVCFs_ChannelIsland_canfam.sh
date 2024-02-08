#!/bin/sh
#SBATCH --job-name=filterGVCF
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/makevcfs.out  #you will need to modify the path 
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/makevcfs.err  #you will need to modify the path
#SBATCH --time=47:00:00 #ten hour run time
#SBATCH -p normal #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --mem-per-cpu=5000MB #this is equivalent to 10G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlyn.mooney@gmail.com 

#load modules that you need

ml biology bcftools/1.16
ml biology samtools/1.16.1


####Each filtering step takes about 13 hours
#Filter out repeat regions and low mappability regions and down get all sites file is biallelic sites that pass filters but may or may not be snps
#Filter canfam3.1
#SECONDS=0
#bcftools view --max-alleles 2 --exclude-types indels /scratch/users/elliea/jazlyn-ellie/grayfox_2023/raw_data/canfam3_hc_gvcfs/canfam3.merged.allsites.vcf.gz | bcftools view -T ^/scratch/users/elliea/jazlyn-ellie/grayfox_2023/genmap_repmask/canfam3_mappability_genmap.1.0.bed -T ^/scratch/users/elliea/jazlyn-ellie/grayfox_2023/genmap_repmask/Canfam3.1_repeatmask.bed | bcftools annotate --rename-chrs /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/canfam3.1_renameChroms.txt -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.gvcf.gz
#echo $SECONDS

#SECONDS=0
#tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.gvcf.gz
#echo $SECONDS

#Filter canfam4
#SECONDS=0 
#bcftools view --max-alleles 2 --exclude-types indels /scratch/users/elliea/jazlyn-ellie/grayfox_2023/raw_data/canfam4_map/canfam4_hc_gvcfs/canfam4.merged.allsites.vcf.gz | bcftools view -T ^/scratch/users/elliea/jazlyn-ellie/grayfox_2023/genmap_repmask/canfam4_mappability_genmap.1.0.bed -T ^/scratch/users/elliea/jazlyn-ellie/grayfox_2023/genmap_repmask/Canfam4_repeatmask.bed | bcftools annotate --rename-chrs /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/canfam4_renameChroms_number.txt -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.gvcf.gz
#echo $SECONDS

#SECONDS=0 
#tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.gvcf.gz
#echo $SECONDS

##Filtering based on depth and allele count
#Canfam3.1
#SECONDS=0
#bcftools filter -e 'AC > 23 && INFO/DP > 204 && INFO/DP < 501' /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.gvcf.gz -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz
#echo $SECONDS

#SECONDS=0
#tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz
#echo $SECONDS

#Canfam4
#SECONDS=0
#bcftools filter -e 'AC > 23 && INFO/DP > 204 && INFO/DP < 501' /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.gvcf.gz -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz
#echo $SECONDS

#SECONDS=0
#tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz
#echo $SECONDS



#####Turn the files into snp only vcfs
#SECONDS=0
#bcftools view -m2 -M2 -v snps /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz 
#echo $SECONDS

#SECONDS=0
#tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz
#echo $SECONDS

#SECONDS=0
#bcftools view -m2 -M2 -v snps /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz
#echo $SECONDS

#SECONDS=0
#tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz
#echo $SECONDS

