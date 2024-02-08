#!/bin/sh
#SBATCH --job-name=gf_af_filterGVCF
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/makevcfs_gf_af.out  #you will need to modify the path 
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/makevcfs_gf_af.err  #you will need to modify the path
#SBATCH --time=47:00:00 #run time
#SBATCH -p normal #the main partition
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
#arctic fox
SECONDS=0
bcftools view --max-alleles 2 --exclude-types indels /scratch/users/elliea/jazlyn-ellie/grayfox_2023/raw_data/arctic_fox/arcticfox_hc_vcf/arcticfox.merged.allsites.vcf.gz | bcftools view -T ^/scratch/users/elliea/jazlyn-ellie/grayfox_2023/genmap_repmask/arcticfox_mappability_genmap.1.0.bed -T ^/scratch/users/elliea/jazlyn-ellie/grayfox_2023/genmap_repmask/ArcticFox_repeatmask.bed | bcftools annotate --rename-chrs /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/arcticfox_renameChroms_number.txt -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/arcticfox_filtered.renameChroms.gvcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/arcticfox_filtered.renameChroms.gvcf.gz
echo $SECONDS

#gray fox
SECONDS=0
bcftools view --max-alleles 2 --exclude-types indels /scratch/users/elliea/jazlyn-ellie/grayfox_2023/raw_data/gray_fox/grayfox_hc_vcf/grayfox.merged.allsites.vcf.gz | bcftools view -T ^/scratch/users/elliea/jazlyn-ellie/grayfox_2023/genmap_repmask/grayfox_mappability_genmap.1.0.bed -T ^/scratch/users/elliea/jazlyn-ellie/grayfox_2023/genmap_repmask/Grayfox_repeatmask.bed | bcftools annotate --rename-chrs /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/grayfox_renameChroms_number.txt -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/grayfox_filtered.renameChroms.gvcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/grayfox_filtered.renameChroms.gvcf.gz
echo $SECONDS



##Filtering based on depth and allele count
#arcticfox
SECONDS=0
bcftools filter -e 'AC > 23 && INFO/DP > 204 && INFO/DP < 501' /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/arcticfox_filtered.renameChroms.gvcf.gz -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/arcticfox_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/arcticfox_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz
echo $SECONDS

#grayfox
SECONDS=0
bcftools filter -e 'AC > 23 && INFO/DP > 204 && INFO/DP < 501' /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/grayfox_filtered.renameChroms.gvcf.gz -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/grayfox_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/grayfox_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz
echo $SECONDS



#####Turn the files into snp only vcfs
SECONDS=0
bcftools view -m2 -M2 -v snps /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/arcticfox_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/arcticfox_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz 
echo $SECONDS

SECONDS=0
tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/arcticfox_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz
echo $SECONDS

SECONDS=0
bcftools view -m2 -M2 -v snps /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/grayfox_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz -Oz -o /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/grayfox_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/grayfox_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz
echo $SECONDS


