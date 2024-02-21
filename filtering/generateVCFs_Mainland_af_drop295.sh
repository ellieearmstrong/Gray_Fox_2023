#!/bin/sh
#SBATCH --job-name=b_af_filterGVCF
#SBATCH --output=/scratch1/jazlynmo/grayfox/makevcfs_af_main.drop295.bcf.out  #you will need to modify the path
#SBATCH --error=/scratch1/jazlynmo/grayfox/makevcfs_af_main.drop295.bcf.err  #you will need to modify the path
#SBATCH --time=48:00:00 #run time
#SBATCH -p qcb #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --mem-per-cpu=5000MB #this is equivalent to 5G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlynmo@usc.edu

#load modules that you need

module load gcc/11.3.0  openblas/0.3.21 bcftools/1.14
module load htslib/1.17

####Each filtering step takes about 13 hours
#Filter out repeat regions and low mappability regions and down get all sites file is biallelic sites that pass filters but may or may not be snps
#arctic fox
#SECONDS=0
#bcftools view -M 2 /scratch1/jazlynmo/grayfox/raw_data/arctic_fox/batch3/sacks_hc_arctic/batch3-mainland-arcticfox.merged.allsites.vcf.gz | bcftools view -T ^/scratch1/jazlynmo/grayfox/genmap_repmask/arcticfox_mappability_genmap.1.0.bed -T ^/scratch1/jazlynmo/grayfox/genmap_repmask/ArcticFox_repeatmask.bed | bcftools annotate --rename-chrs /scratch1/jazlynmo/grayfox/metaData/arcticfox_renameChroms_number.txt -Oz -o /scratch1/jazlynmo/grayfox/vcf/Mainland/arcticfox_filtered.renameChroms.Mainland.gvcf.gz
#echo $SECONDS

#SECONDS=0
#tabix -p vcf /scratch1/jazlynmo/grayfox/vcf/Mainland/arcticfox_filtered.renameChroms.Mainland.gvcf.gz
#echo $SECONDS

###Repeat filtering with arctic fox
#Remove the problematic individual and recode the info column
#delete SRR24465295 from individuals text then retry other filtering
#SECONDS=0
#bcftools view /scratch1/jazlynmo/grayfox/vcf/Mainland/arcticfox_filtered.renameChroms.Mainland.gvcf.gz -s ^SRR24465295 |  bcftools plugin fill-tags -Oz -o /scratch1/jazlynmo/grayfox/vcf/Mainland/arcticfox_filtered.renameChroms.Mainland.drop295.gvcf.gz
#echo $SECONDS

#SECONDS=0
#tabix -p vcf /scratch1/jazlynmo/grayfox/vcf/Mainland/arcticfox_filtered.renameChroms.Mainland.drop295.gvcf.gz
#echo $SECONDS

##Filtering based on depth and allele count
SECONDS=0
bcftools filter -i 'AC > 59 && INFO/DP > 204 && INFO/DP < 501' /scratch1/jazlynmo/grayfox/vcf/Mainland/unfiltered/arcticfox_filtered.renameChroms.Mainland.drop295.gvcf.gz | bcftools view --max-alleles 2 --exclude-types indels -Oz -o /scratch1/jazlynmo/grayfox/vcf/Mainland/arcticfox_filtered.renameChroms.Mainland.drop295.ACgr59_DPgr205lt500.gvcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf /scratch1/jazlynmo/grayfox/vcf/Mainland/arcticfox_filtered.renameChroms.Mainland.drop295.ACgr59_DPgr205lt500.gvcf.gz
echo $SECONDS

#####Turn the files into snp only vcfs
SECONDS=0
bcftools view -m2 -M2 -v snps /scratch1/jazlynmo/grayfox/vcf/Mainland/arcticfox_filtered.renameChroms.Mainland.drop295.ACgr59_DPgr205lt500.gvcf.gz -Oz -o /scratch1/jazlynmo/grayfox/vcf/Mainland/arcticfox_filtered.renameChroms.Mainland.drop295.ACgr59_DPgr205lt500.vcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf /scratch1/jazlynmo/grayfox/vcf/Mainland/arcticfox_filtered.renameChroms.Mainland.drop295.ACgr59_DPgr205lt500.vcf.gz
echo $SECONDS
