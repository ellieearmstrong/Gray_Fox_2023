#!/bin/sh
#SBATCH --job-name=filterGVCF_c3
#SBATCH --output=/scratch1/jazlynmo/grayfox/makevcfs_c3.out  #you will need to modify the path
#SBATCH --error=/scratch1/jazlynmo/grayfox/makevcfs_c3.err  #you will need to modify the path
#SBATCH --time=47:00:00 #ten hour run time
#SBATCH -p qcb #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --mem-per-cpu=5000MB #this is equivalent to 10G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlynmo@usc.edu

#load modules that you need

module load gcc/11.3.0  openblas/0.3.21 bcftools/1.14
module load htslib/1.17


####Each filtering step takes about 13 hours
#Filter out repeat regions and low mappability regions and down get all sites file is biallelic sites that pass filters but may or may not be snps
#Filter canfam3.1
#SECONDS=0
#bcftools view -M 2 /scratch1/jazlynmo/grayfox/raw_data/canfam3_hc_gvcfs/canfam3.merged.allsites.vcf.gz | bcftools view -T ^/scratch1/jazlynmo/grayfox/genmap_repmask/canfam3_mappability_genmap.1.0.bed -T ^/scratch1/jazlynmo/grayfox/genmap_repmask/Canfam3.1_repeatmask.bed | bcftools annotate --rename-chrs /scratch1/jazlynmo/grayfox/metaData/canfam3.1_renameChroms.txt -Oz -o /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/Canfam3.1_filtered.renameChroms.gvcf.gz
#echo $SECONDS

#SECONDS=0
#tabix -p vcf /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/Canfam3.1_filtered.renameChroms.gvcf.gz
#echo $SECONDS

##Filtering based on depth and allele count
#Canfam3.1
SECONDS=0
bcftools filter -i 'AC > 23 && INFO/DP > 164 && INFO/DP < 501' /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/unfiltered/Canfam3.1_filtered.renameChroms.gvcf.gz | bcftools view --max-alleles 2 --exclude-types indels -Oz -o /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/Canfam3.1_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/Canfam3.1_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz
echo $SECONDS


#####Turn the files into snp only vcfs
SECONDS=0
bcftools view -m2 -M2 -v snps /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/Canfam3.1_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz -Oz -o /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/Canfam3.1_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/Canfam3.1_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz
echo $SECONDS
