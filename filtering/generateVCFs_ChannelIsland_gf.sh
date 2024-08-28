#!/bin/sh
#SBATCH --job-name=gf_filterGVCF
#SBATCH --output=/scratch1/jazlynmo/grayfox/makevcf_gf.out  #you will need to modify the path
#SBATCH --error=/scratch1/jazlynmo/grayfox/makevcf_gf.err  #you will need to modify the path
#SBATCH --time=47:00:00 #run time
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
#SECONDS=0
#bcftools view -M 2 /scratch1/jazlynmo/grayfox/raw_data/gray_fox/grayfox_hc_vcf/grayfox.merged.allsites.vcf.gz | bcftools view -T ^/scratch1/jazlynmo/grayfox/genmap_repmask/grayfox_mappability_genmap.1.0.bed -T ^/scratch1/jazlynmo/grayfox/genmap_repmask/Grayfox_repeatmask.bed | bcftools annotate --rename-chrs /scratch1/jazlynmo/grayfox/metaData/grayfox_renameChroms_number.txt -Oz -o /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/grayfox_filtered.renameChroms.gvcf.gz
#echo $SECONDS

#SECONDS=0
#tabix -p vcf /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/grayfox_filtered.renameChroms.gvcf.gz
#echo $SECONDS

#grayfox
SECONDS=0
bcftools filter -e 'AN > 23 && INFO/DP > 204 && INFO/DP < 501' /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/unfiltered/grayfox_filtered.renameChroms.gvcf.gz |  bcftools reheader -s /scratch1/jazlynmo/grayfox/metaData/rename_grayfox | bcftools view --max-alleles 2 --exclude-types indels -Oz -o /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/grayfox_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/grayfox_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz
echo $SECONDS


#####Turn the files into snp only vcfs

SECONDS=0
bcftools view -m2 -M2 -v snps /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/grayfox_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz -Oz -o /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/grayfox_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/grayfox_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz
echo $SECONDS
