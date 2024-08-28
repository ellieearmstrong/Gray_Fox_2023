#!/bin/sh
#SBATCH --job-name=b_gw_filterGVCF
#SBATCH --output=/scratch1/jazlynmo/grayfox/makevcfs_gw_main.drop295.bcf.out  #you will need to modify the path
#SBATCH --error=/scratch1/jazlynmo/grayfox/makevcfs_gw_main.drop295.bcf.err  #you will need to modify the path
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

#Filter out repeat regions and low mappability regions and down get all sites file is biallelic sites that pass filters but may or may not be snps
###Note with this new gvcf individual 295 was never used during joint calling so we do not need to remove
SECONDS=0
bcftools view --max-alleles 2 --exclude-types indels /scratch/users/elliea/jazlyn-ellie/mainlandgf_gwref_merge_allsites_Feb12.g.vcf.gz | bcftools view -T ^/scratch/users/elliea/jazlyn-ellie/graywolf_mappability_genmap.1.0.bed | bcftools annotate --rename-chrs /scratch/users/elliea/jazlyn-ellie/graywolf_renameChroms_number.txt -Oz -o  /scratch/users/elliea/jazlyn-ellie/graywolf_filtered.renameChroms.Mainland.gvcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf  /scratch/users/elliea/jazlyn-ellie/graywolf_filtered.renameChroms.Mainland.gvcf.gz
echo $SECONDS

##Filtering based on depth and allele count
SECONDS=0
bcftools filter -i 'AN > 59 && INFO/DP > 204 && INFO/DP < 501'  /scratch/users/elliea/jazlyn-ellie/graywolf_filtered.renameChroms.Mainland.gvcf.gz -Oz -o  /scratch/users/elliea/jazlyn-ellie/graywolf_filtered.renameChroms.Mainland.drop295.ACgr59_DPgr205lt500.gvcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf  /scratch/users/elliea/jazlyn-ellie/graywolf_filtered.renameChroms.Mainland.drop295.ACgr59_DPgr205lt500.gvcf.gz
echo $SECONDS

#####Turn the files into snp only vcfs
SECONDS=0
bcftools view -m2 -M2 -v snps  /scratch/users/elliea/jazlyn-ellie/graywolf_filtered.renameChroms.Mainland.drop295.ACgr59_DPgr205lt500.gvcf.gz -Oz -o  /scratch/users/elliea/jazlyn-ellie/graywolf_filtered.renameChroms.Mainland.drop295.ACgr59_DPgr205lt500.vcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf  /scratch/users/elliea/jazlyn-ellie/graywolf_filtered.renameChroms.Mainland.drop295.ACgr59_DPgr205lt500.vcf.gz
echo $SECONDS
