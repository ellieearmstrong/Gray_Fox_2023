#Filter out repeat regions and low mappability regions and down get all sites file is biallelic sites that pass filters but may or may not be snps

SECONDS=0
bcftools view --max-alleles 2 --exclude-types indels /scratch/users/elliea/jazlyn-ellie/mainlandgf_gwref_merge_allsites_Feb12.g.vcf.gz | bcftools view -T ^/scratch/users/elliea/jazlyn-ellie/graywolf_mappability_genmap.1.0.bed | bcftools annotate --rename-chrs /scratch/users/elliea/jazlyn-ellie/graywolf_renameChroms_number.txt -Oz -o  /scratch/users/elliea/jazlyn-ellie/graywolf_filtered.renameChroms.Mainland.gvcf.gz
echo $SECONDS

SECONDS=0
tabix -p vcf  /scratch/users/elliea/jazlyn-ellie/graywolf_filtered.renameChroms.Mainland.gvcf.gz
echo $SECONDS

##Filtering based on depth and allele count
SECONDS=0
bcftools filter -i 'AC > 59 && INFO/DP > 204 && INFO/DP < 501'  /scratch/users/elliea/jazlyn-ellie/graywolf_filtered.renameChroms.Mainland.gvcf.gz -Oz -o  /scratch/users/elliea/jazlyn-ellie/graywolf_filtered.renameChroms.Mainland.drop295.ACgr59_DPgr205lt500.gvcf.gz
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
