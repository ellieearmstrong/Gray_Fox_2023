#!/bin/sh
#SBATCH --job-name=filterGVCF
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/FilterMerge.out  #you will need to modify the path
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/FilterMerge.err  #you will need to modify the path
#SBATCH --time=40:00:00 #ten hour run time
#SBATCH -p normal #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --mem-per-cpu=10000MB #this is equivalent to 10G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jaam@ustanford.edu

#load modules that you need

ml biology bcftools/1.16
ml biology samtools/1.16.1

for f in {1..32}
do

#apply mappability filter
#bcftools view -T ^/scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/autosomes_mappability1.0_renameToChroms.bed chrom"$f"_allFoxes.rg.md.haplotypecaller.all.g.renameChroms.gvcf.gz -Oz -o chrom"$f"_allFoxes.rg.md.haplotypecaller.all.g.renameChroms.mappabilityFilter.gvcf.gz

#tabix -p vcf chrom"$f"_allFoxes.rg.md.haplotypecaller.all.g.renameChroms.mappabilityFilter.gvcf.gz

#apply quality filters based on missingness and depth
bcftools filter -i 'AN > 28' chrom"$f"_allFoxes.rg.md.haplotypecaller.all.g.renameChroms.mappabilityFilter.gvcf.gz | bcftools view -e 'QUAL < 20 || INFO/DP > 500 || INFO/DP < 80' -Oz -o chrom"$f"_allFoxes.rg.md.haplotypecaller.all.g.renameChroms.mappabilityFilter.AN.QUAL.DP.gvcf.gz

tabix -p vcf chrom"$f"_allFoxes.rg.md.haplotypecaller.all.g.renameChroms.mappabilityFilter.AN.QUAL.DP.gvcf.gz

echo "chrom $f"



done

