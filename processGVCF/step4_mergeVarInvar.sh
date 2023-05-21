#!/bin/sh
#SBATCH --job-name=mergeGVCF
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/MergeVarInvar.out  #you will need to modify the path
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/MergeVarInvar.err  #you will need to modify the path
#SBATCH --time=20:00:00 #twenty hour run time
#SBATCH -p normal #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --mem-per-cpu=5000MB #this is equivalent to 5G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jaam@stanford.edu

#load modules that you need

ml biology bcftools/1.16
ml biology samtools/1.16.1

for f in {1..32}
do


bcftools concat /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/splitChroms/chrom"$f"_fox_merge_variants_bisnps_autosome.GM.eh2.AN.qualdp.renamedChrom.vcf.gz /scratch/users/elliea/jazlyn-ellie/grayfox_2023/allsites/splitChroms/filter/chrom"$f"_allFoxes.rg.md.haplotypecaller.all.g.renameChroms.mappabilityFilter.AN.QUAL.DP.biallelic.gvcf.gz | bcftools sort - -Oz -o chrom"$f"_allFoxes.rg.md.haplotypecaller.all.g.renameChroms.mappabilityFilter.AN.QUAL.DP.biallelic.varInvar.sorted.vcf.gz

tabix -p vcf chrom"$f"_allFoxes.rg.md.haplotypecaller.all.g.renameChroms.mappabilityFilter.AN.QUAL.DP.biallelic.varInvar.sorted.vcf.gz

echo "chrom $f"



done
