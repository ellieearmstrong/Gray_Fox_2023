#!/bin/sh
#SBATCH --job-name=mergeGVCF
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/mergeIndivsplitchroms.out  #you will need to modify the path
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/mergeIndivsplitchroms.err  #you will need to modify the path
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


bcftools merge chrom"$f"_GFO41F.rg.md.haplotypecaller.all.g.renameChroms.vcf.gz chrom"$f"_SCA16F.rg.md.haplotypecaller.all.g.renameChroms.vcf.gz chrom"$f"_SCLV4F.rg.md.haplotypecaller.all.g.renameChroms.vcf.gz chrom"$f"_SCZ05M.rg.md.haplotypecaller.all.g.renameChroms.vcf.gz chrom"$f"_SMI15F.rg.md.haplotypecaller.all.g.renameChroms.vcf.gz chrom"$f"_SNI05F.rg.md.haplotypecaller.all.g.renameChroms.vcf.gz chrom"$f"_SNI41F.rg.md.haplotypecaller.all.g.renameChroms.vcf.gz chrom"$f"_SRO40F.rg.md.haplotypecaller.all.g.renameChroms.vcf.gz chrom"$f"_SRR7458264.rg.md.haplotypecaller.all.g.renameChroms.vcf.gz chrom"$f"_SRR7458265.rg.md.haplotypecaller.all.g.renameChroms.vcf.gz chrom"$f"_SRR7458266.rg.md.haplotypecaller.all.g.renameChroms.vcf.gz chrom"$f"_SRR7458267.rg.md.haplotypecaller.all.g.renameChroms.vcf.gz chrom"$f"_SRR7458268.rg.md.haplotypecaller.all.g.renameChroms.vcf.gz chrom"$f"_SRR7458269.rg.md.haplotypecaller.all.g.renameChroms.vcf.gz chrom"$f"_SRR7458270.rg.md.haplotypecaller.all.g.renameChroms.vcf.gz chrom"$f"_SRR7458271.rg.md.haplotypecaller.all.g.renameChroms.vcf.gz -Oz -o chrom"$f"_allFoxes.rg.md.haplotypecaller.all.g.renameChroms.gvcf.gz



tabix -p vcf chrom"$f"_allFoxes.rg.md.haplotypecaller.all.g.renameChroms.gvcf.gz

echo "chrom $f" 

done 
