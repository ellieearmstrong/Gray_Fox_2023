#!/bin/sh
#SBATCH --job-name=calls
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/calls.out  #you will need to modify the path
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/calls.err  #you will need to modify the path
#SBATCH --time=10:00:00 #ten hour run time
#SBATCH -p normal #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --mem-per-cpu=500MB #this is equivalent to half G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jaam@stanford.edu


for f in {1..32}
do

zgrep -v "#"  chrom"$f"_allFoxes.rg.md.haplotypecaller.all.g.renameChroms.mappabilityFilter.AN.QUAL.DP.biallelic.varInvar.sorted.vcf.gz | awk '{print $1"\t"$2-1"\t"$2}' > callableSites/chrom"$f"_callableSites.bed

echo "chrom $f"

done
