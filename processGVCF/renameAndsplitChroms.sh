#!/bin/sh
#SBATCH --job-name=processGVCF
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/splitchroms.out  #you will need to modify the path 
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/splitchroms.err  #you will need to modify the path
####SBATCH --time=70:00:00 #ten hour run time
#SBATCH --time=30:00:00 #ten hour run time
#SBATCH -p normal #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --mem-per-cpu=15000MB #this is equivalent to 10G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jaam@ustanford.edu 

#load modules that you need

ml biology bcftools/1.16
ml biology samtools/1.16.1


while read -r infile 
do

#rename chromosomes
bcftools annotate --rename-chrs /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/fox_chrom.map "$infile".vcf -Oz -o "$infile".renameChroms.vcf.gz
tabix -p vcf "$infile".renameChroms.vcf.gz


#split into separate files
for i in {1..32}
do

bcftools view --regions "$i" "$infile".renameChroms.vcf.gz -Oz -o splitChroms/chrom"$i"_"$infile".renameChroms.vcf.gz 
tabix -p vcf splitChroms/chrom"$i"_"$infile".renameChroms.vcf.gz

echo "$i"
done

echo "$infile"

done < allIndivs.txt
