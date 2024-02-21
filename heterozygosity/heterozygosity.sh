#!/bin/sh
#SBATCH --job-name=het
#SBATCH --output=/scratch1/jazlynmo/grayfox/logfiles/het.out  #you will need to modify the path
#SBATCH --error=/scratch1/jazlynmo/grayfox/logfiles/het.err  #you will need to modify the path
#SBATCH --time=06:00:00 #run time
#SBATCH -p qcb #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --mem-per-cpu=2000MB #this is equivalent to 2G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlynmo@usc.edu

#load modules
module load gcc/11.3.0 intel/19.0.4 vcftools

#compute heterozygosity
for f in {Mainland,ChannelIsland}
do

	while read -r p
	do
		vcftools --gzvcf /scratch1/jazlynmo/grayfox/vcf/"$f"/$p.vcf.gz --het --out $p
	done < vcfIn_"$f".txt


done
