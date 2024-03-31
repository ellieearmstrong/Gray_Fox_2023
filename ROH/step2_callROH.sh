#!/bin/sh
#SBATCH --job-name=roh
#SBATCH --output=/scratch1/jazlynmo/grayfox/logfiles/roh.out  #you will need to modify the path
#SBATCH --error=/scratch1/jazlynmo/grayfox/logfiles/roh.err  #you will need to modify the path
#SBATCH --time=1:00:00 #one hour run time
#SBATCH -p qcb #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30000MB #this is equivalent to 30G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlynmo@usc.edu

for ref in {Mainland,ChannelIsland}

do

	while read -r p

	do

	/scratch1/jazlynmo/grayfox/garlic/bin/linux/garlic --tped /scratch1/jazlynmo/grayfox/plink/"$ref"/"$p".tped --tfam /scratch1/jazlynmo/grayfox/plink/"$ref"/"$p".tfam --centromere fakeCentromeres_dog.txt --error 0.001 --gl-type GQ --winsize 100 --auto-overlap-frac --size-bounds 100000 1000000 2000000 3000000 4000000 5000000 --out "$ref"/"$p"_fixWindow100bp

	done < vcfIn_"$ref".txt

done

sleep 180
