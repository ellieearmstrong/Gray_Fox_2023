#!/bin/sh
#SBATCH --job-name=pi32
#SBATCH --output=/scratch1/jazlynmo/grayfox/logfiles/pi_%a_32.out  #you will need to modify the path
#SBATCH --error=/scratch1/jazlynmo/grayfox/logfiles/pi_%a_32.err  #you will need to modify the path
#SBATCH --time=20:00:00 #run time
#SBATCH -p qcb #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --array=25-32
#SBATCH --mem-per-cpu=1000MB #this is equivalent to 1G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlynmo@usc.edu

#load modules
module load gcc/11.3.0 intel/19.0.4 vcftools/0.1.14


#set variables
CHROM=$SLURM_ARRAY_TASK_ID

#grab allele counts per pop and remove any sites with missing data
while read -r infile
do

	for i in {East,West,Hybrid}
	do
		vcftools --gzvcf /project/jazlynmo_738/DataRepository/Canids/Invariant/GrayFox/Mainland/"$infile".gvcf.gz --keep /project/jazlynmo_738/Jazlyn/grayfox/metaData/mainland_"$i"_n12.txt --max-missing 1 --site-pi --chr "chr"$CHROM --out /scratch1/jazlynmo/grayfox/"$infile"_chr"$CHROM"_"$i"_n12
		awk 'NR > 1 {sum += $3; count++} END {if (count > 0) print sum / count}' /scratch1/jazlynmo/grayfox/"$infile"_chr"$CHROM"_"$i"_n12.sites.pi > output_pi_1bp/"$infile"_chr"$CHROM"_"$i"_n12.avg.pi 
	done

	echo "done $infile"

done < gvcfIn_Mainland_chr32.txt
