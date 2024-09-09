#!/bin/sh
#SBATCH --job-name=sfs
#SBATCH --output=/scratch1/jazlynmo/grayfox/logfiles/sfs.out  #you will need to modify the path
#SBATCH --error=/scratch1/jazlynmo/grayfox/logfiles/sfs.err  #you will need to modify the path
#SBATCH --time=20:00:00 #run time
#SBATCH -p qcb #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --mem-per-cpu=2000MB #this is equivalent to 2G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlynmo@usc.edu

#load modules
module load gcc/11.3.0 intel/19.0.4 vcftools/0.1.14

#grab allele counts per pop and remove any sites with missing data
while read -r infile
do

	for i in {East,West,Hybrid}
	do 
		vcftools --gzvcf /project/jazlynmo_738/DataRepository/Canids/Variants/GrayFox/Mainland/"$infile".vcf.gz --keep /project/jazlynmo_738/Jazlyn/grayfox/metaData/mainland_"$i"_n12.txt --max-missing 1 --counts2 --out "$infile"_"$i"_n12
		awk 'NR > 1 {print $5}' "$infile"_"$i"_n12.frq.count | sort | uniq -c | sort -k2,2n | sed -e 1i'countRef\tbin' > "$infile"_"$i"_n12_refAlleleCounts.txt
		
	done 

	
	echo "done $infile"

done < vcfIn_Mainland.txt
