#!/bin/sh
#SBATCH --job-name=grabGT_w
#SBATCH --output=/project/jazlynmo_738/Jazlyn/grayfox/logfiles/grabGT_w.out  #you will need to modify the path
#SBATCH --error=/project/jazlynmo_738/Jazlyn/grayfox/logfiles/grabGT_w.err  #you will need to modify the path
#SBATCH --time=6:00:00 #run time
#SBATCH -p qcb #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --mem-per-cpu=2000MB #this is equivalent to 2G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlynmo@usc.edu

#load modules that you need
module purge
module load gcc/11.3.0 intel/19.0.4 openblas/0.3.21 bcftools/1.14 htslib/1.17 vcftools/0.1.14 python/3.9.12

####Mainland
while read -r infile
do

	vcftools --gzvcf /project/jazlynmo_738/DataRepository/Canids/Invariant/GrayFox/Mainland/"$infile".gvcf.gz --keep /project/jazlynmo_738/Jazlyn/grayfox/metaData/mainland_West_n12.txt --counts --out "$infile"_West_n12

echo $infile

done < gvcfIn_Mainland.txt 
