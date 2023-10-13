#!/bin/sh
#SBATCH --job-name=theta
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/logfiles/theta_%a_C34.out  #you will need to modify the path
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/logfiles/theta_%a_C34.err  #you will need to modify the path
#SBATCH --time=08:00:00 #ten hour run time
#SBATCH -p normal #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --array=1-38
#SBATCH --mem-per-cpu=8000MB #this is equivalent to 10G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jaam@stanford.edu

#set variables
CHROM=$SLURM_ARRAY_TASK_ID
outputDir='/scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/theta'
metaData='/scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData'

#load python
ml python/3.9.0

#loop through populations
for ref in {Canfam3.1,Canfam4} 

do 

	
	outputFile=${outputDir}/'chrom'${CHROM}'_'${ref}'_filtered.renameChroms.ACgr25_DPgr165lt500.segSites'

	python3 CalcKForWattersonsTheta.py /scratch/users/elliea/jazlyn-ellie/grayfox_2023/allsites/splitChroms/chr"$CHROM"_"$ref"_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz ${outputFile}
	
	
echo "done $ref"

done 

#sleep timer
sleep 180
