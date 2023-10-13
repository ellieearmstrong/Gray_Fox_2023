#!/bin/sh
#SBATCH --job-name=thetaAF
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/logfiles/theta_%a_af.out  #you will need to modify the path
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/logfiles/theta_%a_af.err  #you will need to modify the path
#SBATCH --time=05:00:00 #ten hour run time
#SBATCH -p normal #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --array=1-24
#SBATCH --mem-per-cpu=5000MB #this is equivalent to 10G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jaam@stanford.edu

#set variables
CHROM=$SLURM_ARRAY_TASK_ID

outputDir='/scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/theta'

#load python
ml python/3.9.0

#set the reference
ref=arcticfox

	outputFile=${outputDir}/'chrom'${CHROM}'_'${ref}'_filtered.renameChroms.ACgr25_DPgr165lt500.segSites'

        python3 CalcKForWattersonsTheta.py /scratch/users/elliea/jazlyn-ellie/grayfox_2023/allsites/splitChroms/chr"$CHROM"_"$ref"_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz ${outputFile}

echo "done $ref"


#sleep timer
sleep 180
