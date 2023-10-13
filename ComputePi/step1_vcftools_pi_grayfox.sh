#!/bin/sh
#SBATCH --job-name=piGF
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/logfiles/pi_%a_gf.out  #you will need to modify the path
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/logfiles/pi_%a_gf.err  #you will need to modify the path
#SBATCH --time=05:00:00 #ten hour run time
#SBATCH -p normal #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --array=1-32
#SBATCH --mem-per-cpu=500MB #this is equivalent to 10G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jaam@stanford.edu

#set variables
CHROM=$SLURM_ARRAY_TASK_ID

outputDir='/scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/pi'
metaData='/scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData'

#load vcftools
ml biology vcftools/0.1.15

#loop through populations

ref=grayfox


	for pop in {SanMiguel,SanClemente,SanNicolas,SantaCatalina,SantaCruz,SantaRosa,Mainland}
	do

	outputFile=${outputDir}/'chrom'${CHROM}'_'${ref}'_'${pop}'_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz'

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/allsites/splitChroms/chr"$CHROM"_"$ref"_filtered.renameChroms.ACgr25_DPgr165lt500.gvcf.gz --site-pi --keep ${metaData}/${pop}_refGenome.txt --chr chr${CHROM} --out ${outputFile}
	
	done
echo "done $ref"


#sleep timer
sleep 180
