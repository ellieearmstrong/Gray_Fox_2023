#!/bin/sh
#SBATCH --job-name=pi
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/logfiles/pi_%a.out  #you will need to modify the path
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/logfiles/pi_%a.err  #you will need to modify the path
#SBATCH --time=02:00:00 #ten hour run time
#SBATCH -p normal #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --array=1-32
#SBATCH --mem-per-cpu=500MB #this is equivalent to 10G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jaam@stanford.edu


#set up variables
CHROM=$SLURM_ARRAY_TASK_ID
inputDir='/scratch/users/elliea/jazlyn-ellie/grayfox_2023/allsites/splitChroms/filter/mergeVarInvar'
outputDir='/scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/pi'
metaData='/scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData'

#load vcftools
ml biology vcftools/0.1.15


#load master vcf
master_VCF=${inputDir}/'chrom'${CHROM}'_allFoxes.rg.md.haplotypecaller.all.g.renameChroms.mappabilityFilter.AN.QUAL.DP.biallelic.varInvar.sorted.vcf.gz'

#loop through populations
for pop in {SanMiguel,SanClemente,SanNicolas,SantaCatalina,SantaCruz,SantaRosa,Cinereoargenteus}
do

outputFile=${outputDir}/'chrom'${CHROM}'_'${pop}'.rg.md.haplotypecaller.all.g.renameChroms.mappabilityFilter.AN.QUAL.DP.biallelic.varInvar.sorted.pi'

#autosomal pi
vcftools --gzvcf ${master_VCF} --site-pi --keep ${metaData}/${pop}.txt --chr ${CHROM} --out ${outputFile}

done

#sleep timer
#sleep 180
