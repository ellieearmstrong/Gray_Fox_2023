#SBATCH --job-name=slidingwin
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/logfiles/slidingwin.out  #you will need to modify the path
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/logfiles/slidingwin.err  #you will need to modify the path
#SBATCH --time=10:00:00 #run time
#SBATCH -p normal #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --mem-per-cpu=4000MB #this is equivalent to 4G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlyn.mooney@gmail.com

#load python and pysam
ml biology python/2.7.13
ml biology py-pysam/0.14.1_py27

inputDir='/scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/splitChroms'
outputDir='/scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/slidingWindow/'

#finish can fam genomes

SECONDS=0

for ref in {Canfam3.1,Canfam4}
do
        for CHROM in {1..38}
        do
	#command to run script
	master_VCF=${inputDir}/'chr'${CHROM}'_'${ref}'_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz'
	python SlidingWindowHet_Unphased.jar.ab.jam.${ref}.py --vcf ${master_VCF} --window_size 100000 --step_size 10000 --chromNum chr${CHROM} --outpath ${outputDir}

	done

echo "done $ref"

done

echo $SECONDS

#run for gray fox
SECONDS=0
for CHROM in {1..32}
do
	#command to run script
        master_VCF=${inputDir}/'chr'${CHROM}'_grayfox_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz'
        python SlidingWindowHet_Unphased.jar.ab.jam.grayfox.py --vcf ${master_VCF} --window_size 100000 --step_size 10000 --chromNum chr${CHROM} --outpath ${outputDir}

done
echo $SECONDS


#run for arctic fox
SECONDS=0
for CHROM in {1..24}
do
        #command to run script
        master_VCF=${inputDir}/'chr'${CHROM}'_arcticfox_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz'
        python SlidingWindowHet_Unphased.jar.ab.jam.arcticfox.py --vcf ${master_VCF} --window_size 100000 --step_size 10000 --chromNum chr${CHROM} --outpath ${outputDir}

done
echo $SECONDS

