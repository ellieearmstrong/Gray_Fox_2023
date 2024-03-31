#!/bin/sh
#SBATCH --job-name=slideHet
#SBATCH --output=/scratch1/jazlynmo/grayfox/logfiles/slideHet.out  #you will need to modify the path
#SBATCH --error=/scratch1/jazlynmo/grayfox/logfiles/slideHet.err  #you will need to modify the path
#SBATCH --time=10:00:00 #one hour run time
#SBATCH -p qcb #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000MB #this is equivalent to 5G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlynmo@usc.edu

#load python and pysam
module load gcc/11.3.0 openblas/0.3.20 python/3.9.12 py-pysam/0.18.0


inputDir='/scratch1/jazlynmo/grayfox/vcf/Mainland/splitChroms/'
outputDir='/scratch1/jazlynmo/grayfox/Analyses/slidingWindow/output_Mainland/'


#run for gray fox
SECONDS=0
for CHROM in {31..32}
do
        #command to run script
        master_VCF=${inputDir}/'grayfox_filtered.renameChroms.Mainland.drop295.ACgr59_DPgr205lt500_chr'${CHROM}'.vcf.gz'
        python3 SlidingWindowProvideChromSize.py --vcf ${master_VCF} --window_size 100000 --step_size 10000 --chromNum chr${CHROM} --chromLengths grayfox_chromosome_lengths.txt --outpath ${outputDir}

done
echo $SECONDS

#run can fam and graywolf
SECONDS=0
for ref in {Canfam3.1,Canfam4,graywolf}
do
        for CHROM in {37..38}
        do
	#command to run script
	master_VCF=${inputDir}/${ref}'_filtered.renameChroms.Mainland.drop295.ACgr59_DPgr205lt500_chr'${CHROM}'.vcf.gz'
	python3 SlidingWindowProvideChromSize.py --vcf ${master_VCF} --window_size 100000 --step_size 10000 --chromNum chr${CHROM} --chromLengths ${ref}'_chromosome_lengths.txt' --outpath ${outputDir}

	done

echo "done $ref"

done

echo $SECONDS

#run for arctic fox
SECONDS=0
for CHROM in {23..24}
do
        #command to run script
        master_VCF=${inputDir}/'arcticfox_filtered.renameChroms.Mainland.drop295.ACgr59_DPgr205lt500_chr'${CHROM}'.vcf.gz'
        python3 SlidingWindowProvideChromSize.py --vcf ${master_VCF} --window_size 100000 --step_size 10000 --chromNum chr${CHROM} --chromLengths arcticfox_chromosome_lengths.txt --outpath ${outputDir}

done
echo $SECONDS
