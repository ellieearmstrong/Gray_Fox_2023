#!/bin/sh
#SBATCH --job-name=roh
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/logfiles/rohsplitchroms.out  #you will need to modify the path
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/logfiles/rohsplitchroms.err  #you will need to modify the path
#SBATCH --time=02:00:00 #ten hour run time
#SBATCH -p normal #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --mem-per-cpu=40000MB #this is equivalent to 10G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jaam@stanford.edu


####before you run this code you need to: 1) split ROH for each individual 2) make a bedfile with all of your callable sites
####once you have a file with ROH for each individual run this code to reformat, sort, and concatenate
####grep "" *.bed | grep -v "track" | sed -e 's/.bed:/\t/g' | awk '{print $2"\t"$3"\t"$4"\t"$1}' | sort -k 1,1 -k2,2n | sed -e 's/chr//g' > allIndivs_sorted.bed
####awk '{print $1"\t"$2"\t"$3"\t"$4 > "chrom"$1"_allIndivs_sorted.bed"}' allIndivs_sorted.bed

#load modules that you need
ml biology bedtools
for f in {1..32}
do
	echo "$f"
	coverageBed -a chrom"$f"_allIndivs_sorted.bed -b /scratch/users/elliea/jazlyn-ellie/grayfox_2023/allsites/splitChroms/filter/mergeVarInvar/callableSites/chrom"$f"_callableSites.bed -sorted > overlapCallable/chrom"$f"_callableSitesROH.bed
done 
