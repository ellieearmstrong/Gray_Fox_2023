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
#for ref in {arcticfox,grayfox,Canfam3.1,Canfam4};
#do
#grep "" *"$ref".bed | grep -v "track" | sed -e 's/.bed:/\t/g' | awk '{print $2"\t"$3"\t"$4"\t"$1}' | sort -k 1,1 -k2,2n | sed -e 's/chr//g' > allIndivs_"$ref"_sorted.bed
#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4 > "chrom"$1"_allIndivs_'$ref'_sorted.bed"}' allIndivs_"$ref"_sorted.bed
#echo $ref
#done

#load modules that you need
#ml biology bedtools

#arctic fox
#for f in {1..24}
#do
#        coverageBed -a chrom"$f"_allIndivs_arcticfox_sorted.bed -b /scratch/users/elliea/jazlyn-ellie/grayfox_2023/callableSites/arcticfox/chr"$f"_callableSites_arcticfox_merged.bed -sorted > overlapCallableChannelIsland/chrom"$f"_callableSites_arcticfox_ROH.bed
#done
#echo "done arctic fox"

#gray fox
for f in {1..32}
do
	coverageBed -a chrom"$f"_allIndivs_grayfox_sorted.bed -b /scratch/users/elliea/jazlyn-ellie/grayfox_2023/callableSites/grayfox/chr"$f"_callableSites_grayfox_merged.bed -sorted > overlapCallableChannelIsland/chrom"$f"_callableSites_grayfox_ROH.bed
done
echo "done gray fox"

#can fam
#for ref in {Canfam3.1,Canfam4}
#do
#	for f in {1..38}
#	do

#		coverageBed -a chrom"$f"_allIndivs_"$ref"_sorted.bed -b /scratch/users/elliea/jazlyn-ellie/grayfox_2023/callableSites/"$ref"/chr"$f"_callableSites_"$ref"_merged.bed -sorted > overlapCallableChannelIsland/chrom"$f"_callableSites_"$ref"_ROH.bed
#	done
#echo $ref
#done