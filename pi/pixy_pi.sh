#!/bin/sh
#SBATCH --job-name=piAF
#SBATCH --output=/scratch1/jazlynmo/grayfox/logfiles/pi_%a_af.out  #you will need to modify the path
#SBATCH --error=/scratch1/jazlynmo/grayfox/logfiles/pi_%a_af.err  #you will need to modify the path
#SBATCH --time=02:00:00 #ten hour run time
#SBATCH -p qcb #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32 #32 cpus per task
#SBATCH --array=1-8
#SBATCH --mem-per-cpu=1000MB #this is equivalent to 1G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlynmo@usc.edu

#load everything up
module purge
eval "$(conda shell.bash hook)"
module load gcc/8.3.0 python/3.6.8 
conda activate pixy-env

#set variables
CHROM=$SLURM_ARRAY_TASK_ID

####Mainland

while read -r infile
do 

##pixy will not allow the path in the output_prefix you would have to specificy that with a folder
##already know this file include invariant sites so I don't need the check
##compute pi in 50kb windows for each chromosome 

pixy --stats pi --vcf /project/jazlynmo_738/DataRepository/Canids/Invariant/GrayFox/Mainland/"$infile".gvcf.gz --populations /project/jazlynmo_738/Jazlyn/grayfox/metaData/pops_pixy_n12.txt --window_size 50000 --chromosomes "chr"$CHROM --bypass_invariant_check 'yes' --n_cores 32 --output_folder /project/jazlynmo_738/Jazlyn/grayfox/Analyses/pi/output_pi_50Kb/ --output_prefix $infile"_"$CHROM

pixy --stats pi --vcf /project/jazlynmo_738/DataRepository/Canids/Invariant/GrayFox/Mainland/"$infile".gvcf.gz --bed_file /project/jazlynmo_738/Jazlyn/grayfox/callableSites/chr"$CHROM"_callableSites_"$infile"_merged.bed --populations /project/jazlynmo_738/Jazlyn/grayfox/metaData/pops_pixy_n12.txt --chromosomes "chr"$CHROM --bypass_invariant_check 'yes' --n_cores 32 --output_folder /project/jazlynmo_738/Jazlyn/grayfox/Analyses/pi/output_pi_noMiss/ --output_prefix $infile"_"$CHROM

echo $infile

done < gvcfIn_Mainland_2.txt


####Channel Island

#while read -r infile
#do

##pixy will not allow the path in the output_prefix you would have to specificy that with a folder
##already know this file include invariant sites so I don't need the check
##compute weir and cockerham fst in 50kb windows for each chromosome

#pixy --stats pi --vcf /scratch1/jazlynmo/grayfox/vcf/ChannelIsland/"$infile".gvcf.gz --populations /project/jazlynmo_738/Jazlyn/grayfox/metaData/pops_pixy_ChannelIsland.txt --window_size 50000 --chromosomes "chr"$CHROM --bypass_invariant_check 'yes' --n_cores 32 --output_folder /scratch1/jazlynmo/grayfox/Analyses/pi/output_pi_ChannelIsland/ --output_prefix $infile"_"$CHROM

#echo $infile


#done < gvcfIn_ChannelIsland.txt

#sleep 180
