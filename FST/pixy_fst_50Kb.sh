#!/bin/sh
#SBATCH --job-name=fstAF
#SBATCH --output=/project/jazlynmo_738/Jazlyn/grayfox/logfiles/fst_%a_af_50kb.out  #you will need to modify the path
#SBATCH --error=/project/jazlynmo_738/Jazlyn/grayfox/logfiles/fst_%a_af_50kb.err  #you will need to modify the path
#SBATCH --time=1:00:00 #ten hour run time
#SBATCH -p qcb #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32 #one cpu per task
#SBATCH --array=1-38
#SBATCH --mem-per-cpu=1000MB #this is equivalent to 10G
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
##compute weir and cockerham fst in 50kb windows for each chromosome 

pixy --stats fst --vcf /project/jazlynmo_738/DataRepository/Canids/Variants/GrayFox/Mainland/"$infile".vcf.gz --populations /project/jazlynmo_738/Jazlyn/grayfox/metaData/pops_pixy_n12.txt --window_size 50000 --chromosomes "chr"$CHROM --bypass_invariant_check 'yes' --n_cores 32 --output_folder /project/jazlynmo_738/Jazlyn/grayfox/Analyses/FST/output_fst_50Kb/ --output_prefix $infile"_"$CHROM

echo $infile


done < vcfIn_Mainland.txt


####Channel Island

#while read -r infile
#do

##pixy will not allow the path in the output_prefix you would have to specificy that with a folder
##already know this file include invariant sites so I don't need the check
##compute weir and cockerham fst in 100kb windows for each chromosome

#pixy --stats fst --vcf /project/jazlynmo_738/DataRepository/Canids/Variants/GrayFox/ChannelIsland/"$infile".vcf.gz --populations /project/jazlynmo_738/Jazlyn/grayfox/metaData/pops_pixy_ChannelIsland.txt --window_size 50000 --chromosomes "chr"$CHROM --bypass_invariant_check 'yes' --n_cores 32 --output_folder /project/jazlynmo_738/Jazlyn/grayfox/Analyses/FST/output_fst_50Kb_ChannelIsland_50Kb/ --output_prefix $infile"_"$CHROM

#echo $infile


#done < vcfIn_ChannelIsland.txt


sleep 180
