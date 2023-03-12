

##qsub -cwd -V -N slidingWindows -l highp,time=00:20:00,h_data=10G -t 1-38:1 -M eplau -m a runSlidingWindow.sh
##CHROM=${SGE_TASK_ID}

for CHROM in {1..32}; 
do 

inputDir='/scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/splitChroms/'
master_VCF=${inputDir}/'chrom'${CHROM}'_fox_merge_variants_bisnps_autosome.GM.eh2.AN.qualdp.renamedChrom.vcf.gz'


#load python
ml biology python/2.7.13
ml biology py-pysam/0.14.1_py27

#command to run script
#python SlidingWindowHet_Unphased.jar.ab.jam.py --vcf ${master_VCF} --window_size 1000000 --step_size 1000000 --chromNum ${CHROM}
python SlidingWindowHet_Unphased.jar.ab.jam.py --vcf ${master_VCF} --window_size 10000 --step_size 10000 --chromNum ${CHROM}
done


