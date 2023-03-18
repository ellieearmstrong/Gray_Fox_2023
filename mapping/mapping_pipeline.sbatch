#!/bin/bash
##SBATCH --job-name=mapping-variantcall
##SBATCH --time=144:00:00
##SBATCH --mem=240000
##SBATCH --partition=dpetrov
##SBATCH -n 36
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=ellieearmstrong@gmail.com

#load modules
module load biology bwa/0.7.17
module load biology samtools/1.8
module load biology gatk/4.1.4.1

#set reference
REFERENCE=/oak/stanford/groups/dpetrov/armstrong/gray-foxes/other_canid_genomes/GCF_000002285.3_CanFam3.1_genomic.fna
THREADS_BWA=46

#prepare reference for mapping
java -jar /oak/stanford/groups/dpetrov/armstrong/applications/picard.jar CreateSequenceDictionary R=$REFERENCE
samtools faidx $REFERENCE
bwa index -a bwtsw $REFERENCE

##first collect list of chromosomes in genome
cut -f1 "$REFERENCE".fai > intervals.txt
split --number=l/10 intervals.txt

#create intervals list for GATK4
mv xaa intervals1.txt
mv xab intervals2.txt
mv xac intervals3.txt
mv xad intervals4.txt
mv xae intervals5.txt
mv xaf intervals6.txt
mv xag intervals7.txt
mv xah intervals8.txt
mv xai intervals9.txt
mv xaj intervals10.txt

##mapping to genome
for file in *_1.fastq.gz;do
	bn=`basename $file _1.fastq.gz`
	bwa mem -t $THREADS_BWA $REFERENCE ${bn}_1.fastq.gz ${bn}_2.fastq.gz | samtools view -bS - | samtools sort - > ${bn}.sorted.bam
	samtools index ${bn}.sorted.bam
	java -jar /oak/stanford/groups/dpetrov/armstrong/applications/picard.jar AddOrReplaceReadGroups I=${bn}.sorted.bam O=${bn}.sorted.rg.bam SORT_ORDER=coordinate RGID=${bn} RGLB=${bn} RGPL=illumina RGSM=${bn} RGPU=${bn} CREATE_INDEX=True
	java -jar /oak/stanford/groups/dpetrov/armstrong/applications/picard.jar MarkDuplicates I=${bn}.sorted.rg.bam O=${bn}.sorted.rg.md.bam M=${bn}.metrics.md.txt 
	samtools index ${bn}.sorted.rg.md.bam
done

#merge samples-this step is manual for now, then proceed
##merge file if sample has more than one set of reads
##for file in list;do
##/oak/stanford/groups/dpetrov/armstrong/applications/sambamba-0.7.1-linux-static merge -t 20  

#clean folder
#rm *.sorted.bam
#rm *.sorted.bam.bai
#rm *.sorted.rg.bam
#rm *.sorted.rg.bam.bai

#generate depth and coverage stats
#ls *.sorted.rg.md.bam | -j 10 | /oak/stanford/groups/dpetrov/armstrong/applications/mosdepth {} {.}.sorted.rg.md.bam -n 
#done


#run haplotype caller per chromosome
#for file in *.sorted.rg.md.bam;do
#	bn=`basename .sorted.rg.md.bam`
#	gatk HaplotypeCaller -I ${bn}.sorted.rg.md.bam -O ${bn}.rg.md.hcint1.g.vcf -L intervals1.txt -ERC GVCF --emit-ref-confidence BP_RESOLUTION -R $REFERENCE
#	gatk HaplotypeCaller -I ${bn}.sorted.rg.md.bam -O ${bn}.rg.md.hcint2.g.vcf -L intervals2.txt -ERC GVCF --emit-ref-confidence BP_RESOLUTION -R $REFERENCE
#	gatk HaplotypeCaller -I ${bn}.sorted.rg.md.bam -O ${bn}.rg.md.hcint3.g.vcf -L intervals3.txt -ERC GVCF --emit-ref-confidence BP_RESOLUTION -R $REFERENCE
#	gatk HaplotypeCaller -I ${bn}.sorted.rg.md.bam -O ${bn}.rg.md.hcint4.g.vcf -L intervals4.txt -ERC GVCF --emit-ref-confidence BP_RESOLUTION -R $REFERENCE
#	gatk HaplotypeCaller -I ${bn}.sorted.rg.md.bam -O ${bn}.rg.md.hcint5.g.vcf -L intervals5.txt -ERC GVCF --emit-ref-confidence BP_RESOLUTION -R $REFERENCE
#	gatk HaplotypeCaller -I ${bn}.sorted.rg.md.bam -O ${bn}.rg.md.hcint6.g.vcf -L intervals6.txt -ERC GVCF --emit-ref-confidence BP_RESOLUTION -R $REFERENCE
#	gatk HaplotypeCaller -I ${bn}.sorted.rg.md.bam -O ${bn}.rg.md.hcint7.g.vcf -L intervals7.txt -ERC GVCF --emit-ref-confidence BP_RESOLUTION -R $REFERENCE
#	gatk HaplotypeCaller -I ${bn}.sorted.rg.md.bam -O ${bn}.rg.md.hcint8.g.vcf -L intervals8.txt -ERC GVCF --emit-ref-confidence BP_RESOLUTION -R $REFERENCE
#	gatk HaplotypeCaller -I ${bn}.sorted.rg.md.bam -O ${bn}.rg.md.hcint9.g.vcf -L intervals9.txt -ERC GVCF --emit-ref-confidence BP_RESOLUTION -R $REFERENCE
#	gatk HaplotypeCaller -I ${bn}.sorted.rg.md.bam -O ${bn}.rg.md.hcint10.g.vcf -L intervals10.txt -ERC GVCF --emit-ref-confidence BP_RESOLUTION -R $REFERENCE
#done
#

#note-this is only if you need an all sites file, otherwise switch to genomicsDBimport
#for file in *.rg.md.hcint.g.vcf;do
#
