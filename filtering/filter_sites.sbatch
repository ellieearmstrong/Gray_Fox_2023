#!/bin/bash
#SBATCH --job-name=filter-foxes
#SBATCH --time=48:00:00
#SBATCH --mem=240000
#SBATCH --partition=dpetrov,normal,owners,hns
#SBATCH -n 24
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ellieearmstrong@gmail.com


module load biology bedtools
module load biology bcftools
##Sites filtering for gray foxes
#first restrict to only snps
#bcftools view  -m2 -M2 -v snps -O z -o fox_merge_variants_snps.vcf.gz fox_merge_variants.vcf.gz
#tabix -p vcf fox_merge_variants_snps.vcf.gz

####then restrict to only autosomes
#bcftools view -r Scaffold_25__2_contigs__length_69993150,Scaffold_16__1_contigs__length_63979375,Scaffold_10__1_contigs__length_73623094,Scaffold_27__1_contigs__length_47458807,Scaffold_31__1_contigs__length_40660448,Scaffold_22__1_contigs__length_54392966,Scaffold_26__1_contigs__length_49334718,Scaffold_30__1_contigs__length_41889812,Scaffold_21__1_contigs__length_56962569,Scaffold_14__1_contigs__length_66704989,Scaffold_23__1_contigs__length_53374214,Scaffold_15__1_contigs__length_64818041,Scaffold_32__3_contigs__length_92326701,Scaffold_9__1_contigs__length_74450539,Scaffold_6__1_contigs__length_90720172,Scaffold_29__2_contigs__length_79160281,Scaffold_33__2_contigs__length_63093800,Scaffold_24__1_contigs__length_52741366,Scaffold_17__1_contigs__length_59774042,Scaffold_2__1_contigs__length_102603151,Scaffold_20__1_contigs__length_58983324,Scaffold_8__1_contigs__length_79079208,Scaffold_13__1_contigs__length_71581588,Scaffold_18__2_contigs__length_64801261,Scaffold_28__2_contigs__length_62992706,Scaffold_4__1_contigs__length_96410777,Scaffold_5__1_contigs__length_91265922,Scaffold_1__1_contigs__length_117572690,Scaffold_12__1_contigs__length_72908785,Scaffold_3__1_contigs__length_98526141,Scaffold_7__1_contigs__length_82759371,Scaffold_11__1_contigs__length_73040223 -O z -o fox_merge_variants_snps_autosome.vcf.gz fox_merge_variants_snps.vcf.gz
#tabix -p vcf fox_merge_variants_snps_autosome.vcf.gz

#then remove sites with low mappability
bedtools subtract -a fox_merge_variants_snps_autosome.vcf.gz -b grayfox_mappability_genmap.1.0.bed -header | gzip > fox_merge_variants_snps_autosome.GM.vcf.gz
tabix -p vcf fox_merge_variants_snps_autosome.GM.vcf.gz

#rename scaffs so they arent so ugly
#zcat fox_merge_variants_snps_autosome.GM.vcf.gz | sed -e 's/Scaffold_/Scaffold/g' | bgzip -c > fox_merge_variants_snps_autosome.GM.eh.vcf.gz
#zcat fox_merge_variants_snps_autosome.GM.eh.vcf.gz | sed -e 's/__1_contigs_//g' | bgzip -c > fox_merge_variants_snps_autosome.GM.eh2.vcf.gz

#calculate quality stats
bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%AN\t%DP\n' fox_merge_variants_snps_autosome.GM.vcf.gz | bgzip -c > grayfox_quality.tab.gz

#filter for quality
bcftools filter -i 'AN > 28' fox_merge_variants_snps_autosome.GM.vcf.gz | bgzip -c > fox_merge_variants_snps_autosome.GM.AN.vcf.gz
tabix -p vcf fox_merge_variants_snps_autosome.GM.AN.vcf.gz
bcftools view -e 'QUAL < 20 || INFO/DP > 500 || INFO/DP < 80' fox_merge_variants_snps_autosome.GM.AN.vcf.gz -O z -o fox_merge_variants_snps_autosome.GM.AN.qualdp.vcf.gz
tabix -p vcf fox_merge_variants_snps_autosome.GM.AN.qualdp.vcf.gz


