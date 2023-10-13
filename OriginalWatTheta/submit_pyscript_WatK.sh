#!/bin/bash
#$ -cwd 
#$ -V 
#$ -N SegSitesCOCR
#$ -l time=02:00:00,h_data=1G 
#$ -M eplau 
#$ -m a


for i in {CO,CR}
do

for f in {1..22}
do



source /u/local/Modules/default/init/bash
module load python
 
python CalculateK_forWattersonsTheta.py "$i"_Chr"$f"_NewMasterFile_SubsetN30_biallelic_nomiss.recode.vcf "$i"_Chr"$f"_SegSites.txt
#python CalculateK_forWattersonsTheta.py RmIndivsHighAfrAncestry_CO_Chr"$f"_NewMasterFile_Subset.recode.vcf CO_Chr"$f"_SegSites_RmHighAfrAncestry.txt

done
done

