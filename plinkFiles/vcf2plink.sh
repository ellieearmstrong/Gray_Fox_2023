#!/bin/sh

#Load library
ml plink2/2.00a3.7LM

for f in {Mainland,ChannelIsland}
do
	while read -r p

	do

	plink2 --vcf /scratch1/jazlynmo/grayfox/vcf/"$f"/"$p".vcf.gz --allow-extra-chr --chr-set 38 --threads 32 --make-bed --out /scratch1/jazlynmo/grayfox/plink/"$f"/"$p"

	echo $p

	done < vcfIn_"$f".txt

done

#Load new library
#convert to tped
ml gcc/11.3.0  openblas/0.3.21 plink/1.9-beta7

for f in {Mainland,ChannelIsland}
do
	while read -r p

	do

	plink --bfile /scratch1/jazlynmo/grayfox/plink/"$f"/$p  --allow-extra-chr --chr-set 38 --recode transpose --threads 32 --out /scratch1/jazlynmo/grayfox/plink/"$f"/"$p"

	done < vcfIn_"$f".txt

done
