

1) subset YRI from Calres unfiltered file and the 1000 Genomes unfiltered file:

Used script .... to do this

2) print the pos column of the subsetted vcf files
	Printed column using awk, make sure to grep out hashes before 

3) add header to the position column so that I could run my annotation script
	Added header using sed: sed -i 1i'POS' *AllSites*.txt

4) ran my script addPositionAttributes.py to figure out which positions were intronic exonic etc.
	Submitted script as follows: ./runLocal_addPosAtt.sh

Commands in script:
for f in {1..22}; do python addPositionAttributes.py Chr"$f"_AllSitesClare.txt 1 Chr"$f"_allSitesClareAnnot.txt ranges Chr"$f"_exon_1Based.txt Exonic Chr"$f"_intron_1Based_sorted.txt Intronic Chr"$f"_intergenic_UCSCForum_1Based_sorted.txt Intergenic; done

for f in {1..22}; do python addPositionAttributes.py Chr"$f"_AllSites1G.txt 1 Chr"$f"_allSites1GAnnot.txt ranges Chr"$f"_exon_1Based.txt Exonic Chr"$f"_intron_1Based_sorted.txt Intronic Chr"$f"_intergenic_UCSCForum_1Based_sorted.txt Intergenic; done

5) Made files that listed the sites that fell into each region of the genome 

Commands:

for f in {1..22}; do awk '$2==1 {print $1}' Chr"$f"_allSites1GAnnot.txt > Chr"$f"_exonicSites1G.txt; done
for f in {1..22}; do awk '$3==1 {print $1}' Chr"$f"_allSites1GAnnot.txt > Chr"$f"_intronicSites1G.txt; done
for f in {1..22}; do awk '$4==1 {print $1}' Chr"$f"_allSites1GAnnot.txt > Chr"$f"_intergenicSites1G.txt; done


for f in {1..22}; do awk '$2==1 {print $1}' Chr"$f"_allSitesClareAnnot.txt > Chr"$f"_exonicSitesClare.txt; done
for f in {1..22}; do awk '$3==1 {print $1}' Chr"$f"_allSitesClareAnnot.txt > Chr"$f"_intronicSitesClare.txt; done
for f in {1..22}; do awk '$4==1 {print $1}' Chr"$f"_allSitesClareAnnot.txt > Chr"$f"_intergenicSitesClare.txt; done

6)Made all the data in VCFs unphased

Commands:
for f in {1..22}; do python convert_to_unphased.py YRI_Chr"$f"_UnfilteredfromAllPops_Subset.recode.vcf YRI_Chr"$f"_UnfilteredfromAllPops_SubsetUNPHASED.recode.vcf

for f in {1..22}; do python convert_to_unphased.py YRI_Chr"$f"_mergedUnfilteredfromClare_Subset.recode.vcf YRI_Chr"$f"_mergedUnfilteredfromClare_SubsetUNPHASED.recode.vcf; done

7) Used the aforementioned genomic region files to annotate positions in VCFs and calculate K for Wattersons theta

	Script used to calculate Wattersons Theta: CalculateK_forWattersonsTheta_v2_modNonBiallelic.py
	Submitted script as follows: 
					qsub submitWatTheta.sh
					qsub WatTheta2.sh 

8) Pull out sites that were not biallelic

	WHY: subtract non-biallelic segregating sites for each genomic region from numerator since Watterson's theta measures SNP density. HOWEVER, the denominator will remain the same because the denominator is just total sites in each genomic region

	
Commands:
#Find non-bialleic sites
for f in {1..22}; do awk 'NR==FNR{c[$3]++;next};c[$2] > 0' YRI_Chr"$f"_mergedUnfilteredfromClare_WatThetaSummary.out YRI_Chr"$f"_mergedUnfilteredfromClare_WatTheta.out > Chr"$f"_SitesNotBiallelic_mergedClare.txt ;done


for f in {1..22}; do awk 'NR==FNR{c[$3]++;next};c[$2] > 0' YRI_Chr"$f"_UnfilteredfromAllPops_WatThetaSummary.out YRI_Chr"$f"_UnfilteredfromAllPops_WatTheta.out > Chr"$f"_SitesNotBiallelic_UnfilteredAllPops.txt ;done


#Classify non-biallelic sites by genomic region
for f in {1..22}; do awk '$3==1 && $6==1 {print $1"\t"$2}' Chr"$f"_SitesNotBiallelic_UnfilteredAllPops.txt >> SubstractFromTotalSegIntergenicSites_UnfilteredAllPops.txt; done
for f in {1..22}; do awk '$3==1 && $5==1 {print $1"\t"$2}' Chr"$f"_SitesNotBiallelic_UnfilteredAllPops.txt >> SubstractFromTotalSegIntronicSites_UnfilteredAllPops.txt; done 
for f in {1..22}; do awk '$3==1 && $4==1 {print $1"\t"$2}' Chr"$f"_SitesNotBiallelic_UnfilteredAllPops.txt >> SubstractFromTotalSegExonicSites_UnfilteredAllPops.txt; done

for f in {1..22}; do awk '$3==1 && $6==1 {print $1"\t"$2}' Chr"$f"_SitesNotBiallelic_mergedClare.txt >> SubstractFromTotalIntergenicSegSites_mergedClare.txt; done
for f in {1..22}; do awk '$3==1 && $5==1 {print $1"\t"$2}' Chr"$f"_SitesNotBiallelic_mergedClare.txt >> SubstractFromTotalIntronicSegSites_mergedClare.txt; done
for f in {1..22}; do awk '$3==1 && $4==1 {print $1"\t"$2}' Chr"$f"_SitesNotBiallelic_mergedClare.txt >> SubstractFromTotalExonicSegSites_mergedClare.txt; done


#Count non-biallelic sites
for f in {1..22}; do awk '{print $2}' Chr"$f"_SitesNotBiallelic_UnfilteredAllPops.txt | uniq -c | wc -l >> NonBiallelicSitesTotal_UnfilteredAllPops.txt; done

for f in {1..22}; do awk '{print $2}' Chr"$f"_SitesNotBiallelic_mergedClare.txt | uniq -c | wc -l >> NonBiallelicSitesTotal_mergedClare.txt; done

awk '{print $2}' SubstractFromTotalExonicSegSites_mergedClare.txt | uniq -c | wc -l > CountTotalSegExonicSites_mergedClare.txt
awk '{print $2}' SubstractFromTotalIntergenicSegSites_mergedClare.txt | uniq -c | wc -l > CountTotalSegIntergenicSites_mergedClare.txt
awk '{print $2}' SubstractFromTotalIntronicSegSites_mergedClare.txt | uniq -c | wc -l > CountTotalSegIntronicSites_mergedClare.txt

awk '{print $2}' SubstractFromTotalSegIntronicSites_UnfilteredAllPops.txt | uniq -c | wc -l > CountTotalSegIntronicSites_UnfilteredAllPops.txt
awk '{print $2}' SubstractFromTotalSegIntergenicSites_UnfilteredAllPops.txt | uniq -c | wc -l > CountTotalSegIntergenicSites_UnfilteredAllPops.txt
awk '{print $2}' SubstractFromTotalSegExonicSites_UnfilteredAllPops.txt | uniq -c | wc -l > CountTotalSegExonicSites_UnfilteredAllPops.txt

9) Count up segregating sites

Commands:

for f in {1..22}; do grep -w "segregating within an" YRI_Chr"$f"_UnfilteredfromAllPops_WatThetaSummary.out | sed  's/This many sites were segregating within an//g' | sed 's/region//g' |sed 's/Count of segregating sites (K)//g'| sed 's/in your file/total/g' >> YRI_SegregatingSites_UnfilteredAllPops.txt; done

for f in {1..22}; do grep -w "segregating" YRI_Chr"$f"_mergedUnfilteredfromClare_WatThetaSummary.out |sed 's/Count of segregating sites (K)//g'| sed 's/in your file/total/g'  | sed  's/This many sites were segregating within an//g' | sed 's/region//g' >> YRI_SegregatingSites_mergedClare.txt; done

10) Count up Total Sites
for f in {1..22}; do awk '!seen[$1]++' Chr"$f"_exonicSitesClare.txt | wc -l >> TotalExonicSites_UnfilteredClare.txt; done
less TotalExonicSites_UnfilteredClare.txt
for f in {1..22}; do awk '!seen[$1]++' Chr"$f"_intronicSitesClare.txt | wc -l >> TotalIntronicSites_UnfilteredClare.txt; done
for f in {1..22}; do awk '!seen[$1]++' Chr"$f"_intergenicSitesClare.txt | wc -l >> TotalIntergenicSites_UnfilteredClare.txt; done

for f in {1..22}; do awk '!seen[$1]++' Chr"$f"_intergenicSites1G.txt | wc -l >> TotalIntergenicSites_UnfilteredAllPops.txt; done
for f in {1..22}; do awk '!seen[$1]++' Chr"$f"_intronicSites1G.txt | wc -l >> TotalIntronicSites_UnfilteredAllPops.txt; done
for f in {1..22}; do awk '!seen[$1]++' Chr"$f"_exonicSites1G.txt | wc -l >> TotalexonicSites_UnfilteredAllPops.txt; done




