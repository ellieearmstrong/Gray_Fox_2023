ml biology vcftools
ml biology plink/1.90b5.3

#Covert vcf to plink
#vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/fox_merge_variants_snps_autosome.GM.AN.qualdp.vcf.gz --chrom-map fox_chrom.map --plink --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/plinkFiles/fox_merge_variants_snps_autosome.GM.AN.qualdp.vcf.gz 

#Convert ped and map to binaries
#plink --ped /scratch/users/elliea/jazlyn-ellie/grayfox_2023/plinkFiles/fox_merge_variants_snps_autosome.GM.AN.qualdp.vcf.gz.ped --map /scratch/users/elliea/jazlyn-ellie/grayfox_2023/plinkFiles/fox_merge_variants_snps_autosome.GM.AN.qualdp.vcf.gz.map --allow-extra-chr --chr-set 32 --make-bed --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/plinkFiles/fox_merge_variants_bisnps_autosome.GM.eh2.AN.qualdp

#update family id
plink --bfile /scratch/users/elliea/jazlyn-ellie/grayfox_2023/plinkFiles/fox_merge_variants_bisnps_autosome.GM.eh2.AN.qualdp --update-ids recodeFamID.txt --allow-extra-chr --chr-set 32 --make-bed --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/plinkFiles/fox_merge_variants_bisnps_autosome.GM.eh2.AN.qualdp.updatedFID

#recode to tped and tfam per family
plink --bfile /scratch/users/elliea/jazlyn-ellie/grayfox_2023/plinkFiles/fox_merge_variants_bisnps_autosome.GM.eh2.AN.qualdp.updatedFID --chr-set 32 --recode transpose --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/plinkFiles/fox_merge_variants_bisnps_autosome.GM.eh2.AN.qualdp.updatedFID



