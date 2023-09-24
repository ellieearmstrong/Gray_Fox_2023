ml biology plink2.0/a2

#convert vcf to plink
for i in {arcticfox,grayfox}; do plink2 --vcf ../grayfox_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz --allow-extra-chr --chr-set 32 --make-bed --out grayfox_filtered.renameChroms.ACgr25_DPgr165lt500;done
for i in {Canfam3.1,Canfam4}; do plink2 --vcf ../grayfox_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz --allow-extra-chr --chr-set 38 --make-bed --out grayfox_filtered.renameChroms.ACgr25_DPgr165lt500;done


#convert to tped
ml biology plink/1.90b5.3

for i in {arcticfox,grayfox}; do plink --bfile "$i"_filtered.renameChroms.ACgr25_DPgr165lt500 --allow-extra-chr --chr-set 32 --recode transpose --out "$i"_filtered.renameChroms.ACgr25_DPgr165lt500; done

for i in {Canfam3.1,Canfam4}; do plink --bfile "$i"_filtered.renameChroms.ACgr25_DPgr165lt500 --allow-extra-chr --chr-set 38 --recode transpose --out "$i"_filtered.renameChroms.ACgr25_DPgr165lt500; done
