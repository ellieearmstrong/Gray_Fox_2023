
###bypass centromeres
for ref in {Canfam3.1,Canfam4}
do

/scratch/users/elliea/jazlyn-ellie/grayfox_2023/software/garlic/bin/linux/garlic --tped /scratch/users/elliea/jazlyn-ellie/grayfox_2023/plink/"$ref"_filtered.renameChroms.ACgr25_DPgr165lt500.tped --tfam /scratch/users/elliea/jazlyn-ellie/grayfox_2023/plink/"$ref"_filtered.renameChroms.ACgr25_DPgr165lt500.tfam --centromere fakeCentromeres_dog.txt --error 0.001 --winsize 100 --auto-winsize --auto-overlap-frac --out "$ref"

done



for ref in {grayfox,arcticfox}
do

/scratch/users/elliea/jazlyn-ellie/grayfox_2023/software/garlic/bin/linux/garlic --tped /scratch/users/elliea/jazlyn-ellie/grayfox_2023/plink/"$ref"_filtered.renameChroms.ACgr25_DPgr165lt500.tped --tfam /scratch/users/elliea/jazlyn-ellie/grayfox_2023/plink/"$ref"_filtered.renameChroms.ACgr25_DPgr165lt500.tfam --centromere fakeCentromeres_fox.txt --error 0.001 --winsize 100 --auto-winsize --auto-overlap-frac --out "$ref"

done