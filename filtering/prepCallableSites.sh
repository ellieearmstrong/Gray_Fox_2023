
#This code will prep the callable sites

for ref in {arcticfox,grayfox,Canfam3.1,Canfam4} 
do 

#bedtools merge -i callableSites_"$ref".bed > callableSites_"$ref"_merged.bed

#awk '{print $0 > $1"_callableSites_'$ref'_merged.bed"}' callableSites_"$ref"_merged.bed

echo "$ref"

done
