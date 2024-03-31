#!/bin/bash


for i in *.roh.bed
do

	# Read the input file
	while IFS= read -r line; do
    	# Check if the line starts with "track name="
    	if [[ $line == track\ name=* ]]; then
        	# Extract the individual ID from the line
        	ind_id=$(echo "$line" | sed 's/.*Ind:\s*\(\S*\)\s*Pop:.*/\1/p')
    	else
        	# Print the line with the individual ID as the first column
        	echo "$ind_id $line"
    	fi
	done < "$i" > int.out


	awk 'NF > 1  {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' int.out > "indivID_"$i

echo "$i"

done

rm int.out
