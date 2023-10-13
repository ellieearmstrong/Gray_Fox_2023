 



while read -r p
do 

/scratch/users/elliea/jazlyn-ellie/grayfox_2023/software/truffle/truffle --vcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$p" --cpu 4 --segments --missing 1 --maf 0 --nofiltering --out "$p".truffle

echo "$p"

done < vcfIn2.txt
