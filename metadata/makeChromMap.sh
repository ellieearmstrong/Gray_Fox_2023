#pull all scaffold names
zgrep -v "##" fox_merge_variants_bisnps_autosome.GM.eh2.AN.qualdp.vcf.gz | awk '{print $1}' | uniq > chroms.txt

#add a number
awk 'NR >1 {print $0"\t" NR-1}' chroms.txt > fox_chrom.map
