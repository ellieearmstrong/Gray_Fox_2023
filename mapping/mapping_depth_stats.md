# Calculating Mapping and Depth Summary Statistics for Fox Files
```
module load biology samtools
for file in *.sorted.md.rg.bam;do
        bn=`basename $file .sorted.md.rg.bam`
        /oak/stanford/groups/dpetrov/armstrong/applications/mosdepth -x -n -t 10 ${bn} ${bn}.sorted.md.rg.bam
        samtools flagstat ${bn}.sorted.md.rg.bam -O tsv > ${bn}.flagstat.stats
done
```

## Canfam3 

To calculate average across all scaffolds
```
for file in *.mosdepth.summary.txt;do echo $file &&  awk '{ print $4} ;done
```
To calculate average across only autosomes for Canfam3
```
for file in *.mosdepth.summary.txt;do echo $file && head -n 39 $file | awk '{ sum += $4} END { print sum / 38}' ;done
```
