


**Canfam3**

To calculate average across all scaffolds
```
for file in *.mosdepth.summary.txt;do echo $file &&  awk '{ print $4} ;done
```
To calculate average across only autosomes for Canfam3
```
for file in *.mosdepth.summary.txt;do echo $file && head -n 39 $file | awk '{ sum += $4} END { print sum / 38}' ;done
```
