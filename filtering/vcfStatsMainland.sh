#!/bin/sh
#SBATCH --job-name=vcfstats
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/QCMainland.out  #you will need to modify the path
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/QCMainland.err  #you will need to modify the path
#SBATCH --time=24:00:00 #ten hour run time
#SBATCH -p normal #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --mem-per-cpu=1000MB #this is equivalent to 1G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlynmo@usc.edu

#load modules that you need
ml biology bcftools/1.16
ml biology vcftools/0.1.15
#ml python/3.9.0

#get DP and AN per site
bcftools view -s ^SRR24465295 /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/arcticfox_filtered.renameChroms.Mainland.gvcf.gz | bcftools query -f '%DP\n'  > arcticfox_gvcfstats_Mainland_rmSRR24465295.txt
vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/arcticfox_filtered.renameChroms.Mainland.gvcf.gz --remove-indv SRR24465295 --missing-site --out arcticfox_gvcfstats_Mainland_rmSRR24465295
echo "done"

bcftools view -s ^SRR24465295 /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/grayfox_filtered.renameChroms.Mainland.gvcf.gz | bcftools query -f '%DP\n' > grayfox_gvcfstats_Mainland_rmSRR24465295.txt
vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/grayfox_filtered.renameChroms.Mainland.gvcf.gz --remove-indv SRR24465295 --missing-site --out grayfox_gvcfstats_Mainland_rmSRR24465295
echo "done"

bcftools view -s ^SRR24465295 /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.Mainland.gvcf.gz | bcftools query -f '%DP\n' > canFam3.1_gvcfstats_Mainland_rmSRR24465295.txt
vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam3.1_filtered.renameChroms.Mainland.gvcf.gz --remove-indv SRR24465295 --missing-site --out canFam3.1_gvcfstats_Mainland_rmSRR24465295
echo "done"

bcftools view -s ^SRR24465295 /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.Mainland.gvcf.gz | bcftools query -f '%DP\n' > canFam4_gvcfstats_Mainland_rmSRR24465295.txt
vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.Mainland.gvcf.gz --remove-indv SRR24465295 --missing-site --out canFam4_gvcfstats_Mainland_rmSRR24465295
echo "done"

#run python code to concatenate
#for i in {canFam3.1,canFam4,arcticfox,grayfox}; do awk 'NR > 1 {print $5}' "$i"_gvcfstats_Mainland_rmSRR24465295.lmiss | paste "$i"_gvcfstats_Mainland_rmSRR24465295.txt - > "$i"_gvcfstats.DP_NMISS_Mainland_rmSRR24465295.out; echo "$i"; done

#get mean and sd from stats
#python3 calcSumStats_v4.py >> mean_sd_Mainland_rmSRR24465295.txt

#done
