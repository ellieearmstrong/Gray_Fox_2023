#!/bin/sh
#SBATCH --job-name=FST/
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/logfiles/vcftoolsFST.out  #you will need to
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/grayfox_2023/logfiles/vcftoolsFST.err  #you will need to m
#SBATCH --time=30:00:00 #ten hour run time
#SBATCH -p normal #the main partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 #one cpu per task
#SBATCH --mem-per-cpu=1000MB #this is equivalent to 2G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jazlyn.mooney@gmail.com

#load modules that you need
ml biology vcftools

#compute all possible FST/s
for f in {arcticfox,grayfox,Canfam3.1,Canfam4}
	do

	#timer
	SECONDS=0

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanClemente_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanMiguel_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SanClemente_vs_SanMiguel

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanClemente_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanNicolas_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SanClemente_vs_SanNicolas

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanClemente_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaCatalina_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SanClemente_vs_SantaCatalina

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanClemente_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaCruz_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SanClemente_vs_SantaCruz

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanClemente_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaRosa_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SanClemente_vs_SantaRosa

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanClemente_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/Mainland_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SanClemente_vs_Mainland

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanMiguel_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanNicolas_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SanMiguel_vs_SanNicolas

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanMiguel_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaCatalina_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SanMiguel_vs_SantaCatalina

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanMiguel_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaCruz_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SanMiguel_vs_SantaCruz

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanMiguel_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaRosa_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SanMiguel_vs_SantaRosa

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanMiguel_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/Mainland_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SanMiguel_vs_Mainland

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanNicolas_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaCatalina_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SanNicolas_vs_SantaCatalina

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanNicolas_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaCruz_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SanNicolas_vs_SantaCruz

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanNicolas_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaRosa_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SanNicolas_vs_SantaRosa

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SanNicolas_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/Mainland_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SanNicolas_vs_Mainland

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaCatalina_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaCruz_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SantaCatalina_vs_SantaCruz

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaCatalina_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaRosa_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SantaCatalina_vs_SantaRosa

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaCatalina_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/Mainland_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SantaCatalina_vs_Mainland

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaCruz_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaRosa_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SantaCruz_vs_SantaRosa

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaCruz_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/Mainland_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SantaCruz_vs_Mainland

	vcftools --gzvcf /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/SantaRosa_refGenome.txt  --weir-fst-pop /scratch/users/elliea/jazlyn-ellie/grayfox_2023/metaData/Mainland_refGenome.txt --fst-window-size 10000 --fst-window-step 1000  --out /scratch/users/elliea/jazlyn-ellie/grayfox_2023/Analyses/FST/"$f"_filtered.renameChroms.ACgr25_DPgr165lt500_SantaRosa_vs_Mainland

	echo $SECONDS

done

sleep 180
