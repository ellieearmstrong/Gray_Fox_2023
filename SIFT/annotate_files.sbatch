#!/bin/bash
#SBATCH --job-name=annotate-cfam4
#SBATCH --time=10:00:00
#SBATCH --mem=120000
#SBATCH --partition=normal,hns,owners
#SBATCH -n 12
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ellieearmstrong@gmail.com

module load gcc/12.1.0
module load java/18.0.2

java -jar SIFT4G_Annotator.jar -c -i /scratch/users/elliea/jazlyn-ellie/grayfox_2023/vcf/Canfam4_filtered.renameChroms.ACgr25_DPgr165lt500_sift.vcf -d /scratch/users/elliea/jazlyn-ellie/grayfox_2023/SIFT_db_builds/canfam4_test/Canfam4/ -r . -t 16

java -jar SIFT4G_Annotator.jar -c -i Canfam4_filtered.renameChroms.Mainland.ACgr61_DPgr165lt500_sift.vcf -d /scratch/users/elliea/jazlyn-ellie/grayfox_2023/SIFT_db_builds/canfam4_test/Canfam4/ -r . -t 16

java -jar SIFT4G_Annotator.jar -c -i arcticfox_filtered.renameChroms.ACgr25_DPgr165lt500.vcf -d /scratch/users/elliea/jazlyn-ellie/grayfox_2023/SIFT_db_builds/canfam4_test/arcticfox_SIFT -r . -t 16

java -jar SIFT4G_Annotator.jar -c -i arcticfox_filtered.renameChroms.Mainland.ACgr61_DPgr165lt500.vcf -d /scratch/users/elliea/jazlyn-ellie/grayfox_2023/SIFT_db_builds/canfam4_test/arcticfox_SIFT -r . -t 16

java -jar SIFT4G_Annotator.jar -c -i Canfam3.1_filtered.renameChroms.ACgr25_DPgr165lt500.vcf -d CanFam3.1.83 -r . -t 16

java -jar SIFT4G_Annotator.jar -c -i Canfam3.1_filtered.renameChroms.Mainland.ACgr61_DPgr165lt500.vcf -d CanFam3.1.83 -r . -t 16
