#!/bin/bash
#SBATCH --job-name=Pepper_chr_count # Name for job
#SBATCH --partition=shared
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=100G
##SBATCH --core-spec=0 #Uncomment to allow jobs to request all cores on a node
#SBATCH -t 4320
#SBATCH -e /home/ahmccorm/%x-%A_%a.err # Standard err goes to this file
#SBATCH --mail-user ahmccorm@hawaii.edu # email wish to be notified at
#SBATCH --mail-type ALL #events for email. ALL will alert you of job beginning/completion/failure


# Load necessary modules
#module load plink/2.00-20231024-avx2
#module --ignore_cache load bio/PLINK/2.00a3.7-foss-2022a
#module load bio/PLINK/2.00a3.7-foss-2022a
module load bio/BCFtools/1.15.1-GCC-11.3.0

#bcftools view -R keep_chromosomes.txt Pepper_filtered.vcf.gz -Oz -o Pepper_filtered_12chr.vcf.gz
#bcftools index Pepper_filtered_12chr.vcf.gz

bcftools view -r Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12 Pepper_filtered.vcf.gz -Oz -o Pepper_filtered_12chr.vcf.gz
bcftools index Pepper_filtered_12chr.vcf.gz
