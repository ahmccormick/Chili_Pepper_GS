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
module --ignore_cache load bio/PLINK/2.00a3.7-foss-2022a
module load bio/PLINK/2.00a3.7-foss-2022a
module load bio/BCFtools/1.15.1-GCC-11.3.0

#plink2 --vcf Pepper_filtered.vcf.gz --chr-list

#plink2 --vcf Pepper_filtered.vcf.gz --freq --out Pepper_chromosomes
plink2 --vcf Pepper_filtered.vcf.gz --freq --out Pepper_chromosomes --allow-extra-chr

bcftools query -f '%CHROM\n' Pepper_filtered.vcf.gz | sort | uniq > chromosome_list.txt
