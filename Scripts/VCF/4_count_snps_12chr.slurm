#!/bin/bash
#SBATCH --job-name=Pepper_snp_count_12chr # Name for job
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
module load bio/BCFtools/1.15.1-GCC-11.3.0

bcftools view -H Pepper_filtered_12chr.vcf.gz | wc -l > SNP_count.txt
