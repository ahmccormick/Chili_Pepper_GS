#!/bin/bash
#SBATCH --job-name=Pepper_snp_count_12chr # Name for job
#SBATCH --partition=shared
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=900G
##SBATCH --core-spec=0 #Uncomment to allow jobs to request all cores on a node
#SBATCH -t 4320
#SBATCH -e /home/ahmccorm/%x-%A_%a.err # Standard err goes to this file
#SBATCH --mail-user ahmccorm@hawaii.edu # email wish to be notified at
#SBATCH --mail-type ALL #events for email. ALL will alert you of job beginning/completion/failure

# Load necessary modules
#module load plink/2.00-20231024-avx2
#module load tabix/0.2.6

module load bio/PLINK/2.00a3.7-foss-2022a

# Step 1: Perform LD pruning
plink2 --vcf Pepper_filtered_12chr.vcf.gz \
       --indep-pairwise 5000 500 0.2 \
       --out Pepper_filtered_and_ld_pruned \
       --allow-extra-chr

# Step 2: Generate a pruned VCF file based on the pruning results
plink2 --vcf Pepper_filtered_12chr.vcf.gz \
       --extract Pepper_filtered_and_ld_pruned.prune.in \
       --recode vcf bgz \
       --out Pepper_filtered_and_LD_pruned

# Step 3: Index the pruned VCF file
#tabix -p vcf Pepper_filtered_and_LD_pruned.vcf.gz
