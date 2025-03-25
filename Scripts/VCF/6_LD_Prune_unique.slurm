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
module load bio/PLINK/2.00a3.7-foss-2022a


# Step 1: Assign unique variant IDs
plink2 --vcf Pepper_filtered_12chr.vcf.gz \
       --set-all-var-ids '@:#' \
       --make-bed \
       --out Pepper_filtered_unique_ids \
       --allow-extra-chr

# Step 2: Perform LD pruning
plink2 --bfile Pepper_filtered_unique_ids \
       --indep-pairwise 5000 500 0.2 \
       --out Pepper_filtered_and_ld_pruned

# Step 3: Generate a pruned VCF file based on the pruning results
plink2 --bfile Pepper_filtered_unique_ids \
       --extract Pepper_filtered_and_ld_pruned.prune.in \
       --recode vcf bgz \
       --out Pepper_filtered_and_LD_pruned
