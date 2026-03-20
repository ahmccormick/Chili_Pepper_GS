####################################################################################################
# ANNOTATE TABLE S5, S6, S7 with GFF3 gene annotations
# Zhangshugang reference (GCA_030867735.1)
####################################################################################################

library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(dplyr)

#############################################
# 0) Paths
#############################################

gff_path <- "/Users/annamccormick/R/World_Veg_Project/Reference_sequence_cornell/Zhangshugang.gff3"
outdir   <- "~/R/WorldVeg_Capsicum/STI_GWAS/GWAS_tables/"

#############################################
# 1) Load GFF3
#############################################

gff   <- import(gff_path)
genes <- gff[gff$type == "gene"]
mrna  <- gff[gff$type == "mRNA"]
exons <- gff[gff$type == "exon"]
cds   <- gff[gff$type == "CDS"]

#############################################
# 2) Annotation function
#    Expects columns: Chr (integer), Pos (integer)
#############################################

collapse_cl <- function(x) {
  vapply(x, function(z) paste(as.character(z), collapse = "; "), character(1))
}

annotate_snps <- function(snps_df, genes, mrna, exons, cds) {
  
  # format chromosome to match GFF seqnames e.g. Chr01, Chr02 ...
  snps_df <- snps_df %>%
    mutate(Chromosome_gff = sprintf("Chr%02d", as.integer(Chr)))
  
  snp_gr <- GRanges(
    seqnames = snps_df$Chromosome_gff,
    ranges   = IRanges(snps_df$Pos, snps_df$Pos)
  )
  
  # overlaps
  gene_hits <- findOverlaps(snp_gr, genes, ignore.strand = TRUE)
  mrna_hits <- findOverlaps(snp_gr, mrna,  ignore.strand = TRUE)
  exon_hits <- findOverlaps(snp_gr, exons, ignore.strand = TRUE)
  cds_hits  <- findOverlaps(snp_gr, cds,   ignore.strand = TRUE)
  
  # initialise annotation columns
  snps_df$gene_id   <- NA_character_
  snps_df$gene_name <- NA_character_
  snps_df$mrna_id   <- NA_character_
  snps_df$mrna_note <- NA_character_
  snps_df$in_exon   <- FALSE
  snps_df$in_cds    <- FALSE
  
  # fill gene fields
  qg <- queryHits(gene_hits); sg <- subjectHits(gene_hits)
  snps_df$gene_id[qg]   <- mcols(genes)$ID[sg]
  snps_df$gene_name[qg] <- mcols(genes)$Name[sg]
  
  # fill mRNA functional annotation
  qm <- queryHits(mrna_hits); sm <- subjectHits(mrna_hits)
  snps_df$mrna_id[qm]   <- mcols(mrna)$ID[sm]
  snps_df$mrna_note[qm] <- collapse_cl(mcols(mrna)$Note[sm])
  
  # exon / CDS flags
  snps_df$in_exon[unique(queryHits(exon_hits))] <- TRUE
  snps_df$in_cds[unique(queryHits(cds_hits))]   <- TRUE
  
  # classification
  snps_df$annotation <- "Intergenic"
  snps_df$annotation[!is.na(snps_df$gene_id)] <- "Genic"
  snps_df$annotation[snps_df$in_exon]         <- "Exonic"
  snps_df$annotation[snps_df$in_cds]          <- "CDS"
  
  # primary functional label
  snps_df$functional_note <- snps_df$mrna_note
  
  # nearest gene + distance (useful for intergenic SNPs)
  nearest_gene_idx <- nearest(snp_gr, genes, ignore.strand = TRUE)
  snps_df$nearest_gene_id          <- mcols(genes)$ID[nearest_gene_idx]
  snps_df$dist_to_nearest_gene_bp  <- distance(snp_gr, genes[nearest_gene_idx])
  
  # nearest mRNA note
  nearest_mrna_idx      <- nearest(snp_gr, mrna, ignore.strand = TRUE)
  snps_df$nearest_mrna_note <- collapse_cl(mcols(mrna)$Note[nearest_mrna_idx])
  
  # drop the helper column
  snps_df <- snps_df %>% select(-Chromosome_gff)
  
  snps_df
}

#############################################
# 3) Read tables
#############################################

s5 <- read.csv(file.path(outdir, "TableS5_significant_SNPs_2plus_models.csv"),
               stringsAsFactors = FALSE)
s6 <- read.csv(file.path(outdir, "TableS6_multitrait_SNPs_3plus_traits.csv"),
               stringsAsFactors = FALSE)
s7 <- read.csv(file.path(outdir, "TableS7_multitrait_SNPs_3plus_models.csv"),
               stringsAsFactors = FALSE)

# quick check columns are present
stopifnot("Chr" %in% names(s5), "Pos" %in% names(s5))
stopifnot("Chr" %in% names(s6), "Pos" %in% names(s6))
stopifnot("Chr" %in% names(s7), "Pos" %in% names(s7))

#############################################
# 4) Annotate
#############################################

message("Annotating Table S5 (", nrow(s5), " rows)...")
s5_annot <- annotate_snps(s5, genes, mrna, exons, cds)

message("Annotating Table S6 (", nrow(s6), " rows)...")
s6_annot <- annotate_snps(s6, genes, mrna, exons, cds)

message("Annotating Table S7 (", nrow(s7), " rows)...")
s7_annot <- annotate_snps(s7, genes, mrna, exons, cds)

#############################################
# 5) Quick summary
#############################################

for (nm in c("s5_annot", "s6_annot", "s7_annot")) {
  df <- get(nm)
  message("\n--- ", nm, " annotation breakdown ---")
  print(df %>% count(annotation, sort = TRUE))
  message("Intergenic SNPs with nearest gene note available: ",
          sum(!is.na(df$nearest_mrna_note[df$annotation == "Intergenic"])))
}

#############################################
# 6) Write output
#############################################

write.csv(s5_annot,
          file.path(outdir, "TableS5_significant_SNPs_2plus_models_annotated.csv"),
          row.names = FALSE)

write.csv(s6_annot,
          file.path(outdir, "TableS6_multitrait_SNPs_3plus_traits_annotated.csv"),
          row.names = FALSE)

write.csv(s7_annot,
          file.path(outdir, "TableS7_multitrait_SNPs_3plus_models_annotated.csv"),
          row.names = FALSE)

