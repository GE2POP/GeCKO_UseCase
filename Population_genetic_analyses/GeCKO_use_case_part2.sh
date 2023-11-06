# -------------------------------------------------------------#
# 			GeCKO - Domestication wheat - use case			   #
# -------------------------------------------------------------#

set -euo 

# pre-requisite
# --------------------------------------------------------------

# needed files
### --------------------------------------------------------------

# four key output files produced by GeCKO for this use case and available on datagouv.fr
# please uncompress them
# 1. variant_calling_converted.vcf
# 2. 04__Genotype_Locus1_Sample_Locus2_Filtered.vcf (48 759 SNPs)
# 3. zones.bed
# 4. labels_groups.txt

# egglib  files generates by GeCKO_use_case_part1.sh

# --------------------------------------------------------------
### Estimate Diversity Reduction Index (DRI) 

Rscript DRI.R egglib_stats_10perPop.txt > DRI.tsv

#-------------------------------------------------------------------------
### DAPC (does hot handle missing data, imputation is done with beagle)

beagle gt=DurumWheat_20perPop.vcf out=DurumWheat_20perPop_imputed
gzip -d DurumWheat_20perPop_imputed.vcf.gz
Rscript PCA.R DurumWheat_20perPop_imputed.vcf labels_groups.txt

#-------------------------------------------------------------------------
### Fst scan 5A and 4B loci
Rscript Fst_WC.R
