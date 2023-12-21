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
# 3. auto_zones.bed
# 4. labels_groups.txt

# needed envirnoment
### --------------------------------------------------------------
# the two conda environment should be create (once) using
# conda env create --file conda_UC_GECKO.yml
# conda env create --file conda_UC_GECKO_R.yml

cp variant_calling_converted.vcf  DurumWheat_GeCKO_unfiltered_onZavitan.vcf
cp 04__Genotype_Locus1_Sample_Locus2_Filtered.vcf DurumWheat_GeCKO_filtered_onZavitan.vcf
cp auto_zones.bed DW_tarteged_zones.bed


# 1. create filtered dataset based on sample size
# --------------------------------------------------------------


./filterPopStat.sh labels_groups.txt DurumWheat_GeCKO_filtered_onZavitan.vcf 9 DD DC > DurumWheat_10DD_10DC.vcf 
# 23520 SNPs

./filterPopStat.sh labels_groups.txt DurumWheat_GeCKO_filtered_onZavitan.vcf 9 DC DP > DurumWheat_10DC_10DP.vcf 
# 25163 SNPs

./filterPopStat.sh labels_groups.txt DurumWheat_GeCKO_filtered_onZavitan.vcf 9 DP DE > DurumWheat_10DP_10DE.vcf
#25163 SNPs

./filterPopStat.sh labels_groups.txt DurumWheat_GeCKO_filtered_onZavitan.vcf 9 > DurumWheat_10perPop.vcf
#22 225 SNPs

./filterPopStat.sh labels_groups.txt DurumWheat_GeCKO_filtered_onZavitan.vcf 19 > DurumWheat_20perPop.vcf
#6 538 SNPs

# --------------------------------------------------------------
### 2. Compute population genetic stats using egglib 

python egglib_vcf_stats_lseff.py --bed_file DW_tarteged_zones.bed --vcf_file_raw DurumWheat_GeCKO_unfiltered_onZavitan.vcf --vcf_file DurumWheat_10perPop.vcf --labels_file labels_groups.txt
mv egglib_stats.txt egglib_stats_10perPop.txt; mv egglib_stats_pairwise.txt egglib_stats_pairwise_10perPop.txt ; rm egglib_outliers.txt

python egglib_vcf_stats_lseff.py --bed_file DW_tarteged_zones.bed --vcf_file_raw DurumWheat_GeCKO_unfiltered_onZavitan.vcf --vcf_file DurumWheat_10DD_10DC.vcf --labels_file labels_groups.txt
mv egglib_stats.txt egglib_stats_10DD_10DC.txt; rm egglib_stats_pairwise.txt egglib_outliers.txt

python egglib_vcf_stats_lseff.py --bed_file DW_tarteged_zones.bed --vcf_file_raw DurumWheat_GeCKO_unfiltered_onZavitan.vcf --vcf_file DurumWheat_10DC_10DP.vcf --labels_file labels_groups.txt
mv egglib_stats.txt egglib_stats_10DC_10DP.txt; rm egglib_stats_pairwise.txt egglib_outliers.txt

python egglib_vcf_stats_lseff.py --bed_file DW_tarteged_zones.bed --vcf_file_raw DurumWheat_GeCKO_unfiltered_onZavitan.vcf --vcf_file DurumWheat_10DP_10DE.vcf --labels_file labels_groups.txt
mv egglib_stats.txt egglib_stats_10DP_10DE.txt; rm egglib_stats_pairwise.txt egglib_outliers.txt

