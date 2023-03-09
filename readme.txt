# -------------------------------------------------------------#
#          GeCKO - Domestication wheat - use case	       #
# -------------------------------------------------------------#


GeCKO provides four workflows (DataCleaning, ReadsMapping, VariantsCalling, and VcfFiltering), which allow one to easily chain the different analyses needed to transform raw reads into filtered genotypes (VCF file).
To reduce the learning curve, GeCKO workflows use a consistent naming scheme for output files and folders and parameters. Additionally, GeCKO workflows create thorough reports that make it simple to track the progress 
of the analysis and facilitate quality control checks. GeCKO was designed with capture-basedgenotyping in mind so that one can easily chain GeCKO workflows to seamlessly go from raw multiplexed reads, produced by target 
enrichment captures, to a ready-to-use genotype matrix. Each workflow can also be used independently of the others, allowing for a much broader range of applications.

We illustrate the features of GeCKO through a step-by-step analysis that traces the impact of domestication on allele diversity in durum wheat.
We consider 120 accessions representing the three major subspecies involved in the transition steps described above (Figure 2); they represent the wild form T. turgidum dicoccoides (n=30, denoted DD), 
the first domesticated (solid rachis) form (T. turgidum dicoccum, n=30, DC), and the non-hulled cultivated form, (T. turgidum durum, n=60). The latter subspecies has been divided into two subgroups 
depending on whether the varieties originated in the pre- or post-Green Revolution period. The first group consists of lines derived from local varieties called "Landraces" (DP, n=30), and the second 
group consists of "elite" varieties registered in Europe after the Green Revolution (1970-1990; DE, n=30). 

depot dataverse: XXX

# -------------------------------------------------------------------------
### The Diversity Reduction Index (DRI) 

1/ Selection of SNPs according to the number of individuals genotyped per population.
*tools: filterPopStat.sh 
*parameters/filters: Storage of SNPs for which genotyping is available for at least 10 individuals in each of the 4 populations.
*input file(s): 04__Genotype_Locus1_Sample_Locus2_Filtered.vcf ( 48 759 SNPs )
*output file(s): 04__Genotype_Locus1_Sample_Locus2_Filtered_10perpop.vcf (22 225 SNPs)

2/ Estimation of population genetics parameters
*tools: egglib_vcf_stats_lseff.py (EggLib)
module load egglib/3.1
python egglib_vcf_stats_lseff.py --bed_file /home/ardissonm/scratch/DEV/ANALYSES_2022_WORKFLOW/DEV_Cap009_and_Cap010/WORKFLOWS_OUTPUTS/READS_MAPPING/MAPPING_ZAVITAN_REMOVE_DUPLICATES/EXTRACTED_BAMS/REFERENCE_zones/zones.bed --vcf_file_raw /home/ardissonm/scratch/DEV/ANALYSES_2022_WORKFLOW/DEV_Cap009_and_Cap010/WORKFLOWS_OUTPUTS/VARIANTS_CALLING/ZAVITAN_ZONES_119_100_100_REMOVE_DUP/GENOTYPE_GVCFS/variants_calling_converted.vcf --vcf_file ../FILTERS_GENO_POP/04__Genotype_Locus1_Sample_Locus2_Filtered_10perpop.vcf --labels_file labels_groups.txt
*parameters/filters: None
*input file(s): 04__Genotype_Locus1_Sample_Locus2_Filtered_10perpop.vcf
*output file(s): egglib_stats_10DD_10DC_10DP_10DE.txt, egglib_stats_pairwise_10DD_10DC_10DP_10DE.txt 

3/ Calculation of DRI values (The Diversity Reduction Index (DRI) was then calculated for each domestication transition as the ratio between the allelic diversity of the two populations framing this transition; the higher this value, the stronger the loss of diversity)
*tools: DRI.rmd
*parameters/filters: None
*input file(s): egglib_stats_10DD_10DC_10DP_10DE.txt
*output file(s): DRI values (exemple: DRI DD>DC = DD_pi_site/DC_pi_site)

-------------------------------------------------------------------------
### DAPC

1/ Selection of SNPs according to the number of individuals genotyped per population.
*tools: filterPopStat.sh 
*parameters/filters: Storage of SNPs for which genotyping is available for at least 20 individuals in each of the 4 populations.
*input file(s): 04__Genotype_Locus1_Sample_Locus2_Filtered.vcf ( 48 759 SNPs )
*output file(s): 04__Genotype_Locus1_Sample_Locus2_Filtered_20perpop.vcf (6 538 SNPs)

2/ Data imputation
*tools: beagle.22Jul22.46e.jar 
wget https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar
sbatch --partition=agap_normal --job-name=IMPUT --wrap="java -jar beagle.22Jul22.46e.jar gt=04__Genotype_Locus1_Sample_Locus2_Filtered_20perpop.vcf out=04__Genotype_Locus1_Sample_Locus2_Filtered_20perpop_imputed.vcf"
*parameters/filters: default
*input file(s): 04__Genotype_Locus1_Sample_Locus2_Filtered_20perpop.vcf
*output file(s): 04__Genotype_Locus1_Sample_Locus2_Filtered_20perpop_imputed.vcf

3/ DAPC
*tools: DAPC.rmd
*parameters/filters: None
*input file(s): 04__Genotype_Locus1_Sample_Locus2_Filtered_20perpop_imputed.vcf
*output file(s): DAPC.pdf

-------------------------------------------------------------------------
### Fst scan

1/ Selection of SNPs according to the number of individuals genotyped per population.
*tools: filterPopStat.sh 
*parameters/filters: Storage of SNPs for which genotyping is available for at least 10 individuals for each of the two transition populations studied
*input file(s): 04__Genotype_Locus1_Sample_Locus2_Filtered.vcf ( 48 759 SNPs )
*output file(s): 04__Genotype_Locus1_Sample_Locus2_Filtered_10DD_10DC.vcf (23520 SNPs), 04__Genotype_Locus1_Sample_Locus2_Filtered_10DC_10DP.vcf (25163 SNPs), 04__Genotype_Locus1_Sample_Locus2_Filtered_10DP_10DE.vcf (25317 SNPs)

2/ Estimation of population genetics parameters
*tools: egglib_vcf_stats_lseff.py (EggLib)
exemple for DD>DC transition:
module load egglib/3.1
python egglib_vcf_stats_lseff.py --bed_file /home/ardissonm/scratch/DEV/ANALYSES_2022_WORKFLOW/DEV_Cap009_and_Cap010/WORKFLOWS_OUTPUTS/READS_MAPPING/MAPPING_ZAVITAN_REMOVE_DUPLICATES/EXTRACTED_BAMS/REFERENCE_zones/zones.bed --vcf_file_raw /home/ardissonm/scratch/DEV/ANALYSES_2022_WORKFLOW/DEV_Cap009_and_Cap010/WORKFLOWS_OUTPUTS/VARIANTS_CALLING/ZAVITAN_ZONES_119_100_100_REMOVE_DUP/GENOTYPE_GVCFS/variants_calling_converted.vcf --vcf_file ../FILTERS_GENO_POP/04__Genotype_Locus1_Sample_Locus2_Filtered_10DD_10DC.vcf --labels_file labels_groups.txt
*parameters/filters: None
*input file(s): 04__Genotype_Locus1_Sample_Locus2_Filtered_10DD_10DC.vcf (23520 SNPs), 04__Genotype_Locus1_Sample_Locus2_Filtered_10DC_10DP.vcf (25163 SNPs), 04__Genotype_Locus1_Sample_Locus2_Filtered_10DP_10DE.vcf (25317 SNPs)
*output file(s): egglib_stats_10DD_10DC.txt, egglib_stats_10DC_10DP.txt, egglib_stats_10DP_10DE.txt

3/ Fst scan 
*tools: Fst_WC.rmd
*parameters/filters: None
*input file(s): egglib_stats_10DD_10DC.txt, egglib_stats_10DC_10DP.txt, egglib_stats_10DP_10DE.txt
*output file(s): Fst_WC.pdf
