# GeCKO - Domestication wheat - use case


GeCKO provides four workflows (DataCleaning, ReadsMapping, VariantsCalling, and VcfFiltering), which allow one to easily chain the different analyses needed to transform raw reads into filtered genotypes (VCF file). This gitHub project focuses on the use cases detailed in the paper to illustrate the features of GeCKO through a step-by-step analysis that traces the impact of domestication on allele diversity in durum wheat.

Below is a step by step guide to reproduce the use case analyses. For more information on the GeCKO project, please see the [GeCKO gitHub repository](https://github.com/GE2POP/GeCKO).

*This dataset is quite large, if your goal is just to test your GeCKO installation please use the small test dataset provided in the EXAMPLE folder of each workflow (e.g. [data cleaning example](https://github.com/GE2POP/GeCKO)).*

## context
We consider 120 accessions representing the three major subspecies involved in the transition steps described above (Figure 2); they represent the wild form T. turgidum dicoccoides (n=30, denoted DD), 
the first domesticated (solid rachis) form (T. turgidum dicoccum, n=30, DC), and the non-hulled cultivated form, (T. turgidum durum, n=60). The latter subspecies has been divided into two subgroups 
depending on whether the varieties originated in the pre- or post-Green Revolution period. The first group consists of lines derived from local varieties called "Landraces" (DP, n=30), and the second 
group consists of "elite" varieties registered in Europe after the Green Revolution (1970-1990; DE, n=30). 

## reproducing the use case analyse

### setting the environment

1. Installing GeCKO (getting GeCKO workflows) 
GeCKO workflows rely on snakemake, which ensures reproducible and scalable data analysis and make worflows installation straightforward since the only requirement is to have [snakemake](https://snakemake.readthedocs.io/en/stable/) and [conda](https://docs.conda.io/en/latest/) available on your environment (see [GeCKO install procedure](https://github.com/GE2POP/GeCKO#installation) for detailed information.

2. Downloading dataset
The input data needed to reproduce the use case analysis are provided in a [datagouv.fr dedicated repository](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/78MBZY)
Note that this dataset contains input files needed to reproduce GeCKO analyses as well as GeCKO expected output files. The latter can be used to check GeCKO ran as expected or to reproduce the population genetic analysis without the need of runing GeCKO workflows. 

3. clone the GeCKO use case github repository 
By doing so you will get i) GeCKO config files for reproducing the use case, and ii) scripts and conda environment files to reproduce population genetic analysis
```git clone git@github.com:GE2POP/GeCKO_UseCase.git```

### reproducing the GeCKO analysis (from raw reads to filtered VCF)
This should be done on 

1. Copy the fastq.gz files dowloaded from the dataverse in the adequate RAWDATA folder 

```
# for DEV_Cap009.2b
cp dataverse_files/DEV_Cap009.2b_R1.fastq.gz DEV_Cap009.2b/RAWDATA
cp dataverse_files/DEV_Cap009.2b_R2.fastq.gz DEV_Cap009.2b/RAWDATA
# for DEV_Cap0010.2
cp dataverse_files/DEV_Cap010.2_R1.fastq.gz  DEV_Cap010.2/RAWDATA
cp dataverse_files/DEV_Cap010.2_R2.fastq.gz  DEV_Cap010.2/RAWDATA
```

2. Modify config files
The information regarding the fastq files, read index etc. are, by default, retrieved from the config file CONFIG/config_WorkflowName.yml; all this is set with correct values to reproduce the use case analyses. 

The only file you need to adapt is the cluster_config_WorkflowName.yml file used by default to provide specific cluster information (e.g. job queue names) related to this workflow:
- DEV_Cap009.2b/CONFIG/cluster_config_DataCleaning.yml
- DEV_Cap0010.2/CONFIG/cluster_config_DataCleaning.yml
- DEV_Cap009_and_Cap010/CONFIG/cluster_config_ReadMapping.yml
- DEV_Cap009_and_Cap010/CONFIG/cluster_config_VariantCalling.yml
- DEV_Cap009_and_Cap010/CONFIG/cluster_config_VcfFiltering.yml

3. launch Gecko workflow
Assuming you clone the GeCKO repository on your home directory, you shloud have a ~/GECKO folder containing GeCKO workflows and all analyses can be launch using (each command need to be finished before launching the next one):
```
# DataCleaning
cd DEV_Cap009.2b
runGeCKO.sh --workflow DataCleaning --workflow-path ~/GECKO --jobs 50 --job-scheduler SLURM

cd ../DEV_Cap010.2
runGeCKO.sh --workflow DataCleaning --workflow-path ~/GECKO --jobs 50 --job-scheduler SLURM

# ReadMapping
cd ../DEV_Cap009_and_Cap010
runGeCKO.sh --workflow ReadMapping --workflow-path ~/GECKO --jobs 50 --job-scheduler SLURM

# VariantCalling
runGeCKO.sh --workflow VariantCalling --workflow-path ~/GECKO --jobs 50 --job-scheduler SLURM

# VcfFiltering
runGeCKO.sh --workflow VcfFiltering --workflow-path ~/GECKO --jobs 50 --job-scheduler SLURM
```

Note that this launch launches the commands on the master node, it is advise to launch the command as a job on our cluster, for the first datacleaing workflow, this look likes:
``` 
# DataCleaning
cd DEV_Cap009.2b
sbatch --partition=agap_normal --mem=10G --job-name=DC --wrap="runGeCKO.sh --workflow DataCleaning --workflow-path ~/GECKO --jobs 50 --job-scheduler SLURM"
```

### reproducing population genetic analysis

1. copy GeCKO output files in GeCKO_UseCase/Population_genetic_analyses/
The files to be copied can be taken from GeCKO output folders (if you have reproduced the GeCKO analysis) or from the dataverse_file folder:
```
cd GeCKO_UseCase/Population_genetic_analyses/

cp ~/dataverse_files/variants_calling_converted.vcf.gz . 
# DEV_Cap009_and_Cap010/WORKFLOWS_OUTPUTS/VARIANTS_CALLING/GENOTYPE_GVCFS/

cp ~/dataverse_files/04__Genotype_Locus1_Sample_Locus2_Filtered.vcf.gz .
# DEV_Cap009_and_Cap010/WORKFLOWS_OUTPUTS/VCF_FILTERING/

cp ~/dataverse_files/zones.bed .
#

cp ~/dataverse_files/labels_groups.txt .

gzip -d variants_calling_converted.vcf.gz
gzip -d 04__Genotype_Locus1_Sample_Locus2_Filtered.vcf.gz
```

2. create conda environments
```
# for egglib
conda env create --file conda_UC_GECKO.yml
# for R analyses
conda env create --file conda_UC_GECKO_R.yml
chmod +x GeCKO_use_case_*.sh
```

3. compute population statistics
```
conda activate UC-GECKO
./GeCKO_use_case_part1.sh
conda deactivate
```

4. Compute molecular diversity and Generate plots for Fst-scan and DAPC analysis
```
conda activate UC-GECKO
./GeCKO_use_case_part1.sh
conda deactivate
```
Molecular diversity of each population is provided in DRI.tsv. DRI values can be obtained by taking the ratios between the molecular diversity computed for the two populations to compare. 

###
Reference: 

download Triticum dicoccoides reference:
https://www.ncbi.nlm.nih.gov/assembly/GCA_002162155.2
deletion of non-positional scaffolds and simplification of the sequence names of the 14 chromosomes (example: ">chr1A")

> Triticum_dicoccoides.WEWSeq_v.2.0.dna_20210617.fasta






