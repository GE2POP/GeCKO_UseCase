# download reference genome
echo "wheat genome is quite large (> 10Gb) each step is quite long, please be patient"
echo ""
echo "dowloading the wheat genome" 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/162/155/GCF_002162155.1_WEW_v2.0/GCF_002162155.1_WEW_v2.0_genomic.fna.gz
echo "unzipping the wheat genome"
gzip -d GCF_002162155.1_WEW_v2.0_genomic.fna.gz

echo "simplify chromosome names"
cat GCF_002162155.1_WEW_v2.0_genomic.fna | sed 's/ /_/g' | sed 's/,//g' | sed 's/.*chromosome/>chr/g' | sed 's/_WEW_v2.0_whole_genome_shotgun_sequence//g' | sed 's/_//g' > Triticum_dicoccoides.WEWSeq_v.2.0.dna_with_scaffolds_20210617_names.fasta
echo "filtering the reference to keep only chromosomes and remove others contigs/scaffolds"
samtools faidx Triticum_dicoccoides.WEWSeq_v.2.0.dna_with_scaffolds_20210617.fasta
sed -n '/chr/ p' Triticum_dicoccoides.WEWSeq_v.2.0.dna_with_scaffolds_20210617_names.fasta.fai |cut -f 1 >chr.lst
faidx -r chr.lst Triticum_dicoccoides.WEWSeq_v.2.0.dna_with_scaffolds_20210617_names.fasta >Triticum_dicoccoides.WEWSeq_v.2.0.dna_20210617.fasta
rm  Triticum_dicoccoides.WEWSeq_v.2.0.dna_with_scaffolds_20210617_names.fasta chr.lst
echo "done, the wheat reference is ready to be used with GeCKO"
