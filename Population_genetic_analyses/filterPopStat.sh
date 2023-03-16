# filterPopStat.sh

if(( $# < 3 )); then 
	echo " this program allows to filter vcf file by keeping loci having minimum number of individuals per population"
	echo " to keep loci with strictly more than 10 individuals DD and DC:"
	echo " $0 labels_groups.txt 04__Genotype_Locus1_Sample_Locus2_Filtered.vcf 10 DD DC > 04__Genotype_Locus1_Sample_Locus2_Filtered_DD10_DC10.vcf"
	echo " to keep loci with strictly more than 10 individuals in each population:"
	echo " $0 labels_groups.txt 04__Genotype_Locus1_Sample_Locus2_Filtered.vcf 10 > 04__Genotype_Locus1_Sample_Locus2_Filtered_DX10.vcf"
	exit 1; 
fi
pop_file=$1
vcf_file=$2

if(( $# == 3 )); then min=$3; pop1=""; pop2=""; nbPop=0; fi
if(( $# == 4 )); then min=$3; pop1=$4; nbPop=1; fi
if(( $# == 5 )); then min=$3; pop1=$4; pop2=$5; nbPop=2; fi 
awk -v pop1=$pop1 -v pop2=$pop2 -v min=$min -v nbPop=$nbPop '{ 
        if(NR==FNR){indiv2pop[$2]=$1; pops[$1]=$1} 
        else {
            if($1=="#CHROM") {for (c=10; c<=NF; c++){col2pop[c]=indiv2pop[$c]}; };
            if ( $0 !~ /^#/ ) {
                for (pop in pops){count[pop]=0}; 
                for (c=10; c<=NF; c++){if($c !~ /^\./)count[col2pop[c]]++}; 
                minAll=NF; minPop=""; for (pop in pops){if(count[pop] < minAll) {minAll=count[pop]; minPop=pop} };
                if (nbPop==0){ if (minAll > min){print $0}} # print minAll" "minPop}
                if (nbPop==1){ if (count[pop1]>min ){print $0}} #pop1" "count[pop1]"\t"pop2" "count[pop2]" "$0} }
                if (nbPop==2){ if (count[pop1]>min && count[pop2]>min){print $0 }} #{print pop1" "count[pop1]"\t"pop2" "count[pop2]" "$0} }
                #tot=0;for (pop in pops){printf "%s:%s_",pop,count[pop]; tot+=count[pop]};printf("tot_"tot"\t"$0"\n");
                #tot=0;for (pop in pops){printf "%s\t%s\t",pop,count[pop]; tot+=count[pop]};printf("\n")
            }
            else { print $0}
        }
    }' ${pop_file} ${vcf_file} 