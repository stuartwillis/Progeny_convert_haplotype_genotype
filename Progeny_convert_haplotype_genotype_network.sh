#!/bin/bash

#take strings from one file, search for them in another file, output matches to file with output_prefix.txt

if [ -z "$4" ]
then
echo -e "\n"
echo "Not enough arguments."
echo "Correct usage: SCRIPT.sh [PLINK || HAPLOVIEW || HAPLOTYPE || NETWORK] genotype-file marker-file output_prefix [optional: individual-filter-proportion]"
echo "CHOOSE one mode: reformat from Progeny to PLINK (Plink1.9 ped/map) or HAPLOVIEW, HAPLOTYPE using Shapeit (output is haps/sample combined with your phenotype data), OR haplotype then make a NETWORK"
echo "MAKE SURE YOUR INPUT FILES ARE UNIX-DELIMTED WITH ONE(1) HANGING LINE!!!"
echo "Optional proportion to filter individuals by missing data must be between 0 and 1 (default is 0.2)"
echo -e "\n\n"
echo "Example usage: bash Progeny_haplotype_genotype.sh REFORMAT progeny_genotypes.txt marker_file.txt Omykiss_2018 0.1"
echo -e "\n"
exit 1
#else
fi

#check filtering value is between 0.1 and 1

if [ "$1" == "PLINK" ]; then
	echo -e "\n"
    echo "I will reformat the files into "$1" format."
	echo -e "\n"
elif [ "$1" == "HAPLOVIEW" ]; then
	echo -e "\n"
    echo "I will reformat the files into "$1" format."
	echo -e "\n"
elif [ "$1" == "HAPLOTYPE" ]; then
	echo -e "\n"
    echo "I will reformat the files and "$1" them for you too. I'm such a nice script."
	echo -e "\n"
elif [ "$1" == "NETWORK" ]; then
	echo -e "\n"
    echo "I will reformat the files, haplotype tyem, and make a haplotype "$1" for you too. I'm such a nice script."
	echo -e "\n"
else
	echo -e "\n"
    echo "Incorrect Usage"
	echo "Correct usage: SCRIPT.sh [PLINK || HAPLOVIEW || HAPLOTYPE || NETWORK] genotype-file marker-file output_prefix [optional: filter-proportion]"
	echo -e "\n"
	echo "Example usage: bash Progeny_haplotype_genotype.sh REFORMAT progeny_genotypes.txt marker_file.txt Omykiss_2018 0.9"
	echo -e "\n"
    exit 1
fi

if [ -z "$5" ]
then
	FILT=0.2
else
	filtint=`echo "scale=0; $5*100" | bc | sed 's/\.[0-9]*//g'`
	if [ "$filtint" -gt 100 ] 
	then 
		echo "filter proportion must be between 0 and 1"
		exit 1
	elif [ "$filtint" -lt 0 ] 
	then
		echo "filter proportion must be between 0.1 and 1"
		exit 1
	fi
FILT=$5
fi


#QUERYb=( `cut -f1 $1 `)
LENb=( `wc -l $3 `)
LENb=$(($LENb - 1))

echo "There appear to be "$LENb" markers to haplotype (minus Sex marker if included)."

rm markers.txt &> /dev/null 
sed 's/-/\./g' $3 > markers.txt

LENa=( `wc -l $2 `)
LENa=$(($LENa - 1))

echo "There appear to be "$LENa" individuals to haplotype."

rm genotypes.txt &> /dev/null 
head -n1 $2 | sed 's/-/\./g' > genotypes.txt
grep -v "Individual" $2 >> genotypes.txt

echo "Expected output files will be "$4".map, "$4".ped, and "$4".hap (if applicable). Is this ok? [yes or no]"

read Indcorrect

if [ "$Indcorrect" == "no" ]; then
        echo "Please rename files or choose a different output prefix."
        echo "Please also remember to include a hanging line on both genotype and marker files."
        exit 1
elif [ "$Indcorrect" == "yes" ]; then
            echo "Proceeding with file conversion"
else
        echo "Incorrect Input"
        exit 1
fi

####write R code to a file (reformat_Progeny_PLINK.R) to re-format Progeny file to PLINK .ped and .map


echo "######------------read in export from Progeny and marker file with positions, export in PLINK format (.ped and .map)" > reformat_Progeny_PLINK.R

echo "#extract genotypes using allele names which must match the \"assay\" field with \"-\" translated to \".\"" >> reformat_Progeny_PLINK.R
echo "#extract genetic sex as \"X\" or \"Y\" in the appropriate field...which field? Code unknowns appropriately!!" >> reformat_Progeny_PLINK.R
echo "progeny_genotypes=read.delim(\"genotypes.txt\",colClasses = c(\"factor\"))" >> reformat_Progeny_PLINK.R

echo "#read in markers and marker classes (headings must be \"Assay\",\"Chromosome\",\"Position\",\"Presumed Type\"; all others will be ignored; order does not matter)" >> reformat_Progeny_PLINK.R
echo "#this includes real or dummy chromosomes and positions, but only a SINGLE chromosome will be analyzed" >> reformat_Progeny_PLINK.R
echo "#optionally include \"Sex Marker\" (others are ignored) in the column \"Presumed Type\" for sex assignment; otherwise all will be marked unknown" >> reformat_Progeny_PLINK.R
echo "#in this case include the sex marker FIRST in the marker file (chromosome will be ignored)" >> reformat_Progeny_PLINK.R
echo "#allele names are made from assay assuming the marker names are a simple extension of that (i.e. Assay-A1 or Assay.A1)" >> reformat_Progeny_PLINK.R
echo "markers<-read.delim(\"markers.txt\",colClasses = c(\"factor\"))" >> reformat_Progeny_PLINK.R
echo "hap_markers<-droplevels(markers[!markers\$Presumed.Type==\"Sex Marker\",c(\"Assay\",\"Chromosome\",\"Position\")])" >> reformat_Progeny_PLINK.R
echo "hap_markers<-hap_markers[hap_markers\$Chromosome==levels(hap_markers\$Chromosome)[1],]" >> reformat_Progeny_PLINK.R
echo "hap_markers<-hap_markers[order(hap_markers\$Position),]" >> reformat_Progeny_PLINK.R

echo "blanks<-as.integer(rep(\"0\",times=nrow(progeny_genotypes)))" >> reformat_Progeny_PLINK.R

echo "if (\"Sex Marker\" %in% levels(markers\$Presumed.Type)) {" >> reformat_Progeny_PLINK.R
echo "  sex_marker<-markers[markers\$Presumed.Type==\"Sex Marker\",\"Assay\"]" >> reformat_Progeny_PLINK.R
echo "  sex<-as.character(progeny_genotypes[,paste(sex_marker,\".A2\",sep=\"\")])" >> reformat_Progeny_PLINK.R
echo "  sex<-replace(sex,sex==\"X\",\"2\")" >> reformat_Progeny_PLINK.R
echo "  sex<-replace(sex,sex==\"Y\",\"1\")" >> reformat_Progeny_PLINK.R
echo "  sex<-as.integer(sex)" >> reformat_Progeny_PLINK.R
echo "} else {" >> reformat_Progeny_PLINK.R
echo "  sex<-blanks" >> reformat_Progeny_PLINK.R
echo "}" >> reformat_Progeny_PLINK.R

echo "marker_cols<-NULL" >> reformat_Progeny_PLINK.R
echo "for (m in 1:nrow(hap_markers)) {" >> reformat_Progeny_PLINK.R
echo "  marker_cols<-c(marker_cols,paste(hap_markers\$Assay[m],\".A1\",sep=\"\"))" >> reformat_Progeny_PLINK.R
echo "  marker_cols<-c(marker_cols,paste(hap_markers\$Assay[m],\".A2\",sep=\"\"))" >> reformat_Progeny_PLINK.R
echo "}" >> reformat_Progeny_PLINK.R
echo "PLINK_genotypes<-cbind(as.character(progeny_genotypes\$Individual.Name),as.character(progeny_genotypes\$Individual.Name),blanks,blanks,sex,blanks,progeny_genotypes[,marker_cols])" >> reformat_Progeny_PLINK.R
echo "PLINK_map<-cbind(as.integer(rep(\"1\",times=nrow(hap_markers))),as.character(hap_markers\$Assay),as.integer(rep(\"0\",times=nrow(hap_markers))),as.character(hap_markers\$Position))" >> reformat_Progeny_PLINK.R

echo "#filter individuals with missing data above 0.X" >> reformat_Progeny_PLINK.R
echo "prop_miss<-NULL" >> reformat_Progeny_PLINK.R
echo "for (n in 1:nrow(PLINK_genotypes)) {" >> reformat_Progeny_PLINK.R
echo "  freqs<-as.data.frame(table( unlist(PLINK_genotypes[n,marker_cols]) ) )" >> reformat_Progeny_PLINK.R
echo "  if (\"0\" %in% freqs\$Var1 ) {" >> reformat_Progeny_PLINK.R
echo "      miss<-(freqs[freqs\$Var1==\"0\",\"Freq\"]/length(marker_cols))" >> reformat_Progeny_PLINK.R
echo "    } else {" >> reformat_Progeny_PLINK.R
echo "      miss<-0" >> reformat_Progeny_PLINK.R
echo "    }" >> reformat_Progeny_PLINK.R
echo "    prop_miss<-c(prop_miss,miss)" >> reformat_Progeny_PLINK.R
echo "}" >> reformat_Progeny_PLINK.R

echo "PLINK_genotypes_filtered<-PLINK_genotypes[prop_miss<"$FILT",]" >> reformat_Progeny_PLINK.R

echo "write.table(PLINK_genotypes_filtered,file=\"genotypes.ped\",append=FALSE,quote=FALSE,sep=\"\t\",col.names=F,row.names=F)" >> reformat_Progeny_PLINK.R
echo "write.table(PLINK_map,file=\"genotypes.map\",append=FALSE,quote=FALSE,sep=\"\t\",col.names=F,row.names=F)" >> reformat_Progeny_PLINK.R
echo "quit(save = \"no\", status = 0, runLast = FALSE);" >> reformat_Progeny_PLINK.R

###execute reformat R file (reformat_Progeny_PLINK.R)
R CMD BATCH reformat_Progeny_PLINK.R

cp genotypes.ped $4.ped
cp genotypes.map $4.map

if [ "$1" == "PLINK" ]
then
	#rm markers.txt &> /dev/null 
	#rm genotypes.txt &> /dev/null 
	#rm reformat_Progeny_PLINK.R &> /dev/null 
	#rm genotypes.ped &> /dev/null 
	#rm genotypes.map &> /dev/null 
	echo "PLINK files created. Haplotyping option not chosen."
	exit 1
elif [ "$1" == "HAPLOVIEW" ]
then
	cut -f2,4 $4.map > temp.txt
	mv temp.txt $4.map
	exit 1
fi

echo "PLINK files created. Proceeding with haplotyping."

SHAPE=( ` which shapeit ` )

if [ -z "$SHAPE" ]
then
	echo "I couldn't find the shapeit executable; please provide the path including the executable name"
	echo "e.g. ~/miniconda3/envs/plink_env/bin/shapeit"
	read SHAPEIT
	TEST=( ` echo $SHAPEIT | grep "shapeit" ` )
	if [ -z "$TEST" ]; then
		echo "Nope, still couldn't read it. You done messed up. \#epicfail"
		exit 1
	fi
else
	SHAPEIT=$SHAPE
fi

echo "I understand the shapeit executable to be at " $SHAPEIT

###run ShapeIt on input formatted by the R script
$SHAPEIT --input-ped genotypes.ped genotypes.map -O genotypes.phased --force

###write R file (parse_shapeit_output.R) to convert ShapeIt output to a more usable format

echo "#######-----------read in output from ShapeIt (space-delimited), reformat and number haplotypes, output in XXX format" > parse_shapeit_output.R

echo "shapeit_haplotypes=read.delim(\"genotypes.phased.haps\",sep=\" \",header=FALSE)" >> parse_shapeit_output.R
echo "npos<-nrow(shapeit_haplotypes)" >> parse_shapeit_output.R
echo "nindivH<-(ncol(shapeit_haplotypes)-5)/2" >> parse_shapeit_output.R
echo "shapeit_samples=read.delim(\"genotypes.phased.sample\",sep=\" \",header=FALSE)" >> parse_shapeit_output.R
echo "colnames(shapeit_samples)<-as.character(unlist(shapeit_samples[1,]))" >> parse_shapeit_output.R
echo "shapeit_samples<-shapeit_samples[3:nrow(shapeit_samples),]" >> parse_shapeit_output.R
echo "nindivS<-nrow(shapeit_samples)" >> parse_shapeit_output.R
  
echo "#format check" >> parse_shapeit_output.R
echo "if (!nindivH==nindivS) {" >> parse_shapeit_output.R
echo "  print(\"FORMAT ERROR: number of individuals is different in input files!!!!\")" >> parse_shapeit_output.R
echo "} else {" >> parse_shapeit_output.R

echo "shapeit_haplotypes_nucleotides<-shapeit_haplotypes" >> parse_shapeit_output.R
echo "for (m in 1:nrow(shapeit_haplotypes_nucleotides)) {" >> parse_shapeit_output.R
echo "  major_allele=as.character(shapeit_haplotypes_nucleotides[m,4])" >> parse_shapeit_output.R
echo "  minor_allele=as.character(shapeit_haplotypes_nucleotides[m,5])" >> parse_shapeit_output.R
echo "  shapeit_haplotypes_nucleotides[m,]<-replace(shapeit_haplotypes_nucleotides[m,],shapeit_haplotypes_nucleotides[m,]==\"0\",major_allele)" >> parse_shapeit_output.R
echo "  shapeit_haplotypes_nucleotides[m,]<-replace(shapeit_haplotypes_nucleotides[m,],shapeit_haplotypes_nucleotides[m,]==\"1\",minor_allele)" >> parse_shapeit_output.R
echo "}" >> parse_shapeit_output.R

echo "shapeit_haplotypes_nucleotidesVERT<-as.data.frame(t(shapeit_haplotypes_nucleotides))" >> parse_shapeit_output.R
echo "colnames(shapeit_haplotypes_nucleotidesVERT)<-unlist(shapeit_haplotypes_nucleotidesVERT[2,])" >> parse_shapeit_output.R
echo "shapeit_haplotypes_nucleotidesVERTmeta<-shapeit_haplotypes_nucleotidesVERT[1:5,]" >> parse_shapeit_output.R
echo "shapeit_haplotypes_nucleotidesVERT<-shapeit_haplotypes_nucleotidesVERT[6:nrow(shapeit_haplotypes_nucleotidesVERT),]" >> parse_shapeit_output.R
echo "shapeit_haplotypes_nucleotidesVERT[,1:npos]<-lapply(shapeit_haplotypes_nucleotidesVERT[,1:npos],as.character)" >> parse_shapeit_output.R
echo "####create multi-locus hapltoypes by pasting nucleotides together" >> parse_shapeit_output.R
echo "shapeit_haplotypes_nucleotidesVERT\$MLhap<-apply( shapeit_haplotypes_nucleotidesVERT[ , 1:npos] , 1 , paste , collapse = \"\" )" >> parse_shapeit_output.R
echo "shapeit_haplotypes_nucleotidesVERT\$indiv<-rep(shapeit_samples\$ID_2,each=2)" >> parse_shapeit_output.R
echo "nucleotide_haplotypes<-shapeit_haplotypes_nucleotidesVERT[,c(\"indiv\",\"MLhap\")]" >> parse_shapeit_output.R

echo "####make a table of the number of occurences (frequency) of unique multi-locus haplotypes" >> parse_shapeit_output.R
echo "unique_haplotypes<-as.data.frame(table(nucleotide_haplotypes\$MLhap))" >> parse_shapeit_output.R
echo "colnames(unique_haplotypes)<-c(\"MLhap\",\"Freq\")" >> parse_shapeit_output.R
echo "####number the haplotypes in descending order by frequency (most common first)" >> parse_shapeit_output.R
echo "unique_haplotypes\$hap_num<-as.integer((max(rank(unique_haplotypes\$Freq,ties.method=\"last\"))+1)-rank(unique_haplotypes\$Freq,ties.method=\"random\"))" >> parse_shapeit_output.R
echo "####create vector of haplotype numbers starting with zero (0) for VCF" >> parse_shapeit_output.R
echo "unique_haplotypes\$hap_num_ref0<-unique_haplotypes\$hap_num-1" >> parse_shapeit_output.R
echo "####determine which haplotypes (by number) each individual has" >> parse_shapeit_output.R
echo "haps_vector<-NULL" >> parse_shapeit_output.R
echo "haps_vector_ref0<-NULL" >> parse_shapeit_output.R
echo "for (i in 1:nrow(nucleotide_haplotypes)) {" >> parse_shapeit_output.R
echo "  haps_vector<-c(haps_vector,as.character(unique_haplotypes[unique_haplotypes[,1] %in% nucleotide_haplotypes[i,\"MLhap\"],\"hap_num\"]))" >> parse_shapeit_output.R
echo "  haps_vector_ref0<-c(haps_vector_ref0,as.character(unique_haplotypes[unique_haplotypes[,1] %in% nucleotide_haplotypes[i,\"MLhap\"],\"hap_num_ref0\"]))" >> parse_shapeit_output.R
echo "}" >> parse_shapeit_output.R
echo "nucleotide_haplotypes\$hap_num<-haps_vector" >> parse_shapeit_output.R
echo "nucleotide_haplotypes\$hap_num_ref0<-haps_vector_ref0" >> parse_shapeit_output.R

echo "#create genotypes for each individual based on assigned haplotypes" >> parse_shapeit_output.R
echo "nuc_genotypes<-NULL" >> parse_shapeit_output.R
echo "num_genotypes<-NULL" >> parse_shapeit_output.R
echo "num_genotypes_ref0<-NULL" >> parse_shapeit_output.R
echo "for (i in 1:nindivS) {" >> parse_shapeit_output.R
echo "  df<-nucleotide_haplotypes[nucleotide_haplotypes\$indiv %in% shapeit_samples\$ID_2[i],]" >> parse_shapeit_output.R
echo "  nuc_genotypes<-c(nuc_genotypes,paste(df[,\"MLhap\"],collapse=\"/\"))" >> parse_shapeit_output.R
echo "  num_genotypes<-c(num_genotypes,paste(df[,\"hap_num\"],collapse=\"/\"))" >> parse_shapeit_output.R
echo "  num_genotypes_ref0<-c(num_genotypes_ref0,paste(df[,\"hap_num_ref0\"],collapse=\"/\"))" >> parse_shapeit_output.R
echo "}" >> parse_shapeit_output.R
echo "genotypes<-as.data.frame(cbind(as.character(shapeit_samples\$ID_2),nuc_genotypes,num_genotypes,num_genotypes_ref0))" >> parse_shapeit_output.R
echo "colnames(genotypes)<-c(\"indiv\",\"MLgeno\",\"NUMgeno\",\"NUMgenoREF0\")" >> parse_shapeit_output.R

echo "write.table(genotypes,file=\"inferred_genotypes.txt\",append=FALSE,quote=FALSE,sep=\"\t\",col.names=T,row.names=F)" >> parse_shapeit_output.R
echo "write.table(nucleotide_haplotypes,file=\"inferred_haplotypes.txt\",append=FALSE,quote=FALSE,sep=\"\t\",col.names=T,row.names=F)" >> parse_shapeit_output.R
echo "write.table(unique_haplotypes,file=\"haplotypes_table.txt\",append=FALSE,quote=FALSE,sep=\"\t\",col.names=T,row.names=F)" >> parse_shapeit_output.R

echo "#create a VCF file with the new haplotypes as a multi-allelic marker" >> parse_shapeit_output.R
echo "unique_haplotypes_ordered<-unique_haplotypes[order(unique_haplotypes\$hap_num,decreasing = FALSE),]" >> parse_shapeit_output.R
echo "format=\"##fileformat=VCFv4.2\"" >> parse_shapeit_output.R
echo "write.table(format,file=\"inferred_haplotypes.vcf\",append=FALSE,quote=FALSE,sep=\"\t\",col.names=F,row.names=F)" >> parse_shapeit_output.R
echo "header<-paste(c(\"#CHROM\",\"POS\",\"ID\",\"REF\",\"ALT\",\"QUAL\",\"FILTER\",\"INFO\",\"FORMAT\",as.character(genotypes\$indiv)),collapse=\"\t\")" >> parse_shapeit_output.R
echo "write.table(header,file=\"inferred_haplotypes.vcf\",append=TRUE,quote=FALSE,sep=\"\t\",col.names=F,row.names=F)" >> parse_shapeit_output.R
echo "vcf_genotypes<-paste(c(\"1\",as.character(shapeit_haplotypes[1,3]),\"MLhaplotype\",as.character(unique_haplotypes_ordered[1,1]),paste(unique_haplotypes_ordered[2:nrow(unique_haplotypes_ordered),1],collapse=\",\"),\".\",\"PASS\",\".\",\"GT\",as.character(num_genotypes_ref0)),collapse=\"\t\")" >> parse_shapeit_output.R
echo "write.table(vcf_genotypes,file=\"inferred_haplotypes.vcf\",append=TRUE,quote=FALSE,sep=\"\t\",col.names=F,row.names=F)" >> parse_shapeit_output.R

echo "#create a hapmap output for Gapit" >> parse_shapeit_output.R
echo "nucleotide_haplotypes_hapmap<-nucleotide_haplotypes[,c(1,rep(4,times=(nrow(unique_haplotypes)-1)))]" >> parse_shapeit_output.R
echo "num_genotypes_ref0_nosep<-NULL" >> parse_shapeit_output.R
echo "alternate_allele_genotypes<-data.frame(matrix(nrow = length(shapeit_samples\$ID_2), ncol = 0))" >> parse_shapeit_output.R
echo "rownames(alternate_allele_genotypes)<-shapeit_samples\$ID_2" >> parse_shapeit_output.R
echo "for (c in 1:(ncol(nucleotide_haplotypes_hapmap)-1)) {" >> parse_shapeit_output.R
echo "  for (i in 1:nrow(shapeit_samples)) {" >> parse_shapeit_output.R
echo "    df<-nucleotide_haplotypes_hapmap[nucleotide_haplotypes_hapmap\$indiv %in% shapeit_samples\$ID_2[i],c+1]" >> parse_shapeit_output.R
echo "    df[df==0]<-\"A\"" >> parse_shapeit_output.R
echo "    df[df==(c)]<-\"T\"" >> parse_shapeit_output.R
echo "    df[!(df==\"A\"|df==\"T\")]<-\"A\"" >> parse_shapeit_output.R
echo "    num_genotypes_ref0_nosep<-c(num_genotypes_ref0_nosep,paste(df,collapse=\"\"))" >> parse_shapeit_output.R
echo "  }" >> parse_shapeit_output.R
echo "  alternate_allele_genotypes[,c]<-num_genotypes_ref0_nosep" >> parse_shapeit_output.R
echo "  num_genotypes_ref0_nosep<-NULL" >> parse_shapeit_output.R
echo "}" >> parse_shapeit_output.R
echo "alternate_allele_genotypesVERT<-NULL" >> parse_shapeit_output.R
echo "alternate_allele_genotypesVERT\$rs<-paste(c(rep(\"haplotype_\",times=ncol(alternate_allele_genotypes))),c(2:(ncol(alternate_allele_genotypes)+1)),sep = \"\")" >> parse_shapeit_output.R
echo "alternate_allele_genotypesVERT\$alleles<-rep(\"NA\",times=ncol(alternate_allele_genotypes))" >> parse_shapeit_output.R
echo "alternate_allele_genotypesVERT\$chrom<-rep(\"Chr\",times=ncol(alternate_allele_genotypes))" >> parse_shapeit_output.R
echo "alternate_allele_genotypesVERT\$pos<-c(rep(as.integer(shapeit_haplotypes[1,3]),times=ncol(alternate_allele_genotypes)))+0:(ncol(alternate_allele_genotypes)-1)" >> parse_shapeit_output.R
echo "alternate_allele_genotypesVERT\$strand<-rep(\"NA\",times=ncol(alternate_allele_genotypes))" >> parse_shapeit_output.R
echo "alternate_allele_genotypesVERT\$assembly<-rep(\"NA\",times=ncol(alternate_allele_genotypes))" >> parse_shapeit_output.R
echo "alternate_allele_genotypesVERT\$center<-rep(\"NA\",times=ncol(alternate_allele_genotypes))" >> parse_shapeit_output.R
echo "alternate_allele_genotypesVERT\$protLSID<-rep(\"NA\",times=ncol(alternate_allele_genotypes))" >> parse_shapeit_output.R
echo "alternate_allele_genotypesVERT\$assayLSID<-rep(\"NA\",times=ncol(alternate_allele_genotypes))" >> parse_shapeit_output.R
echo "alternate_allele_genotypesVERT\$panel<-rep(\"NA\",times=ncol(alternate_allele_genotypes))" >> parse_shapeit_output.R
echo "alternate_allele_genotypesVERT\$QCcode<-rep(\"NA\",times=ncol(alternate_allele_genotypes))" >> parse_shapeit_output.R
echo "alternate_allele_genotypesVERT<-cbind(alternate_allele_genotypesVERT,as.data.frame(t(alternate_allele_genotypes)))" >> parse_shapeit_output.R
echo "colnames(alternate_allele_genotypesVERT)=c(\"rs\",\"alleles\",\"chrom\",\"pos\",\"strand\",\"assembly\",\"center\",\"protLSID\",\"assayLSID\",\"panel\",\"QCcode\",rownames(alternate_allele_genotypes))" >> parse_shapeit_output.R

echo "write.table(alternate_allele_genotypesVERT,file=\"inferred_haplotypes_genotypes_hapmap.txt\",append=FALSE,quote=FALSE,sep=\"\t\",col.names=T,row.names=F,na=\"NA\")" >> parse_shapeit_output.R

echo "}" >> parse_shapeit_output.R
echo "quit(save = \"no\", status = 0, runLast = FALSE);" >> parse_shapeit_output.R

###execute R file to parse ShapeIt output(parse_shapeit_output.R)
R CMD BATCH parse_shapeit_output.R


mv inferred_genotypes.txt $4.inferred_genotypes.txt
mv inferred_haplotypes.txt $4.inferred_haplotypes.txt
mv haplotypes_table.txt $4.haplotypes_table.txt
mv genotypes.phased.haps $4.phased.haps
mv genotypes.phased.sample $4.phased.sample
mv inferred_haplotypes.vcf $4.phased.vcf
mv inferred_haplotypes_genotypes_hapmap.txt $4.phased.hapmap_genotypes.txt


if [ "$1" == "HAPLOTYPE" ]
then
	echo "Haplotyping finished, exiting."
	exit 1
elif [ "$1" == "NETWORK" ]
then
	echo "Haplotyping finished, proceeding to network inference."
	cp $4.inferred_haplotypes.txt haplotypes.txt
fi

###write R file (make_fasta_network.R) to convert ShapeIt output to a more usable format

echo "######writing fasta" > make_fasta_network.R
echo "nucleotide_haplotypes<-read.delim(\"haplotypes.txt\",colClasses = c(\"factor\"))" >> make_fasta_network.R
echo "is.odd <- function(x) x %% 2 != 0" >> make_fasta_network.R
echo "fasta<-NULL" >> make_fasta_network.R
echo "for (i in 1:nrow(nucleotide_haplotypes)) {" >> make_fasta_network.R
echo "  if (is.odd(i)) {" >> make_fasta_network.R
echo "    fasta<-rbind(fasta,paste(\">\",as.character(nucleotide_haplotypes\$indiv[i]),\"A\",sep=\"\"))  " >> make_fasta_network.R
echo "  } else {" >> make_fasta_network.R
echo "    fasta<-rbind(fasta,paste(\">\",as.character(nucleotide_haplotypes\$indiv[i]),\"B\",sep=\"\"))  " >> make_fasta_network.R
echo "  }" >> make_fasta_network.R
echo "  fasta<-rbind(fasta,as.character(nucleotide_haplotypes\$MLhap[i]))" >> make_fasta_network.R
echo "}" >> make_fasta_network.R
echo "write.table(fasta,file=\"haplotypes.fasta\",append=FALSE,quote=FALSE,sep=\"\t\",col.names=F,row.names=F)" >> make_fasta_network.R

echo "######making haplotype network" >> make_fasta_network.R
echo "list.of.packages <- c(\"pegas\")" >> make_fasta_network.R
echo "new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,\"Package\"])]" >> make_fasta_network.R
echo "if(length(new.packages)) print(\"RALERT: Installing dependencies for first time use....\")" >> make_fasta_network.R
echo "if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')" >> make_fasta_network.R
echo "require(pegas)" >> make_fasta_network.R
echo "haps<-read.dna(\"haplotypes.fasta\", format=\"fasta\")" >> make_fasta_network.R
echo "NET<-haploNet(haplotype(haps))" >> make_fasta_network.R

echo "##making PDFs of various network figures" >> make_fasta_network.R
echo "pdf(\"rplot.pdf\",width = 10, height = 5) " >> make_fasta_network.R
echo "plot(NET, size = 1, fast = FALSE, labels = TRUE,cex=0.5, show.mutation = 3 )" >> make_fasta_network.R
echo "dev.off()" >> make_fasta_network.R
echo "pdf(\"rplot2.pdf\",width = 10, height = 5) " >> make_fasta_network.R
echo "plot(NET, size = 1, fast = FALSE, labels = TRUE,cex=0.5, show.mutation = 0 )" >> make_fasta_network.R
echo "dev.off()" >> make_fasta_network.R
echo "pdf(\"rplot3.pdf\",width = 10, height = 5) " >> make_fasta_network.R
echo "plot(NET, size = 1, fast = FALSE, labels = TRUE,cex=0.5, show.mutation = 0 ,threshold=c(1,2))" >> make_fasta_network.R
echo "dev.off()" >> make_fasta_network.R
echo "pdf(\"rplot4.pdf\",width = 10, height = 5) " >> make_fasta_network.R
echo "plot(NET, size = attr(NET, \"freq\"), fast = FALSE, scale.ratio=100,labels = FALSE,cex=0.5, show.mutation = 0,threshold=0)" >> make_fasta_network.R
echo "dev.off()" >> make_fasta_network.R
echo "pdf(\"rplot5.pdf\",width = 10, height = 5) " >> make_fasta_network.R
echo "plot(NET, size = attr(NET, \"freq\"), fast = FALSE, scale.ratio=100,labels = TRUE,cex=0.5, show.mutation = 0,threshold=0)" >> make_fasta_network.R
echo "dev.off()" >> make_fasta_network.R
echo "quit(save = \"no\", status = 0, runLast = FALSE);" >> make_fasta_network.R


###execute R file to parse ShapeIt output(parse_shapeit_output.R)
R CMD BATCH make_fasta_network.R

mv haplotypes.fasta $4.haplotypes.fasta
mv rplot.pdf $4.Rplot.pdf
mv rplot2.pdf $4.Rplot2.pdf
mv rplot3.pdf $4.Rplot3.pdf
mv rplot4.pdf $4.Rplot4.pdf
mv rplot5.pdf $4.Rplot5.pdf

rm markers.txt &> /dev/null 
rm genotypes.txt &> /dev/null 
rm reformat_Progeny_PLINK.R &> /dev/null 
rm genotypes.ped &> /dev/null 
rm genotypes.map &> /dev/null 
rm genotypes.phased.haps &> /dev/null 
rm genotypes.phased.sample &> /dev/null 
rm parse_shapeit_output.R &> /dev/null 
rm make_fasta_network.R &> /dev/null 
rm *.Rout &> /dev/null 
rm inferred_genotypes.txt &> /dev/null 
rm inferred_haplotypes.txt &> /dev/null 
rm haplotypes_table.txt &> /dev/null 
rm haplotypes.txt &> /dev/null 
rm shapeit_* &> /dev/null 

