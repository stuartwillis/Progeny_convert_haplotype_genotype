#Progeny_convert_haplotype_genotype.sh

A bash script to convert Progeny exported files to PLINK or Haploview formats, or haplotype the input genotype file and output inferred haplotypes and the multi-locus genotypes of the individuals.

Two inputs are required for conversion or haplotyping:

1. A genotype file exported from Progeny
2. A file of marker names, position, type, etc.

The genotype file may include genotypes other than those to be haplotyped/tested, and they will be ignored (along with any other columns).

Which markers to haplotype/import are specified by the marker file.

Marker column headings must include "Assay","Chromosome","Position","Presumed Type"; all others will be ignored and order does not matter.

This includes real or dummy chromosomes and positions, but only a SINGLE chromosome by name will be analyzed (the first encountered ignoring the sex marker), and markers will be haplotyped in the order of positions.

Optionally include a "Sex Marker" (others are ignored) in the column "Presumed Type" for sex assignment; otherwise all individuals will be marked unknown

If included, the sex marker must be FIRST in the marker file (chromosome will be ignored)

Allele names are made from the assay field assuming the marker names are a simple extension of that (i.e. Assay-A1 or Assay.A1), so be sure these match genotypes file.

Correct usage: SCRIPT.sh [PLINK || HAPLOVIEW || HAPLOTYPE] genotype-file marker-file output_prefix [optional: individual-filter-proportion]

CHOOSE one mode: REFORMAT from Progeny to PLINK (Plink1.9 ped/map) or HAPLOTYPE using Shapeit (output is haps/sample combined with your phenotype data)

MAKE SURE YOUR INPUT FILES ARE UNIX-DELIMTED WITH ONE(1) HANGING LINE

Optional proportion to filter individuals by missing data must be between 0 and 1 (default is 0.2). More missing data results in lower confidence haplotypes.


##Dependencies

Shapeit2 - https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download

R base (available in the path)


##Output

0utput will be preceded with output_prefix.

Output from PLINK/HAPLOVIEW conversion includes:

output_prefix.map
output_prefix.ped

with the map in the appropriate format.

Output from haplotyping includes:

output_prefix.phased.haps

(raw haplotype output from Shapeit2)

output_prefix.phased.sample

(sample information for each individual from Shapeit2; matches the order of the .haps file)

output_prefix.haplotypes_table.txt

(the inferred haplotypes with their frequency and numerical code/rank by frequency)

output_prefix.inferred_genotypes.txt

(the multi-allelic genotypes of each individual, as both nucleotides and rank codes, delimited by '/')

output_prefix.inferred_haplotypes.txt

(the nucleotide and rank code haplotypes for each individual, one line per haplotype, two lines per individual)

output_prefix.phased.vcf

(multi-allelic vcf with individuals scored for haplotypes; no filter for minimum haplotype frequency is applied)

output_prefix.phased.hapmap_genotypes.txt

(genotypes of each alternate allele relative to the reference coded in hapmap format for Gapit. Each individual is scored for T=alternate or A=reference OR any other allele, for each of all inferred alternate haplotypes. This was, each haplotype can be scored against the phenotypic variation observed with all other alleles. For simultaneous testing, see the R package haplo.stats)

output_prefix.haplotypes.fasta

(fasta format file of nucleotide hapltoypes, two lines per individual suffixed A and B)



###Example usage: bash Progeny_convert_haplotype_genotype.sh HAPLOTYPE progeny_genotypes.txt marker_file.txt output_prefx 0.1