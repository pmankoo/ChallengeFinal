Code for Annotation and parsing of variant calls to address questions posed by the Challenge.
=============================================================================================

For this challenge, each variant must be annotated with the following pieces of information:
------------
1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) If there are multiple possibilities, annotate with the most deleterious possibility.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from Broad Institute ExAC Project API
(API documentation is available here: http://exac.hms.harvard.edu/)
6. Additional optional information from ExAC that you feel might be relevant.

Usage:
------------
Assuming you already have snpeff installed on your machine, please modify the path in "ToolsConfig.py" leading to the /snpEff directory.
- For convenience, the output of snpEff annotated file "Challenge_data.snpeff.vcf" is provided so the examiner dosen't need to re-run SnpEff.
- [OPTIONAL]: If there is a desire to run "SnpEff", please uncomment the following line in the "variantannot_Mankoo.py" code after modifying the path to /snpEff directory.
- #run_snpeff('Challenge_data.vcf') 

### Starting from scratch
	python3 variantannot_Mankoo.py

OUTPUT COLUMNS APPENDED TO INPUT VCF FILE:
------------------------------------------------------------
HEADERS				| Description |
--------			| ----------- |
NumAltAlleles			| The total number of alternate alleles per variant.
MostDeleteriousPerVariant	| The most deleterious effect per allele per variant as defined by SnpEff.
EXAC_AF_PERALLELE		| The allele frequency per allele per variant from EXAC
normalResults			| For normal sample, reporting "total read depth, allelic counts per allele, %Reads supporting an allele per variant for all alleles"
vaf5Results			| For vaf5 sample, reporting "total read depth, allelic counts per allele, %Reads supporting an allele per variant for all alleles"


Extra Information:
--------------------
- The code took 25.372075080871582 seconds to run on my MacBook Pro wth pre-calculated snpEff and about 87 seconds with snpEff run.
- The Output vcf formatted file provided is "Challenge_data_snpeff_modified.vcf".
