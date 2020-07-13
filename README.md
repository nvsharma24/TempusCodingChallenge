# TempusCodingChallenge - Nitya Sharma

This program parses and annotates a VCF file with locus specific information and extracts external annotations from the ExAC database. Results from these analyses can be used to identify putatively deleterious causal variants. 

Download the .tar file from https://github.com/nvsharma24/TempusCodingChallenge/AnnotateVCFwithExAC/R/ and install AnnotateVCFwithExAC_0.1.0.tar from your local directory 

Necessary packages to run the scripts are: 
- VariantAnnotation (Bioconductor - use BiocManager::install("VariantAnnotation") #NOTE:must be used for latest version of R (4.0)#
- httr
- rlang

annotatevcf() takes in a .vcf files and outputs a dataframe of 13 columns that can be written to a tab delimited file and includes information gathered from both the VCF file and the ExAC database. 

The columns are as follows:

- "chr" - chromosome
- "start" - start position
- "end" - end position
- "ref_allele" - reference allele
- "alt_allele" - alternate allele
- "read_depth" - read depth at locus
- "perc_alt.reads" - percent of reads supporting alternate alelle
- "perc_ref.reads" - percent of reads supporting reference allele
- "var_type" - type of variant (snp, indel, complex, etc)
- "Allele_Freq", - frequency of alelle in population (from ExAC)
- "Mutation_type" - major consequence of variant (from ExAC)
- "SIFT_Score" - SIFT score (from ExAC)
- "HGNC" - HGNC gene symbol (from EXAC)


AnnotatebyExac() can be used independently to query the ExAC database and generates a dataframe with specific information from that database only. This function is used within annotatevcf() 

query.exac() simply queries the ExAC database and outputs a list of lists with locus specific annotations. This function is used within AnnotatebyExac()
Refer to AnnotateVCF_v3.R script for more detailed information regarding inputs for AnnotatebyExac and query.exac().
