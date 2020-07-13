# TempusCodingChallenge

This program parses and annotates a VCF file with locus specific information and extracts external annotations from the ExAC database. Results from these analyses can be used to identify putatively deleterious causal variants. 

Download the .tar file from https://github.com/nvsharma24/TempusCodingChallenge/AnnotateVCFwithExAC/R/ and install AnnotateVCFwithExAC_0.1.0.tar from your local directory 

Necessary packages to run the scripts are: 
- VariantAnnotation (Bioconductor - use BiocManager::install("VariantAnnotation") #NOTE:must be used for latest version of R (4.0)#
- httr
- rlang

annotatevcf() takes in a .vcf files and outputs a dataframe of 13 columns that can be written to a file and includes information gathered from both the VCF file and the ExAC database. 
AnnotatebyExac() can be used independently to query the ExAC database and generates a dataframe with specific information from that database only. This function is used within annotatevcf() 
query.exac() simply queries the ExAC database and outputs a list of lists with locus specific annotations. This function is used within AnnotatebyExac()
Refer to script for more detailed information regarding inputs for AnnotatebyExac and query.exac()
