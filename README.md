# TempusCodingChallenge

This program parses and annotates a VCF file with locus specific information and extracts external annotations from the ExAC database.

The output of annotatevcf() is a dataframe of 13 columns that can be written to a file and included information gathered from both the VCF file and the ExAC database. 
AnnotatebyExac() can be used independently to query the ExAC database and annotate a dataframe with specific information from that database only. 
