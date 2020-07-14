#' @name AnnotateVCFwithExAC
#' @title Annotate VCF file and add annotations from ExAC database
#' @description Parses an input VCF file to extract specific variant annotations from the file and externally from ExAC database
#' @author Nitya Sharma
#' @param filename input VCF file 
#' @return Returns a dataframe with 13 columns
#' "chr" - chromosome
#' "start" - start position
#' "stop" - stop position
#' "ref_allele" - reference allele
#' "alt_allele" - alternate allele
#' "read_depth" - read depth at locus
#' "perc_alt.readss" - percent of reads supporting alternate alelle
#' "perc_ref.reads" - percent of reads supporting reference allele
#' "var_type" - type of variant (snp, indel, complex, etc)
#' "Allele_Freq", - frequency of alelle in population (from ExAC)
#' "Mutation_type" - major consequence of variant (from ExAC)
#' "SIFT_Score" - SIFT score (from ExAC)
#' "HGNC" - HGNC gene symbol (from EXAC)
#' @references http://exac.hms.harvard.edu
#' @examples
#' annotated_output <- annotatevcf(vcf_inputfile)

annotatevcf <- function(fn){
  
  # ReadVCF and separate multiallelic variant sites to individual rows
  vcf <- VariantAnnotation::readVcf(fn)
  vcf.ex <- VariantAnnotation::expand(vcf)
  #rm(vcf)
  
  # make GRanges object into DF 
  chr <- data.frame(DelayedArray::rowRanges(vcf.ex))[,"seqnames"]
  start <- data.frame(DelayedArray::rowRanges(vcf.ex))[,"start"]
  end <- data.frame(DelayedArray::rowRanges(vcf.ex))[,"end"]
  alt <- data.frame(VariantAnnotation::alt(vcf.ex))  
  ref <- data.frame(VariantAnnotation::ref(vcf.ex))
  
  # Extract Info data
  myparams <- c("DP", "AO", "RO", "TYPE" )
  myinfo <- VariantAnnotation::info(vcf.ex)[,myparams]
  numaltreads <- myinfo[,"AO"]
  perc_alt_reads <- myinfo[,"AO"]/myinfo[,"DP"]*100
  perc_ref_reads <- myinfo[,"RO"]/myinfo[,"DP"]*100
  annodata <- cbind(chr, start, end, ref, alt, myinfo[,"DP"], numaltreads, perc_alt_reads, perc_ref_reads,myinfo[,"TYPE"])
  exacq <- paste(annodata[,1],"-", annodata[,2], "-", annodata[,4], "-", annodata[,5], sep = "")
  annodata <- cbind(annodata, exacq)
  colnames(annodata) <- c("chr", "start", "end", "ref_allele", "alt_allele", "read_depth", "num_alt_reads", "perc_alt-reads", "perc_ref-reads", "var_type", "exacq")
  rownames(annodata) <- exacq
  
  # Extract annotation data from ExAC db
  exac.out <- query.exac(annodata$exacq)
  
  # Order keys of exac output to that in the annodata dataframe
  ordered <- names(exac.out)[order(match(names(exac.out), annodata$exacq))]
  
  # Make df from exac data with relevant info
  exacdf <- AnnotatebyExac(exac.out,ordered)
  rownames(exacdf) <- exacq

  anno.out <- cbind(data.frame(annodata), exacdf)
  anno.out <- anno.out[!colnames(anno.out) == "exacq"]
  
  return(anno.out)
  
}

#' @name query.exac
#' @title Query ExAC database
#' @description Queries ExAC database and extracts annotations for variant IDs
#' @param vector of query IDs
#' @examples
#' query.exac(query.vec)
#' @return Returns a list of lists of annotations extracted from ExAC database
query.exac <- function(query.vec){
  exacq <- httr::POST(url = "http://exac.hms.harvard.edu/rest/bulk/variant",
                      body = jsonlite::toJSON(query.vec), encode = "json"
  )
  json.res <- httr::content(exacq)
  return(json.res)
}

#' @name AnnotatebyExac
#' @title Parses ExAC database annotation 
#' @description This function parses the annotation info outputted by ExAC and organizes it into a dataframe with specific annotations -- Allele Frequency, Type of Mutation (Major Consequence of mutation), SIFT score to predict impact and HGNC symbol
#' @param exacquery_out
#' @param queryIDs
#' @examples
#' AnnotatebyExac(exac.res,ids)
#' @return Returns a list of lists of annotations extracted from ExAC database
AnnotatebyExac <- function(exac.res, ids){
  
  exacdf <- data.frame()
  for(i in ids){
    qry <- i
    
    alfreq <- exac.res[[i]]$variant$allele_freq
    alfreq <- if (!is.null(alfreq) ) alfreq else "NA"
    
    # Not every variant has annotation data
    if (!is.null(exac.res[[i]]$variant$vep_annotations) && !rlang::is_empty(exac.res[[i]]$variant$vep_annotations)){
      muttype <- unlist(exac.res[[i]]$variant$vep_annotations[[1]])["major_consequence"]
      
      #### if is.null OR is empty string
      muttype <- if (!is.null(muttype) && !muttype == "" ) muttype else "NA"
      
      sift <- unlist(exac.res[[i]]$variant$vep_annotations[[1]])["SIFT"]
      sift <- if (!is.null(sift) && !sift == "" ) sift else "NA" 
      
      hgnc <-  unlist(exac.res[[i]]$variant$vep_annotations[[1]])["SYMBOL"]
      hgnc <- if (!is.null(hgnc) && !hgnc== "") hgnc else "NA" 
    }
    else{ 
      muttype <- "NA"
      sift <- "NA"
      hgnc <- "NA"
    }
    exacdf <- rbind(exacdf, cbind(qry, alfreq, muttype, sift, hgnc))
  }
  
  colnames(exacdf) <- c("exacq", "Allele_Freq", "Mutation_type", "SIFT_Score", "HGNC")
  return(exacdf)
}

