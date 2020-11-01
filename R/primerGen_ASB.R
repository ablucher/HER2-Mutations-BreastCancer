#rm(list=ls())
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg19")

########################################################################
#Script for Samuel Tsang - Mutation Library/ Primer Generation
#Script original author: JZhang@MDAnderson for Samuel Tsang
#Script edits: ABlucher@OHSU for Samuel Tsang
  #10-28-20; update BioCManager and packages; update file paths to here() for clear project directory mangagement
########################################################################

#first un-comment and run the install command, then run the library() command
#BiocManager::install(version = "3.12")
library(BiocManager)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19") #this takes a little while
library(BSgenome.Hsapiens.UCSC.hg19)
#install.packages("here")
library(here) #this tells R where our main project folder is
              #if you use this package to specify filenames, then they are compatable across PC/Macs :) see bottom of script

primerGen <- function(variantFile, orfFile, outName,
                      baseToMatch = 25, misMatch = 3, genome = "BSgenome.Hsapiens.UCSC.hg19",flanking = 30,
                      chrom = "Chromosome", start = "Start_Position", end = "End_Position", entrez = "Entrez_Gene_ID", 
                      hugo = "Hugo_Symbol", ref = "Reference_Allele", tumorAllele1 = "Tumor_Seq_Allele1", 
                      tumorAllele2 = "Tumor_Seq_Allele2", strand = "Transcript_Strand", 
                      type = "Variant_Type"){

  variants <- read.delim(variantFile, sep = "\t", header = TRUE, as.is = TRUE)
  variants <- cbind(variants, index = paste("indel", 1:nrow(variants)))
  orfSeq <- read.delim(orfFile, sep = "\t", header = TRUE, as.is = TRUE)
  primers <- apply(variants, 1, findPrimers, orfSeq = orfSeq, 
                   misMatch = misMatch, baseToMatch = baseToMatch, flanking = flanking, 
                   genome = genome, chrom = chrom, start = start, end = end, entrez = entrez, 
                   hugo = hugo, ref = ref, tumorAllele1 = tumorAllele1, 
                   tumorAllele2 = tumorAllele2, strand = strand, type = type, index = "index")
  if(class(primers) == "matrix"){
    
    primers <- t(primers)
    
    colnames(primers) <- c("symbol", "chrom", "start", "end", "reference", "strand", "mutation",
                           
                           "cloneID", "mappedPrimer", "AA_Change", "Entrez_Gene_ID", "Variant_Type", "Variant_Classification")
  }else{
    primers <- do.call("rbind", primers)
  }
  write.table(primers, outName, sep = "\t", col.names = TRUE, row.names = FALSE, 
              quote = FALSE)
  
  return(primers)
}


#### new set of code base on discussion July 20 2012 #######
#### orfSeq = newOrfCollection.tsv
#### apply(mutations, 1, findPrimers)
findPrimers <- function(mut, orfSeq, flanking = 30, misMatch = 0, 
                        baseToMatch = 20, organism = "Hsapiens", 
                        genome = "BSgenome.Hsapiens.UCSC.hg19", orf = "ORF.sequence", 
                        clone = "Clone.ID", chrom = "Chromosome", start = "Start_Position", 
                        end = "End_Position", entrez = "Entrez_Gene_ID", hugo = "Hugo_Symbol", 
                        ref = "Reference_Allele", tumorAllele1 = "Tumor_Seq_Allele1", 
                        tumorAllele2 = "Tumor_Seq_Allele2", strand = "Transcript_Strand", 
                        type = "Variant_Type", AA_Change = "Amino_Acid_Change", Entrez_Gene_ID = "Entrez_Gene_ID",
                        Variant_Type = "Variant_Type", Variant_Classification = "Variant_Classification", index = "index"){
  
  require(genome, character.only = TRUE) 
  orfs <- orfSeq[gsub(" ", "", orfSeq[, "LocusLink.ID"]) %in% 
                   gsub(" ", "", mut[entrez])  |  gsub(" ", "", orfSeq[, "Gene.Symbol"])
                 %in% gsub(" ", "", mut[hugo]), , drop = FALSE]
  
  if(nrow(orfs) < 1){
    return(NA)
  }
  orfs <- try(makeStringSet(orfs[, orf], orfs[, clone]), silent = TRUE)
  if(class(orfs) == "try-error"){
    
    print(paste("Failed to get orf string set for", mut[index]))
    return(NA)
  }    
  mutAllele <- getMutAllele(mut[ref], mut[c(tumorAllele1, tumorAllele2)])
  mutFlanking <- try(getFlanking(paste("chr", gsub(" ", "", mut[chrom]), sep = ""), 
                                 as.numeric(mut[start]), as.numeric(mut[end]), what = mut[type],
                                 flanking = flanking, organism = organism, genome = genome), silent = TRUE)
  if(class(mutFlanking) == "try-error"){
    
    print(paste("Failed to get flanking for mutation", mut[index]))
    return(NA)
  }
  matched <- try(matchORFs(mutFlanking, orfs, mutBase = gsub(" ", "", mutAllele), 
                           refBase = gsub(" ", "", mut[ref]),
                           flanking = flanking, misMatch = misMatch, 
                           strand = ifelse(gsub(" ", "", mut[strand]) == "-", "negative", "positive"), 
                           baseToMatch = baseToMatch), silent = TRUE)
  if(class(matched) == "try-error"){
    print(paste("Failed to get matched orfs for mutation", mut[index]))
    return(NA)
    
  }
  if(!all(is.na(matched))){
    
    matched <- cbind(symbol = mut[hugo], chrom = mut[chrom], start = mut[start], end = mut[end], reference = mut[ref],
                      
                     strand = mut[strand], mutation = mutAllele, cloneID = names(matched), mappedPrimer = unlist(matched),   
 
                     AA_Change = mut[AA_Change], Entrez_Gene_ID = mut[Entrez_Gene_ID], 
                     
                     Variant_Type = mut[Variant_Type], Variant_Classification = mut[Variant_Classification])
  }
  
  return(matched)    
}

### mutFlanking - object returned by running getFlanking ###
### orfs - DNAStringSet object based on orf sequencess using 
###          the makeStringSet function 
matchORFs <- function(mutFlanking, orfs, mutBase, refBase,
                      flanking = 30, misMatch = 2, strand = c("positive", "negative"), 
                      baseToMatch = 15){
  if(baseToMatch < flanking){
    cutL <- (length(mutFlanking[["left"]]) + 1 - baseToMatch):length(mutFlanking[["left"]]) 
    cutR <- 1:baseToMatch
  }else{
    cutL <- 1:length(mutFlanking[["left"]]) 
    cutR <- 1:length(mutFlanking[["right"]])
  }
  if(strand == "negative"){ 
    
    ml <- vmatchPattern(reverseComplement(mutFlanking[["left"]][cutL]), 
                        
                        orfs, max.mismatch = misMatch)
    mr <- vmatchPattern(reverseComplement(mutFlanking[["right"]][cutR]), 
                        orfs, max.mismatch = misMatch)
  }else{
    
    ml <- vmatchPattern(mutFlanking[["left"]][cutL], orfs, max.mismatch = misMatch)
    mr <- vmatchPattern(mutFlanking[["right"]][cutR], orfs, max.mismatch = misMatch) 
  }
  matched <- cbind(l = sapply(ml, function(x) length(x) != 0), 
                   r = sapply(mr, function(x) length(x) != 0))
  sums <- apply(matched, 1, sum)
  matchedOrfs <- list()
  if(any(sums == 2)){
    for(id in names(which(sums == 2))){    
      
      matchedOrfs[[id]] <- getORFPrimer(ml[[id]], mr[[id]], 
                                        seq = orfs[[id]], mutBase = mutBase, refBase = refBase, 
                                        flanking = flanking, strand = strand)
    }  
    
  }else{
    if(any(sums == 1)){
      for(id in names(which(sums == 1))){
        matchedOrfs[[id]] <- getORFPrimer(ml[[id]], mr[[id]], 
                                          seq = orfs[[id]], mutBase = mutBase, refBase = refBase, 
                                          flanking = flanking, strand = strand)
      }
    }else{
      return(NA)
    }
  }
  return(matchedOrfs)
}


### called by matchORFs
### left - DNAString object of the left flanking
### right - DNAString object of the right flanking
### refBase - reference base
getORFPrimer <- function(left, right, seq, mutBase, refBase,
                         flanking = 30, strand){
  if(strand == "negative"){
    temp <- left
    left <- right
    right <- temp
    mutBase <- as.character(reverseComplement(DNAString(mutBase)))
  }
  seqLength <- length(seq)
  if(length(left) != 0 && length(right) != 0){
    return(paste(tolower(as.character(seq[getBoundry(end(left), seqLength, 
                                                     flanking, what = "left")])),  toupper(mutBase),  
                 tolower(as.character(seq[getBoundry(start(right), seqLength, flanking, 
                                                     what = "right")])),  sep = ""))
  }
  if(refBase == "-"){
    adjust <- 1
  }else{
    adjust <- nchar(refBase) + 1
  }
  if(length(left) == 0){
    return(paste(tolower(as.character(seq[getBoundry((start(right) - adjust), 
                                                     seqLength, flanking, "left")])), 
                 toupper(mutBase),  tolower(as.character(seq[getBoundry(start(right), seqLength, 
                                                                        flanking, "right")])), sep = ""))
  }
  if(length(right) == 0){
    return(paste(tolower(as.character(seq[getBoundry(end(left), seqLength, flanking, "left")])), 
                 toupper(mutBase),  tolower(as.character(seq[getBoundry((end(left) + 
                                                                           adjust), seqLength, flanking, "right")])),  sep = ""))
  }
  return(NA)
}

getBoundry <- function(loc, segLength, flanking, 
                       what = c("left", "right")){
  what <- match.arg(what)
  if(what == "left"){
    left <- (loc - (flanking - 1)):loc
    return(left[left > 0 & left <= segLength])
  }else{
    right <- loc:(loc + (flanking - 1))
    return(right[right <= segLength])
  }
}

### returned left sequence ends before loc and right sequence
### starts baseMuted base pairs after loc 
getFlanking <- function(chrom, start, end, what, 
                        flanking = 30, organism = "Hsapiens", 
                        genome = "BSgenome.Hsapiens.UCSC.hg19", 
                        baseMuted = 1){
  
  require(genome, character.only = TRUE)
  seq <- as(get(organism, pos = grep(genome, search()))[[chrom]],
            "DNAString")
  if(toupper(what) %in% c("DEL", "SNV")){
    return(c(left = seq[(start - flanking):(start - 1)], 
             right = seq[(end + 1):(end + flanking)]))
  }
  if(toupper(what) %in% "INS"){
    return(c(left = seq[(start - flanking + 1):start], 
             right = seq[end:(end + flanking - 1)]))   
  }
  stop(paste("variant type", what, "is not supported"))
}

makeStringSet <- function(CDNASeqs, seqNames){
  seqs <- DNAStringSet(CDNASeqs)
  names(seqs) <- seqNames
  return(seqs)
}


getMutAllele <- function(ref, mutAlleles){
  return(unique(unlist(mutAlleles[!mutAlleles %in% ref])))
}


checkORFAlignment <- function(orf, genome = "BSgenome.Hsapiens.UCSC.hg19", 
                              seqCol = "ORF.sequence", symbolCol = "Gene.Symbol", misMatch = 5){
  require(genome, character.only = TRUE)
  chrom <- unique(names(geneToLoc(as.character(orf[symbolCol]), what = "symbol")))
  seq <- as(get(organism, pos = grep(genome, 
                                     search()))[[paste("chr", chrom, sep = "")]], "DNAString")
  ml <- matchPattern(as.character(orf[seqCol]), seq, max.mismatch = misMatch)
  
}

#old file paths
#variantFile <- "C:/Users/samue/OneDrive/Work/KS Lab/Protocols/PrimerGen/Request/Request.txt"
#orfFile <- "C:/Users/samue/OneDrive/Work/KS Lab/Protocols/PrimerGen/ORF_collection_v3.txt"
#outName <- "C:/Users/samue/OneDrive/Work/KS Lab/Protocols/PrimerGen/Output_primer/Output.txt"
#primerGen(variantFile, orfFile, outName)

#TESTING, looks good
#new file paths; if we specify using this way, then file names are compatable across platforms
#variantFile <- here("data", "test_GeneratePrimerSequences", "Request.txt")
#orfFile <- here("data","test_GeneratePrimerSequences",  "ORF_collection_v3.txt")
#outName <- here("output", "Output_test.tsv")

#RUN FOR ERBB2 mutations
variantFile <- here("data", "resources_erbb2", "Request_ERBB2.tsv")
orfFile <- here("data","resources_ORF",  "ORF_collection_v3.txt")
outName <- here("output", "Output_ERBB2_Test2.tsv")

#samuel ->run this line to create primer sequences and output file :)
#run 10/29/20
primerGen(variantFile, orfFile, outName)





