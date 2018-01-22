
# R script for making CRISPR/Cas9 library
# Meer Mustafa, 11.2.16


# need to install the genome string once


# perhaps do a check: check to see if packages are already installed
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)


# input: "chr#", chrstart, chrend
# output is one file/table containing the sgRNA sequences, their location, strand, and on-target score

create_Cas9_library = function (PAM, chr, chrstart, chrend) {
  # load library
  # * may need to install BSgenome
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  # set sequence for entire chromosome
  x = getSeq(Hsapiens, chr)
  
  # get seq for region of interest
  y = getSeq(Hsapiens,chr, chrstart, chrend)
  
  # whole genome
  #z = getSeq(Hsapiens)
  
  # assign PAM pattern
  patternNGG = DNAString(PAM)
  
  # find locations of PAM sequence on the region of interest
  NGGforwardmatchfixed = matchPattern(pattern = patternNGG, subject = y, max.mismatch = 0, fixed = FALSE)
  # and - strand
  CCNforwardmatchfixed = matchPattern(pattern = reverseComplement(patternNGG), subject = y, max.mismatch = 0, fixed = FALSE)
  
  #assign the start and end columns of the NGG and CCN DNAStrings to vectors
  NGGstart = start(NGGforwardmatchfixed)+chrstart
  NGGstop = end(NGGforwardmatchfixed)+chrstart
  CCNstart = start(CCNforwardmatchfixed)+chrstart
  CCNend = end(CCNforwardmatchfixed)+chrstart
  
  # verify that pattern starts and stops are the same length, they should be
  if (length(NGGstart) != length(NGGstop)) {stop}
  
  # get the sgRNA positions from the PAM positions
  ######## CONFIRM THESE GUIDE START AND STOP POSITIONS
  sgRNAstart_NGG = NGGstart - 21
  sgRNAstop_NGG = NGGstart - 2
  negative_strand_sgRNAstart_NGG = CCNstart + 2
  negative_strand_sgRNAend_NGG = CCNend + 19
  # NOTE THAT NEGATIVE STRAND SGRNA NEEDS TO BE REVERSE COMPLEMENTED WHEN GETSEQ IS DONE
  
  
  ##### get sgRNA sequences from 5-3' NGG on + strand
  topstrand_sgRNA_NGG_library = DNAStringSet(getSeq(Hsapiens, chr, sgRNAstart_NGG, sgRNAstop_NGG), use.names = F)
  topstrand_sgRNA_NGG_library[0:10]
  
  # and sgRNAs on the - strand
  bottomstrand_sgRNA_NGG_library = DNAStringSet(getSeq(Hsapiens, chr, negative_strand_sgRNAstart_NGG, negative_strand_sgRNAend_NGG), use.names = F)
  bottomstrand_sgRNA_NGG_library[0:10]
  # TAKE THE REVERSE COMPLEMENT OF THESE - strand, CCN SGRNAS
  bottomstrand_sgRNA_NGG_library = reverseComplement(bottomstrand_sgRNA_NGG_library)
  bottomstrand_sgRNA_NGG_library[0:10]
  
  
  # Put the library of sgRNA sequences, strands, positions into a df and write to file
  NGG_based_library_with_positions = data.frame(chr = rep(chr, 
                                                          length(topstrand_sgRNA_NGG_library)),  
                                                strand = rep("+", times = length(topstrand_sgRNA_NGG_library)), 
                                                sgRNA_start = sgRNAstart_NGG, 
                                                sgRNA_end = sgRNAstop_NGG, 
                                                sgRNA_sequence = topstrand_sgRNA_NGG_library)
  head(NGG_based_library_with_positions)
  nrow(NGG_based_library_with_positions)
  
  # and - strand
  CCN_based_library_with_positions = data.frame(chr = rep(chr, times = length(bottomstrand_sgRNA_NGG_library)), 
                                                strand = rep("-", times = length(bottomstrand_sgRNA_NGG_library)), 
                                                sgRNA_start = negative_strand_sgRNAstart_NGG, 
                                                sgRNA_end = negative_strand_sgRNAend_NGG, 
                                                sgRNA_sequence = bottomstrand_sgRNA_NGG_library)
  head(CCN_based_library_with_positions)
  nrow(CCN_based_library_with_positions)
  
  # write the separated bed
  #write.table(NGG_based_library_with_positions, file = "NGG_based_library.bed", sep = "\t", row.names = F, col.names = T, quote = F)
  #write.table(CCN_based_library_with_positions, file = "CCN_based_library.bed", sep = "\t", row.names = F, col.names = T, quote = F)
  
  
  # can try to integrate the two libraries to become one, while maintaining info about strandedness
  combined_top_and_bottom_strand_library = rbind(NGG_based_library_with_positions, CCN_based_library_with_positions)
  # and sort by sgRNA start
  combined_sorted_top_and_bottom_library = combined_top_and_bottom_strand_library[order(combined_top_and_bottom_strand_library$sgRNA_start),]
  
  # assign it to variable to be saved outside this function
  assign(x = 'combined_top_and_bottom_strand_library',
        combined_sorted_top_and_bottom_library,
        envir = globalenv())
  
  # save output
  #write.table(combined_sorted_top_and_bottom_library, file = "combined_top_bottom_strand_library.txt", sep = "\t", row.names = F, col.names = T, quote = F)
  
  # save BED file for viewing
  #write.table(combined_sorted_top_and_bottom_library[,c(1,3,4)], file = "combined_top_bottom_strand_library.bed", sep = "\t", row.names = F, col.names = F, quote = F)
  

  
}


# run the function which will write the files
create_Cas9_library('NGG',chr = 'chr8', 127900000, 130900000)

str(combined_top_and_bottom_strand_library)



# expansion ideas: 

# different CRISPR (or even other family) systems
# can specific SpCas9, SaCas9, LbCpf1, AsCpf1, etc

# different genomes
# if statement flags to check for genomes (e.g. hg19, mm10).
# if they don't exist, install and library those genomes














