
##### ----- CREATE FILES FOR ON-TARGET SCORING
# need 30 bp seq described below, compared to 20bp seq
### do this manually, step-by-step

# get seq for region of interest
y = getSeq(Hsapiens, 'chr8', 127735000, 130750000)

# whole genome
#z = getSeq(Hsapiens)

# assign PAM pattern
patternNGG = DNAString('NGG')

# find locations of PAM sequence on the region of interest
NGGforwardmatchfixed = matchPattern(pattern = patternNGG, subject = y, max.mismatch = 0, fixed = FALSE)
# and - strand
CCNforwardmatchfixed = matchPattern(pattern = reverseComplement(patternNGG), subject = y, max.mismatch = 0, fixed = FALSE)

#assign the start and end columns of the NGG and CCN DNAStrings to vectors
NGGstart = start(NGGforwardmatchfixed)+127735000  # 12769999
NGGstop = end(NGGforwardmatchfixed)+127735000
CCNstart = start(CCNforwardmatchfixed)+127735000
CCNend = end(CCNforwardmatchfixed)+127735000


#### FOR ON-TARGET SCORES COMPUTED BY "compute_on_target_scores_MM.py", NEED A 30MER SEQUENCE WHICH IS NNNN20MERNGGNNN
#######  ON TARGET SCORES
# NEED 4 BASES 5' THE 20MER AND 3 BASES 3' OF THE LAST G IN NGG
######## CONFIRM THESE GUIDE START AND STOP POSITIONS
ONTARGETsgRNAstart_NGG = NGGstart - 25
ONTARGETsgRNAstop_NGG = NGGstart + 4
ONTARGETnegative_strand_sgRNAstart_NGG = CCNstart - 4
ONTARGETnegative_strand_sgRNAend_NGG = CCNend + 23
# NOTE THAT NEGATIVE STRAND SGRNA NEEDS TO BE REVERSE COMPLEMENTED WHEN GETSEQ IS DONE

# create new libraries with 30mers
##### actually get sgRNA sequences from 5-3' NGG on + strand
ONTARGETtopstrand_sgRNA_NGG_library = DNAStringSet(getSeq(Hsapiens, 'chr8', ONTARGETsgRNAstart_NGG, ONTARGETsgRNAstop_NGG), use.names = F)
head(ONTARGETtopstrand_sgRNA_NGG_library)
length(ONTARGETtopstrand_sgRNA_NGG_library)

ONTARGETbottomstrand_sgRNA_NGG_library = DNAStringSet(getSeq(Hsapiens, 'chr8', ONTARGETnegative_strand_sgRNAstart_NGG, ONTARGETnegative_strand_sgRNAend_NGG), use.names = F)
# TAKE THE REVERSE COMPLEMENT OF THESE CCN SGRNAS
ONTARGETbottomstrand_sgRNA_NGG_library = reverseComplement(ONTARGETbottomstrand_sgRNA_NGG_library)
length(ONTARGETbottomstrand_sgRNA_NGG_library)


# save the file so that Doench 2016 algorithm can compute on-target score for each
# WRITE THE NGG BASED SGRNA LIBRARY TO A DF
ONTARGETNGG_based_library_with_positions = data.frame(chr = rep("chr8", length(ONTARGETtopstrand_sgRNA_NGG_library)),  strand = rep("+", times = length(ONTARGETtopstrand_sgRNA_NGG_library)), sgRNA_start = ONTARGETsgRNAstart_NGG, sgRNA_end = ONTARGETsgRNAstop_NGG, sgRNA_sequence = ONTARGETtopstrand_sgRNA_NGG_library)
head(ONTARGETNGG_based_library_with_positions)
nrow(ONTARGETNGG_based_library_with_positions)

# WRITE THE CCN BASED SGRNA LIBRARY
ONTARGETCCN_based_library_with_positions = data.frame(chr = rep("chr8", times = length(ONTARGETbottomstrand_sgRNA_NGG_library)), strand = rep("-", times = length(ONTARGETbottomstrand_sgRNA_NGG_library)), sgRNA_start = ONTARGETnegative_strand_sgRNAstart_NGG, sgRNA_end = ONTARGETnegative_strand_sgRNAend_NGG, sgRNA_sequence = ONTARGETbottomstrand_sgRNA_NGG_library)
head(ONTARGETCCN_based_library_with_positions)
nrow(ONTARGETCCN_based_library_with_positions)


# write the top strand 30mers LIBRARIES to file
on_target_top_strand_only_sequence = data.frame(ONTARGETtopstrand_sgRNA_NGG_library)
head(on_target_top_strand_only_sequence)
write.table(on_target_top_strand_only_sequence, file = "on_target_scores_NGG_based_library.bed", sep = "\t", row.names = F, col.names = F, quote = F)

# write the bottom strand 30mers to file
on_target_bottom_strand_only_sequence = data.frame(ONTARGETbottomstrand_sgRNA_NGG_library)
head(on_target_bottom_strand_only_sequence)
write.table(on_target_bottom_strand_only_sequence, file = "on_target_scores_CCN_based_library.bed", sep = "\t", row.names = F, col.names = F, quote = F)




### WORK ON MAKING THIS JUMP TO PYTHON FOR ON-TARGET SCORING SEAMLESS ------
############ FILTER THE RAW SGRNAS, AND SCORE THEM
# compute the on-target score of sgRNA using Doench python script
# RUN MY PYTHON SCRIPT "compute_on_target_scores_MM.py" to create "on_target_NGG_top_strand_seq_and_score.txt" file
cat('About to start computing on-target scores!')
# maybe separate function for reading in 30mers files, executing the python script, and filtering

# takes two inputs, the 30mer seq '\t' on-target scores for two files - 1. top and 2. bottom strand
# outputs a filtered df and file

### RUN THIS COMMAND ON COMMAND LINE
# -t for top strand input file, -b for bottom strand input file, -x for top strand 30mer seq and score, -y for bottom strand 30mer seq and score
system(paste('python /Users/meer/Dropbox/SLab histone/CRISPR_MYC_library/Code/Rule_Set_2_scoring_v1/TEST/compute_on_target_scores_MM.py 
             -t on_target_scores_NGG_based_library.bed -b on_target_scores_CCN_based_library.bed 
             -x on_target_NGG_top_strand_seq_and_score.txt 
             -y on_target_CCN_bottom_strand_seq_and_score.txt'))













# move onto the next R script once you have scored the guides --------






#### ARCHIVED AFTER THIS POINT ------ 






# hard code that's put into the function up top, kept here for archive
##### ----- CREATE RAW SGRNA LIBRARY
######## start of HARD CODE

library(BSgenome.Hsapiens.UCSC.hg19)
x = getSeq(Hsapiens, "chr8")
y = getSeq(Hsapiens, "chr8", 127735000, 130750000)

# whole genome
z = getSeq(Hsapiens)

x[127700015:127700025]
y[15:25]



# assign pattern
patternNGG = DNAString('NGG')

# assign
NGGforwardmatchfixed = matchPattern(pattern = patternNGG, subject = y, max.mismatch = 0, fixed = FALSE)
CCNforwardmatchfixed = matchPattern(pattern = reverseComplement(patternNGG), subject = y, max.mismatch = 0, fixed = FALSE)
NGGforwardmatchfixed
CCNforwardmatchfixed




#assign the start and end columns of the NGG and CCN DNAStrings to vectors
NGGstart = start(NGGforwardmatchfixed)+127735000  # 12769999
NGGstop = end(NGGforwardmatchfixed)+127735000
CCNstart = start(CCNforwardmatchfixed)+127735000
CCNend = end(CCNforwardmatchfixed)+127735000

head(NGGstart)
length(NGGstart)
head(NGGstop)
length(NGGstop)
# 147788 NGG on positive strand

head(CCNstart)
length(CCNstart)
head(CCNend)
length(CCNend)

# sgRNA start is 18 positions upstream of the N inside NGG
######## CONFIRM THESE GUIDE START AND STOP POSITIONS
sgRNAstart_NGG = NGGstart - 21
sgRNAstop_NGG = NGGstart - 2
negative_strand_sgRNAstart_NGG = CCNstart + 2
negative_strand_sgRNAend_NGG = CCNend + 19
# NOTE THAT NEGATIVE STRAND SGRNA NEEDS TO BE REVERSE COMPLEMENTED WHEN GETSEQ IS DONE\





head(sgRNAstart_NGG)
length(sgRNAstart_NGG)
head(sgRNAstop_NGG)
length(sgRNAstop_NGG)

head(negative_strand_sgRNAstart_NGG)
head(negative_strand_sgRNAend_NGG)


##### actually get sgRNA sequences from 5-3' NGG on + strand
topstrand_sgRNA_NGG_library = DNAStringSet(getSeq(Hsapiens, 'chr8', sgRNAstart_NGG, sgRNAstop_NGG), use.names = F)
topstrand_sgRNA_NGG_library[0:10]


bottomstrand_sgRNA_NGG_library = DNAStringSet(getSeq(Hsapiens, 'chr8', negative_strand_sgRNAstart_NGG, negative_strand_sgRNAend_NGG), use.names = F)
bottomstrand_sgRNA_NGG_library[0:10]
# TAKE THE REVERSE COMPLEMENT OF THESE CCN SGRNAS
bottomstrand_sgRNA_NGG_library = reverseComplement(bottomstrand_sgRNA_NGG_library)
bottomstrand_sgRNA_NGG_library[0:20]




# WRITE THE NGG BASED SGRNA LIBRARY
NGG_based_library_with_positions = data.frame(chr = rep("chr8", length(topstrand_sgRNA_NGG_library)),  strand = rep("+", times = length(topstrand_sgRNA_NGG_library)), sgRNA_start = sgRNAstart_NGG, sgRNA_end = sgRNAstop_NGG, sgRNA_sequence = topstrand_sgRNA_NGG_library)
head(NGG_based_library_with_positions)
nrow(NGG_based_library_with_positions)

# WRITE THE CCN BASED SGRNA LIBRARY
CCN_based_library_with_positions = data.frame(chr = rep("chr8", times = length(bottomstrand_sgRNA_NGG_library)), strand = rep("-", times = length(bottomstrand_sgRNA_NGG_library)), sgRNA_start = negative_strand_sgRNAstart_NGG, sgRNA_end = negative_strand_sgRNAend_NGG, sgRNA_sequence = bottomstrand_sgRNA_NGG_library)
CCN_based_library_with_positions[1:20,]
nrow(CCN_based_library_with_positions)


write.table(NGG_based_library_with_positions, file = "NGG_based_library.bed", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(CCN_based_library_with_positions, file = "CCN_based_library.bed", sep = "\t", row.names = F, col.names = T, quote = F)




# can try to integrate the two
combined_top_and_bottom_strand_library = rbind(NGG_based_library_with_positions, CCN_based_library_with_positions)
# and sort by sgRNA start
combined_sorted_top_and_bottom_library = combined_top_and_bottom_strand_library[order(combined_top_and_bottom_strand_library$sgRNA_start),]
# save output
write.table(combined_sorted_top_and_bottom_library, file = "combined_top_bottom_strand_library.bed", sep = "\t", row.names = F, col.names = T, quote = F)









# THIS WAS DONE EARLIER
#######  ON TARGET SCORES
# NEED 4 BASES 4' THE 20MER AND 3 BASES 3' OF THE LAST G IN NGG
######## CONFIRM THESE GUIDE START AND STOP POSITIONS
ONTARGETsgRNAstart_NGG = NGGstart - 25
ONTARGETsgRNAstop_NGG = NGGstart + 4
ONTARGETnegative_strand_sgRNAstart_NGG = CCNstart - 4
ONTARGETnegative_strand_sgRNAend_NGG = CCNend + 23
# NOTE THAT NEGATIVE STRAND SGRNA NEEDS TO BE REVERSE COMPLEMENTED WHEN GETSEQ IS DONE




head(ONTARGETsgRNAstart_NGG)
length(ONTARGETsgRNAstart_NGG)
head(ONTARGETsgRNAstop_NGG)
length(ONTARGETsgRNAstop_NGG)

head(ONTARGETnegative_strand_sgRNAstart_NGG)
head(ONTARGETnegative_strand_sgRNAend_NGG)

x[ONTARGETsgRNAstart_NGG[1]:ONTARGETsgRNAstop_NGG[1]]
reverseComplement(x[ONTARGETnegative_strand_sgRNAstart_NGG[1]:ONTARGETnegative_strand_sgRNAend_NGG[1]])






# READ IN THE 30MER SCORED SGRNAS AND ATTACH ONTO THE EXISTING SGRNA LIBRARY FILE ------
# try again with unfiltered
on_target_scores_top = read.table('on_target_NGG_top_strand_seq_and_score.txt', sep =  '\t', col.names = c('sgRNA_30mer', 'on_target_score'))
head(on_target_scores_top)
nrow(on_target_scores_top)
nrow(NGG_based_library_with_positions)
head(NGG_based_library_with_positions)
tail(NGG_based_library_with_positions)
full_top_w_score = cbind(NGG_based_library_with_positions[1:nrow(NGG_based_library_with_positions)-1, ] , on_target_score = on_target_scores_top$on_target_score)
head(full_top_w_score)
tail(full_top_w_score)


# bottom strand
on_target_scores_bottom = read.table('on_target_CCN_bottom_strand_seq_and_score.txt', sep =  '\t', col.names = c('sgRNA_30mer', 'on_target_score'))
head(on_target_scores_bottom)
tail(on_target_scores_bottom)
nrow(on_target_scores_bottom)
nrow(CCN_based_library_with_positions)
head(CCN_based_library_with_positions)
tail(CCN_based_library_with_positions)
full_bottom_w_score = cbind(CCN_based_library_with_positions[1:nrow(CCN_based_library_with_positions)-1, ] , on_target_score = on_target_scores_bottom$on_target_score)
head(full_bottom_w_score)
tail(full_bottom_w_score)

# concatenate the two score containing dfs
combined_sgRNA_and_on_target_scores = rbind(full_top_w_score, full_bottom_w_score)
head(combined_sgRNA_and_on_target_scores)

# and sort by sgRNA start
combined_sgRNA_and_on_target_scores = combined_sgRNA_and_on_target_scores[order(combined_sgRNA_and_on_target_scores$sgRNA_start),]
head(combined_sgRNA_and_on_target_scores)







# TAKE OUT NEGATIVE VALUES FROM THE -18 SUBTRACTION ABOVE^^^

# check to see ifthey merged in order, VERY IMPORTANT
# may NEED TO SORT THE VECTORS IF THEY DIDN'T



# verify by looking at both strands of DNA 
# do this again but with REVERSE COMPLEMENT OF THE GETSEQ


















##### ARCHIVE-----



############################ filtering sgRNA MYC library
# Meer Mustafa Nov. 1, 2016



# for removing sgRNA with >5 Ns and >4 Ts
letterFrequency(DNAstringset, "T")





# read in table of sgRNAs under each peak
sgRNAlibraryPEAKS = read.table("sgRNAs_MYC_domain_peak_targets_Database.txt", sep = "\t", header = TRUE)
sgRNAlibraryNONPEAKS = read.table("sgRNAs_MYC_domain_nonpeak_targets_Database.txt", sep = "\t", header = TRUE)

# how many guides in each category
nrow(sgRNAlibraryPEAKS)
nrow(sgRNAlibraryNONPEAKS)

# total sgRNA library size
nrow(sgRNAlibraryPEAKS) + nrow(sgRNAlibraryNONPEAKS)

# percent of library that targets PEAKS
(nrow(sgRNAlibraryPEAKS) / (nrow(sgRNAlibraryPEAKS) + nrow(sgRNAlibraryNONPEAKS))) * 100



# FILTER RAW PEAK SGRNAS BY PILEUP HEIGHT THRESHOLD
threshold10_pileup = subset(sgRNAlibraryPEAKS, sgRNAlibraryPEAKS$peak.pileup.at.peak.summit > 10)
nrow(nonpeak_threshold10_pileup)
# gives 71K guides

threshold20_pileup = subset(sgRNAlibraryPEAKS, sgRNAlibraryPEAKS$peak.pileup.at.peak.summit > 20)
nrow(nonpeak_threshold20_pileup)
# gives 62K guides

nonpeak_threshold100_pileup = subset(sgRNAlibraryPEAKS, sgRNAlibraryPEAKS$peak.pileup.at.peak.summit > 100)
nrow(nonpeak_threshold100_pileup)
# gives 18K guides













######## METHOD 1
# HOW TO MAKE THE DNASTRINGSET WITH THE CORRECT COORDINATES

sgRNAseq = NULL
i = 0
#USE DNASTRING SET TO HOLD MUTLIPLE DNASTRING
for (i in range(i:length(newsgRNAstart))) {
  # take the ith vector of the newsgRNA start and stops
  sgRNAseq[[i]] = y[newsgRNAstart:newsgRNAstop]
  
  i = i + 1
}




# MTHOD 2 DNASTRINGSET TO GET LIBRARY, THEN DICTIONARY IT
sgRNANGGlibrary = DNAStringSet(y[newsgRNAstart:newsgRNAstop])


# TAKES TOO LONG TO PROCESS
# readjust coordinates to add the positions to give MYC domain coordinates
sgRNANGGlibrary = DNAStringSet(getSeq(Hsapiens, 'chr8', sgRNAstart+127700000, sgRNAstop+130700000))




for (i in newsgRNAstart) {
  DNAStringSet(y[c(newsgRNAstart) : c(newsgRNAstop)])
  i+1
}








##### finding overlap of sgRNA and peaks

peaksgRNAmatchup = NULL

# counter for bed files
i = 1

#**** confirm that BED files are in the working directory in R
#**** confirm that BED files are MYC domain only (chr8)
for (f in list.files (, pattern="(.*).bed") ) {
  
  #TEST INDIVIDUAL FILE
  #f = "A549dhs_peaks_myc.bed"
  
  # read in all matched BED files as bed. tab delimited. no header on BED file.
  bed = read.table(f, sep = "\t", header = FALSE)
  
  
  for (startposition in newsgRNAstart) {if (startposition %in% range(bed$V2:bed$V3)) {cat("found match!")}}
  
  
  for (line in bed) {
    # for each line in the currently loaded BED file:
    for (startposition in newsgRNAstart) {
      
      # if (start is inside the range of bed$V2:bed$V3) {
      
      if (startposition %in% range(bed$V2:bed$V3)) {
        # create new vector that holds the vector of pileups or 0s
        peaksgRNAmatchup[[i]] = bed$V6   # which is the pileup value
        #else (start in range(bed$V2:bed$V3))
      } else {
        peaksgRNAmatchup[[i]] = 0
      }
      
      
      #else (start in range(bed$V2:bed$V3)) {
      #peaksgRNAmatchup[[i]] = 0
      # else (start is not inside range of bed$V2:bed$V3) {
      # peaksgRNAmatchup[[i]] = 0
      
      # add 1 to start position to go onto the next sgRNA position
      i = i + 1
    }
    
    # underneath the if/else, should be appending either pileup or 0 onto a vector
    
  }
  
  
  
  # so the output desired variable is peaksgRNAmatchup, a vector
  
  
  # add one to the counter to go to the next BED file
  i = i + 1
  
}

