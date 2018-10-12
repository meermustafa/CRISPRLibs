
############### CREATE FILES FOR ON-TARGET SCORING ###############
# To use Doench 2014 algorithm, need 30 bp seq described below, compared to 20bp seq

# get seq for region of interest
y = getSeq(Hsapiens, chr, start, end)

# whole genome
#z = getSeq(Hsapiens)

# assign PAM pattern
patternNGG = DNAString('NGG')

# find locations of PAM sequence on the region of interest
NGGforwardmatchfixed = matchPattern(pattern = patternNGG, subject = y, max.mismatch = 0, fixed = FALSE)
# and - strand
CCNforwardmatchfixed = matchPattern(pattern = reverseComplement(patternNGG), subject = y, max.mismatch = 0, fixed = FALSE)

#assign the start and end columns of the NGG and CCN DNAStrings to vectors
NGGstart = start(NGGforwardmatchfixed)
NGGstop = end(NGGforwardmatchfixed)
CCNstart = start(CCNforwardmatchfixed)
CCNend = end(CCNforwardmatchfixed)


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
system(paste('python compute_on_target_scores_MM.py 
             -t on_target_scores_NGG_based_library.bed -b on_target_scores_CCN_based_library.bed 
             -x on_target_NGG_top_strand_seq_and_score.txt 
             -y on_target_CCN_bottom_strand_seq_and_score.txt'))
