#Q1: Import and align DNA sequences
#load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

#install Bioconductor packages if needed
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
if (!requireNamespace("msa", quietly = TRUE))        BiocManager::install("msa")

#load libraries into the session
library(Biostrings)
library(msa)

#import sequences
dna <- readDNAStringSet("sequences.fasta")

length(dna)
dna  #shows a preview of headers + sequence lengths

#align sequences
#align with one method, I used ClustalW (ClustalW or Muscle)
aln <- msa(dna, method = "ClustalW")

#view the alignment object
aln

#convert the alignment to DNAStringSet to later do the consensus and GC count
aligned_set <- as(aln, "DNAStringSet")

#confirm alignment width so that all of the sequences should now have same length
width(aligned_set)

#Q2: Measure alignment quality
#convert alignment to a character matrix
aln_mat <- as.matrix(aligned_set)

#look a the overall gap proportion
total_cells <- length(aln_mat)
gap_cells   <- sum(aln_mat == "-")
gap_prop    <- gap_cells / total_cells
gap_prop

#look at the gaps by column
gap_by_col      <- colSums(aln_mat == "-")
gap_prop_by_col <- gap_by_col / nrow(aln_mat)

#sees how many columns are gappy
summary(gap_prop_by_col)

#tells us how many columns have >20% gaps
cols_over_20pct_gaps <- sum(gap_prop_by_col > 0.20)
cols_over_20pct_gaps

#average pairwise identity
#an identity distance of 0 means identical and larger means more different
library(seqinr)
aln_seqinr <- msaConvert(aln, type = "seqinr::alignment")
dist_mat   <- dist.alignment(aln_seqinr, "identity")

D <- as.matrix(dist_mat)
diag(D) <- NA  #ignores self-comparisons

avg_dist <- mean(D, na.rm = TRUE) #average identity distance across pairs
avg_id   <- 1 - avg_dist #converts to average identity proportion
avg_dist
avg_id

#Answer: The overall gap proportion is very low at 0.000078% so that means there
#was practically nothing added in order to properly align the sequences. The alignment
#columns had nothing more than >20% gaps meaning that there were no badly aligned sections.
#Also the average pairwise identity score was 0.981/98% so based on all of that the alignment is very good

#Q3: Calculate the consensus sequence for the alignment
#converts the alignment to matrix
aln_mat <- as.matrix(aligned_set)

#consensus per column:
#ignore gaps ("-")
#ignore anything not A/C/G/T (like N or other symbols)
cons_vec <- apply(aln_mat, 2, function(col) {
  col <- toupper(col)
  col <- col[col != "-"]
  col <- col[col %in% c("A","C","G","T")]
  
  if (length(col) == 0) return("N")  #sees if the column has no valid bases then it will mark it as N
  names(which.max(table(col)))
})

consensus_no_gaps <- paste0(cons_vec, collapse = "")
consensus_no_gaps

#confirms the length being 642
nchar(consensus_no_gaps)

#writes consensus to FASTA
writeLines(
  c(">consensus_sequence", consensus_no_gaps),
  "consensus_sequence.fasta"
)

#Q4: Calculate GC content for the overall alignment
#makes alignment matrix into one long vector of bases
all_chars <- toupper(as.vector(aln_mat))

#keeps only valid DNA bases (ignore gaps and any other symbols)
valid_bases <- all_chars[all_chars %in% c("A","C","G","T")]

#counts the bases
A_count <- sum(valid_bases == "A")
C_count <- sum(valid_bases == "C")
G_count <- sum(valid_bases == "G")
T_count <- sum(valid_bases == "T")

total_bases <- A_count + C_count + G_count + T_count

#calculates the GC percentage
gc_percent <- 100 * (G_count + C_count) / total_bases

gc_percent

#Answer: GC percent is 51% and 

#Q5: Compares individuals finds any outliers and describes mutations
#which sequences are most different?
#farthest pair/largest identity distance
farthest_idx <- which(D == max(D, na.rm = TRUE), arr.ind = TRUE)[1, ]
sample1 <- rownames(D)[farthest_idx[1]]
sample2 <- colnames(D)[farthest_idx[2]]
max_dist <- D[farthest_idx[1], farthest_idx[2]]

sample1; sample2; max_dist

#most different overall/what has the highest mean distance to all others
mean_dist <- rowMeans(D, na.rm = TRUE)
most_diff_sample <- names(which.max(mean_dist))
most_diff_sample
mean_dist[most_diff_sample]

#what kinds of mutations does the most different individual have?
#compares "most_diff_sample" to the alignment consensus (cons_vec that was already made in Q3)
sample_row <- toupper(aln_mat[most_diff_sample, ])
ref_row    <- cons_vec

#positions where sample differs from consensus and ignores unknown/ref N
diff_pos <- which(sample_row != ref_row & ref_row != "N")

#classifies the differences
is_transition <- function(ref, alt) {
  (ref == "A" & alt == "G") | (ref == "G" & alt == "A") |
    (ref == "C" & alt == "T") | (ref == "T" & alt == "C")
}

mut_type <- character(length(diff_pos))

for (k in seq_along(diff_pos)) {
  i <- diff_pos[k]
  ref <- ref_row[i]
  alt <- sample_row[i]
  
#ignores any unexpected symbols
  if (!(ref %in% c("A","C","G","T","-"))) { mut_type[k] <- "Other/ambiguous"; next }
  if (!(alt %in% c("A","C","G","T","-"))) { mut_type[k] <- "Other/ambiguous"; next }
  
  if (ref == "-" & alt %in% c("A","C","G","T")) {
    mut_type[k] <- "Insertion"
  } else if (ref %in% c("A","C","G","T") & alt == "-") {
    mut_type[k] <- "Deletion"
  } else if (ref %in% c("A","C","G","T") & alt %in% c("A","C","G","T")) {
    mut_type[k] <- if (is_transition(ref, alt)) "SNP (transition)" else "SNP (transversion)"
  } else {
    mut_type[k] <- "Other/ambiguous"
  }
}

mutations <- data.frame(
  alignment_position = diff_pos,
  consensus_base = ref_row[diff_pos],
  sample_base = sample_row[diff_pos],
  mutation_type = mut_type,
  row.names = NULL
)

mutations
table(mutations$mutation_type)

#Answer: There is 1 deletion and 7 substitions.

# Q6: Identify the gene by comparing sequence to a database
#for BLAST uses a gapless DNA sequence (consensus is already gapless)
blast_query <- consensus_no_gaps

# Write to FASTA for BLAST
writeLines(
  c(">consensus_query_for_BLAST", blast_query),
  "Q6_consensus_query_for_BLAST.fasta"
)

nchar(blast_query)

#Q6 answer
#I BLASTed the consensus DNA sequence against the NCBI nucleotide collection (nt)
#the best match was the Homo sapiens HBB gene hemoglobin beta
#the accession number of the top hit was LC121775.1
#the match showed 100% query coverage, 100% identity and an E-value of 0.0
#showing a perfect alignment to the human beta-globin gene.


#Q7: Finds the most different individual, translates the sequence, and then writes protein FASTA
#I did not re-load packages here because Biostrings was already loaded in Q1.
#and im not re-calculating distances because D and mean_dist were already created in Q5

#most_diff_sample already exists from Q5 but this prints it again just to keep things clear
most_diff_sample
mean_dist[most_diff_sample]

#extract the most different individual's aligned DNA sequence from the alignment matrix
most_diff_aligned <- paste0(aln_mat[most_diff_sample, ], collapse = "")

#remove gaps to get the unaligned DNA sequence
most_diff_no_gaps <- gsub("-", "", most_diff_aligned)

#clean sequence to only A/C/G/T (prevents translation errors if any weird characters exist)
most_diff_clean <- toupper(gsub("[^ACGT]", "", most_diff_no_gaps))

nchar(most_diff_clean)

#trims to a multiple of 3 so that the translation stays in frame
trim_len <- nchar(most_diff_clean) - (nchar(most_diff_clean) %% 3)
coding_seq <- substr(most_diff_clean, 1, trim_len)

nchar(coding_seq)

#translate DNA -> protein (if any not well defined codons somehow make it here they become "X")
protein_seq <- as.character(
  Biostrings::translate(DNAString(coding_seq), if.fuzzy.codon = "X")
)

#quick preview
substr(protein_seq, 1, 80)

#write the protein sequence to a FASTA file
writeLines(
  c(paste0(">", most_diff_sample, "_translated_protein"), protein_seq),
  "Q7_most_different_protein.fasta"
)

#Question 8 answer: The top match was to hemoglobin subunit beta HBB and the hit
#with the lowest E value had an accession number of KAI2558340.1.

#Question 9 answer: Per information listed on OMIM, the common diseases shown were
#Sickle cell disease and Beta-thalassemia both of which cause mutations in the HBB gene.
#Even though the information found in question 7 showed some SNP differences, the BLAST
#result was able to show that it matched a normal human beta-globin so based on all of
#that there is nothing to show that this person from question 7 has the disease.