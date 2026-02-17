#Lab 6

#install and load
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msa")
library(msa)

#other packages needed for the lab
BiocManager::install("Biostrings")
install.packages(c("seqinr", "phangorn"))

library(Biostrings)
library(seqinr)
library(phangorn)

#steps 7 and 8: read each FASTA separately and then combine

files <- sprintf("Sequence_%d.fasta", 1:5)
stopifnot(all(file.exists(files)))

dna_list <- lapply(files, readDNAStringSet)
dna_all  <- do.call(c, dna_list)

#prints sequences and then renames them to easier labels (S1–S5)
dna_all
names(dna_all) <- paste0("S", seq_along(dna_all))
dna_all


#step 9: MUSCLE

aln <- msa(dna_all, method = "Muscle")

#step 10: gaps

print(aln, show = "complete")

#makes the  alignment a matrix so bases can be counted
aligned_set <- as(aln, "DNAStringSet")
aln_mat <- as.matrix(aligned_set)

gap_count <- sum(aln_mat == "-")
gap_count


#step 11: alignment length


aln_length <- width(aligned_set)[1]
aln_length


#step 12: GC content but without exclude gaps

tab <- table(aln_mat)
G <- if ("G" %in% names(tab)) tab[["G"]] else 0
C <- if ("C" %in% names(tab)) tab[["C"]] else 0
A <- if ("A" %in% names(tab)) tab[["A"]] else 0
T <- if ("T" %in% names(tab)) tab[["T"]] else 0

#figures out number of gc bases excluding gaps
total_no_gaps <- A + T + G + C
gc_percent <- 100 * (G + C) / total_no_gaps
gc_percent


## Step 13: seqinr format + distance matrix


aln_seqinr <- msaConvert(aln, type = "seqinr::alignment")
dist_mat <- dist.alignment(aln_seqinr, "identity")
dist_mat


## Step 14: closest + farthest


##converts the distance object to matrix and removes self-comparisons so that self-distances (0) don’t affect results
D <- as.matrix(dist_mat)
diag(D) <- NA

closest_idx  <- which(D == min(D, na.rm = TRUE), arr.ind = TRUE)[1,]
farthest_idx <- which(D == max(D, na.rm = TRUE), arr.ind = TRUE)[1,]

closest  <- c(rownames(D)[closest_idx[1]],  colnames(D)[closest_idx[2]],  D[closest_idx[1],  closest_idx[2]])
farthest <- c(rownames(D)[farthest_idx[1]], colnames(D)[farthest_idx[2]], D[farthest_idx[1], farthest_idx[2]])

closest
farthest


#step 15: translates one sequence


one_seq <- aligned_set[[1]]

#remove gaps
seq_raw <- gsub("-", "", as.character(one_seq))

#replace anything not A/C/G/T with N
seq_clean <- toupper(gsub("[^ACGT]", "N", seq_raw))

#trim to multiple of 3 (i had to do this because it kept giving me errors until I did this)
trim_len <- nchar(seq_clean) - (nchar(seq_clean) %% 3)
coding_seq <- substr(seq_clean, 1, trim_len)

#translate - if theres any ambiguous codons, translate puts an X
aa_seq <- as.character(Biostrings::translate(DNAString(coding_seq), if.fuzzy.codon = "X"))

substr(aa_seq, 1, 80)





## Step 16: write alignment through phangorn

phy <- msaConvert(aln, type = "phangorn::phyDat")
write.phyDat(phy, file = "alignment.fasta", format = "fasta")
