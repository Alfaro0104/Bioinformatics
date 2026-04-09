#Installs BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

#Installs Biostrings if needed
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings", ask = FALSE, update = FALSE)
}

#Loads in the biostrings package
library(Biostrings)

#Reads the gene alignment
dna <- readDNAStringSet("metazoa_alignment.gene.fasta")

#Shows the sequence names
names(dna)

#Finds the Homo sapiens sequence
human_idx <- grep("Homo_sapiens|Homo sapiens", names(dna), ignore.case = TRUE)

#Extracts it
human_dna <- dna[human_idx]

#Removes any alignment gaps
human_seq_clean <- gsub("-", "", as.character(human_dna[[1]]))

#Trims to a multiple of 3 because the translation reads DNA in codons
#and each codon is 3 nucleotides long
trim_n <- nchar(human_seq_clean) %% 3
if (trim_n != 0) {
  human_seq_clean <- substr(human_seq_clean, 1, nchar(human_seq_clean) - trim_n)
}

#Translates to protein
human_protein <- translate(DNAString(human_seq_clean))

#Saves as fasta so that I can BLAST it for later
human_protein_set <- AAStringSet(human_protein)
names(human_protein_set) <- names(human_dna)
writeXStringSet(human_protein_set, filepath = "Homo_sapiens_protein.fasta", format = "fasta")

#Prints result
human_protein_set
