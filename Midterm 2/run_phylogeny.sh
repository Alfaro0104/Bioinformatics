#!/bin/bash

# Make an output folder for the phylogeny results
mkdir -p tree_results

# Run a maximum likelihood phylogenetic analysis with 100 bootstrap replicates
# RAxML-NG is appropriate for this dataset because it is a nucleotide alignment
# and maximum likelihood is a standard method for estimating phylogenies from DNA sequence data
raxml-ng \
  --all \
  --msa metazoa_alignment.5k.fasta \
  --model GTR+G \
  --bs-trees 100 \
  --threads auto \
  --prefix tree_results/metazoa_5k
