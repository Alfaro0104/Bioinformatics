#!/bin/bash

#Make an output folder for the phylogeny results
mkdir -p tree_results

#Ran a maximum likelihood phylogenetic analysis with 100 bootstrap replicates
#RAxML-NG is good for this data because it is a concatenated nucleotide alignment and not partitioned
raxml-ng \
  --all \
  --msa metazoa_alignment.5k.fasta \
  --model GTR+G \
  --bs-trees 100 \
  --threads auto \
  --prefix tree_results/metazoa_5k
