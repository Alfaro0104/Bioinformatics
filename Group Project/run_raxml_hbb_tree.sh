#!/bin/bash

# RAxML-NG maximum-likelihood tree for hemoglobin beta protein sequences
# Input: aligned FASTA file generated from R MSA script
# Output: ML tree with bootstrap support values

raxml-ng \
  --all \
  --msa hemoglobin_beta_alignment_6_species.fasta \
  --msa-format FASTA \
  --model LG+G \
  --prefix HBB_tree_bootstrap \
  --bs-trees 100 \
  --threads 2 \
  --seed 12345
