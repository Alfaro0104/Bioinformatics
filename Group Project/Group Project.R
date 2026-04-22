#Alignment part 1

#Install packages once if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
if (!requireNamespace("msa", quietly = TRUE)) BiocManager::install("msa")

library(Biostrings)
library(msa)

#Read the hemoglobin beta FASTA sequences
human <- readAAStringSet("P68871 HBB_HUMAN.fasta")
chicken <- readAAStringSet("P02112 HBB_CHICK.fasta")
bar_goose <- readAAStringSet("P02118 HBB_ANSIN.fasta")
greylag_goose <- readAAStringSet("P02117 HBB_ANSAN.fasta")
andean_goose <- readAAStringSet("P07036 HBB_CHLME.fasta")
vulture <- readAAStringSet("P68061 HBB_AEGMO.fasta")

#Combine all sequences into one AAStringSet
seqs <- c(
  human,
  chicken,
  bar_goose,
  greylag_goose,
  andean_goose,
  vulture
)

#Rename sequences so the alignment is easier to read
names(seqs) <- c(
  "Human_HBB",
  "Chicken_HBB",
  "Bar_headed_goose_HBB",
  "Greylag_goose_HBB",
  "Andean_goose_HBB",
  "Cinereous_vulture_HBB"
)

#Check that sequences loaded correctly
seqs
names(seqs)

#Run multiple sequence alignment using ClustalW
alignment <- msa(seqs, method = "ClustalW")

#Print alignment in the console
print(alignment)

#Save a nice PDF version of the alignment
msaPrettyPrint(
  alignment,
  output = "asis",
  file = "hemoglobin_beta_alignment_6_species.pdf",
  showNames = "left",
  showLogo = "top"
)

#Save the alignment as a FASTA file
writeXStringSet(
  AAStringSet(unmasked(alignment)),
  filepath = "hemoglobin_beta_alignment_6_species.fasta"
)




#Substitution detection part 2
# Convert the MSA to a character matrix
aligned <- as.matrix(alignment)

# Remove any columns that are only gaps, just in case
aligned <- aligned[, colSums(aligned != "-") > 0]

# Identify positions where at least one species differs
variable_positions <- which(
  apply(aligned, 2, function(column) length(unique(column)) > 1)
)

# Create a table showing each variable amino acid position
variable_table <- data.frame(
  Position = variable_positions,
  t(aligned[, variable_positions]),
  row.names = NULL
)

# Rename columns clearly
colnames(variable_table) <- c("Position", rownames(aligned))

# Print table in console
print(variable_table)

# Save table as CSV
write.csv(
  variable_table,
  "variable_amino_acid_substitutions_6_species.csv",
  row.names = FALSE
)



#GO analysis part 3

if (!requireNamespace("UniprotR", quietly = TRUE)) install.packages("UniprotR")
library(UniprotR)

accessions <- c(
  "P68871",
  "P02112",
  "P02118",
  "P02117",
  "P07036",
  "P68061"
)

go_results <- GetProteinGOInfo(accessions)

write.csv(
  go_results,
  "HBB_GO_terms_6_species.csv",
  row.names = FALSE
)

PlotGoInfo(go_results)

png(
  filename = "HBB_GO_terms_plot.png",
  width = 1200,
  height = 900,
  res = 150
)

PlotGoInfo(go_results)

dev.off()

go_results <- GetProteinGOInfo(accessions)

View(go_results)
dim(go_results)
colnames(go_results)


go_summary <- go_results[, c(
  "Gene.Ontology..biological.process.",
  "Gene.Ontology..molecular.function.",
  "Gene.Ontology..cellular.component."
)]

write.csv(
  go_summary,
  "HBB_GO_summary_6_species.csv",
  row.names = TRUE
)




# Make a simple GO annotation count plot from the GO results table

go_categories <- c(
  Biological_Process = "Gene.Ontology..biological.process.",
  Molecular_Function = "Gene.Ontology..molecular.function.",
  Cellular_Component = "Gene.Ontology..cellular.component."
)

go_counts <- sapply(go_categories, function(colname) {
  sum(!is.na(go_results[[colname]]) & go_results[[colname]] != "")
})

go_counts

png(
  filename = "HBB_GO_category_counts.png",
  width = 1200,
  height = 900,
  res = 150
)

barplot(
  go_counts,
  main = "GO annotation categories for hemoglobin beta proteins",
  ylab = "Number of proteins with GO annotations",
  ylim = c(0, max(go_counts) + 1)
)

dev.off()
