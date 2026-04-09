#UniProt accession for human POLG that I got from the BLASTp result being NP_001119603.1 and then
#went to UniProt and the equivalent value was P54098
acc <- "P54098"

#Builds the UniProt REST query
query <- URLencode(paste0("accession:", acc), reserved = TRUE)
u <- paste0(
  "https://rest.uniprot.org/uniprotkb/search?query=",
  query,
  "&fields=accession,gene_names,go_p,go_f,go_c&format=tsv"
)

#Reads the table from UniProt
go_tab <- read.delim(u, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

#Pulls any GO columns
bp_raw <- go_tab[[3]]
mf_raw <- go_tab[[4]]
cc_raw <- go_tab[[5]]

#Splits GO terms
bp_terms <- if (!is.na(bp_raw) && bp_raw != "") unlist(strsplit(bp_raw, "; ")) else character(0)
mf_terms <- if (!is.na(mf_raw) && mf_raw != "") unlist(strsplit(mf_raw, "; ")) else character(0)
cc_terms <- if (!is.na(cc_raw) && cc_raw != "") unlist(strsplit(cc_raw, "; ")) else character(0)

#Builds a GO table
go_df <- data.frame(
  Ontology = c(rep("Biological Process", length(bp_terms)),
               rep("Molecular Function", length(mf_terms)),
               rep("Cellular Component", length(cc_terms))),
  Term = c(bp_terms, mf_terms, cc_terms),
  stringsAsFactors = FALSE
)

#Saves GO terms
write.csv(go_df, "uniprot_go_terms.csv", row.names = FALSE)

#Counts GO terms by ontology
go_counts <- table(go_df$Ontology)

#Plots a GO summary
pdf("go_info_plot.pdf", width = 8, height = 5)
barplot(go_counts,
        main = "GO sub-ontologies for POLG",
        ylab = "Number of GO terms")
dev.off()

#Prints one example from each ontology
cat("Biological Process example:", if (length(bp_terms) > 0) bp_terms[1] else "None found", "\n")
cat("Molecular Function example:", if (length(mf_terms) > 0) mf_terms[1] else "None found", "\n")
cat("Cellular Component example:", if (length(cc_terms) > 0) cc_terms[1] else "None found", "\n")
