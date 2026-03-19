#The first couple lines of code here are to install the packages needed for this lab
#and to set the working directory to my correct folder
setwd("C:\Users\Alfar\Documents\GitHub\Bioinformatics\lab10")
BiocManager::install("Biostrings")
BiocManager::install("GenomicAlignments")
install.packages("UniprotR")
install.packages("protti")
install.packages("r3dmol")

#These lines load the packages into the R session
library(Biostrings)
library(UniprotR)
library(protti)


#This makes an amino acid sequence object from the translated sequence I got from Lab 6
#and the stop codons were taken out to make the BLAST more efficient which is why there are "*" present
aa_seq <- AAStringSet(c(MyProtein = "MTHQSHAYHIVKPSP*PLTGALSALLMTSGLAM*FHFHSITLLILGLLTNTLTIYQ*WRDVTRESTYQGHHTPPVQKGLRYGIILFITSEVFFFAGFF*AFYHSSLAPTPQLGGHWPPTGITPLNPLEVPLLNTSVLLASGVSIT*AHHSLIENNRNQIIQALLITILLGLYFTLLQASEYFESPFTISDGIYGSTFFVATGFHGLHVIIGSTFLTICFIRQLIFHFTSKHHFGFEAAA*YWHFVDVV*LFLYVSIY"))

#This turns the amino acid sequence into a FASTA file for the UniProt BLAST
writeXStringSet(aa_seq, filepath = "protein_sequence.fasta", format = "fasta")

#This line reads the UniProt accession numbers I got from a text file
#so that the file should have one accession per line
acc_df <- read.table("uniprot_accessions.txt", header = FALSE, stringsAsFactors = FALSE)

#This part actually turns the accession column into a character vector
accessions <- as.character(acc_df$V1)

#Retrieves the Gene Ontology information for the accession numbers i got from earlier
#and saves the GO results to an object called go_info
go_info <- GetProteinGOInfo(accessions)

#Prints the GO results to the console so that i can actually see them
go_info

#Saves the default GO plot to a PNG image file in my GitHub folder
png("GO_plot_basic.png", width = 1200, height = 900)

#Plots the GO information
PlotGoInfo(go_info)

#Closes the graphics file so that the image file is actually written correctly
dev.off()

#Takes my GO results and turns them into a clean table to be saved as a csv file
go_df <- as.data.frame(go_info)
write.csv(go_df, "GO_terms_results.csv", row.names = FALSE)

#Finds any pathology or biotechnology related information associated with the protein
pathology_info <- GetPathology_Biotech(accessions)

#Finds any disease information associated with the protein
disease_info <- Get.diseases(accessions)

#Saves both outputs frome earlier to a CSV files
write.csv(as.data.frame(pathology_info), "pathology_results.csv", row.names = FALSE)
write.csv(as.data.frame(disease_info), "disease_results.csv", row.names = FALSE)

#Uses the protti package to get any UniProt information for the accession numbers
#and also has the column xref_pdb which has any PDB IDs I can use
uniprot_info <- fetch_uniprot(accessions)

#Saves the UniProt information table
write.csv(uniprot_info, "uniprot_info_results.csv", row.names = FALSE)

#Pulls out the PDB IDs from the xref_pdb column
#and unlist() is used because the PDB IDs may be stored as lists
pdb_ids <- unique(unlist(uniprot_info$xref_pdb))

#takes out any missing or empty values
pdb_ids <- pdb_ids[!is.na(pdb_ids) & pdb_ids != ""]
#and removes the semicolon so that the next line works
pdb_ids <- gsub(";", "", pdb_ids)

#Gathers any structural information from the Protein Data Bank
pdb_info <- fetch_pdb(pdb_ids)

#Saves the PDB results to a CSV file
write.csv(pdb_info, "pdb_results.csv", row.names = FALSE)

#Gets me the AlphaFold structural prediction information for the accession numbers
af_info <- fetch_alphafold_prediction(accessions)


#For question 8:My GO results showed me that this gene is used in cell respiration
#and the mitochondria electron transport chain specifically referring to cytochrome c.
#It also mentioned the cytochrome-c oxidase activity and that it is located in the inner
#mitochondrial membrane.

#For question 9: The disease that I found was called Leber hereditary neuropathy which is
#a mitochondrial disorder that leads to vision loss from a lack of energy production that happens
#because of a disorder of the mitochondrial electron transport chain specifically at
#complex IV.
