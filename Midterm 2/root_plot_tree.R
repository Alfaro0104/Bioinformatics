#Installs ape if it is not already installed
if (!requireNamespace("ape", quietly = TRUE)) {
  install.packages("ape")
}

#Loads the ape package so the tree can be read and rooted
library(ape)

#Reads in the tree from RAxML-NG
tree <- read.tree("tree_results/metazoa_5k.raxml.support")

#Prints the tip labels so that the taxon names can be checked
tree$tip.label

#Roots the tree on the branch to Plakina jani and Grantia compressa
rooted_tree <- root(tree, outgroup = c("Plakina_jani", "Grantia_compressa"), resolve.root = TRUE)

#Save the rooted tree in a Newick format
write.tree(rooted_tree, file = "tree_results/metazoa_5k_rooted.newick")

#Saves the rooted tree as a PDF figure with the actual bootstrap support values
pdf("tree_results/metazoa_5k_rooted_tree.pdf", width = 12, height = 8)
plot(rooted_tree, cex = 0.7)
nodelabels(text = rooted_tree$node.label, frame = "rect", bg = "lightblue", cex = 0.8)
dev.off()
