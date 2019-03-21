#!/usr/bin/env Rscript

library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
library("colorspace")

# Usage: Rscript ./ggtree.R <raxml.tre> <popmap.txt>

# Get user arguments
args = commandArgs(trailingOnly=TRUE)
# Newick tree file output by RAxML as bipartitions
intree <- args[1]
# Population map file in format "sample\tgroup"
popmap <- args[2]

# Define output file names
outname1 <- paste(intree,"_nobrln_tree.svg",sep="") 
outname2 <- paste(intree,"_brln_tree.svg",sep="") 
# Read in population map
pops <- read.csv(popmap,header = F, sep = "\t")
# Get number of populations
numpop <- length(unique(pops$V2))

# Define taxonomic units to group by
poplist <- split(pops$V1, pops$V2)

# Read in tree file
raxml <- read.tree(intree)

# Assign taxonomic unit groups to tree
mytree <- groupOTU(raxml, poplist)

# Plot tree without branch lengths showing only bootstraps over 70
# Adjust width to fit content if required
svg(outname1, width = 10)
ggtree(mytree, aes(color=group),branch.length="none", layout="rectangular") + geom_tiplab(size=1) +
  scale_color_manual(values=c("black", rainbow_hcl(numpop))) + theme(legend.position="right")+ 
  geom_text2(size = 2,aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  scale_size_manual(values=c(1, .1))  + ggplot2::xlim(0, 60)
dev.off()
# Plot tree with branch lengths showing only bootstraps over 70
svg(outname2, width = 10)
ggtree(mytree, aes(color=group), layout="rectangular") + geom_tiplab(size=1) +
  scale_color_manual(values=c("black", rainbow_hcl(numpop)))  + theme(legend.position="right")+ 
  geom_text2(size = 2,aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70))+
  scale_size_manual(values=c(1, .1)) 
dev.off()
