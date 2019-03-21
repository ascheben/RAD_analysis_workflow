#!/usr/bin/env Rscript

library(gdsfmt)
library(SNPRelate)
library(gridExtra)
library(ggrepel)

# Usage: Rscript scripts/pca.R <input.vcf> <popmap.txt> <out_prefix>

# Get user arguments
args <- commandArgs(trailingOnly = TRUE)
# Input VCF file
invcf <- args[1]
# Read tab-separated popmap ("sample\tgroup") matching names in VCF
popfile <- read.csv(args[2], header = FALSE, sep = "\t")
# Add popmap header
colnames(popfile) <- c("sample", "group")
# Set output file names
outname <- args[3] 
gdsfile <- paste(outname, "gds", sep=".")
pdffile1 <- paste(outname, "_label_pca.pdf", sep="")
pdffile2 <- paste(outname, "_pca.pdf", sep="")
# Convert VCF to GDS
snpgdsVCF2GDS(invcf,gdsfile, method="biallelic.only")
# Read in genotypes
genofile <- snpgdsOpen(gdsfile, allow.duplicate=TRUE)
# LD pruning option (turned off by default)
set.seed(1000)
snpset <- snpgdsLDpruning(genofile,autosome.only=FALSE, slide.max.bp = 1)
# PCA 
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile,autosome.only=FALSE, snp.id=snpset.id,num.thread=2)
pc.percent <- pca$varprop*100
round.pc.percent <- round(pc.percent, 2)
# Get % explained and store in vectors below
ev1pc <- round.pc.percent[1]
ev2pc <- round.pc.percent[2]
ev1pc <- paste("PC1 (", ev1pc, "%)", sep ="")
ev2pc <- paste("PC2 (", ev2pc, "%)", sep ="")
# Prepare PCA results
tab <- data.frame(sample.id = pca$sample.id,
                 EV1 = pca$eigenvect[,1], # the first eigenvector
                 EV2 = pca$eigenvect[,2], # the second eigenvector
                 stringsAsFactors = FALSE)
tab$sample.id <- sub(".*__", "", tab$sample.id)
tabpops <- merge(tab,popfile, by.x = "sample.id" , by.y = "sample" )
# Set plot theme
theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             strip.background=element_blank(),axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"))
# Plot with ggrepel labels
p <- ggplot(tabpops,aes(x=EV1,y=EV2,color=group, label=sample.id))
p <- p + geom_point() + theme + guides(fill=guide_legend(title="New Legend Title")) + geom_text_repel() + xlab(ev1pc) + ylab(ev2pc)
# Plot without labels
p_nolab <- ggplot(tabpops,aes(x=EV1,y=EV2,color=group))
p_nolab <- p_nolab + geom_point() + theme + xlab(ev1pc) + ylab(ev2pc)
# Save plots as pdf
pdf(pdffile1,width = 10,height = 10)
yy <- grid.arrange(p,nrow=1)
op <- par(no.readonly=TRUE)
par(op)
dev.off()
pdf(pdffile2,width = 10,height = 10)
yy <- grid.arrange(p_nolab,nrow=1)
op <- par(no.readonly=TRUE)
par(op)
dev.off()
