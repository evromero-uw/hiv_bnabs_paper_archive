## Written by Abigail Clyde. Last updated 8.18.25

### Make a tree of all participant sequences
## Align all participants to HXB2 in one .fasta file

library(ggplot2)
library(readxl)
library(stringr)
library(stringdist)
library(seqinr)
library(ape)
library(phytools)
library(RColorBrewer)

## functions
read.fa <- function(fileName, format="fasta"){ 
  
  myal <- read.alignment(fileName, format="fasta") 
  outal <- as.data.frame(matrix(NA, nrow=length(myal$seq), ncol=3))
  colnames(outal) <- c('name','seq', 'longName')
  outal[,1] <- myal$nam
  outal[,3] <- myal$nam
  outal[,2] <- unlist(myal$seq)
  return(outal)
  
}

write.fa <- function(names,dna,fileName,addBracket=FALSE){
  if(addBracket|any(grep('^[^>]',names)))names<-paste('>',names,sep='')
  dna <- sub(' +$','',dna,perl=TRUE)
  output <- paste(names,dna,sep="\n")
  writeLines(output,sep="\n",con=fileName)
}

mywrite.fa <- function(myalignment, myfile){ write.fa(myalignment[,1],myalignment[,2],myfile,addBracket=TRUE) }

## make tree
setwd("~/10-1074_full_tree")
#### input/output directories
myNEWICKdir <- "~/10-1074_full_tree"
alfile <- "~/10-1074_total.fasta"
display.brewer.all()
brewer.pal(9, "Set1")
utimecols <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928")
mytimepts <- c("1HB3", "1HC2", "1HC3", '1HD1', '1HD4K', '1HD5K', '1HD6K', '1HD7K', '1HD9K', '1HD10K', '1HD11K')
assignCols <- cbind.data.frame(mytimepts, utimecols[1:length(mytimepts)])
names(assignCols) <- c("timePT", "Color")
assignCols
ids <- "TEN_total"

myal <- read.fa(alfile)
### create newick file with IqTree
mywrite.fa(myal[-1,], "temp.fasta")
setwd("~/iqtree-2.3.4-macOS")
system("bin/iqtree2 -s ~/10-1074_full_tree/temp.fasta -m TN -T AUTO -redo --prefix ~/10-1074_full_tree/10-1074_total")
setwd("~/10-1074_full_tree")

treefile <- "~/10-1074_total.treefile"
mytree <- read.tree(file=treefile)
mytree <- midpoint_root(mytree)
mytree <- ladderize(mytree)
part <- toupper(unlist(lapply(strsplit(mytree$tip.label, split="_"), function(x) x[2])))
part
unique(part)
mycols <- assignCols$Color[match(part, assignCols$timePT)] 
utimepts <- unique(part)
legcol <- assignCols$Color[match(utimepts, assignCols$timePT)]
figure <- "TEN_full_tree.pdf"
pdf(file=figure, pointsize=14, width=28, height=28, useDingbats=F) 
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

plot(mytree, cex=2, tip.color=mycols, show.tip.label = F, cex.main=4, type="fan", edge.width = 4) 
tiplabels(pch=15, col=mycols, lwd=4, cex=2)
add.scale.bar("bottomright", lwd=4)
#legend("topright", legend=utimepts, fill=legcol, bty='n', cex=2)
dev.off()

