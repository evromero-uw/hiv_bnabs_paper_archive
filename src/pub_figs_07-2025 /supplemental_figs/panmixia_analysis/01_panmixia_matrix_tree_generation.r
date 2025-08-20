## Panmixia analysis figure

## First I generate trees from fasta files of each sample from aligned fasta files
library(ggplot2)
library(stringr)
library(seqinr)
library(ape)
library(phytools)

setwd("~/final_panmixia/")

##functions
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

myNTdir <- "~/fastas/"
myNEWICKdir <- "~/trees/"
myTREEdir <- "~/pdf_trees/"
utimecols <- c("chartreuse3", "orange", "lightblue2", "dodgerblue", "navy",   "thistle1")
mytimepts <- c("D0", "W1", "W4", "W8", "W12", "W16")
assignCols <- cbind.data.frame(mytimepts, utimecols[1:length(mytimepts)])
names(assignCols) <- c("timePT", "Color")
assignCols
NTfiles <- list.files(myNTdir, pattern=".fasta")
if("HXB2.fasta" %in% NTfiles){ NTfiles <- NTfiles[-which(NTfiles=="HXB2.fasta")] }
filenames <- gsub(".fasta", "", NTfiles)
filenames

### generate trees for intermixed samples using iqtree, 10000 bootstraps
for(i in 1:length(filenames)){
  print(filenames[i])
  ### read in NT alignment
  alfile <- paste(myNTdir, NTfiles[i], sep="")
  myal <- read.fa(alfile)
  mywrite.fa(myal[-1,], "temp.fasta")
  setwd("~/iqtree-2.3.4-macOS")
  system(paste("bin/iqtree2 -s ~/temp.fasta -m TN -B 10000 -T AUTO --prefix ~/trees/", 
               filenames[i], "_NT", sep=""))
  setwd("~/new_sept_2024")
}

treefiles <- list.files(myNEWICKdir, pattern=".treefile")
treefiles
id2<- gsub(".treefile", "", treefiles)
utimecols <- c("black", "black","black","red", "red")
mytimepts <- c("MC", "83", "88", "KX", "KY")
assignCols <- cbind.data.frame(mytimepts, utimecols[1:length(mytimepts)])
names(assignCols) <- c("timePT", "Color")
assignCols

## generate matrices
for (i in 1:length(treefiles)) {
  tree<- paste(myNEWICKdir, treefiles[i], sep="")
  tree <- read.tree(file=tree)
  mytree <- midpoint_root(tree)
  mytree <- ladderize(mytree)
  matrix <- cophenetic.phylo(mytree)
  matrix
  write.csv(matrix, file=paste(myNEWICKdir, paste(id2[i],"_10000matrix.csv", "" ), ""))
  ### these .csvs will be used in next py script to generate p-values for panmixia

  ##plot trees of intermixed samples
  treesource <-substring(mytree$tip.label,1,2)
  treesource
  utree <- unique(treesource)
  mycols <- assignCols$Color[match(treesource, assignCols$timePT)] 
  legcol <- assignCols$Color[match(utree, assignCols$timePT)]
  figure <- paste(myNEWICKdir, id2[i], ".10000bootTree.pdf",sep="")
  pdf(file=figure, pointsize=14, width=16, height=28, useDingbats=F) 
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  
  plot(mytree, cex=2, tip.color=mycols, main=id2[i], show.tip.label = F, cex.main=4) 
  tiplabels(pch=15, col=mycols, lwd=2, cex=2)
  add.scale.bar("bottomright")
  legend("topright", legend=c("SMRT-UMI", "SGA"), fill=c("black", "red"), bty='n', cex=2)
  dev.off()
}


####example plot in figure

id<-"1HD1_d0"
utimecols <- c("black", "black","black","red", "red")
mytimepts <- c("MC", "83", "88", "KX", "KY")
assignCols <- cbind.data.frame(mytimepts, utimecols[1:length(mytimepts)])
names(assignCols) <- c("timePT", "Color")
assignCols
tree<- ""
tree <- read.tree(file=tree)
mytree <- midpoint_root(tree)
mytree <- ladderize(mytree)
matrix <- cophenetic.phylo(mytree)
matrix
treesource <-substring(mytree$tip.label,1,2)
treesource
utree <- unique(treesource)
mycols <- assignCols$Color[match(treesource, assignCols$timePT)] 
legcol <- assignCols$Color[match(utree, assignCols$timePT)]
figure <- paste("1hd1_d0.expanded.tree.pdf")
pdf(file=figure, pointsize=14, width=18, height=46, useDingbats=F) 
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

plot(mytree, cex=2, tip.color=mycols, show.tip.label = F, cex.main=4, edge.width = 3) 
tiplabels(pch=15, col=mycols, lwd=2, cex=4)
add.scale.bar("bottomright")
legend("topright", legend=c("SMRT-UMI", "SGA"), fill=c("black", "red"), bty='n', cex=2)
dev.off()

