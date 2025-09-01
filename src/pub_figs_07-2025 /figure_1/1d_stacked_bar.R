#Written by Abigail Clyde, last updated 8.17.25

## This R script counts sequences sequenced by SMRT-UMI vs viral load by pulling from .fasta files and puts these into a dataframe
## Then makes the scatter plots comparing SMRT-UMI to viral load

library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(tidyverse)
library(scales)
library(svglite)

##functions

#read a fasta file
#fileName:name of file
#returns: dataframe with columns name (name between > and the first ' '), seq (sequence), and longName (the whole > line)
read.fa <- function(fileName, format="fasta"){ 
  
  myal <- read.alignment(fileName, format="fasta") 
  outal <- as.data.frame(matrix(NA, nrow=length(myal$seq), ncol=3))
  colnames(outal) <- c('name','seq', 'longName')
  outal[,1] <- myal$nam
  outal[,3] <- myal$nam
  outal[,2] <- unlist(myal$seq)
  return(outal)
  
}

##make sequences list of those sequenced by SMRT-UMI
##There is one .fasta file per participant with all timepoints sequenced from both SMRT-UMI and previously sequenced by SGA aligned
## each fasta file is named: study_participantID.fasta 
setwd("~/dir")
seqfiles <- list.files(pattern=".fasta")
##isolate participant IDs from fasta names
ids <- unlist(lapply(strsplit(seqfiles, split="_"), function(x) x[2]))
filenames <- gsub(".fasta", "", seqfiles)

## empty dataframe for sequence counts
Total <- as.data.frame(matrix(NA, ncol=3, nrow=length(utimepts)))
## Read in NT alignments and collate sequence data
for(i in 1:length(ids)){
  myal <- read.fa(pasta(seqfiles[i], sep=""))
  myal[,2] <- sapply(myal[,2], toupper)
  myal[,2] <- gsub("U", "T", myal[,2])
  myal[,1] <- gsub("-", "_", myal[,1])
  sequences <- myal[,1]
  sequences
  ## isolate sequences sequenced by SMRT-UMI (fasta seq names = study_participant_timept)
  ##SGA sequences start with KX or KY (NCBI accession numbers), while SMRT-UMI seqs start with study number (835 or 885)
  UMI <- substr(sequences, 1,2)
  sequences <- sequences[which(UMI=="83")] ##
  sequences
  timepts <- toupper(unlist(lapply(strsplit(sequences, split="_"), function(x) x[3])))
  utimepts <- unique(timepts)
  utimepts <- unique(timepts[-1]) ### exclude HXB2 reference sequence (first seq in alignment)
  numval <- as.numeric(substring(utimepts,2,nchar(utimepts)))
  utimepts <- utimepts[order(numval)]
  utimepts 
 
  seq_count <- as.data.frame(matrix(NA, ncol=3, nrow=length(utimepts)))
  timept_v <- vector("list", length(utimepts))
  names(seq_count) <- c("ID", "TimePT", "Count")
  for(j in 1:length(utimepts)){
    seq_count$ID <- ids[i]
    seq_count$TimePT[j] <- utimepts[j]
    seq_count$Count[j] <- length(sequences[which(timepts==utimepts[j])])
  }
  seq_count
  Total <- rbind(seq_count, Total)
}
Total$sample <- paste(Total$ID, Total$TimePT, sep="_") ## gives ID_timepoint for easy ID
write.csv(Total, file="seq_count.csv")

## make graphs comparing SMRT-UMI counts to viral load
## Manually added viral load numbers for each timepoint/participant to above .csv "seq_count.csv" (same .csv with  5th column "VL")
seq_count <- read.csv("seq_count.csv")
theme_plot <- function (){
  font <- "Arial"
  theme_bw() %+replace%
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size=6, family=font),
      axis.title = element_text(size=6, family=font),
      #legend.title=element_text(size=16, family=font, face="bold"),
      #legend.position="right",
      plot.title = element_text(size=6, family=font)
    )
}

corr <- cor.test(x=seq_count$VL, y=seq_count$Count, method = 'spearman')
ggplot(seq_count, aes(x=VL, y=Count)) +
  geom_point(size=0.5) +
  xlab("Viral Load (c/mL)") +
  ylab("Sequences Recovered by SMRT-UMI") +
  scale_x_log10(labels = comma, breaks = c(1000,10000,100000,1000000)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(0,5,10,25,50,100,200,300)) +
  theme_plot() +
  coord_cartesian(
    ylim = c(1, 350),
    expand = TRUE,
    default = FALSE,
    clip = "on"
  ) +
  annotate("label", x=100000, y=1, label=paste("r=", corr$estimate, sep=""), label.size=NA)
ggsave("SMRT_UMI_vs_VL.svg", width = 5, height = 4.25, units="cm", dpi=600)
