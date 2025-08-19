### Written by Abigail Clyde, last updated 8.17.25

## This R script makes the stacked bar plots of sequences recovered by participant and timepoint (Figure 1c)
## Input for this script is a .csv of sequences from each study that were recovered for each participant at each timepoint by SMRT-UMI
## the .csv is organized as: columns: | ID | Timept _weeks | Sequence_Count |
## Each row of .csv is a sample we sequenced by SMRT-UMI (i.e. there are 4 rows of participant1, one for each timepoint sequenced)

## To convert this .csv into a stacked bar plot" 

library(ggplot2)
library(ragg)
library(RColorBrewer)
library(svglite)

setwd("~/dir")
theme_plot <- function (){
  font <- "Arial"
  theme_bw() %+replace%
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size=6, family=font),
      axis.title = element_text(size=6, family=font),
      plot.title = element_text(size=6, family=font)
    )
}
data <- read.csv("smrtumi_sequences.csv")
data$TimePT <- factor(data$TimePT, levels=c("W16", "W12", "W8", "W4", "W1", "D0")) ### if the weeks are not already ordered by timepoint, this does so
blues <- c(brewer.pal(6, "Blues"))
ggplot(data, aes(fill=TimePT, y=Count, x=ID)) +
  geom_bar(position ="stack", stat="identity") + 
  theme_plot() + 
  labs(title = "10-1074", x= "Participant", y="Sequences Recovered by SMRT-UMI", fill="Time Points") +
  scale_fill_manual(values = blues) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("smrtumi_sequences.svg", width=6, height =4.5, units="cm", dpi=600)