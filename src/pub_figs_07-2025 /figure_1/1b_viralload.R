### Written by Abigail Clyde, last updated 8.17.25

#### This is an R script to generate viral load curves for the 10-1074 and 3BNC117 bnAb studies (Figure 1b)
## Input for this script is a .csv for each study from Caskey et al 2016 and Caskey et al 2017
## Rows are days post bnAb infusion (timepts), Columns are participant IDs and second set of columns that are participant1_2 and have same viral load values but only in timepoint rows that were sampled by SMRT-UMI
## This second set of columns will be used to make points on VL curve to indicate which timepoints were sampled for each participant
## Viral loads are reported in log10 
## The same code was used for the 10-1074 viral load figure panel and the 3BNC117 viral load panel

library(ggplot2)
library(RColorBrewer)
library(svglite)
setwd("~/dir")
viralload <- read.csv("viralload.csv")
## this returns a data frame with columns: timepts | participant1 | participant2 | etc. | participant1_2 | participant 2_2 | etc.
## and rows: 0, 1, 2, etc. (days)

## set theme for figures
theme_plot <- function (){
  font <- "Arial"
  theme_bw() %+replace%
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border=element_rect(color="black", size =0.5, fill=NA),
      axis.text = element_text(size=6, family=font),
      axis.title = element_text(size=6, family=font)
    )
}
display.brewer.all(colorblindFriendly = TRUE)
colors <- brewer.pal(n = 11, name = 'Paired')
colors[11] <- "#B15928"
colors
## the colors vector returns this list of colors: "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#B15928"

ggplot(data=viralload, aes(x=timepts)) + 
  geom_line(aes(y=participant1), size=0.5, color=colors[1]) +
  geom_line(aes(y=participant2), size=0.5, color=colors[2])  +
  geom_line(aes(y=participant3), size=0.5, color=colors[3])  +
  geom_line(aes(y=participant4), size=0.5, color=colors[4]) +
  geom_line(aes(y=participant5), size=0.5, color=colors[5]) +
  geom_point(aes(y=participant1_2), size=0.5, color=colors[1]) +
  geom_point(aes(y=participant2_2), size=0.5, color=colors[2]) +
  geom_point(aes(y=participant3_2), size=0.5, color=colors[3]) +
  geom_point(aes(y=participant4_2), size=0.5, color=colors[4]) +
  geom_point(aes(y=participant5_2), size=0.5, color=colors[5]) +
  xlab("Weeks since infusion") + ylab("HIV-1 RNA (log10 copies/mL)") + 
  ggtitle("Viral Load") +
  theme_AC() +
  scale_y_continuous(breaks=c(1,2,3,4,5), labels = c("1", "2", "3", "4", "5"), limits=c(0.5, 5.5)) + 
  scale_x_continuous(breaks=c(0,7,28,56,84,112), labels =c("0", "1", "4", "8", "12", "16")) ### this converts timepoint labels from days to weeks
ggsave("viralload.svg", width=5, height =4.2, units="cm", dpi=600)

### saving as .svg allows later manipulation of size in Illustrator without changing proportions of chart
### Plot legend was added in Adobe Illustrator

