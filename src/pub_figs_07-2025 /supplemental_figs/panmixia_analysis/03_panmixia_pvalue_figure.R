## Written by Abigail Clyde. Last updated 8.20.25

## Supp. fig 2
## This is an R script for plotting p-values generated from matrices comparing SGA and SMRT-UMI sequences from the same sample 

library(ggplot2)
library(svglite)
theme_plot <- function (){
  font <- "Arial"
  theme_bw() %+replace%
    theme(
      panel.grid.major.x = element_line(color="lightgrey", size=0.2),
      panel.grid.major.y = element_line(color="lightgrey", size=0.2),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size=6, family=font),
      axis.title = element_text(size=6, family=font),
      legend.title=element_text(size=6, family=font),
      legend.text = element_text(size=6, family=font),
      legend.position="right",
      plot.title = element_text(size=6, family=font)
    )
}

## matrix p-value .csv exported from 02_SGA_SMRT-UMI_panmixia.ipynb was binned by p-values and converted to .csv below
pval<-read.csv("~/panmixia/binned_pvalues.csv")
pval
ggplot(pval, aes(x=category, y=X., fill=p.10..4)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values=c("grey", "black"))+
  scale_x_discrete(labels=c("p>1E-300", "p>1E-100", "p>1E-30", "p>1E-10", "p>1E-4", "p>1E-3", "p>1E-2", "p>1E-1", "p>0.5")) +
  scale_y_continuous(breaks=c(2, 4, 6, 8, 10))+
  xlab("P-value") +
  ylab("Number of participants") + 
  theme_AC() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  guides(fill = guide_legend(title = "Compartmentalized?"))
ggsave("panmixia_pval.svg", height=4, width=6, units="cm")
