## Written by Abigail Clyde. Last updated 8.20.25

### This is an R script for generating correlations between time to rebound and viral load in 10-1074 samples

## First, read in .csv containing data
## This includes rebound timepoint (first timepoint post-nadir with a 0.5log10 increase in VL, confirmed by subsequent appointment), percentage pre-existing mutations, and VL at Day 0 and Week 1


library(ggplot2)
library(ggpmisc)
library(svglite)

df <- read.csv("~/time_to_rebound_MCA0885.csv")
colnames(df) <- c("id", "weeks", "days", "twentyfive", "thirtytwo", "thirtyfour", "totalres", "vl_log", "vl", "vl_log_w1", "vl_w1")
df[ , colSums(is.na(df))==0]

## now exclude 1HD9K which is an outlier (100% resistance at day 0)
g1 <- subset(df, id=="1HD9K")
df <- df[-9,]
df
theme_plot <- function (){
  font <- "Arial"
  theme_bw() %+replace%
    theme(
      panel.grid.major.x = element_line(color="lightgrey", size=0.2),
      panel.grid.major.y = element_line(color="lightgrey", size=0.2),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size=6, family=font),
      axis.title = element_text(size=6, family=font),
      #legend.title=element_text(size=16, family=font, face="bold"),
      #legend.position="right",
      plot.title = element_text(size=6, family=font)
    )
}

### In panel A, we compare days to rebound with % pre-existing resistance mutations
corr <- cor.test(x=df$days, y=df$totalres, method = 'pearson')
corr
ggplot(data=df, aes(days, totalres)) +
  geom_point() +
  ylab("Percent resistance mutations at day 0") + 
  xlab("Rebound timepoint(weeks)") +
  scale_x_continuous(breaks=c(0,7,14,28,56), labels =c("0", "1", "2", "4", "8"))+
  scale_y_continuous(limits=c(0, 1)) +
  geom_smooth(method=lm, se=FALSE, color="black", size = 0.25) +
  theme_plot() +
  geom_count()+
  scale_size_continuous(breaks = round)+
  geom_point(data=g1, color="red")+
  stat_cor(method="pearson", label.x=45, label.y=0.99, size=6/.pt)+
  geom_smooth(method=lm, color="black", se=FALSE, size=0.5)
ggsave("time_to_rebound_res.svg", width=8, height=5, units="cm", dpi=300)


## Panel B compares viral load at day 0 to % pre-existing resistance mutations
corr.b <- cor.test(x=df$vl, y=df$totalres, method = 'pearson')
corr.b
ggplot(data=df, aes(vl, totalres)) +
  geom_point(size=0.5) +
  ylab("Percent resistance mutations at day 0") + 
  xlab("Viral load at day 0 (thousand copies/mL)") +
  scale_x_log10(breaks=c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000), labels=c("0", "10", "20", "30", "40", "50", "60", "70", "80")) +
  scale_y_continuous(limits=c(0, 1), labels=c(0, 25, 50, 75, 100)) +
  geom_smooth(method=lm, se=FALSE, color="black", size = 0.25) +
  theme_plot() +
  geom_point(data=g1, size=1, color="red")+
  annotate("label", x=75000, y=1, label=paste("r=", corr.b$estimate, sep=""), label.size=NA, size=6/.pt) +
  geom_smooth(method=lm, color="black", se=FALSE, size=0.5)
ggsave("vl_res.svg", width=8, height=5, units="cm", dpi=300)


## Panel C compares viral laod at week 1 to % pres-existing resistance mutations
corr.c <- cor.test(x=df$vl_w1, y=df$totalres, method = 'pearson')
corr.c
ggplot(data=df, aes(vl_w1, totalres)) +
  geom_point(size=0.5) +
  ylab("Percent resistance mutations at Day 0") + 
  xlab("Viral load at Week 1 (copies/mL)") +
  geom_smooth(method=lm, se=FALSE, color="black", size = 0.25) +
  #scale_x_continuous(trans="log10", breaks=c(1, 100, 1000, 2000, 5000, 10000)) +
  scale_x_log10(breaks=c(0, 500, 1000, 2000, 5000, 10000), labels =c("0", "500", "1000", "2000", "5000", "10000")) + 
  scale_y_continuous(limits=c(0, 1), labels=c(0, 25, 50, 75, 100)) +
  theme_plot() +
  geom_point(size=1) +
  ylab("Percent resistance mutations at Day 0") + 
  xlab("Viral load at Week 1 (copies/mL)") +
  scale_x_log10(breaks=c(0, 500, 1000, 2000, 5000, 10000), labels =c("0", "500", "1000", "2000", "5000", "10000")) + 
  scale_y_continuous(limits=c(0, 1), labels=c(0, 25, 50, 75, 100)) +
  theme_plot() +
  geom_point(data=g1, size=1, color="red") +
  annotate("label", x=75000, y=1, label=paste("r=", corr.b$estimate, sep=""), label.size=NA, size=6/.pt) +
  geom_smooth(method=lm, color="black", se=FALSE, size=0.5)
ggsave("vl_w1_res.svg", width=8, height=5, units="cm", dpi=300)
