# Supplementary Figure 2. Inference of demographic history using SMC++ for *Schistosoma mansoni* subpopulations
### Load libraries
```
library("ggplot2")
library("scales")
```
### Load data
```
smc2 <- read.csv("plot2.csv", header=TRUE)
smc2 <- subset(smc, label!='SEN6' & label!='MAY3' & label!='MAY5' &
                  label!='SEN1' & label!='SEN2' )
```
### Plot data
```
ggplot(data=subset(smc2)) + 
  geom_line(aes(x=(x),y=log10(y), color=label), size=1) +
  scale_x_log10(expand=c(0,0), breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values=c(
    "KOC1"="#EE6677",
    "KOC2"="#EE6677",
    "KOC3"="#EE6677",
    "LAB"="#4477AA",
    "LAB1"="#4477AA",
    "MAY1"="#66CCEE",
    "MAY2"="#66CCEE",
    "MAY3"="#66CCEE",
    "MAY4"="#66CCEE",
    "MAY5"="#66CCEE",
    "GUAD"="#AA3377",
    "PR"="#c4a484",
    "KEN"="#57d9b3",
    "CAM"="pink",
    "TZ1"="#228833", 
    "TZ2"="#228833", 
    "ISL1"="#CCBB44",
    "ISL2"="#CCBB44",
    "ISL3"="#CCBB44",
    "SEN1"="#ff9d00",
    "SEN2"="#ff9d00",
    "SEN3"="#ff9d00",
    "SEN4"="#ff9d00",
    "SEN5"="#ff9d00",
    na.value="grey75")) +

  xlab("Years before present") + 
  labs(y=expression(bold("Effective population size "*bolditalic((N[e]))))) +
  coord_cartesian(xlim = c(500, 10^5)) +
  theme_bw() + 
  theme(legend.key.size = unit(0.5,"cm"), 
        legend.title =element_blank(),
        axis.text=element_text(size=5, face="bold", color="black"),
        axis.title=element_text(size=7, face="bold", color="black"),
        panel.grid = element_blank())  +
  annotation_logticks(sides="lb",outside = FALSE,
                      short = unit(0.05, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.1,"cm"))
```

