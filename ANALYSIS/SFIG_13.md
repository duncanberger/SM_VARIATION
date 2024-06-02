# Supplmentary figure 13: Coefficient of variation scores
### Load packages
```
library("ggplot2")
```
### Read in data
```
cv <- read.table("cv_scores.csv", header=FALSE, sep=",")
```
### Plot data
```
ggplot(data=cv) + 
  geom_point(aes(x=V1, y=V2)) +
  scale_y_continuous(expand=c(0,0), limits=c(0.45,0.65)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,21), breaks=c(0,2,4,6,8,10,12,14,16,18,20)) +
  xlab("Number of clusters (K)") + 
  ylab("Coefficient of Variation Error") +
  theme( legend.position="none",panel.grid = element_blank(), 
         axis.text.y=element_text(face="bold", color="black", size=8),
         axis.title.y = element_text(face="bold", color="black", size=10),
         axis.text.x=element_text(face="bold", color="black", size=8),
         axis.title.x = element_text(face="bold", color="black", size=10),  
         panel.border = element_rect(color="black",fill = NA),
         panel.background = element_blank())
```
