# Supplementary Figure 5: Functional profiling of *Sm*.TRPM<sub>PZQ</sub> variants
### Load packages
```{r}
library("dplyr")
library("ggplot2")
library("reshape2")
```
### Load data
```
pzq <- read.csv("pzq.csv", header=TRUE)
pzq_group <- pzq %>% group_by(Conc, Sample) %>% summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))
```
### Plot dose-response curves
```
A <- ggplot(data=subset(pzq_group, Sample=="Wild type "), aes(x=(Conc), y=Response_mean)) + theme_bw() + 
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE, size=0.6, color="black") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(side = "b", outside = TRUE, short = unit(.25,"mm"), mid = unit(0.5,"mm"), long = unit(1.5,"mm")) +
  xlab("Concentration (M)") +ylab("Response (% wild type)") + 
  scale_y_continuous(limits=c(-40,200), expand=c(0,0), breaks=c(-40,-20,0,20,40,60,80,100,120,140,160,180,200)) +
  ggtitle(label = "Wild type") +
  geom_point(data=subset(pzq, Sample=="Wild type "), aes(x=(Conc), y=Response), size=1.2, shape=1, alpha=0.75, color="grey50") +
  theme( panel.grid = element_blank(), legend.title = element_blank(), 
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         plot.title = element_text(face="bold", color="black", size=8, ,hjust = 1),
         axis.title.x = element_text(face="bold", color="black", size=8)) + coord_cartesian(clip = "off")

B <- ggplot(data=subset(pzq_group, Sample=="L102V"), aes(x=(Conc), y=Response_mean)) + theme_bw() + 
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE, size=0.6, color="black") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(side = "b", outside = TRUE, short = unit(.25,"mm"), mid = unit(0.5,"mm"), long = unit(1.5,"mm")) +
  xlab("Concentration (M)") +ylab("Response (% wild type)") + 
  scale_y_continuous(limits=c(-40,200), expand=c(0,0), breaks=c(-40,-20,0,20,40,60,80,100,120,140,160,180,200)) +
  ggtitle(label = "L102V") +
  geom_point(data=subset(pzq, Sample=="L102V"), aes(x=(Conc), y=Response), size=1.2, shape=1, alpha=0.75, color="grey50") +
  theme( panel.grid = element_blank(), legend.title = element_blank(), 
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         plot.title = element_text(face="bold", color="black", size=8, ,hjust = 1),
         axis.title.x = element_text(face="bold", color="black", size=8)) + coord_cartesian(clip = "off")

C <- ggplot(data=subset(pzq_group, Sample=="T105N "), aes(x=(Conc), y=Response_mean)) + theme_bw() + 
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE, size=0.6, color="black") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(side = "b", outside = TRUE, short = unit(.25,"mm"), mid = unit(0.5,"mm"), long = unit(1.5,"mm")) +
  xlab("Concentration (M)") +ylab("Response (% wild type)") + 
  scale_y_continuous(limits=c(-40,200), expand=c(0,0), breaks=c(-40,-20,0,20,40,60,80,100,120,140,160,180,200)) +
  ggtitle(label = "T105N ") +
  geom_point(data=subset(pzq, Sample=="T105N "), aes(x=(Conc), y=Response), size=1.2, shape=1, alpha=0.75, color="grey50") +
  theme( panel.grid = element_blank(), legend.title = element_blank(), 
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         plot.title = element_text(face="bold", color="black", size=8, ,hjust = 1),
         axis.title.x = element_text(face="bold", color="black", size=8)) + coord_cartesian(clip = "off")

D <- ggplot(data=subset(pzq_group, Sample=="T1005I "), aes(x=(Conc), y=Response_mean)) + theme_bw() + 
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE, size=0.6, color="black") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(side = "b", outside = TRUE, short = unit(.25,"mm"), mid = unit(0.5,"mm"), long = unit(1.5,"mm")) +
  xlab("Concentration (M)") +ylab("Response (% wild type)") + 
  scale_y_continuous(limits=c(-40,200), expand=c(0,0), breaks=c(-40,-20,0,20,40,60,80,100,120,140,160,180,200)) +
  ggtitle(label = "T1005I ") +
  geom_point(data=subset(pzq, Sample=="T1005I "), aes(x=(Conc), y=Response), size=1.2, shape=1, alpha=0.75, color="grey50") +
  theme( panel.grid = element_blank(), legend.title = element_blank(), 
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         plot.title = element_text(face="bold", color="black", size=8, ,hjust = 1),
         axis.title.x = element_text(face="bold", color="black", size=8)) + coord_cartesian(clip = "off")

E <- ggplot(data=subset(pzq_group, Sample=="M1068I "), aes(x=(Conc), y=Response_mean)) + theme_bw() + 
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE, size=0.6, color="black") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(side = "b", outside = TRUE, short = unit(.25,"mm"), mid = unit(0.5,"mm"), long = unit(1.5,"mm")) +
  xlab("Concentration (M)") +ylab("Response (% wild type)") + 
  scale_y_continuous(limits=c(-40,200), expand=c(0,0), breaks=c(-40,-20,0,20,40,60,80,100,120,140,160,180,200)) +
  ggtitle(label = "M1068I ") +
  geom_point(data=subset(pzq, Sample=="M1068I "), aes(x=(Conc), y=Response), size=1.2, shape=1, alpha=0.75, color="grey50") +
  theme( panel.grid = element_blank(), legend.title = element_blank(), 
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         plot.title = element_text(face="bold", color="black", size=8, ,hjust = 1),
         axis.title.x = element_text(face="bold", color="black", size=8)) + coord_cartesian(clip = "off")

F <- ggplot(data=subset(pzq_group, Sample=="S1341N "), aes(x=(Conc), y=Response_mean)) + theme_bw() + 
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE, size=0.6, color="black") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(side = "b", outside = TRUE, short = unit(.25,"mm"), mid = unit(0.5,"mm"), long = unit(1.5,"mm")) +
  xlab("Concentration (M)") +ylab("Response (% wild type)") + 
  scale_y_continuous(limits=c(-40,200), expand=c(0,0), breaks=c(-40,-20,0,20,40,60,80,100,120,140,160,180,200)) +
  ggtitle(label = "S1341N ") +
  geom_point(data=subset(pzq, Sample=="S1341N "), aes(x=(Conc), y=Response), size=1.2, shape=1, alpha=0.75, color="grey50") +
  theme( panel.grid = element_blank(), legend.title = element_blank(), 
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         plot.title = element_text(face="bold", color="black", size=8, ,hjust = 1),
         axis.title.x = element_text(face="bold", color="black", size=8)) + coord_cartesian(clip = "off")

G <- ggplot(data=subset(pzq_group, Sample=="F1399Y"), aes(x=(Conc), y=Response_mean)) + theme_bw() + 
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE, size=0.6, color="black") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(side = "b", outside = TRUE, short = unit(.25,"mm"), mid = unit(0.5,"mm"), long = unit(1.5,"mm")) +
  xlab("Concentration (M)") +ylab("Response (% wild type)") + 
  scale_y_continuous(limits=c(-40,200), expand=c(0,0), breaks=c(-40,-20,0,20,40,60,80,100,120,140,160,180,200)) +
  ggtitle(label = "F1399Y") +
  geom_point(data=subset(pzq, Sample=="F1399Y"), aes(x=(Conc), y=Response), size=1.2, shape=1, alpha=0.75, color="grey50") +
  theme( panel.grid = element_blank(), legend.title = element_blank(), 
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         plot.title = element_text(face="bold", color="black", size=8, ,hjust = 1),
         axis.title.x = element_text(face="bold", color="black", size=8)) + coord_cartesian(clip = "off")

H <- ggplot(data=subset(pzq_group, Sample=="V1405I "), aes(x=(Conc), y=Response_mean)) + theme_bw() + 
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE, size=0.6, color="black") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(side = "b", outside = TRUE, short = unit(.25,"mm"), mid = unit(0.5,"mm"), long = unit(1.5,"mm")) +
  xlab("Concentration (M)") +ylab("Response (% wild type)") + 
  scale_y_continuous(limits=c(-40,200), expand=c(0,0), breaks=c(-40,-20,0,20,40,60,80,100,120,140,160,180,200)) +
  ggtitle(label = "V1405I ") +
  geom_point(data=subset(pzq, Sample=="V1405I "), aes(x=(Conc), y=Response), size=1.2, shape=1, alpha=0.75, color="grey50") +
  theme( panel.grid = element_blank(), legend.title = element_blank(), 
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         plot.title = element_text(face="bold", color="black", size=8, ,hjust = 1),
         axis.title.x = element_text(face="bold", color="black", size=8)) + coord_cartesian(clip = "off")

I <- ggplot(data=subset(pzq_group, Sample=="L1476I "), aes(x=(Conc), y=Response_mean)) + theme_bw() + 
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE, size=0.6, color="black") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(side = "b", outside = TRUE, short = unit(.25,"mm"), mid = unit(0.5,"mm"), long = unit(1.5,"mm")) +
  xlab("Concentration (M)") +ylab("Response (% wild type)") + 
  scale_y_continuous(limits=c(-40,200), expand=c(0,0), breaks=c(-40,-20,0,20,40,60,80,100,120,140,160,180,200)) +
  ggtitle(label = "L1476I ") +
  geom_point(data=subset(pzq, Sample=="L1476I "), aes(x=(Conc), y=Response), size=1.2, shape=1, alpha=0.75, color="grey50") +
  theme( panel.grid = element_blank(), legend.title = element_blank(), 
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         plot.title = element_text(face="bold", color="black", size=8, ,hjust = 1),
         axis.title.x = element_text(face="bold", color="black", size=8)) + coord_cartesian(clip = "off")

J <- ggplot(data=subset(pzq_group, Sample=="T1624K "), aes(x=(Conc), y=Response_mean)) + theme_bw() + 
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE, size=0.6, color="black") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(side = "b", outside = TRUE, short = unit(.25,"mm"), mid = unit(0.5,"mm"), long = unit(1.5,"mm")) +
  xlab("Concentration (M)") +ylab("Response (% wild type)") + 
  scale_y_continuous(limits=c(-40,200), expand=c(0,0), breaks=c(-40,-20,0,20,40,60,80,100,120,140,160,180,200)) +
  ggtitle(label = "T1624K ") +
  geom_point(data=subset(pzq, Sample=="T1624K "), aes(x=(Conc), y=Response), size=1.2, shape=1, alpha=0.75, color="grey50") +
  theme( panel.grid = element_blank(), legend.title = element_blank(), 
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         plot.title = element_text(face="bold", color="black", size=8, ,hjust = 1),
         axis.title.x = element_text(face="bold", color="black", size=8)) + coord_cartesian(clip = "off")

K <- ggplot(data=subset(pzq_group, Sample=="R1843Q "), aes(x=(Conc), y=Response_mean)) + theme_bw() + 
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE, size=0.6, color="black") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(side = "b", outside = TRUE, short = unit(.25,"mm"), mid = unit(0.5,"mm"), long = unit(1.5,"mm")) +
  xlab("Concentration (M)") +ylab("Response (% wild type)") + 
  scale_y_continuous(limits=c(-40,200), expand=c(0,0), breaks=c(-40,-20,0,20,40,60,80,100,120,140,160,180,200)) +
  ggtitle(label = "R1843Q ") +
  geom_point(data=subset(pzq, Sample=="R1843Q "), aes(x=(Conc), y=Response), size=1.2, shape=1, alpha=0.75, color="grey50") +
  theme( panel.grid = element_blank(), legend.title = element_blank(), 
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         plot.title = element_text(face="bold", color="black", size=8, ,hjust = 1),
         axis.title.x = element_text(face="bold", color="black", size=8)) + coord_cartesian(clip = "off")

L <- ggplot(data=subset(pzq_group, Sample=="Y1544C"), aes(x=(Conc), y=Response_mean)) + theme_bw() + 
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE, size=0.6, color="black") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(side = "b", outside = TRUE, short = unit(.25,"mm"), mid = unit(0.5,"mm"), long = unit(1.5,"mm")) +
  xlab("Concentration (M)") +ylab("Response (% wild type)") + 
  scale_y_continuous(limits=c(-40,200), expand=c(0,0), breaks=c(-40,-20,0,20,40,60,80,100,120,140,160,180,200)) +
  ggtitle(label = "Y1544C") +
  geom_point(data=subset(pzq, Sample=="Y1544C"), aes(x=(Conc), y=Response), size=1.2, shape=1, alpha=0.75, color="grey50") +
  theme( panel.grid = element_blank(), legend.title = element_blank(), 
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         plot.title = element_text(face="bold", color="black", size=8, ,hjust = 1),
         axis.title.x = element_text(face="bold", color="black", size=8)) + coord_cartesian(clip = "off")

M <- ggplot(data=subset(pzq_group, Sample=="Q1670K "), aes(x=(Conc), y=Response_mean)) + theme_bw() + 
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE, size=0.6, color="black") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(side = "b", outside = TRUE, short = unit(.25,"mm"), mid = unit(0.5,"mm"), long = unit(1.5,"mm")) +
  xlab("Concentration (M)") +ylab("Response (% wild type)") + 
  scale_y_continuous(limits=c(-40,200), expand=c(0,0), breaks=c(-40,-20,0,20,40,60,80,100,120,140,160,180,200)) +
  ggtitle(label = "Q1670K ") +
  geom_point(data=subset(pzq, Sample=="Q1670K "), aes(x=(Conc), y=Response), size=1.2, shape=1, alpha=0.75, color="grey50") +
  theme( panel.grid = element_blank(), legend.title = element_blank(), 
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         plot.title = element_text(face="bold", color="black", size=8, ,hjust = 1),
         axis.title.x = element_text(face="bold", color="black", size=8)) + coord_cartesian(clip = "off")
```
### Merge plots
```
((A|B)/(C|D)/(E|F)/(G|H)) 

(A|I)/(J|K)/(L|M)/(J|K)
```
