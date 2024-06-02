# Figure 3: Functional profiling of *Sm*.TRPM<sub>PZQ</sub> variants
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
PZQ_A <- ggplot(data=subset(pzq_group, Sample!="Q1670K " & Sample!="Y1544C" & Sample!="R1843Q " & Sample!="T1624K " ), 
                aes(x=(Conc), y=Response_mean, color=Sample, fill=Sample)) + theme_bw() + 
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE, size=0.005, alpha=0.8) +
  geom_errorbar(aes(ymin=Response_mean-Response_se, ymax=Response_mean+Response_se, x=(Conc)), width=.045) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(side = "b", outside = TRUE, short = unit(.25,"mm"), mid = unit(0.5,"mm"), long = unit(1.5,"mm")) +
  xlab("Concentration (M)") +ylab("Response (% wild type)") + 
  scale_y_continuous(limits=c(-20,160), expand=c(0,0), breaks=c(-20,0,20,40,60,80,100,120,140,160)) +
  geom_point(data=subset(pzq_group, Sample!="Q1670K " & Sample!="Y1544C" & Sample!="R1843Q " & Sample!="T1624K " ), 
             aes(x=(Conc), y=Response_mean, shape=Sample), size=1.4) +
  scale_color_manual(values=c("#008400","#008400","#008400","#008400","#008400","#008400","#008400","#008400","black")) +
  scale_fill_manual(values=c("#008400","#008400","#008400","#008400","#008400","#008400","#008400","#008400","black")) +
  scale_shape_manual(values=c(2,21,5,23,1,25,24,6,22)) +
  theme( panel.grid = element_blank(), legend.title = element_blank(), 
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         axis.line = element_line(colour = "black"),
         panel.border = element_blank(),
         plot.title = element_text(face="bold", color="black", size=8, ,hjust = 1),
         axis.title.x = element_text(face="bold", color="black", size=8)) + coord_cartesian(clip = "off")

PZQ_B <- ggplot(data=subset(pzq_group, Sample=="Q1670K " | Sample=="Y1544C" | Sample=="R1843Q " | Sample=="T1624K " ), 
                aes(x=(Conc), y=Response_mean, color=Sample, fill=Sample)) + theme_bw() + 
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE, size=0.005, alpha=0.8) +
  geom_errorbar(aes(ymin=Response_mean-Response_se, ymax=Response_mean+Response_se, x=(Conc)), width=.045) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(side = "b", outside = TRUE, short = unit(.25,"mm"), mid = unit(0.5,"mm"), long = unit(1.5,"mm")) +
  xlab("Concentration (M)") +ylab("Response (% wild type)") + 
  scale_y_continuous(limits=c(-20,160), expand=c(0,0), breaks=c(-20,0,20,40,60,80,100,120,140,160)) +
  geom_point(data=subset(pzq_group, Sample=="Q1670K " | Sample=="Y1544C" | Sample=="R1843Q " | Sample=="T1624K " ), 
             aes(x=(Conc), y=Response_mean, shape=Sample), size=1.8) +
  scale_color_manual(values=c("#ff0000","#ff8400","#ff8400","#ff0000")) +
  scale_fill_manual(values=c("#ff0000","#ff8400","#ff8400","#ff0000")) +
  scale_shape_manual(values=c(22,0,17,6)) +
  theme( panel.grid = element_blank(), legend.title = element_blank(), 
         axis.text.x=element_text(color="black", size=6),axis.text.y=element_text(color="black", size=6),
         axis.title.y = element_text(face="bold", color="black", size=8),
         axis.line = element_line(colour = "black"),
         panel.border = element_blank(),
         plot.title = element_text(face="bold", color="black", size=8,hjust = 1),
         axis.title.x = element_text(face="bold", color="black", size=8)) + coord_cartesian(clip = "off")
```
### Merge plots
```
PZQ_A|PZQ_B
```
