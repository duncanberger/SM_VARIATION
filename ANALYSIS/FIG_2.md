# Figure 2
### Figure 2a: f3
```
# Import f3 values (from scikit-allel results)
f3_subset <- read.table("f3.csv", header=TRUE, sep=",")

# Plot f3 values for each population 
f3 <- ggplot(data=(f3_subset)) + 
  geom_linerange(aes(y=AVG,x=POPs, color=POPs, fill=POPs,ymin=AVG-SE, ymax=AVG+SE), width=.2 )+
  geom_point(aes(y=AVG,x=POPs, color=POPs, fill=POPs),outlier.alpha = 0, width=0.68) +
  geom_hline(yintercept = 0, color='grey50', alpha=0.5) +
  scale_y_continuous(expand=c(0,0), limits=c(-0.3,0.5), breaks=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.5)) +
  xlab("") + scale_color_manual(
  values=c("ADMIXED_TORO"="#EE6677","OG_TORO"="#EE6677","OG_MKAD"="#66CCEE",
  "ADMIXED_MKAD"="#66CCEE","ADMIXED_EA"='#CCBB44',"ADMIXED_GN"="#228833")) + 
  scale_fill_manual(values=c(
  "ADMIXED_TORO"="#EE6677","OG_TORO"="#EE6677",
  "OG_MKAD"="#66CCEE", "ADMIXED_MKAD"="#66CCEE",
  "ADMIXED_EA"='#CCBB44',"ADMIXED_GN"="#228833")) +
  theme_bw() + ylab("")+ 
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold", size=7),
        axis.text=element_text(face="bold", size=6),
        axis.text.y=element_blank(),
        axis.ticks = element_line(size = 0.25),
        axis.text.x=element_blank(),
        legend.position = "none",
        panel.grid = element_blank()) +
  coord_flip()
```
### FIG 2b: Patterson's D
```
# Import Patterson's D values (from scikit-allel results)
d_subset <- read.table("d.csv", header=TRUE, sep=",")

# Plot
d <- ggplot(data=(d_subset)) + 
  geom_linerange(aes(y=AVG,x=POPs, color=POPs, fill=POPs,ymin=AVG-SE, ymax=AVG+SE), width=.2 )+
  geom_point(aes(y=AVG,x=POPs, color=POPs, fill=POPs),outlier.alpha = 0, width=0.68) +
  geom_hline(yintercept = 0, color='grey50', alpha=0.5) +
  scale_y_continuous(expand=c(0,0), limits=c(-0.5,0.5), breaks=c(-0.5,-0.25,0,0.25,0.5)) +
  xlab("") + scale_color_manual(
  values=c("ADMIXED_TORO"="#EE6677","OG_TORO"="#EE6677","OG_MKAD"="#66CCEE",
  "ADMIXED_MKAD"="#66CCEE","ADMIXED_EA"='#CCBB44',"ADMIXED_GN"="#228833")) + 
  scale_fill_manual(values=c(
  "ADMIXED_TORO"="#EE6677","OG_TORO"="#EE6677",
  "OG_MKAD"="#66CCEE", "ADMIXED_MKAD"="#66CCEE",
  "ADMIXED_EA"='#CCBB44',"ADMIXED_GN"="#228833")) +
  theme_bw() + ylab("")+ 
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold", size=7),
        axis.text=element_text(face="bold", size=6),
        axis.ticks = element_line(size = 0.25),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        legend.position = "none",
        panel.grid = element_blank()) +
  coord_flip()
```
### FIG 2c: Ancestry_hmm
```
# Import ancestry_hmm results (taken from the log files of each run for each population and merged into a *.txt file)
timing_ahmm <- read.table("ahmm.summary.txt", header=FALSE)

# Plot the highest 10%
test1 <- ggplot() + 
  geom_half_violin(data=subset(timing_sprime, V7>=236192 & V8=="ea" & gen>=1),aes(x=V8,y=(gen)),fill="#66CCEE",color="white", alpha=1,side="r", trim = TRUE) +
  geom_half_violin(data=subset(timing_sprime, V7>=218993 & V8=="mkad" & gen>=1),aes(x=V8,y=(gen)),fill="#228833",color="white", alpha=1,side="r", trim = TRUE)  +
  geom_half_violin(data=subset(timing_sprime, V7>=991196 & V8=="gn" & gen>=1),aes(x=V8,y=(gen)),fill="#CCBB44",color="white", alpha=1,side="r", trim = TRUE) +
  geom_half_boxplot(data=subset(timing_sprime, V7>991196 & V8=="gn" & gen>=1),aes(x=V8,y=gen),color="grey50",  side="l",outlier.alpha = 0, width=0.32, color="grey25", alpha=0) +
  geom_half_boxplot(data=subset(timing_sprime, V7>=236192 & V8=="ea" & gen>=1),aes(x=V8,y=gen),color="grey50",  side="l",outlier.alpha = 0, width=0.32, color="grey25", alpha=0) +
  geom_half_boxplot(data=subset(timing_sprime, V7>218993 & V8=="mkad" & gen>=1),aes(x=V8,y=gen),color="grey50",side="l",outlier.alpha = 0, width=0.32, color="grey25", alpha=0) +
  coord_flip() + theme_bw() + ylab("")+ xlab("") +
  scale_y_continuous(expand=c(0,0), limits=c(0,500)) +
  theme(strip.background = element_blank(),
        axis.title = element_text(face="bold", size=7),
        axis.text=element_text(face="bold", size=6),
        #       axis.text.y=element_blank(),
        #       axis.text.x=element_blank(),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none",
        panel.grid = element_blank())

### FIG 2d: Ancestry_hmm

timing_sprime <- read.table("all.sprime.re.txt", sep='\t', header=FALSE)
timing_sprime$gen <- 1/((timing_sprime$V2/100)*0.99)

test2 <- ggplot() + 
  geom_half_violin(data=subset(timing_ahmm, V2=="gn"),aes(y=V1,x=V2),fill="#228833",color="white", alpha=1,  side="r", trim = TRUE) +
  geom_half_violin(data=subset(timing_ahmm, V2=="ea"),aes(y=V1,x=V2),fill="#CCBB44",color="white", alpha=1,   side="r", trim = TRUE) +
  geom_half_violin(data=subset(timing_ahmm, V2=="mkad"),aes(y=V1,x=V2), fill="#66CCEE",color="white", alpha=1, side="r", trim = TRUE) +
  geom_half_boxplot(data=subset(timing_ahmm, V2=="gn"),aes(y=V1,x=V2),color="grey50",  side="l",outlier.alpha = 0, width=0.32, color="grey25", alpha=0) +
  geom_half_boxplot(data=subset(timing_ahmm, V2=="ea"),aes(y=V1,x=V2),color="grey50",  side="l",outlier.alpha = 0, width=0.32, color="grey25", alpha=0) +
  geom_half_boxplot(data=subset(timing_ahmm, V2=="mkad"),aes(y=V1,x=V2),color="grey50",  side="l", outlier.alpha = 0, width=0.32, color="grey25", alpha=0) +
  coord_flip() +  theme_bw() + ylab("")+ xlab("") +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme(strip.background = element_blank(),
        axis.title = element_text(face="bold", size=7),
        axis.text=element_text(face="bold", size=6),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none",
        panel.grid = element_blank())




### FIG 2e: Pie charts

pie_chart <- read.table("counts.csv", sep=",")

Kocoge <- ggplot(data=subset(pie_chart, V4=="Kocoge")) + geom_bar(aes(y=V2, fill=V1,x=""),width = 1, stat = "identity", color="black") + coord_polar("y", start=0) + 
  theme_minimal() + scale_fill_manual(values=c("#edf8fb","#8c96c6","#6e016b")) + theme(axis.title.x = element_blank(), axis.text=element_blank(),axis.title.y = element_blank(),panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),legend.position = "none") 
Bugoto <- ggplot(data=subset(pie_chart, V4=="Bugoto")) + geom_bar(aes(y=V2, fill=V1,x=""),width = 1, stat = "identity", color="black") + coord_polar("y", start=0) + 
  theme_minimal() + scale_fill_manual(values=c("#edf8fb","#8c96c6","#6e016b")) + theme(axis.title.x = element_blank(), axis.text=element_blank(),axis.title.y = element_blank(),panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),legend.position = "none") 
Bwondha <- ggplot(data=subset(pie_chart, V4=="Bwondha")) + geom_bar(aes(y=V2, fill=V1,x=""),width = 1, stat = "identity", color="black") + coord_polar("y", start=0) + 
  theme_minimal() + scale_fill_manual(values=c("#edf8fb","#8c96c6","#6e016b")) + theme(axis.title.x = element_blank(), axis.text=element_blank(),axis.title.y = element_blank(),panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),legend.position = "none") 
Chole_beach <- ggplot(data=subset(pie_chart, V4=="Chole_beach")) + geom_bar(aes(y=V2, fill=V1,x=""),width = 1, stat = "identity", color="black") + coord_polar("y", start=0) + 
  theme_minimal() + scale_fill_manual(values=c("#edf8fb","#8c96c6","#6e016b")) + theme(axis.title.x = element_blank(), axis.text=element_blank(),axis.title.y = element_blank(),panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),legend.position = "none") 
IS_TOP <- ggplot(data=subset(pie_chart, V4=="IS_TOP")) + geom_bar(aes(y=V2, fill=V1,x=""),width = 1, stat = "identity", color="black") + coord_polar("y", start=0) + 
  theme_minimal() + scale_fill_manual(values=c("#edf8fb","#8c96c6","#6e016b")) + theme(axis.title.x = element_blank(), axis.text=element_blank(),axis.title.y = element_blank(),panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),legend.position = "none") 
IS_BOT <- ggplot(data=subset(pie_chart, V4=="IS_BOT")) + geom_bar(aes(y=V2, fill=V1,x=""),width = 1, stat = "identity", color="black") + coord_polar("y", start=0) + 
  theme_minimal() + scale_fill_manual(values=c("#edf8fb","#8c96c6","#6e016b")) + theme(axis.title.x = element_blank(), axis.text=element_blank(),axis.title.y = element_blank(),panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),legend.position = "none") 
Musubi <- ggplot(data=subset(pie_chart, V4=="Musubi")) + geom_bar(aes(y=V2, fill=V1,x=""),width = 1, stat = "identity", color="black") + coord_polar("y", start=0) + 
  theme_minimal() + scale_fill_manual(values=c("#edf8fb","#8c96c6","#6e016b")) + theme(axis.title.x = element_blank(), axis.text=element_blank(),axis.title.y = element_blank(),panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),legend.position = "none") 
Luhama <- ggplot(data=subset(pie_chart, V4=="Luhama")) + geom_bar(aes(y=V2, fill=V1,x=""),width = 1, stat = "identity", color="black") + coord_polar("y", start=0) + 
  theme_minimal() + scale_fill_manual(values=c("#edf8fb","#8c96c6","#6e016b")) + theme(axis.title.x = element_blank(), axis.text=element_blank(),axis.title.y = element_blank(),panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),legend.position = "none") 
Nyamatongo <- ggplot(data=subset(pie_chart, V4=="Nyamatongo")) + geom_bar(aes(y=V2, fill=V1,x=""),width = 1, stat = "identity", color="black") + coord_polar("y", start=0) + 
  theme_minimal() + scale_fill_manual(values=c("#edf8fb","#8c96c6","#6e016b")) + theme(axis.title.x = element_blank(), axis.text=element_blank(),axis.title.y = element_blank(),panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),legend.position = "none") 
Sweya <- ggplot(data=subset(pie_chart, V4=="Sweya")) + geom_bar(aes(y=V2, fill=V1,x=""),width = 1, stat = "identity", color="black") + coord_polar("y", start=0) + 
  theme_minimal() + scale_fill_manual(values=c("#edf8fb","#8c96c6","#6e016b")) + theme(axis.title.x = element_blank(), axis.text=element_blank(),axis.title.y = element_blank(),panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),legend.position = "none") 
Bukagabo_Beach <- ggplot(data=subset(pie_chart, V4=="Bukagabo_Beach")) + geom_bar(aes(y=V2, fill=V1,x=""),width = 1, stat = "identity", color="black") + coord_polar("y", start=0) + 
  theme_minimal() + scale_fill_manual(values=c("#edf8fb","#8c96c6","#6e016b")) + theme(axis.title.x = element_blank(), axis.text=element_blank(),axis.title.y = element_blank(),panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),legend.position = "none") 


plot_grid(Kocoge,Bugoto, IS_TOP,IS_BOT,Musubi,Luhama,Nyamatongo,Sweya,Chole_beach, Bwondha,Bukagabo_Beach, nrow=4)



