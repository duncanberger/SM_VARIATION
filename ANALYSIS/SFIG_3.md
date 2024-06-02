# Supplementary Figure 3: Genome-wide integrated haplotype scores (iHS) were calculated independently using unrelated accessions from four populations
### Load libraries
```
# Load libraries
library("ggplot2")
library("reshape2")
library("dplyr")
```
### Load data
```
ihs <- read.table("all2.ihs.nosibs.txt", header=FALSE)
```
### Summarize iHS scores
```
ihs$START <- (RoundTo(ihs$V3, multiple = 5000, FUN = floor)+1)
ihs$STOP <- (RoundTo(ihs$V3, multiple = 5000, FUN = ceiling))
ihs_median <- aggregate((abs(V8))~V1+START+V10+STOP, ihs, median)
ihs_length <- aggregate((abs(V8))~V1+START+V10+STOP, ihs, FUN=length)
ihs_summary0 <- merge(ihs_median, ihs_length, by.x = c("V1","START","V10","STOP"), by.y = c("V1","START","V10","STOP"))
names(ihs_summary0) <- c("CHROM", "START","POP", "STOP","IHS","COUNT")

ihs_summary0_EA <- subset(ihs_summary0, POP=="EA")
ihs_summary0_EAt <- (subset(ihs_summary0_EA, IHS>=2) %>% arrange(CHROM, START) %>% 
  group_by(CHROM, grp = cumsum(c(1, diff(START) > 50000))) %>% 
  filter(n() >= 2) %>% select(CHROM,START))
ihs_summary0_EAta <- merge(ihs_summary0_EA, ihs_summary0_EAt, by = c("CHROM","START"), all = TRUE)
ihs_summary0_EAta[is.na(ihs_summary0_EAta)] <- 0
ihs_summary0_EAta$COL <- ifelse(ihs_summary0_EAta$grp>0, "red",
                                ifelse(ihs_summary0_EAta$grp==0,"grey","grey"))
ihs_summary0_EAta_min <- aggregate(ihs_summary0_EAta$START, by = list(ihs_summary0_EAta$grp, ihs_summary0_EAta$CHROM), min)
ihs_summary0_EAta_max <- aggregate(ihs_summary0_EAta$START, by = list(ihs_summary0_EAta$grp, ihs_summary0_EAta$CHROM), max)
ihs_summary0_EAta_windows <- merge(ihs_summary0_EAta_min,ihs_summary0_EAta_max, by=c("Group.1","Group.2"),all = TRUE)
#write.table(ihs_summary0_EAta_windows,"ihs_summary0_EAta_windows.csv",row.names = FALSE, sep = ",")
```
### Plot iHS scores (per-location)
```
ihs_EA <- ggplot(data=subset(ihs_summary0_EAta, COUNT>10  & CHROM!="PAR1" & CHROM!="PAR2")) + 
  geom_point(aes(x=((START)/1000000),y=(IHS), color=COL), size=0.01) +
  scale_color_manual(values = rep(c("grey50","red"))) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,6), breaks=c(0.0,3.0,6.0)) +
  theme_bw() + ylab("iHS score") + xlab("Position (Mb)") +
  theme(legend.position="none",panel.grid = element_blank(),panel.spacing = unit(0.2, "lines"),
    strip.background = element_blank(),strip.text.x=element_text(face="bold", color="black", size=6),
    axis.text.x=element_text(face="bold", color="black", size=4), axis.text.y=element_text(face="bold", color="black", size=6), axis.title.x = element_text(face="bold", color="black", size=7),
    axis.title.y = element_text(face="bold", color="black", size=7)) +
  facet_grid(.~CHROM, space = "free_x", scales="free_x")

ihs_summary0_MKAD <- subset(ihs_summary0, POP=="MKAD")
ihs_summary0_MKADt <- (subset(ihs_summary0_MKAD, IHS>=2) %>% arrange(CHROM, START) %>% 
                       group_by(CHROM, grp = cumsum(c(1, diff(START) > 50000))) %>% 
                       filter(n() >= 2) %>% select(CHROM,START))
ihs_summary0_MKADta <- merge(ihs_summary0_MKAD, ihs_summary0_MKADt, by = c("CHROM","START"), all = TRUE)
ihs_summary0_MKADta[is.na(ihs_summary0_MKADta)] <- 0
ihs_summary0_MKADta$COL <- ifelse(ihs_summary0_MKADta$grp>0, "red",
                                ifelse(ihs_summary0_MKADta$grp==0,"grey","grey"))
ihs_summary0_MKADta_min <- aggregate(ihs_summary0_MKADta$START, by = list(ihs_summary0_MKADta$grp, ihs_summary0_MKADta$CHROM), min)
ihs_summary0_MKADta_max <- aggregate(ihs_summary0_MKADta$START, by = list(ihs_summary0_MKADta$grp, ihs_summary0_MKADta$CHROM), max)
ihs_summary0_MKADta_windows <- merge(ihs_summary0_MKADta_min,ihs_summary0_MKADta_max, by=c("Group.1","Group.2"),all = TRUE)

ihs_MKAD <- ggplot(data=subset(ihs_summary0_MKADta, COUNT>10  & CHROM!="PAR1" & CHROM!="PAR2")) + 
  geom_point(aes(x=((START)/1000000),y=(IHS), color=COL), size=0.01) +
  scale_color_manual(values = rep(c("grey50","red"))) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,6), breaks=c(0.0,3.0,6.0)) +
  theme_bw() + ylab("iHS score") + xlab("Position (Mb)") +
  theme(legend.position="none",panel.grid = element_blank(),panel.spacing = unit(0.2, "lines"),
        strip.background = element_blank(),strip.text.x=element_text(face="bold", color="black", size=6),
        axis.text.x=element_text(face="bold", color="black", size=4), axis.text.y=element_text(face="bold", color="black", size=6), axis.title.x = element_text(face="bold", color="black", size=7),
        axis.title.y = element_text(face="bold", color="black", size=7)) +
  facet_grid(.~CHROM, space = "free_x", scales="free_x")

ihs_summary0_GN <- subset(ihs_summary0, POP=="GN")
ihs_summary0_GNt <- (subset(ihs_summary0_GN, IHS>=2) %>% arrange(CHROM, START) %>% 
                         group_by(CHROM, grp = cumsum(c(1, diff(START) > 50000))) %>% 
                         filter(n() >= 2) %>% select(CHROM,START))
ihs_summary0_GNta <- merge(ihs_summary0_GN, ihs_summary0_GNt, by = c("CHROM","START"), all = TRUE)
ihs_summary0_GNta[is.na(ihs_summary0_GNta)] <- 0
ihs_summary0_GNta$COL <- ifelse(ihs_summary0_GNta$grp>0, "red",
                                  ifelse(ihs_summary0_GNta$grp==0,"grey","grey"))
ihs_summary0_GNta_min <- aggregate(ihs_summary0_GNta$START, by = list(ihs_summary0_GNta$grp, ihs_summary0_GNta$CHROM), min)
ihs_summary0_GNta_max <- aggregate(ihs_summary0_GNta$START, by = list(ihs_summary0_GNta$grp, ihs_summary0_GNta$CHROM), max)
ihs_summary0_GNta_windows <- merge(ihs_summary0_GNta_min,ihs_summary0_GNta_max, by=c("Group.1","Group.2"),all = TRUE)

ihs_GN <- ggplot(data=subset(ihs_summary0_GNta, COUNT>10  & CHROM!="PAR1" & CHROM!="PAR2")) + 
  geom_point(aes(x=((START)/1000000),y=(IHS), color=COL), size=0.01) +
  scale_color_manual(values = rep(c("grey50","red"))) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,6), breaks=c(0.0,3.0,6.0)) +
  theme_bw() + ylab("iHS score") + xlab("Position (Mb)") +
  theme(legend.position="none",panel.grid = element_blank(),panel.spacing = unit(0.2, "lines"),
        strip.background = element_blank(),strip.text.x=element_text(face="bold", color="black", size=6),
        axis.text.x=element_text(face="bold", color="black", size=4), axis.text.y=element_text(face="bold", color="black", size=6), axis.title.x = element_text(face="bold", color="black", size=7),
        axis.title.y = element_text(face="bold", color="black", size=7)) +
  facet_grid(.~CHROM, space = "free_x", scales="free_x")

ihs_summary0_TR <- subset(ihs_summary0, POP=="TORO")
ihs_summary0_TRt <- (subset(ihs_summary0_TR, IHS>=2) %>% arrange(CHROM, START) %>% 
                       group_by(CHROM, grp = cumsum(c(1, diff(START) > 50000))) %>% 
                       filter(n() >= 2) %>% select(CHROM,START))
ihs_summary0_TRta <- merge(ihs_summary0_TR, ihs_summary0_TRt, by = c("CHROM","START"), all = TRUE)
ihs_summary0_TRta[is.na(ihs_summary0_TRta)] <- 0
ihs_summary0_TRta$COL <- ifelse(ihs_summary0_TRta$grp>0, "red",
                                ifelse(ihs_summary0_TRta$grp==0,"grey","grey"))
ihs_summary0_TRta_min <- aggregate(ihs_summary0_TRta$START, by = list(ihs_summary0_TRta$grp, ihs_summary0_TRta$CHROM), min)
ihs_summary0_TRta_max <- aggregate(ihs_summary0_TRta$START, by = list(ihs_summary0_TRta$grp, ihs_summary0_TRta$CHROM), max)
ihs_summary0_TRta_windows <- merge(ihs_summary0_TRta_min,ihs_summary0_TRta_max, by=c("Group.1","Group.2"),all = TRUE)

ihs_TORO <- ggplot(data=subset(ihs_summary0_TRta, COUNT>10  & CHROM!="PAR1" & CHROM!="PAR2")) + 
  geom_point(aes(x=((START)/1000000),y=(IHS), color=COL), size=0.01) +
  scale_color_manual(values = rep(c("grey50","red"))) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,6), breaks=c(0.0,3.0,6.0)) +
  theme_bw() + ylab("iHS score") + xlab("Position (Mb)") +
  theme(legend.position="none",panel.grid = element_blank(),panel.spacing = unit(0.2, "lines"),
        strip.background = element_blank(),strip.text.x=element_text(face="bold", color="black", size=6),
        axis.text.x=element_text(face="bold", color="black", size=4), axis.text.y=element_text(face="bold", color="black", size=6), axis.title.x = element_text(face="bold", color="black", size=7),
        axis.title.y = element_text(face="bold", color="black", size=7)) +
  facet_grid(.~CHROM, space = "free_x", scales="free_x")
```
### Merge plots
```
ihs_plots<- plot_grid(ihs_MKAD,ihs_EA,ihs_GN,ihs_TORO, ncol=1, labels = c("a","b","c","d"), label_size = 11, align = "v")
ggsave(filename = "sfig_ihs.png",ihs_plots, units = c("cm"), width = 18.5, height=17.5)
```
