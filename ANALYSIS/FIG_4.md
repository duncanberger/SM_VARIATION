# Figure 4: Relatedness between *Schistosoma mansoni* accessions from endemic regions
### Load libraries
```
# Load libraries
library("ggtern")
library("ggplot2")
library("reshape2")
library("dplyr")
library("ggtree")
```
### FIG 4a: Ternary plot
```
ngs <- read.table("ngsrelated.new2.txt", header=TRUE)
ngs2 <- read.table("both.kin.txt", header=TRUE)
host <- read.table("host3.list", header=FALSE)
ngs_q <- (merge(host,ngs2, all=TRUE, by.x = "V2", by.y='ida'))
ngs_q10 <- (merge(host,ngs_q, all=TRUE, by.x = "V2", by.y='idb'))
ngs_q10$match <-ifelse(ngs_q2$V3.y == ngs_q2$V3.x, "match", 
                      ifelse(ngs_q2$V3.y != ngs_q2$V3.x, "fail","else"))

ngs_q10$KIN <-ifelse(ngs_q10$J9 <0.6 & ngs_q10$J7>0.05 & ngs_q10$J8>0.2, "FS", 
                     ifelse(ngs_q10$J9 <0.6 & ngs_q10$J7<=0.2 & ngs_q10$J8>0.2, "HS",
                            ifelse(ngs_q10$J9 <0.8 & ngs_q10$J7<=0.2 & ngs_q10$J8>0.2, "FC", "X")))

ngs_q10$THET <-ifelse(ngs_q10$theta >0.354, "MZ",
                      ifelse(ngs_q10$J7>0.7, "MZ", 
                             ifelse(ngs_q10$theta >0.177 & ngs_q10$J8>0.05, "FD", 
                                    ifelse(ngs_q10$theta >0.0884 & ngs_q10$theta<0.177 & ngs_q10$J8>0.05 , "SD",
                                           ifelse(ngs_q10$theta >0.0442 & ngs_q10$theta <0.0884 & ngs_q10$J8>0.05, "TTD", "X")))))

ggtern(data=subset(ngs_q10),aes(J9,J8,J7)) + 
  scale_color_manual(values=c("#974181","#7daf5c","#d57946","#4e9ab7","grey50","grey50")) +
  theme_bw(base_size=8) +
  theme_showarrows() +
  geom_mask() +
  geom_point(size=1.25, aes(fill=THET,color=THET, shape=as.factor(match))) +
  Tlab(expression(K[1]), labelarrow = "K_1") + 
  Llab(expression(K[0]), labelarrow = "K_2") +
  Rlab(expression(K[2]), labelarrow = "K_0") 
```
### FIG 4c: Pre- vs Post-treatment nucleotide diversity
```
infra_pi<- read.table('auto_pi.pos2.txt', header=TRUE)  

infra_ND <- ggplot(data=subset(infra_pi, pop2!="pop2"), aes(x=as.factor(pop), y=log10(as.numeric(avg_pi)), fill=pop2, color=pop2)) + 
  stat_halfeye(adjust=1, justification=-0.5,width=0.4, alpha=0.5) +
  theme_bw() + 
  scale_y_continuous(expand=c(0,0), breaks=c(-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1), limits=c(-4,10)) +
  coord_cartesian(ylim = c(-4.5, -1)) +
  xlab("TEST") +
  scale_color_manual(values=c("#5C5992FF","#F0D77BFF")) +
  scale_fill_manual(values=c("#5C5992FF","#F0D77BFF")) +
  geom_boxplot(aes(fill=pop2, color=pop2),width=0.25, outlier.alpha = 0, alpha=0) +
  labs(y=expression(bold(-log[10]*(bolditalic("\U03C0"))))) +
  theme(panel.grid=element_blank(), 
        legend.text = element_text(face="bold"), 
        legend.position = "none", 
        axis.text = element_text(face = "bold", size=7, color="black"),
        axis.title = element_text(face = "bold", size=9, color="black"))
```
