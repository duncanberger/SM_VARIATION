# Supplementary Figure 14
### Load libraries
```
library("ggplot2")
library("reshape2")
library("dplyr")
```
### FIG 14a: Principal component analysis (by sequencing batch)
```
# Load data
key <- read.table("890_meta.csv", header=TRUE, sep=",",comment.char = "")
df<- read.delim("pruned.eigenvec", header=TRUE, sep="\t")

# Merge and get Eigens
df_2 <- (merge(key, df, all=TRUE, by.y = "IID", by.x='sample_ID_old')) 
df_3 <- df_2[!is.na(df_2$PC1), ]
eigens<-read.delim("pruned.eigenval",header=F)
sum_eigs<-sum(eigens$V1)
sum_eigs<-lapply(eigens$V1,function(x){
  rt<-(x/sum_eigs)*100
  rt<-round(rt)
  return(rt)
})

# Plot by sequencing batch
by_batch <- ggplot(data=subset(df_3),aes((PC1),(PC2))) +
  geom_point(data=subset(df_3, Simple2!="Tanzania" ), aes(fill=study2, shape=as.factor(study2)),
             size = 2, alpha=0.75, color="black") +
  geom_point(data=subset(df_3, Simple2=="Tanzania" ), aes(fill=study2, shape=as.factor(study2)),
             size = 2, alpha=0.75, color="black") +
  xlab(paste0("PC1 (",sum_eigs[[1]],"%)")) + 
  ylab(paste0("PC2 (",sum_eigs[[2]],"%)")) +
  theme_bw() +
  scale_shape_manual(values=c(23,22,24,21,25)) +
  scale_x_continuous(limits=c(-0.05,0.4), expand=c(0,0),
                     breaks=c(-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4)) +
  scale_y_continuous(limits=c(-0.2,0.30), expand=c(0,0), 
                     breaks=c(-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3)) +
  theme(panel.grid=element_blank(), 
        legend.text = element_text(face="bold"), 
        axis.text = element_text(face = "bold", size=7, color="black"),
        axis.title = element_text(face = "bold", size=9, color="black"))
```
### FIG 14b: Principal component analysis (by sequencing batch) and location sampled
```
pcaz <- ggplot(data=subset(df_3),aes((PC1),(PC2))) +
  geom_point(data=subset(df_3, Simple2!="Tanzania"), aes(fill=Simple2, shape=as.factor(study_ID)),
             size = 2, alpha=0.75, color="black") +
  geom_point(data=subset(df_3, Simple2=="Tanzania"), aes(fill=Simple2, shape=as.factor(study_ID)),
             size = 2, alpha=0.75, color="black") +
  xlab(paste0("PC1 (",sum_eigs[[1]],"%)")) + 
  ylab(paste0("PC2 (",sum_eigs[[2]],"%)")) +
  theme_bw() +
  scale_shape_manual(values=c(23,22,24,21)) +
  scale_fill_manual(values=c(
    "Tororo"="#EE6677",
    "Lake Albert"="#4477AA",
    "Mayuge"="#66CCEE",
    "Rodhaini"="black",
    "SA"="#AA3377",
    "Kenya"="#57d9b3",
    "Tanzania"="#228833", 
    "Islands"="#CCBB44",
    "WA"="#ff9d00")) +
  scale_x_continuous(limits=c(-0.05,0.4), expand=c(0,0),
                     breaks=c(-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4)) +
  scale_y_continuous(limits=c(-0.2,0.30), expand=c(0,0), 
                     breaks=c(-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3)) +
  theme(panel.grid=element_blank(), 
        legend.text = element_text(face="bold"), 
        axis.text = element_text(face = "bold", size=7, color="black"),
        axis.title = element_text(face = "bold", size=9, color="black"))
```
### FIG 14c: Coefficient of inbreeding by sequencing batch
```
ibc <- read.table("ibc.csv", header=TRUE, sep=",")
ibc_2 <- (merge(key, ibc, all=FALSE, by.y = "INDV.", by.x='sample_ID_old')) 

c_INB <- ggplot(data=subset(ibc_2, species=="Schistosoma mansoni")) + 
  geom_point(aes(y=F, x=study2, color=Simple4),position = position_jitter(width = .12), size = .5) +
  stat_halfeye(adjust=1, justification=-0.5,width=0.5,aes(fill=Simple4)) +
  theme_bw() + 
  scale_y_continuous(expand=c(0,0), breaks=c(-1,-0.5,0,0.5,1)) +
  coord_cartesian(ylim = c(-1, 1)) +
  xlab("TEST") +
  ylab("Coefficient of inbreeding (F)") +
  geom_boxplot(aes(y=F, x=study2, color=Simple4), width=0.12, outlier.alpha = 0, color="black", alpha=0, fill="black") +
  scale_color_manual(values=c("2_TORO"="#EE6677",
                              "1_LAB"="#4477AA",
                              "3_MAY"="#66CCEE",
                              "8_PR"="#6b3702",
                              "5_KEN"="#57d9b3",
                              "6_CAM"="#e783eb",
                              "4_TZ"="#228833", 
                              "4_ISL"="#CCBB44",
                              "7_SEN"="#aa3377",
                              "8_GUAD"="#ff9d00"),
                     na.value="grey85") +
  theme(legend.position = "none",
        axis.text=element_text(size=7, face="bold", color="black"),
        axis.title=element_text(size=9, face="bold", color="black"),
        panel.grid = element_blank()) 
```
### Merge plots
```
a_grid_sfb <- plot_grid(by_batch,by_batch, pcaz,pcaz, nrow=2, rel_widths = c(1,0.15))
b_grid_sfb <- plot_grid(a_grid_sfb, c_INB, nrow=2, rel_heights = c(1,0.6))
ggsave(filename = "bc_grid_sfb.svg",b_grid_sfb, units = c("cm"), width = 17, height=20, device = "svg")
```
