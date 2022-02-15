# Figure 1
```
# Load libraries
library("phytools")
library("phangorn")
library("ape")
library("ggplot2")
library("reshape2")
library("dplyr")
library("ggtree")
```
### Figure 1b: Phylogeny
```
# Import metadata and IQTREE newick file
key <- read.table("890_meta.csv", header=TRUE, sep=",",comment.char = "")
tree <- midpoint(read.newick("iqtree.nwk")) 

# Plot tree (colors may be different than final figure)
ggtree::ggtree((tree), layout = "ape",aes(color=(Simple2)), size=0.2) %<+% subset(key) +
  scale_color_manual(values=c(
    "Tororo"="#EE6677",
    "Lake Albert"="#4477AA",
    "Mayuge"="#66CCEE",
    "Rodhaini"="black",
    "SA"="#AA3377",
    "Kenya"="#57d9b3",
    "Tanzania"="#228833", 
    "Islands"="#CCBB44",
    "WA"="#ff9d00"),
    na.value="grey85") +
  theme(legend.position = "none") +
  geom_treescale(linesize = 0.5,width=0.05, color="grey85", x = 0.5)
```
### Figure 1c: Principal component analysis 
```
# Import metadata, import eigenvectors and eigenvalues from PLINK results
key <- read.table("890_meta.csv", header=TRUE, sep=",",comment.char = "")
df<- read.delim("pruned.eigenvec", header=TRUE, sep="\t")
df_2 <- (merge(key, df, all=TRUE, by.y = "IID", by.x='sample_ID_old')) 
df_3 <- df_2[!is.na(df_2$PC1), ]
eigens<-read.delim("pruned.eigenval",header=F)
sum_eigs<-sum(eigens$V1)
sum_eigs<-lapply(eigens$V1,function(x){
  rt<-(x/sum_eigs)*100
  rt<-round(rt)
  return(rt)
})

# Plot PCA
pcaz <- ggplot(data=subset(df_3),aes((PC1),(PC2))) +
  geom_point(data=subset(df_3, Simple2!="Tanzania"), aes(fill=Simple2, shape=as.factor(study_ID)),
             size = 2, alpha=1, color="black") +
  geom_point(data=subset(df_3, Simple2=="Tanzania"), aes(fill=Simple2, shape=as.factor(study_ID)),
             size = 2, alpha=1, color="black") +
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
        legend.position = "none", 
        axis.text = element_text(face = "bold", size=7, color="black"),
        axis.title = element_text(face = "bold", size=9, color="black"))
```
#### FIG 1e: Nucleotide diversity
```
# Import PIXY results
pi_all <- read.table("combined.pixy.5000.pop.CHR_pi.txt", header=TRUE)

# Subsample to 5000 windows per-population (for plotting clarity)
pi_all_subset2 <- pi_all %>% group_by(pop)  %>%  sample_n(size=5000)
pi_all_subset2$V1 <- factor(pi_all_subset2$pop, levels=c("LAB", "TORO", "MKAD","EA","GN"))
labz <- c("Lake Albert",  "Eastern Uganda", "Southern Uganda","Koome Islands", "Northern Tanzania")
level_order <- c('LAB', 'TORO', 'MKAD',"EA","GN")

# Plot for all populations
piz <- ggplot(data=pi_all,aes(x=factor(pop, level = level_order),y=log10(avg_pi))) + 
  geom_point(data=pi_all_subset2, aes(color=pop),position = position_jitter(width = .12), size = .05) +
  stat_halfeye(adjust=1, justification=-0.5,width=0.5,aes(fill=pop)) +
  theme_bw() + 
  scale_y_continuous(expand=c(0,0), breaks=c(-7,-6,-5,-4,-3,-2,-1)) +
  coord_cartesian(ylim = c(-6, -1)) +
  xlab("TEST") +
  geom_boxplot(width=0.12, outlier.alpha = 0, color="black", alpha=0, fill="black") +
  scale_x_discrete(labels= labz) +
  labs(y=expression(bold(-log[10]*(bolditalic("\U03C0"))))) +
  theme(legend.position = "none",
        axis.text=element_text(size=7, face="bold", color="black"),
        axis.title=element_text(size=9, face="bold", color="black"),
        panel.grid = element_blank()) +
  scale_fill_manual(values=c(
    "EA"="#CCBB44",
    "LAB"="#4477AA",
    "MKAD"="#66CCEE",
    "Rodhaini"="black",
    "WA"="#AA3377",
    "GN"="#228833", 
    "TORO"="#EE6677",
    "SA"="AA3377",
    na.value="grey75")) +
  scale_color_manual(values=c(
    "EA"="#CCBB44",
    "LAB"="#4477AA",
    "MKAD"="#66CCEE",
    "Rodhaini"="black",
    "WA"="#AA3377",
    "GN"="#228833", 
    "TORO"="#EE6677",
    "SA"="AA3377",
    na.value="grey75"))
```
#### Figure 1e: Fixation index
```
# Import PIXY results (renamed files)
fst_popgen <- read.table("combined.pixy.5000.pop.CHR_fst.txt", header=FALSE)

# Subset for each pairs of populations and correct values <0 to 0
fst_popgen_2 <- subset(fst_popgen, V6!="NaN" & V6!="NA")
fst_popgen_2[fst_popgen_2 < 0] <- 0
fst_popgen_2_MKAD_EA <- subset(fst_popgen_2, V1=="MKAD" & V2=="EA")
fst_popgen_2_MKAD_GN <- subset(fst_popgen_2, V1=="MKAD" & V2=="GN")
fst_popgen_2_MKAD_TORO <- subset(fst_popgen_2, V1=="MKAD" & V2=="TORO")
fst_popgen_2_MKAD_LAB <- subset(fst_popgen_2, V1=="MKAD" & V2=="LAB")
fst_popgen_2_EA_GN <- subset(fst_popgen_2, V1=="EA" & V2=="GN")
fst_popgen_2_EA_TORO <- subset(fst_popgen_2, V1=="EA" & V2=="TORO")
fst_popgen_2_EA_LAB <- subset(fst_popgen_2, V1=="EA" & V2=="LAB")
fst_popgen_2_GN_TORO <- subset(fst_popgen_2, V1=="GN" & V2=="TORO")
fst_popgen_2_GN_LAB <- subset(fst_popgen_2, V1=="GN" & V2=="LAB")
fst_popgen_2_LAB_TORO <- subset(fst_popgen_2, V1=="LAB" & V2=="TORO")

# Calculate mean and median (example for one comparison only)
mean(fst_popgen_2_MKAD_EA$V6) # 0.002411672
mean(fst_popgen_2_MKAD_GN$V6) # 0.004077169
mean(fst_popgen_2_MKAD_TORO$V6) # 0.02636811
mean(fst_popgen_2_MKAD_LAB$V6) # 0.03472619
mean(fst_popgen_2_EA_GN$V6) # 0.003823846
mean(fst_popgen_2_EA_TORO$V6) # 0.02490525
mean(fst_popgen_2_EA_LAB$V6) # 0.03127308
mean(fst_popgen_2_GN_TORO$V6) # 0.03576755
mean(fst_popgen_2_GN_LAB$V6) # 0.04589715
mean(fst_popgen_2_LAB_TORO$V6) # 0.02703682
bstrap_means <- c()
bstrap_medians <- c()

# Calculate bootstrap means and medians
for (i in 1:100) { 
  bstrap_medians <- c(bstrap_medians,mean(sample(fst_popgen_2_LAB_TORO$V6,size=length(fst_popgen_2_LAB_TORO$V6),replace=TRUE)))
  bstrap_means <- c(bstrap_means,mean(sample(fst_popgen_2_LAB_TORO$V6,size=length(fst_popgen_2_LAB_TORO$V6),replace=TRUE)))
}
quantile(bstrap_medians,c(0.025,0.975))

# Manually copy FST values into csv file and import
fst_lot <- read.table("fst.csv", sep=",")

# Make FST plot
fst_plot <- ggplot(fst_lot, aes(V1, (V2), fill= V3)) + 
  geom_tile() + theme_bw()+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  scale_fill_viridis(option="magma", begin=0.5, end=0.05) +
  labs(y=expression(bold(-log[10]*(bolditalic("\U03C0"))))) +
  theme(
    axis.text=element_text(size=7, face="bold", color="black"),
    axis.title=element_text(size=9, face="bold", color="black"),
    panel.grid = element_blank())
```
### Figure 1f: ADMIXTURE
```
# Import metadata
key <- read.table("890_meta.csv", header=TRUE, sep=",",comment.char = "")

# Import theme
Admixture_theme <- theme(panel.grid=element_blank(), 
                         axis.ticks.x = element_blank(),
                         axis.text.y = element_text(face="bold", color="black"),
                         axis.text.x = element_blank(),
                         strip.text=element_text(face="bold"),
                         axis.title.y=element_text(face="bold",size=9),
                         panel.background = element_blank(),
                         strip.background = element_blank(),
                         panel.border = element_rect(color="black",fill=NA))

# Import ADMIXTURE results (manually re-ordered in excel 
admix <- read.table("admix.5.txt", header=FALSE, sep='\t')
admix_1a <- melt(admix,id.vars = c("V1","V8","V9"))
admix_2a <- merge(admix_1a, key, by.x="V1", by.y="sample_ID")
admix_2a$V1 <- reorder(admix_2a$V1, admix_2a$V9) 

# Plot ADMIXTURE results
ad_5 <- ggplot(data=subset(admix_2a)) +
  geom_col(aes(x=as.factor(V1), y=as.numeric(value), fill=variable),
           position="fill",width=1,show.legend = F, alpha=1)  + 
  facet_grid(.~Simple3,scales="free", space = "free_x") + 
  xlab("") + ylab("") +
  scale_fill_manual(values=c("V2"="black","V3"="#CCBB44","V4"="#66CCEE","V5"="#4477AA","V6"="#AA3377",
                             "V7"="purple","V8"="#AA3377","V9"="grey")) +
  scale_color_manual(values=c("V2"="black","V3"="#CCBB44","V4"="#66CCEE","V5"="#4477AA","V6"="#AA3377",
                              "V7"="purple","V8"="#AA3377","V9"="grey")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() + Admixture_theme + 
  theme(strip.text = element_blank(),
        axis.text.y=element_text(size=7, face="bold", color="black"),
        panel.border = element_blank(),
        panel.spacing = unit(0.3, "lines"))
```
