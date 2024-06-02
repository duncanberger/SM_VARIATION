# Supplementary Figure 1: Diversity and differentiation of Lake Victoria *Schistosoma mansoni* populations
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
#### Supplementary figure 1a: Nucleotide diversity
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
#### Supplementary figure 1b: Fixation index
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
