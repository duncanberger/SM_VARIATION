# Figure 5
```
# Load libraries
library("ggplot2")
```
### FIG 5a: Lollipop plots
```
# Import files


prots <- read.table("Smp_246790.domains2.txt", header=FALSE, sep="|")
trp <- read.table("new_allele_counts_trp.txt", header=FALSE, check.names = FALSE)
freq2 <- read.table("het_freq2.txt", header=FALSE)

X <- ggplot(subset(prots)) +
  geom_gene_arrow(fill = "white", aes(xmin=0,xmax=V3,y=2), 
                  arrow_body_height = unit(4, "mm")) +
  geom_subgene_arrow(aes(xsubmin = V2, xsubmax = V3, xmin = 0, xmax = 2268, y =2, fill = V1),
                     arrow_body_height = unit(4, "mm"), alpha=0.5) +
  xlim(-10,2300) + 
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title=element_blank(),
        #      axis.text.x = element_blank(),
        axis.ticks = element_blank())

Y <- ggplot() +
  geom_linerange(data=subset(subset(trp, V2=="1/1")),aes(x=as.numeric(V3),ymax=V1/550,ymin=0, color=V4), alpha=0.75) + 
  geom_point(data=subset(subset(trp, V2=="1/1" )), aes(x=as.numeric(V3),y=V1/550, color=V4), size=0.75, alpha=0.75) +
  xlim(-10,2300)  + scale_fill_manual(values=c(
    "conservative_inframe_deletion"="#1b9e77",
    "conservative_inframe_insertion"="#d95f02",
    "disruptive_inframe_deletion"="#7570b3",
    "disruptive_inframe_insertion"="#e7298a",
    "frameshift_variant"="#66a61e",
    "missense_variant"="#e6ab02",
    "splice_region_variant"="#a6761d",
    "synonymous_variant"="grey")) +
  scale_color_manual(values=c(
    "conservative_inframe_deletion"="#1b9e77",
    "conservative_inframe_insertion"="#d95f02",
    "disruptive_inframe_deletion"="#7570b3",
    "disruptive_inframe_insertion"="#e7298a",
    "frameshift_variant"="#66a61e",
    "missense_variant"="#e6ab02",
    "splice_region_variant"="#a6761d",
    "synonymous_variant"="grey")) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid = element_blank(),
        axis.title=element_blank(),
        axis.text =element_text(size=5, color="black", face="bold")) + 
  scale_y_reverse(expand=c(0,0)) +
  ylim(0,1) 

Z <- ggplot() +
  geom_linerange(data=subset(subset(trp, V2=="1/0" | V2=="0/1")),aes(x=V3,ymax=V1/550,ymin=0, color=V4), alpha=0.75) + 
  geom_point(data=subset(subset(trp, V2=="1/0" | V2=="0/1")), aes(x=V3,y=V1/550, color=V4), size=0.75, alpha=0.75) +
  xlim(-10,2300)  + 
  scale_fill_manual(values=c(
    "conservative_inframe_deletion"="#1b9e77",
    "conservative_inframe_insertion"="#d95f02",
    "disruptive_inframe_deletion"="#7570b3",
    "disruptive_inframe_insertion"="#e7298a",
    "frameshift_variant"="#66a61e",
    "missense_variant"="#e6ab02",
    "splice_region_variant"="#a6761d",
    "synonymous_variant"="grey")) +
  scale_color_manual(values=c(
    "conservative_inframe_deletion"="#1b9e77",
    "conservative_inframe_insertion"="#d95f02",
    "disruptive_inframe_deletion"="#7570b3",
    "disruptive_inframe_insertion"="#e7298a",
    "frameshift_variant"="#66a61e",
    "missense_variant"="#e6ab02",
    "splice_region_variant"="#a6761d",
    "synonymous_variant"="grey")) +
  theme(legend.position = "none",panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid = element_blank(),axis.title=element_blank(),
        axis.text =element_text(size=5, color="black", face="bold"))+
  ylim(0,1) +
  scale_y_reverse(limits=c(1,0)) 

plot_grid(Y,X,Z, ncol=1, align="v", axis="l", rel_heights = c(4,1,4))
```

