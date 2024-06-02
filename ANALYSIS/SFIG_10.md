# Supplementary Figure 10
### Load libraries
```
# Load libraries
library("ggplot2")
library("reshape2")
library("dplyr")
```
### Load and summarize data
```
sub_cov <- read.table("chr3.cov")
mean_covv <- read.table("3.cov.all")
struc_ex <- read.table("struc_246790.txt")
all_cov_win <- merge(mean_covv, sub_cov, by.x = c("V1"), by.y = c("V1"))
all_cov_win$diff <- all_cov_win$V5/all_cov_win$V2.x
```
### Plot structural variants in specific regions
```
MK0130 <- ggplot() + 
  geom_hline(yintercept = 1, color="grey50", alpha=0.8) +
  geom_rect(aes(xmin=2.627353,xmax=2.774807, ymin=0,ymax=10), fill="lightblue",alpha=0.25) +
  geom_rect(aes(xmin=2.949450,xmax=3.027139, ymin=0,ymax=10), fill="lightblue",alpha=0.25) +
  geom_rect(aes(xmin=3.132321,xmax=3.339629, ymin=0,ymax=10), fill="plum1",alpha=0.25) +
  geom_label(aes(x=2.70108, y=8.5), label="Smp_246790.5",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_label(aes(x=2.988295, y=8.5), label="Smp_345310.1",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_label(aes(x=3.235975, y=8.5), label="207 kb deletion",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_line(data=subset(all_cov_win, V1=="schisto_infrapop7718860" ),
            aes(x=((V4)/1000000),y=(V5/V2.x)), size=0.3, color="black")  + 
  scale_x_continuous(expand=c(0,0), limits=c(2.5,3.4)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,10)) +
  theme_bw() + xlab("Position (Mb)") + ylab("Normalized depth") +
  theme(panel.grid=element_blank(), 
        legend.position = "none", 
        axis.text = element_text(face = "bold", size=6, color="black"),
        axis.title = element_text(face = "bold", size=7, color="black"))


MK0025 <- ggplot() + 
  geom_hline(yintercept = 1, color="grey50", alpha=0.8) +
  geom_rect(aes(xmin=2.627353,xmax=2.774807, ymin=0,ymax=10), fill="lightblue",alpha=0.25) +
  geom_rect(aes(xmin=2.949450,xmax=3.027139, ymin=0,ymax=10), fill="lightblue",alpha=0.25) +
  geom_rect(aes(xmin=3.132321,xmax=3.339629, ymin=0,ymax=10), fill="plum1",alpha=0.25) +
  geom_label(aes(x=2.70108, y=8.5), label="Smp_246790.5",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_label(aes(x=2.988295, y=8.5), label="Smp_345310.1",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_label(aes(x=3.235975, y=8.5), label="207 kb deletion",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_line(data=subset(all_cov_win, V1=="schisto_infrapop7718693" ),
            aes(x=((V4)/1000000),y=(V5/V2.x)), size=0.3, color="black")  + 
  scale_x_continuous(expand=c(0,0), limits=c(2.5,3.4)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,10)) +
  theme_bw() + xlab("Position (Mb)") + ylab("Normalized depth") +
  theme(panel.grid=element_blank(), 
        legend.position = "none", 
        axis.text = element_text(face = "bold", size=6, color="black"),
        axis.title = element_text(face = "bold", size=7, color="black"))

GN0001 <- ggplot() + 
  geom_hline(yintercept = 1, color="grey50", alpha=0.8) +
  geom_rect(aes(xmin=2.627353,xmax=2.774807, ymin=0,ymax=10), fill="lightblue",alpha=0.25) +
  geom_rect(aes(xmin=2.949450,xmax=3.027139, ymin=0,ymax=10), fill="lightblue",alpha=0.25) +
  geom_rect(aes(xmin=3.132321,xmax=3.339629, ymin=0,ymax=10), fill="plum1",alpha=0.25) +
  geom_label(aes(x=2.70108, y=8.5), label="Smp_246790.5",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_label(aes(x=2.988295, y=8.5), label="Smp_345310.1",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_label(aes(x=3.235975, y=8.5), label="207 kb deletion",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_line(data=subset(all_cov_win, V1=="6356STDY9699564" ),
            aes(x=((V4)/1000000),y=(V5/V2.x)), size=0.3, color="black")  + 
  scale_x_continuous(expand=c(0,0), limits=c(2.5,3.4)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,10)) +
  theme_bw() + xlab("Position (Mb)") + ylab("Normalized depth") +
  theme(panel.grid=element_blank(), 
        legend.position = "none", 
        axis.text = element_text(face = "bold", size=6, color="black"),
        axis.title = element_text(face = "bold", size=7, color="black"))

EA0084 <- ggplot() + 
  geom_hline(yintercept = 1, color="grey50", alpha=0.8) +
  geom_rect(aes(xmin=2.627353,xmax=2.774807, ymin=0,ymax=10), fill="lightblue",alpha=0.25) +
  geom_rect(aes(xmin=2.949450,xmax=3.027139, ymin=0,ymax=10), fill="lightblue",alpha=0.25) +
  geom_rect(aes(xmin=3.132321,xmax=3.339629, ymin=0,ymax=10), fill="plum1",alpha=0.25) +
  geom_label(aes(x=2.70108, y=8.5), label="Smp_246790.5",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_label(aes(x=2.988295, y=8.5), label="Smp_345310.1",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_label(aes(x=3.235975, y=8.5), label="207 kb deletion",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_line(data=subset(all_cov_win, V1=="5582STDY7771053" ),
            aes(x=((V4)/1000000),y=(V5/V2.x)), size=0.3, color="black")  + 
  scale_x_continuous(expand=c(0,0), limits=c(2.5,3.4)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,10)) +
  theme_bw() + xlab("Position (Mb)") + ylab("Normalized depth") +
  theme(panel.grid=element_blank(), 
        legend.position = "none", 
        axis.text = element_text(face = "bold", size=6, color="black"),
        axis.title = element_text(face = "bold", size=7, color="black"))

MK0063 <- ggplot() + 
  geom_hline(yintercept = 1, color="grey50", alpha=0.8) +
  geom_rect(aes(xmin=2.627353,xmax=2.774807, ymin=0,ymax=10), fill="lightblue",alpha=0.25) +
  geom_rect(aes(xmin=2.949450,xmax=3.027139, ymin=0,ymax=10), fill="lightblue",alpha=0.25) +
  geom_rect(aes(xmin=3.132321,xmax=3.339629, ymin=0,ymax=10), fill="plum1",alpha=0.25) +
  geom_label(aes(x=2.70108, y=8.5), label="Smp_246790.5",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_label(aes(x=2.988295, y=8.5), label="Smp_345310.1",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_label(aes(x=3.235975, y=8.5), label="207 kb deletion",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_line(data=subset(all_cov_win, V1=="schisto_infrapop7718731" ),
            aes(x=((V4)/1000000),y=(V5/V2.x)), size=0.3, color="black")  + 
  scale_x_continuous(expand=c(0,0), limits=c(2.5,3.4)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,10)) +
  theme_bw() + xlab("Position (Mb)") + ylab("Normalized depth") +
  theme(panel.grid=element_blank(), 
        legend.position = "none", 
        axis.text = element_text(face = "bold", size=6, color="black"),
        axis.title = element_text(face = "bold", size=7, color="black"))

EA0094 <- ggplot() + 
  geom_hline(yintercept = 1, color="grey50", alpha=0.8) +
  geom_rect(aes(xmin=2.627353,xmax=2.774807, ymin=0,ymax=10), fill="lightblue",alpha=0.25) +
  geom_rect(aes(xmin=2.949450,xmax=3.027139, ymin=0,ymax=10), fill="lightblue",alpha=0.25) +
  geom_rect(aes(xmin=3.132321,xmax=3.339629, ymin=0,ymax=10), fill="plum1",alpha=0.25) +
  geom_label(aes(x=2.70108, y=8.5), label="Smp_246790.5",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_label(aes(x=2.988295, y=8.5), label="Smp_345310.1",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_label(aes(x=3.235975, y=8.5), label="207 kb deletion",color = "black",alpha=0, label.size = 0,size=2.5, fontface = "bold") +
  geom_line(data=subset(all_cov_win, V1=="5582STDY7771043" ),
            aes(x=((V4)/1000000),y=(V5/V2.x)), size=0.3, color="black")  + 
  scale_x_continuous(expand=c(0,0), limits=c(2.5,3.4)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,10)) +
  theme_bw() + xlab("Position (Mb)") + ylab("Normalized depth") +
  theme(panel.grid=element_blank(), 
        legend.position = "none", 
        axis.text = element_text(face = "bold", size=6, color="black"),
        axis.title = element_text(face = "bold", size=7, color="black"))

```
### Merge plots
```
sv_cov_all <- plot_grid(EA0094,MK0063,EA0084,GN0001,MK0130,MK0025, ncol=1, align="v")
ggsave(filename = "sv_cov_all.png",sv_cov_all, units = c("cm"), width = 18.5, height=22.5)
```
