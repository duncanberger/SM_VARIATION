# Supplementary figures 7-9
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
### Load tree
```
tree <- midpoint(read.newick("iqtree.nwk")) 

key <- read.table("890_meta.csv", header=TRUE, sep=",",comment.char = "")
variants_l <- read.table("phy_variants.csv", sep=",", header=FALSE)

key2 <- merge(key,variants_l, by.x = c("sample_ID"), by.y = c("V1"), all.y = TRUE)
key2$Simple5 <- key2$Simple2
key2$Simple6 <- ifelse(key2$V2!="1.00", key2$Simple5, NA)
key3 <- subset(key2, select = -c(sample_ID))
```
### Plot trees
```
key2 <- merge(key,variants_l, by.x = c("sample_ID"), by.y = c("V1"), all.y = TRUE, all.x=TRUE)
key3 <- subset(key2, select = -c(sample_ID))

HOM_1843 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Arg1843Gln" & V2=="HOM") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HET_1843 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Arg1843Gln" & V2=="HET") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HOM_1670 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Gln1670Lys" & V2=="HOM") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HET_1670 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Gln1670Lys" & V2=="HET") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HOM_1915 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Gln1915fs" & V2=="HOM") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HET_1915 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Gln1915fs" & V2=="HET") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HOM_1020 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Ile1020Ile" & V2=="HOM") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HET_1020 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Ile1020Ile" & V2=="HET") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HOM_102 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Leu102Val" & V2=="HOM") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HET_102 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Leu102Val" & V2=="HET") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HOM_1476 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Leu1476Ile" & V2=="HOM") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HET_1476 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Leu1476Ile" & V2=="HET") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HOM_1068 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Met1068Ile" & V2=="HOM") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HET_1068 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Met1068Ile" & V2=="HET") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HOM_1 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Met1fs" & V2=="HOM") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HET_1 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Met1fs" & V2=="HET") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HOM_1399 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Phe1399Tyr" & V2=="HOM") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HET_1399 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Phe1399Tyr" & V2=="HET") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HOM_1341 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Ser1341Asn" & V2=="HOM") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HET_1341 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Ser1341Asn" & V2=="HET") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HOM_1005 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Thr1005Ile" & V2=="HOM") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HET_1005 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Thr1005Ile" & V2=="HET") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HOM_105 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Thr105Asn" & V2=="HOM") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HET_105 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Thr105Asn" & V2=="HET") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HOM_1624 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Thr1624Lys" & V2=="HOM") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HET_1624 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Thr1624Lys" & V2=="HET") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HOM_1554 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Tyr1554Cys" & V2=="HOM") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HET_1554 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Tyr1554Cys" & V2=="HET") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HOM_1405 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Val1405Ile" & V2=="HOM") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HET_1405 <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="p.Val1405Ile" & V2=="HET") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HOM_DEL <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="DEL" & V2=="HOM") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")

HET_DEL <- ggtree(tree, layout="ape", aes(color=Simple5), size=0.2) %<+% subset(key3, V3=="DEL" & V2=="HET") +
  geom_tippoint(size=0.25)+
  scale_color_manual(values=c("Tororo"="#EE6677","Lake Albert"="#4477AA","Mayuge"="#66CCEE","Rodhaini"="black","SA"="#AA3377","Kenya"="#57d9b3","Tanzania"="#228833",  "Islands"="#CCBB44", "Guadeloupe"="#ff9d00","Senegal"="#aa3377","Puerto Rico"="#6b3702",  "Cameroon"="#e783eb"), na.value="grey90") +
  theme(legend.position = "none")
```
### Merge plots
```
test_b1 <- plot_grid(HET_102,HOM_102,
          HET_105,HET_105,
          HET_1005,HOM_1005,
          HET_1068,HOM_1068,ncol=2,
          labels=c("p.L102V","p.L102V",
                   "p.T105N","p.T105NXXXXXX",
                   "p.T1005I","p.T1005I",
                   "p.M1068I","p.M1068I"), label_size = 8)

test_b2 <- plot_grid(HET_1341,HOM_1341,
                     HET_1399,HOM_1399,
                     HET_1405,HOM_1405,
                     HET_1476,HOM_1476,
                     ncol=2,
                     labels=c(
                       "p.S1341N","p.S1341N",
                       "p.F1399Y","p.F1399Y",
                       "p.V1405I","p.V1405I",
                       "p.L1476I","p.L1476I"), label_size = 8)

test_b3 <-plot_grid(HET_1554,HOM_1554,
          HET_1624,HOM_1624,
          HET_1670,HOM_1670,
          HET_1843,HOM_1843,
          ncol=2,
          labels=c("p.Y1554C","p.Y1554CXXXX",
                   "p.T1624K","p.T1624KXXXX",
                   "p.Q1670Lys","p.Q1670KXXXX",
                   "p.R1843Q","p.R1843Q"), label_size = 8)

test_b4 <- plot_grid(HET_1,HOM_1,
          HET_1915,HOM_1915,
          HET_DEL,HOM_DEL,
          HET_DEL,HOM_DEL,
          ncol=2,
          labels=c("p.M1fs","p.M1fs",
                   "p.Q1915fs","p.Q1915fs",
                   "~150 kb deletion","~150 kb deletion","REM","REM"), label_size = 8)

```




