# This script analyses data to address the question: 
# How have responses of flowering time and fitness to germination timing evolved across the clade?
# Do species with greater compensatory flowering plasticity have greater fitness stability in the face of variable germination timing?
# Is there evidence of climate adaptation?

# Users, please specify your own path to load data in this script.

# load libraries
library(ape) # version 5.7.1
library(phytools) # version 1.9.16
library(ggnewscale)  # version 0.5.0
library(caper) # version 1.0.3
library(tidyverse) # version 2.0.0
library(ggpubr) # version 0.6.0
library(ggrepel) # version 0.9.3
library(cowplot) # version 1.1.1
library(car) # version 3.1.2

remotes::install_github("clauswilke/colorblindr")
library(colorblindr) # version 0.1.0

BiocManager::install("ggtree", force = TRUE) # version 3.6.2
library(ggtree) # https://yulab-smu.top/treedata-book/chapter7.html


#### Testing for phylogenetic signal in slopes ####
# Blomberg's K
# phylosig from phytools including standard error

# read in phylo
all.phylo <- read.tree("./Germination.Fitness/Raw.Data/tree_pruned.new")

# species we need
sp.list=c("Caulanthus_anceps","Caulanthus_coulteri","Caulanthus_inflatus",
          "Streptanthus_breweri","Streptanthus_diversifolius","Streptanthus_drepanoides",
          "Streptanthus_glandulosus","Streptanthus_insignis","Streptanthus_polygaloides",
          "Streptanthus_tortuosus")

# prune phylogeny to only includes our species
phylo=keep.tip(all.phylo, sp.list)

# testing for phylogenetic signal in slopes

# read in slopes dataframe, also includes standard error of slopes
model.slopes=read.csv("./pheno.slopes.csv")

# days.2.bud and sept.1.bud are direct opposites of each other, difference in p is due to randomization

phylosig(phylo, model.slopes$days.2.bud, method = "K", test = TRUE, se = model.slopes$days.2.bud.se, nsim=1000)
# K = 1.18966, p = 0.136
phylosig(phylo, model.slopes$sept.1.bud, method = "K", test = TRUE, se = model.slopes$sept.1.bud.se,nsim=1000)
# K = 1.18966, p = 0.151

# revised slopes
phylosig(phylo, model.slopes$days.2.bud.final, method = "K", test = TRUE, se = model.slopes$days.2.bud.final.se, nsim=1000)
# K = 1.18857, p = 0.143
phylosig(phylo, model.slopes$sept.1.bud.final, method = "K", test = TRUE, se = model.slopes$sept.1.bud.final.se,nsim=1000)
# K = 1.18857, p = 0.149

# days.2.flower and sept.1.flower are direct opposites of each other, difference in p is due to randomization

phylosig(phylo, model.slopes$days.2.flower, method = "K", test = TRUE, se = model.slopes$days.2.flower.se, nsim=1000)
# K = 1.22415, p = 0.145
phylosig(phylo, model.slopes$sept.1.flower, method = "K", test = TRUE, se = model.slopes$sept.1.flower.se,nsim=1000)
# K = 1.22415, p = 0.156

# phylogenetic signal of slopes for combinations of populations

phylosig(phylo, model.slopes$days.2.bud.CAAN1.CAIN3, method = "K", test = TRUE, se = model.slopes$days.2.bud.CAAN1.CAIN3.se, nsim=1000)
# K = 1.14555, p = 0.195
phylosig(phylo, model.slopes$sept.1.bud.CAAN1.CAIN3, method = "K", test = TRUE, se = model.slopes$sept.1.bud.CAAN1.CAIN3.se,nsim=1000)
# K = 1.14555, p = 0.199
phylosig(phylo, model.slopes$days.2.bud.CAAN1.CAIN4, method = "K", test = TRUE, se = model.slopes$days.2.bud.CAAN1.CAIN4.se, nsim=1000)
# K = 1.16093, p = 0.144
phylosig(phylo, model.slopes$sept.1.bud.CAAN1.CAIN4, method = "K", test = TRUE, se = model.slopes$sept.1.bud.CAAN1.CAIN4.se,nsim=1000)
# K = 1.16093, p = 0.175
phylosig(phylo, model.slopes$days.2.bud.CAAN2.CAIN3, method = "K", test = TRUE, se = model.slopes$days.2.bud.CAAN2.CAIN3.se, nsim=1000)
# K = 1.16645, p = 0.150
phylosig(phylo, model.slopes$sept.1.bud.CAAN2.CAIN3, method = "K", test = TRUE, se = model.slopes$sept.1.bud.CAAN2.CAIN3.se,nsim=1000)
# K = 1.16645, p = 0.167
phylosig(phylo, model.slopes$days.2.bud.CAAN2.CAIN4, method = "K", test = TRUE, se = model.slopes$days.2.bud.CAAN2.CAIN4.se, nsim=1000)
# K = 1.18198, p = 0.125
phylosig(phylo, model.slopes$sept.1.bud.CAAN2.CAIN4, method = "K", test = TRUE, se = model.slopes$sept.1.bud.CAAN2.CAIN4.se,nsim=1000)
# K = 1.18198, p = 0.123

phylosig(phylo, model.slopes$days.2.flower.CAAN1.CAIN3, method = "K", test = TRUE, se = model.slopes$days.2.flower.CAAN1.CAIN3.se, nsim=1000)
# K = 1.16982, p = 0.212
phylosig(phylo, model.slopes$sept.1.flower.CAAN1.CAIN3, method = "K", test = TRUE, se = model.slopes$sept.1.flower.CAAN1.CAIN3.se,nsim=1000)
# K = 1.14555, p = 0.191
phylosig(phylo, model.slopes$days.2.flower.CAAN1.CAIN4, method = "K", test = TRUE, se = model.slopes$days.2.flower.CAAN1.CAIN4.se, nsim=1000)
# K = 1.17559, p = 0.160
phylosig(phylo, model.slopes$sept.1.flower.CAAN1.CAIN4, method = "K", test = TRUE, se = model.slopes$sept.1.flower.CAAN1.CAIN4.se,nsim=1000)
# K = 1.17559, p = 0.159
phylosig(phylo, model.slopes$days.2.flower.CAAN2.CAIN3, method = "K", test = TRUE, se = model.slopes$days.2.flower.CAAN2.CAIN3.se, nsim=1000)
# K = 1.19008, p = 0.149
phylosig(phylo, model.slopes$sept.1.flower.CAAN2.CAIN3, method = "K", test = TRUE, se = model.slopes$sept.1.flower.CAAN2.CAIN3.se,nsim=1000)
# K = 1.19008, p = 0.141
phylosig(phylo, model.slopes$days.2.flower.CAAN2.CAIN4, method = "K", test = TRUE, se = model.slopes$days.2.flower.CAAN2.CAIN4.se, nsim=1000)
# K = 1.19488, p = 0.108
phylosig(phylo, model.slopes$sept.1.flower.CAAN2.CAIN4, method = "K", test = TRUE, se = model.slopes$sept.1.flower.CAAN2.CAIN4.se,nsim=1000)
# K = 1.19488, p = 0.120

# fitness slopes
fitness.slopes=read.csv("./fitness.slopes.2.csv")

phylosig(phylo, fitness.slopes$pflwr, method = "K", test = TRUE, se = fitness.slopes$pflwr.se,nsim=1000)
# K = 0.973273, p = 0.908
phylosig(phylo, fitness.slopes$seed_num_nb, method = "K", test = TRUE, se = fitness.slopes$seed_num_nb.se,nsim=1000)
# K = 1.06775, p = 0.842
phylosig(phylo, fitness.slopes$year1fit, method = "K", test = TRUE, se = fitness.slopes$year1fit.se,nsim=1000)
# K = 0.682429, p = 0.9997
phylosig(phylo, fitness.slopes$seed_mass, method = "K", se = fitness.slopes$seed_mass.se,test = TRUE, nsim=1000)
# K = 0.771932, p = 0.442
phylosig(phylo, fitness.slopes$pflwr.weighted, method = "K", test = TRUE, se = fitness.slopes$pflwr.weighted.se,nsim=1000)
# K = 0.950963, p = 0.905

# revised slopes

phylosig(phylo, fitness.slopes$pflwr.final, method = "K", test = TRUE, se = fitness.slopes$pflwr.final.se,nsim=1000)
# K = 0.951959, p = 0.896
phylosig(phylo, fitness.slopes$seed_num_nb.final, method = "K", test = TRUE, se = fitness.slopes$seed_num_nb.final.se,nsim=1000)
# K = 1.06552, p = 0.863
phylosig(phylo, fitness.slopes$year1fit.final, method = "K", test = TRUE, se = fitness.slopes$year1fit.final.se,nsim=1000)
# K = 0.681723, p = 0.998
phylosig(phylo, fitness.slopes$seed_mass.final, method = "K", se = fitness.slopes$seed_mass.final.se,test = TRUE, nsim=1000)
# K = 0.773631, p = 0.42

# phylogenetic signal of slopes for combinations of populations

phylosig(phylo, fitness.slopes$pflwr.CAAN1.CAIN3, method = "K", test = TRUE, se = fitness.slopes$pflwr.CAAN1.CAIN3.se, nsim=1000)
# K = 0.937076, p = 0.877
phylosig(phylo, fitness.slopes$pflwr.CAAN1.CAIN4, method = "K", test = TRUE, se = fitness.slopes$pflwr.CAAN1.CAIN4.se, nsim=1000)
# K = 0.935965, p = 0.870
phylosig(phylo, fitness.slopes$pflwr.CAAN2.CAIN3, method = "K", test = TRUE, se = fitness.slopes$pflwr.CAAN2.CAIN3.se, nsim=1000)
# K = 0.952129, p = 0.943
phylosig(phylo, fitness.slopes$pflwr.CAAN2.CAIN4, method = "K", test = TRUE, se = fitness.slopes$pflwr.CAAN2.CAIN4.se, nsim=1000)
# K = 0.951017, p = 0.920

phylosig(phylo, fitness.slopes$seed_num_nb.CAAN1.CAIN3, method = "K", test = TRUE, se = fitness.slopes$seed_num_nb.CAAN1.CAIN3.se,nsim=1000)
# K = 1.0327, p = 0.799
phylosig(phylo, fitness.slopes$seed_num_nb.CAAN1.CAIN4, method = "K", test = TRUE, se = fitness.slopes$seed_num_nb.CAAN1.CAIN4.se,nsim=1000)
# K = 0.971969, p = 0.897
phylosig(phylo, fitness.slopes$seed_num_nb.CAAN2.CAIN3, method = "K", test = TRUE, se = fitness.slopes$seed_num_nb.CAAN2.CAIN3.se,nsim=1000)
# K = 1.04099, p = 0.757
phylosig(phylo, fitness.slopes$seed_num_nb.CAAN2.CAIN4, method = "K", test = TRUE, se = fitness.slopes$seed_num_nb.CAAN2.CAIN4.se,nsim=1000)
# K = 0.977187, p = 0.885

phylosig(phylo, fitness.slopes$year1fit.CAAN1.CAIN3, method = "K", test = TRUE, se = fitness.slopes$year1fit.CAAN1.CAIN3.se, nsim=1000)
# K = 0.702078, p = 0.992
phylosig(phylo, fitness.slopes$year1fit.CAAN1.CAIN4, method = "K", test = TRUE, se = fitness.slopes$year1fit.CAAN1.CAIN4.se, nsim=1000)
# K = 0.669168, p = 0.993
phylosig(phylo, fitness.slopes$year1fit.CAAN2.CAIN3, method = "K", test = TRUE, se = fitness.slopes$year1fit.CAAN2.CAIN3.se, nsim=1000)
# K = 0.705075, p = 0.994
phylosig(phylo, fitness.slopes$year1fit.CAAN2.CAIN4, method = "K", test = TRUE, se = fitness.slopes$year1fit.CAAN2.CAIN4.se, nsim=1000)
# K = 0.678227, p = 0.6782270.992

phylosig(phylo, fitness.slopes$seed_mass.CAAN1.CAIN3, method = "K", test = TRUE, se = fitness.slopes$seed_mass.CAAN1.CAIN3.se,nsim=1000)
# K = 0.774558, p = 0.409
phylosig(phylo, fitness.slopes$seed_mass.CAAN1.CAIN4, method = "K", test = TRUE, se = fitness.slopes$seed_mass.CAAN1.CAIN4.se,nsim=1000)
# K = 0.747243, p = 0.421
phylosig(phylo, fitness.slopes$seed_mass.CAAN2.CAIN3, method = "K", test = TRUE, se = fitness.slopes$seed_mass.CAAN2.CAIN3.se,nsim=1000)
# K = 0.780504, p = 0.424
phylosig(phylo, fitness.slopes$seed_mass.CAAN2.CAIN4, method = "K", test = TRUE, se = fitness.slopes$seed_mass.CAAN2.CAIN4.se,nsim=1000)
# K = 0.750859, p = 0.408

#### Figure 5 #####

all.phylo <- ggtree::read.tree("./tree_pruned.new")

sp.list=c("Caulanthus_anceps","Caulanthus_coulteri","Caulanthus_inflatus",
          "Streptanthus_breweri","Streptanthus_diversifolius","Streptanthus_drepanoides",
          "Streptanthus_glandulosus","Streptanthus_insignis","Streptanthus_polygaloides",
          "Streptanthus_tortuosus")

# prune phylogeny to only includes our species
phylo=keep.tip(all.phylo, sp.list)

phylo$tip.label = c("Caulanthus inflatus","Caulanthus coulteri","Caulanthus anceps",
                    "Streptanthus glandulosus","Streptanthus insignis","Streptanthus polygaloides",
                    "Streptanthus diversifolius","Streptanthus tortuosus","Streptanthus breweri",
                    "Streptanthus drepanoides")

ggplot(phylo) + geom_tree() + theme_tree()

tree = ggtree(phylo, size =1)
tree

# read in slopes
all.slopes=read.csv("./all.slopes.csv")

pheno.slopes$sp = c("Caulanthus inflatus","Caulanthus coulteri","Caulanthus anceps",
                    "Streptanthus glandulosus","Streptanthus insignis","Streptanthus polygaloides",
                    "Streptanthus diversifolius","Streptanthus tortuosus","Streptanthus breweri",
                    "Streptanthus drepanoides")

# climate pc1
all.slopes$PC1 = c(-3.4434154,-1.6218771,-1.6381408,0.2230000,-0.9859651,1.3163704,0.1026840,4.1567596,0.7433119,1.1472724)

# add tip labels
lb = phylo$tip.label
tips = data.frame(label = lb, Species = paste(all.slopes$sp))

tree.2 = tree %<+% tips +
  geom_tiplab(aes(label = Species),as_ylab = TRUE)
tree.2

# get slopes we are interested in 

main.ms.slopes.fig = all.slopes[,c(3,11)] %>%
  dplyr::rename(T2B = days.2.bud)%>%
  dplyr::rename(Seeds = seed_num_nb)
row.names(main.ms.slopes.fig)=all.slopes$sp

# each slope separate

T2B = data.frame(main.ms.slopes.fig$T2B)
row.names(T2B) = all.slopes$sp
colnames(T2B) = "Time to First Bud"

Seeds = data.frame(main.ms.slopes.fig$Seeds)
row.names(Seeds) = all.slopes$sp
colnames(Seeds) = "Number of Seeds"

p = gheatmap(tree.2, T2B, width = 0.5,legend_title = "Time to First Bud",
             colnames_position = "top")+
  scale_fill_viridis_c(option="A", name="Time to First Bud")

p2 = p + new_scale_fill()

p3 = gheatmap(p2, Seeds, legend_title = "Number of Seeds", offset = 0.0008, width = 0.5,
              colnames_position = "top")+
  scale_fill_viridis_c(option="A", name="Number of Seeds")

p4 = p3 + theme(text = element_text (size = 22), legend.title = element_text(size = 14),
                legend.text = element_text(size = 14)) 

p5 = p4 %<+% all.slopes + geom_tippoint(aes(color = PC1), size = 10)
p5

#ggsave("./Germination.Fitness/Results/color.T2B.Seed.slopes.ontree.pc1.pdf", height = 10, width = 12)

phylo_plot.2 = edit_colors(p5, desaturate) # https://ggplot2-book.org/scales-colour
phylo_plot.3 = ggdraw(phylo_plot.2)
phylo_plot.3

# figure with colored pc1 was manually merged with figure of gray saturated heat map

#ggsave("./Germination.Fitness/Results/color.gray.T2B.Seed.slopes.ontree.pc1.pdf", height = 10, width = 12)

#### Figure 5 REVISED #####

all.phylo <- ggtree::read.tree("./tree_pruned.new")

sp.list=c("Caulanthus_anceps","Caulanthus_coulteri","Caulanthus_inflatus",
          "Streptanthus_breweri","Streptanthus_diversifolius","Streptanthus_drepanoides",
          "Streptanthus_glandulosus","Streptanthus_insignis","Streptanthus_polygaloides",
          "Streptanthus_tortuosus")

# prune phylogeny to only includes our species
phylo=keep.tip(all.phylo, sp.list)

phylo$tip.label = c("Caulanthus inflatus","Caulanthus coulteri","Caulanthus anceps",
                    "Streptanthus glandulosus","Streptanthus insignis","Streptanthus polygaloides",
                    "Streptanthus diversifolius","Streptanthus tortuosus","Streptanthus breweri",
                    "Streptanthus drepanoides")

ggplot(phylo) + geom_tree() + theme_tree()

tree = ggtree(phylo, size =1)
tree

# read in slopes
all.slopes=read.csv("./all.slopes.csv")

all.slopes$sp = c("Caulanthus inflatus","Caulanthus coulteri","Caulanthus anceps",
                    "Streptanthus glandulosus","Streptanthus insignis","Streptanthus polygaloides",
                    "Streptanthus diversifolius","Streptanthus tortuosus","Streptanthus breweri",
                    "Streptanthus drepanoides")

# climate pc1
all.slopes$PC1 = c(-3.4434154,-1.6218771,-1.6381408,0.2230000,-0.9859651,1.3163704,0.1026840,4.1567596,0.7433119,1.1472724)

# add tip labels
lb = phylo$tip.label
tips = data.frame(label = lb, Species = paste(all.slopes$sp))

tree.2 = tree %<+% tips +
  geom_tiplab(aes(label = Species),as_ylab = TRUE)
tree.2

# get slopes we are interested in 

main.ms.slopes.fig = all.slopes[,c(33,41)] %>%
  dplyr::rename(T2B = days.2.bud.final)%>%
  dplyr::rename(Seeds = seed_num_nb.final)
row.names(main.ms.slopes.fig)=all.slopes$sp

# each slope separate

T2B = data.frame(main.ms.slopes.fig$T2B)
row.names(T2B) = all.slopes$sp
colnames(T2B) = "Time to First Bud"

Seeds = data.frame(main.ms.slopes.fig$Seeds)
row.names(Seeds) = all.slopes$sp
colnames(Seeds) = "Number of Seeds"

p = gheatmap(tree.2, T2B, width = 0.5,legend_title = "Time to First Bud",
             colnames_position = "top")+
  scale_fill_viridis_c(option="A", name="Time to First Bud")

p2 = p + new_scale_fill()

p3 = gheatmap(p2, Seeds, legend_title = "Number of Seeds", offset = 0.0008, width = 0.5,
              colnames_position = "top")+
  scale_fill_viridis_c(option="A", name="Number of Seeds")

p4 = p3 + theme(text = element_text (size = 22), legend.title = element_text(size = 14),
                legend.text = element_text(size = 14)) 

p5 = p4 %<+% all.slopes + geom_tippoint(aes(color = PC1), size = 10)
p5

# ggsave("./Germination.Fitness/Results/color.T2B.Seed.slopes.ontree.pc1.final.pdf", height = 10, width = 12)

phylo_plot.2 = edit_colors(p5, desaturate) # https://ggplot2-book.org/scales-colour
phylo_plot.3 = ggdraw(phylo_plot.2)
phylo_plot.3

# figure with colored pc1 was manually merged with figure of gray saturated heat map

# ggsave("./Germination.Fitness/Results/color.gray.T2B.Seed.slopes.ontree.pc1.final.pdf", height = 10, width = 12)

#### Figure S8 ####

all.phylo <- ggtree::read.tree("./tree_pruned.new")

sp.list=c("Caulanthus_anceps","Caulanthus_coulteri","Caulanthus_inflatus",
          "Streptanthus_breweri","Streptanthus_diversifolius","Streptanthus_drepanoides",
          "Streptanthus_glandulosus","Streptanthus_insignis","Streptanthus_polygaloides",
          "Streptanthus_tortuosus")

# prune phylogeny to only includes our species
phylo=keep.tip(all.phylo, sp.list)

phylo$tip.label = c("Caulanthus inflatus","Caulanthus coulteri","Caulanthus anceps",
                    "Streptanthus glandulosus","Streptanthus insignis","Streptanthus polygaloides",
                    "Streptanthus diversifolius","Streptanthus tortuosus","Streptanthus breweri",
                    "Streptanthus drepanoides")

ggplot(phylo) + geom_tree() + theme_tree()

tree = ggtree(phylo, size =1)
tree

# read in slopes
all.slopes=read.csv("Germination.Fitness/Results/all.slopes.csv")

pheno.slopes$sp = c("Caulanthus inflatus","Caulanthus coulteri","Caulanthus anceps",
                    "Streptanthus glandulosus","Streptanthus insignis","Streptanthus polygaloides",
                    "Streptanthus diversifolius","Streptanthus tortuosus","Streptanthus breweri",
                    "Streptanthus drepanoides")

# add tip labels
lb = phylo$tip.label
tips = data.frame(label = lb, Species = paste(all.slopes$sp))

tree.2 = tree %<+% tips +
  geom_tiplab(aes(label = Species),as_ylab = TRUE)
tree.2

# get slopes we are interested in 

supp.ms.slopes.fig = all.slopes[,c(6,9,13,12)] %>%
  dplyr::rename(FBD = sept.1.bud)%>%
  dplyr::rename(Pflwr = pflwr)%>%
  dplyr::rename(SeedMass = seed_mass)%>%
  dplyr::rename(yr1fit = year1fit)

row.names(supp.ms.slopes.fig)=all.slopes$sp

# each slope separate

FBD = data.frame(supp.ms.slopes.fig$FBD)
row.names(FBD) = all.slopes$sp
colnames(FBD) = "First Bud Date"

Pflwr = data.frame(supp.ms.slopes.fig$Pflwr)
row.names(Pflwr) = all.slopes$sp
colnames(Pflwr) = "Prob. of Flowering"

SeedMass = data.frame(supp.ms.slopes.fig$SeedMass)
row.names(SeedMass) = all.slopes$sp
colnames(SeedMass) = "Total Seed Mass"

Year1fit = data.frame(supp.ms.slopes.fig$yr1fit)
row.names(Year1fit) = all.slopes$sp
colnames(Year1fit) = "Year 1 Fitness"

p.supp = gheatmap(tree.2, FBD, width = 0.5,legend_title = "First Bud Date",
                  colnames_position = "top")+
  scale_fill_viridis_c(option="A", name="First Bud Date")

p2.supp = p.supp + new_scale_fill()

p3.supp = gheatmap(p2.supp, Pflwr, legend_title = "Prob. of Flowering", offset = 0.0008, width = 0.5,
                   colnames_position = "bottom")+
  scale_fill_viridis_c(option="A", name="Prob. of Flowering")

p4.supp = p3.supp + new_scale_fill()

p5.supp = gheatmap(p4.supp, SeedMass, legend_title = "Total Seed Mass", offset = 0.0016, width = 0.5,
                   colnames_position = "top")+
  scale_fill_viridis_c(option="A", name="Total Seed Mass")

p6.supp = p5.supp + new_scale_fill()

p7.supp = gheatmap(p6.supp, Year1fit, legend_title = "Year 1 Fitness", offset = 0.0024, width = 0.5,
                   colnames_position = "bottom")+
  scale_fill_viridis_c(option="A", name="Year 1 Fitness")

p8.supp = p7.supp + theme(text = element_text (size = 22), legend.title = element_text(size = 12),
                          legend.text = element_text(size = 12)) 
p8.supp

phylo_plot.2.supp = edit_colors(p8.supp, desaturate) # https://ggplot2-book.org/scales-colour
phylo_plot.3.supp = ggdraw(phylo_plot.2.supp)
phylo_plot.3.supp

# ggsave("./Germination.Fitness/Results/color.gray.supp.slopes.ontree.pdf", height = 10, width = 12)

#### Figure S8 REVISED ####

all.phylo <- ggtree::read.tree("./tree_pruned.new")

sp.list=c("Caulanthus_anceps","Caulanthus_coulteri","Caulanthus_inflatus",
          "Streptanthus_breweri","Streptanthus_diversifolius","Streptanthus_drepanoides",
          "Streptanthus_glandulosus","Streptanthus_insignis","Streptanthus_polygaloides",
          "Streptanthus_tortuosus")

# prune phylogeny to only includes our species
phylo=keep.tip(all.phylo, sp.list)

phylo$tip.label = c("Caulanthus inflatus","Caulanthus coulteri","Caulanthus anceps",
                    "Streptanthus glandulosus","Streptanthus insignis","Streptanthus polygaloides",
                    "Streptanthus diversifolius","Streptanthus tortuosus","Streptanthus breweri",
                    "Streptanthus drepanoides")

ggplot(phylo) + geom_tree() + theme_tree()

tree = ggtree(phylo, size =1)
tree

# read in slopes
all.slopes=read.csv("./all.slopes.csv")

all.slopes$sp = c("Caulanthus inflatus","Caulanthus coulteri","Caulanthus anceps",
                    "Streptanthus glandulosus","Streptanthus insignis","Streptanthus polygaloides",
                    "Streptanthus diversifolius","Streptanthus tortuosus","Streptanthus breweri",
                    "Streptanthus drepanoides")

# add tip labels
lb = phylo$tip.label
tips = data.frame(label = lb, Species = paste(all.slopes$sp))

tree.2 = tree %<+% tips +
  geom_tiplab(aes(label = Species),as_ylab = TRUE)
tree.2

# get slopes we are interested in 

supp.ms.slopes.fig = all.slopes[,c(35,37,45,43)] %>%
  dplyr::rename(FBD = sept.1.bud.final)%>%
  dplyr::rename(Pflwr = pflwr.final)%>%
  dplyr::rename(SeedMass = seed_mass.final)%>%
  dplyr::rename(yr1fit = year1fit.final)

row.names(supp.ms.slopes.fig)=all.slopes$sp

# each slope separate

FBD = data.frame(supp.ms.slopes.fig$FBD)
row.names(FBD) = all.slopes$sp
colnames(FBD) = "First Bud Date"

Pflwr = data.frame(supp.ms.slopes.fig$Pflwr)
row.names(Pflwr) = all.slopes$sp
colnames(Pflwr) = "Prob. of Flowering"

SeedMass = data.frame(supp.ms.slopes.fig$SeedMass)
row.names(SeedMass) = all.slopes$sp
colnames(SeedMass) = "Total Seed Mass"

Year1fit = data.frame(supp.ms.slopes.fig$yr1fit)
row.names(Year1fit) = all.slopes$sp
colnames(Year1fit) = "Year 1 Fitness"

p.supp = gheatmap(tree.2, FBD, width = 0.5,legend_title = "First Bud Date",
                  colnames_position = "top")+
  scale_fill_viridis_c(option="A", name="First Bud Date")

p2.supp = p.supp + new_scale_fill()

p3.supp = gheatmap(p2.supp, Pflwr, legend_title = "Prob. of Flowering", offset = 0.0008, width = 0.5,
                   colnames_position = "bottom")+
  scale_fill_viridis_c(option="A", name="Prob. of Flowering")

p4.supp = p3.supp + new_scale_fill()

p5.supp = gheatmap(p4.supp, SeedMass, legend_title = "Total Seed Mass", offset = 0.0016, width = 0.5,
                   colnames_position = "top")+
  scale_fill_viridis_c(option="A", name="Total Seed Mass")

p6.supp = p5.supp + new_scale_fill()

p7.supp = gheatmap(p6.supp, Year1fit, legend_title = "Year 1 Fitness", offset = 0.0024, width = 0.5,
                   colnames_position = "bottom")+
  scale_fill_viridis_c(option="A", name="Year 1 Fitness")

p8.supp = p7.supp + theme(text = element_text (size = 22), legend.title = element_text(size = 12),
                          legend.text = element_text(size = 12)) 
p8.supp

phylo_plot.2.supp = edit_colors(p8.supp, desaturate) # https://ggplot2-book.org/scales-colour
phylo_plot.3.supp = ggdraw(phylo_plot.2.supp)
phylo_plot.3.supp

ggsave("./Germination.Fitness/Results/color.gray.supp.slopes.ontree.final.pdf", height = 10, width = 12)



#### PGLS to relate phenology and fitness slopes ####
# read in tree
all.phylo <- read.tree("./tree_pruned.new")

sp.list=c("Caulanthus_anceps","Caulanthus_coulteri","Caulanthus_inflatus",
          "Streptanthus_breweri","Streptanthus_diversifolius","Streptanthus_drepanoides",
          "Streptanthus_glandulosus","Streptanthus_insignis","Streptanthus_polygaloides",
          "Streptanthus_tortuosus")

# prune phylogeny to only includes our species
phylo=keep.tip(all.phylo, sp.list)

all.slopes=read.csv("./all.slopes.csv")

# pgls using caper package

comp.data.1<-comparative.data(phylo, all.slopes, names.col="sp", vcv.dim=2, warn.dropped=TRUE, vcv = TRUE)

# d2bud and seed_num_nb

d2bud.seed.num.pgls = pgls(days.2.bud~seed_num_nb, data=comp.data.1)
summary(d2bud.seed.num.pgls) # not-significant, p = 0.07

# revised slopes

d2bud.seed.num.pgls.final = pgls(days.2.bud.final~seed_num_nb.final, data=comp.data.1)
summary(d2bud.seed.num.pgls.final) # not-significant, p = 0.07

# population sensitivity

d2bud.seed.num.pgls.pop = pgls(days.2.bud.CAAN1.CAIN3~seed_num_nb.CAAN1.CAIN3, data=comp.data.1)
summary(d2bud.seed.num.pgls.pop) # not-significant, p = 0.06

d2bud.seed.num.pgls.pop = pgls(days.2.bud.CAAN1.CAIN4~seed_num_nb.CAAN1.CAIN4, data=comp.data.1)
summary(d2bud.seed.num.pgls.pop) # not-significant, p = 0.09

d2bud.seed.num.pgls.pop = pgls(days.2.bud.CAAN2.CAIN3~seed_num_nb.CAAN2.CAIN3, data=comp.data.1)
summary(d2bud.seed.num.pgls.pop) # not-significant, p = 0.07

d2bud.seed.num.pgls.pop = pgls(days.2.bud.CAAN2.CAIN4~seed_num_nb.CAAN2.CAIN4, data=comp.data.1)
summary(d2bud.seed.num.pgls.pop) # not-significant, p = 0.10

# d2bud and pflwr

d2bud.pflwr.pgls = pgls(days.2.bud~pflwr, data=comp.data.1)
summary(d2bud.pflwr.pgls) # not significant p = 0.78

d2bud.pflwr.weighted.pgls = pgls(days.2.bud~pflwr.weighted, data=comp.data.1)
summary(d2bud.pflwr.weighted.pgls) # not significant p = 0.75

# revised slopes

d2bud.pflwr.pgls.final = pgls(days.2.bud.final~pflwr.final, data=comp.data.1)
summary(d2bud.pflwr.pgls.final) # not significant p = 0.78

# population sensitivity

d2bud.pflwr.pgls.pop = pgls(days.2.bud.CAAN1.CAIN3~pflwr.CAAN1.CAIN3, data=comp.data.1)
summary(d2bud.pflwr.pgls.pop) # not-significant, p = 0.74

d2bud.pflwr.pgls.pop = pgls(days.2.bud.CAAN1.CAIN4~pflwr.CAAN1.CAIN4, data=comp.data.1)
summary(d2bud.pflwr.pgls.pop) # not-significant, p = 0.77

d2bud.pflwr.pgls.pop = pgls(days.2.bud.CAAN2.CAIN3~pflwr.CAAN2.CAIN3, data=comp.data.1)
summary(d2bud.pflwr.pgls.pop) # not-significant, p = 0.72

d2bud.pflwr.pgls.pop = pgls(days.2.bud.CAAN2.CAIN4~pflwr.CAAN2.CAIN4, data=comp.data.1)
summary(d2bud.pflwr.pgls.pop) # not-significant, p = 0.75

# d2bud and yr1fit 

d2bud.yr1fit.pgls = pgls(days.2.bud~year1fit, data=comp.data.1)
summary(d2bud.yr1fit.pgls) # not significant p = 0.49

# revised slopes

d2bud.yr1fit.pgls.final = pgls(days.2.bud.final~year1fit.final, data=comp.data.1)
summary(d2bud.yr1fit.pgls.final) # not significant p = 0.50

# population sensitivity

d2bud.yr1fit.pgls.pop = pgls(days.2.bud.CAAN1.CAIN3~year1fit.CAAN1.CAIN3, data=comp.data.1)
summary(d2bud.yr1fit.pgls.pop) # not-significant, p = 0.41

d2bud.yr1fit.pgls.pop = pgls(days.2.bud.CAAN1.CAIN4~year1fit.CAAN1.CAIN4, data=comp.data.1)
summary(d2bud.yr1fit.pgls.pop) # not-significant, p = 0.58

d2bud.yr1fit.pgls.pop = pgls(days.2.bud.CAAN2.CAIN3~year1fit.CAAN2.CAIN3, data=comp.data.1)
summary(d2bud.yr1fit.pgls.pop) # not-significant, p = 0.44

d2bud.yr1fit.pgls.pop = pgls(days.2.bud.CAAN2.CAIN4~year1fit.CAAN2.CAIN4, data=comp.data.1)
summary(d2bud.yr1fit.pgls.pop) # not-significant, p = 0.61

#### Figure S10 ####

d2bud.seed.num.pgls.plot= ggplot(all.slopes, aes(y = days.2.bud, x = seed_num_nb))+
  geom_point(size=2)+
  theme_classic(base_size=22)+
  labs(y = "Time to First Bud (days)",x = "Number of Seeds")
d2bud.seed.num.pgls.plot

#ggsave("Germination.Fitness/Results/days2bud.seed.num.pgls.pdf", height = 7, width = 7)

d2bud.pflwr.pgls.plot= ggplot(all.slopes, aes(y = pflwr, x = year1fit))+
  geom_point(size=2)+
  theme_classic(base_size=22)+
  labs(y = "Time to First Bud (days)",x = "Probability of Flowering")
d2bud.pflwr.pgls.plot

#ggsave("Germination.Fitness/Results/days2bud.pflwr.pgls.pdf", height = 7, width = 7)

d2bud.yr1fit.pgls.plot= ggplot(all.slopes, aes(y = days.2.bud, x = year1fit))+
  geom_point(size=2)+
  theme_classic(base_size=22)+
  labs(y = "Time to First Bud (days)",x = "Year 1 fitness")
d2bud.yr1fit.pgls.plot
#ggsave("Germination.Fitness/Results/days2bud.yr1fit.pgls.pdf", height = 7, width = 7)

#### Figure S10 REVISED ####

d2bud.seed.num.pgls.plot= ggplot(all.slopes, aes(y = days.2.bud.final, x = seed_num_nb.final))+
  geom_point(size=4)+
  theme_classic(base_size=22)+
  labs(y = "Time to First Bud (days)",x = "Number of Seeds")
d2bud.seed.num.pgls.plot

#ggsave("Germination.Fitness/Results/days2bud.seed.num.pgls.final.pdf", height = 7, width = 7)

d2bud.pflwr.pgls.plot= ggplot(all.slopes, aes(y = pflwr.final, x = year1fit.final))+
  geom_point(size=4)+
  theme_classic(base_size=22)+
  labs(y = "Time to First Bud (days)",x = "Probability of Flowering")
d2bud.pflwr.pgls.plot

#ggsave("Germination.Fitness/Results/days2bud.pflwr.pgls.final.pdf", height = 7, width = 7)

d2bud.yr1fit.pgls.plot= ggplot(all.slopes, aes(y = days.2.bud.final, x = year1fit.final))+
  geom_point(size=4)+
  theme_classic(base_size=22)+
  labs(y = "Time to First Bud (days)",x = "Year 1 fitness")
d2bud.yr1fit.pgls.plot

#ggsave("Germination.Fitness/Results/days2bud.yr1fit.pgls.final.pdf", height = 7, width = 7)

#### PGLS to relate germination slopes to fitness slopes ####
# Germination timing slopes come from Worthy et al. 2024 Ecology

all.slopes=read.csv("./all.slopes.csv")

# add slopes from germ pheno: germination_fraction~cohort
# had to average two CAIN and two CAAN slopes
all.slopes$pheno.slopes = c(-0.45898,-0.1238,-0.021865,-0.59355,-0.33908,-0.25232,-0.40204,-0.2401,
                            -0.36986,-0.36455)
all.slopes$pheno.slopes.CAAN1.CAIN3 = c(-0.4941,-0.1238,-0.07615,-0.59355,-0.33908,-0.25232,-0.40204,-0.2401,
                            -0.36986,-0.36455)
all.slopes$pheno.slopes.CAAN1.CAIN4 = c(-0.42386,-0.1238,-0.07615,-0.59355,-0.33908,-0.25232,-0.40204,-0.2401,
                                        -0.36986,-0.36455)
all.slopes$pheno.slopes.CAAN2.CAIN3 = c(-0.4941,-0.1238,0.032425,-0.59355,-0.33908,-0.25232,-0.40204,-0.2401,
                                        -0.36986,-0.36455)
all.slopes$pheno.slopes.CAAN2.CAIN4 = c(-0.42386,-0.1238,0.03242,-0.59355,-0.33908,-0.25232,-0.40204,-0.2401,
                                        -0.36986,-0.36455)

all.phylo <- read.tree("./tree_pruned.new")

sp.list=c("Caulanthus_anceps","Caulanthus_coulteri","Caulanthus_inflatus",
          "Streptanthus_breweri","Streptanthus_diversifolius","Streptanthus_drepanoides",
          "Streptanthus_glandulosus","Streptanthus_insignis","Streptanthus_polygaloides",
          "Streptanthus_tortuosus")

# prune phylogeny to only includes our species
phylo=keep.tip(all.phylo, sp.list)

comp.data<-comparative.data(phylo, all.slopes, names.col="sp", vcv.dim=2, warn.dropped=TRUE, vcv = TRUE)

# seed_num_nb and pheno slopes
pheno.slopes.seed.num.pgls = pgls(seed_num_nb~pheno.slopes, data=comp.data)
summary(pheno.slopes.seed.num.pgls) # not-significant, p = 0.50

# revised slopes
pheno.slopes.seed.num.pgls.final = pgls(seed_num_nb.final~pheno.slopes, data=comp.data)
summary(pheno.slopes.seed.num.pgls.final) # not-significant, p = 0.51

# population sensitivity

pheno.slopes.seed.num.pgls.pop = pgls(seed_num_nb.CAAN1.CAIN3~pheno.slopes.CAAN1.CAIN3, data=comp.data)
summary(pheno.slopes.seed.num.pgls.pop) # not-significant, p = 0.31

pheno.slopes.seed.num.pgls.pop = pgls(seed_num_nb.CAAN1.CAIN4~pheno.slopes.CAAN1.CAIN4, data=comp.data)
summary(pheno.slopes.seed.num.pgls.pop) # not-significant, p = 0.64

pheno.slopes.seed.num.pgls.pop = pgls(seed_num_nb.CAAN2.CAIN3~pheno.slopes.CAAN2.CAIN3, data=comp.data)
summary(pheno.slopes.seed.num.pgls.pop) # not-significant, p = 0.41

pheno.slopes.seed.num.pgls.pop = pgls(seed_num_nb.CAAN2.CAIN4~pheno.slopes.CAAN2.CAIN4, data=comp.data)
summary(pheno.slopes.seed.num.pgls.pop) # not-significant, p = 0.77

# year1fit and pheno slopes

pheno.slopes.year1fit.pgls = pgls(year1fit~pheno.slopes, data=comp.data)
summary(pheno.slopes.year1fit.pgls) # not-significant, p = 0.71

# revised slopes

pheno.slopes.year1fit.pgls.final = pgls(year1fit.final~pheno.slopes, data=comp.data)
summary(pheno.slopes.year1fit.pgls.final) # not-significant, p = 0.70

# population sensitivity

pheno.slopes.year1fit.pgls.pop = pgls(year1fit.CAAN1.CAIN3~pheno.slopes.CAAN1.CAIN3, data=comp.data)
summary(pheno.slopes.year1fit.pgls.pop) # not-significant, p = 0.97

pheno.slopes.year1fit.pgls.pop = pgls(year1fit.CAAN1.CAIN4~pheno.slopes.CAAN1.CAIN4, data=comp.data)
summary(pheno.slopes.year1fit.pgls.pop) # not-significant, p = 0.57

pheno.slopes.year1fit.pgls.pop = pgls(year1fit.CAAN2.CAIN3~pheno.slopes.CAAN2.CAIN3, data=comp.data)
summary(pheno.slopes.year1fit.pgls.pop) # not-significant, p = 0.93

pheno.slopes.year1fit.pgls.pop = pgls(year1fit.CAAN2.CAIN4~pheno.slopes.CAAN2.CAIN4, data=comp.data)
summary(pheno.slopes.year1fit.pgls.pop) # not-significant, p = 0.50

#### Figure S9 ####

# no line b/c not significant
pheno.slopes.seed.num.pgls.plot= ggplot(all.slopes, aes(y = seed_num_nb, x = pheno.slopes))+
  geom_point(size=2)+
  theme_classic(base_size=22)+
  labs(y = "Number of Seeds", x = "Germination Specialization")
pheno.slopes.seed.num.pgls.plot

#ggsave("Germination.Fitness/Results/seed.num.pheno.slopes.pgls.pdf", height = 7, width = 7)

pheno.slopes.year1fit.pgls.plot= ggplot(all.slopes, aes(y = year1fit, x = pheno.slopes))+
  geom_point(size = 2)+
  theme_classic(base_size=22)+
  labs(y = "Year 1 fitness\n (p(flower)*seed number)", x = "Germination Specialization")
pheno.slopes.year1fit.pgls.plot

#ggsave("Germination.Fitness/Results/yr1fit.pheno.slopes.pgls.pdf", height = 7, width = 7)

#### Figure S9 REVISED ####

# no line b/c not significant
pheno.slopes.seed.num.pgls.plot= ggplot(all.slopes, aes(y = seed_num_nb.final, x = pheno.slopes))+
  geom_point(size=4)+
  theme_classic(base_size=22)+
  labs(y = "Number of Seeds", x = "Germination Specialization")
pheno.slopes.seed.num.pgls.plot

#ggsave("Germination.Fitness/Results/seed.num.pheno.slopes.pgls.final.pdf", height = 7, width = 7)

pheno.slopes.year1fit.pgls.plot= ggplot(all.slopes, aes(y = year1fit.final, x = pheno.slopes))+
  geom_point(size = 4)+
  theme_classic(base_size=22)+
  labs(y = "Year 1 fitness\n (p(flower)*seed number)", x = "Germination Specialization")
pheno.slopes.year1fit.pgls.plot

#ggsave("Germination.Fitness/Results/yr1fit.pheno.slopes.pgls.final.pdf", height = 7, width = 7)

#### PGLS to relate phenology and fitness slopes to climate ####

all.phylo <- read.tree("./tree_pruned.new")

sp.list=c("Caulanthus_anceps","Caulanthus_coulteri","Caulanthus_inflatus",
          "Streptanthus_breweri","Streptanthus_diversifolius","Streptanthus_drepanoides",
          "Streptanthus_glandulosus","Streptanthus_insignis","Streptanthus_polygaloides",
          "Streptanthus_tortuosus")

# prune phylogeny to only includes our species
phylo=keep.tip(all.phylo, sp.list)

plot(phylo, no.margin = TRUE, font = 3, cex = .75)

# merge climate data with slopes data
# read in slopes  dataframe

all.slopes=read.csv("./all.slopes.csv")

# Get PC1

locs = read.csv("./georeferencing_clean.csv")
climate = read.csv("./all_herbarium_climate.csv") %>%
  mutate(type = "herbarium")

climsummaries = climate %>% 
  filter(clim_year > 1990) %>% 
  filter(clim_year < 2016) %>% 
  group_by(specimen,type) %>% 
  dplyr::summarize(Cwd = mean(cwd), PPT_CV = (sd(ppt_mm)/mean(ppt_mm)),PPT = mean(ppt_mm), 
                   Tmin_SD=sd(tmin), Tmax_SD=sd(tmax),Tmin = mean(tmin), Tmax = mean(tmax)) %>%
  left_join(., locs) %>%
  dplyr::rename(id = specimen)
# data is summarized over years 1991-2016

climsummaries.2 = climsummaries %>%
  group_by(folder) %>%
  dplyr::summarize(Cwd = mean(Cwd), PPT_CV = mean(PPT_CV),PPT = mean(PPT), 
                   Tmin_SD=sd(Tmin_SD), Tmax_SD=sd(Tmax_SD),Tmin = mean(Tmin), Tmax = mean(Tmax))

# subset herbarium data by each species
climsummaries.3 = subset(climsummaries.2, climsummaries.2$folder %in% c("c_anceps","c_coulteri","c_inflatus",
                                                                        "s_breweri","s_diversifolius",
                                                                        "s_drepanoides","s_glandulosus",
                                                                        "s_insignis","s_polygaloides",
                                                                        "s_tortuosus"))

# PC of herbarium data
all.data.4pc.year = climsummaries.3 %>%
  ungroup() %>%
  dplyr:: select(Cwd,PPT,PPT_CV,Tmin,Tmax,Tmax_SD,Tmin_SD)

all.data.pc.year = prcomp(all.data.4pc.year, scale  = TRUE, center = TRUE)
all.data.pc.year
biplot(all.data.pc.year)
all.data.year.pc.dat = data.frame(all.data.pc.year$x)
all.data.year.2 = cbind(climsummaries.3, all.data.year.pc.dat)

all.data.year.2$folder = c("CAAN","CACO1","CAIN","STBR3","STDI","STDR2","STGL1","STIN","STPO","STTO_TM2")
colnames(all.data.year.2)[1] = "Pop"
all.data.year.final = left_join(all.slopes,all.data.year.2, by = "Pop")

comp.data<-comparative.data(phylo, all.data.year.final, names.col="sp", vcv.dim=2, warn.dropped=TRUE)

# days.2.bud~PC1
days2bud.pc1.pgls = pgls(days.2.bud~PC1, data=comp.data)
summary(days2bud.pc1.pgls) # significant R2 = 0.63
# PC1 is associated with all variables except temp_SD

# revised slopes

days2bud.pc1.pgls.final = pgls(days.2.bud.final~PC1, data=comp.data)
summary(days2bud.pc1.pgls.final) # significant R2 = 0.67
# PC1 is associated with all variables except temp_SD

# population sensitivity

days2bud.pc1.pgls.pop = pgls(days.2.bud.CAAN1.CAIN3~PC1, data=comp.data)
summary(days2bud.pc1.pgls.pop) # significant R2 = 0.65

days2bud.pc1.pgls.pop = pgls(days.2.bud.CAAN1.CAIN4~PC1, data=comp.data)
summary(days2bud.pc1.pgls.pop) # significant R2 = 0.66

days2bud.pc1.pgls.pop = pgls(days.2.bud.CAAN2.CAIN3~PC1, data=comp.data)
summary(days2bud.pc1.pgls.pop) # significant R2 = 0.67

days2bud.pc1.pgls.pop = pgls(days.2.bud.CAAN2.CAIN4~PC1, data=comp.data)
summary(days2bud.pc1.pgls.pop) # significant R2 = 0.68

# pflwr~PC1

flwprob.pc1.pgls = pgls(pflwr~PC1, data=comp.data)
summary(flwprob.pc1.pgls) # not significant
# PC1 is associated with all variables except temp_SD

# revised slope

flwprob.pc1.pgls.final = pgls(pflwr.final~PC1, data=comp.data)
summary(flwprob.pc1.pgls.final) # not significant
# PC1 is associated with all variables except temp_SD

flwprob.pc1.pgls.weighted = pgls(pflwr.weighted~PC1, data=comp.data)
summary(flwprob.pc1.pgls.weighted) # not significant

# population sensitivity

flwprob.pc1.pgls.pop = pgls(pflwr.CAAN1.CAIN3~PC1, data=comp.data)
summary(flwprob.pc1.pgls.pop) # not significant

flwprob.pc1.pgls.pop = pgls(pflwr.CAAN1.CAIN4~PC1, data=comp.data)
summary(flwprob.pc1.pgls.pop) # not significant

flwprob.pc1.pgls.pop = pgls(pflwr.CAAN2.CAIN3~PC1, data=comp.data)
summary(flwprob.pc1.pgls.pop) # not significant

flwprob.pc1.pgls.pop = pgls(pflwr.CAAN2.CAIN4~PC1, data=comp.data)
summary(flwprob.pc1.pgls.pop) # not significant

# seed_num~PC1

seed.num.pc1.pgls = pgls(seed_num_nb~PC1, data=comp.data)
summary(seed.num.pc1.pgls) # significant p = 0.02
# PC1 is associated with all variables except temp_SD

# revised slopes

seed.num.pc1.pgls.final = pgls(seed_num_nb.final~PC1, data=comp.data)
summary(seed.num.pc1.pgls.final) # significant p = 0.02
# PC1 is associated with all variables except temp_SD

# population sensitivity

seed.num.pc1.pgls.pop = pgls(seed_num_nb.CAAN1.CAIN3~PC1, data=comp.data)
summary(seed.num.pc1.pgls.pop) # significant R2 = 0.61

seed.num.pc1.pgls.pop = pgls(seed_num_nb.CAAN1.CAIN4~PC1, data=comp.data)
summary(seed.num.pc1.pgls.pop) # not-significant, p = 0.055 R2 = 0.39

seed.num.pc1.pgls.pop = pgls(seed_num_nb.CAAN2.CAIN3~PC1, data=comp.data)
summary(seed.num.pc1.pgls.pop) # significant R2 = 0.63

seed.num.pc1.pgls.pop = pgls(seed_num_nb.CAAN2.CAIN4~PC1, data=comp.data)
summary(seed.num.pc1.pgls.pop) # not-significant, p = 0.05 R2 = 0.40

# year1fit~PC1
fit.pc1.pgls = pgls(year1fit~PC1, data=comp.data)
summary(fit.pc1.pgls) # not significant

# revised slopes
fit.pc1.pgls.final = pgls(year1fit.final~PC1, data=comp.data)
summary(fit.pc1.pgls.final) # not significant

# population sensitivity

fit.pc1.pgls.pop = pgls(year1fit.CAAN1.CAIN3~PC1, data=comp.data)
summary(fit.pc1.pgls.pop) # not significant

fit.pc1.pgls.pop = pgls(year1fit.CAAN1.CAIN4~PC1, data=comp.data)
summary(fit.pc1.pgls.pop) # not-significant

fit.pc1.pgls.pop = pgls(year1fit.CAAN2.CAIN3~PC1, data=comp.data)
summary(fit.pc1.pgls.pop) # not-significant

fit.pc1.pgls.pop = pgls(year1fit.CAAN2.CAIN4~PC1, data=comp.data)
summary(fit.pc1.pgls.pop) # not-significant

#### Figure 6 ####

all.data.year.final$Pop.2 = c("CAIN","CACO","CAAN","STGL","STIN","STPO","STDI",
                              "STTO","STBR","STDR")
all.data.year.final$pred.days2bud.pc1 = predict(days2bud.pc1.pgls)

#define annotations to axis labels for PCA axes
text_x_low = text_grob("Hot & dry", size = 12,col="gray35")
text_x_high = text_grob("Cool & wet", size = 12,col="gray35")

days2bud.pc1.pgls.plot= ggplot(all.data.year.final, aes(y = days.2.bud, x = PC1, label = Pop.2))+
  geom_line(aes(y = pred.days2bud.pc1))+
  geom_point(size = 4)+
  geom_label_repel(size = 5)+
  theme_classic(base_size=22)+
  labs(y = "Time to First Bud (days)")+
  annotation_custom(text_x_high,xmin=3.8,xmax= 4,ymin =-1.07,ymax = -1.07) +
  annotation_custom(text_x_low,xmin=-2.8,xmax=-3,ymin =-1.07,ymax = -1.07)+
  coord_cartesian( clip="off") #this keeps it from clipping off the stuff outside the plot
days2bud.pc1.pgls.plot

#ggsave("Germination.Fitness/Results/days2bud.pc1.pgls.pdf", height = 7, width = 7)

all.data.year.final$pred.flwprob.pc1 = predict(flwprob.pc1.pgls)

#define annotations to axis labels for PCA axes
text_x_low = text_grob("Hot & dry", size = 12,col="gray35")
text_x_high = text_grob("Cool & wet", size = 12,col="gray35")

# no line because not significant
pflwr.pc1.pgls.plot= ggplot(all.data.year.final, aes(y = pflwr, x = PC1, label = Pop.2))+
  geom_point(size = 4)+
  geom_label_repel(size = 5)+
  theme_classic(base_size=22)+
  labs(y = "Probability of Flowering")+
  annotation_custom(text_x_high,xmin=3.8,xmax= 4,ymin =-0.04,ymax = -0.04) +
  annotation_custom(text_x_low,xmin=-2.8,xmax=-3,ymin =-0.04,ymax = -0.04)+
  coord_cartesian( clip="off") #this keeps it from clipping off the stuff outside the plot
pflwr.pc1.pgls.plot

#ggsave("Germination.Fitness/Results/pflwr.pc1.pgls.pdf", height = 7, width = 7)

all.data.year.final$pred.seed.num.pc1 = predict(seed.num.pc1.pgls)

#define annotations to axis labels for PCA axes
text_x_low = text_grob("Hot & dry", size = 12,col="gray35")
text_x_high = text_grob("Cool & wet", size = 12,col="gray35")

# no line because not significant
seed.num.pc1.pgls.plot= ggplot(all.data.year.final, aes(y = seed_num_nb, x = PC1, label = Pop.2))+
  geom_line(aes(y = pred.seed.num.pc1))+
  geom_point(size=4)+
  geom_label_repel(size=5)+
  theme_classic(base_size=22)+
  labs(y = "Number of Seeds")+
  annotation_custom(text_x_high,xmin=3.8,xmax= 4,ymin =-0.033,ymax = -0.033) +
  annotation_custom(text_x_low,xmin=-2.8,xmax=-3,ymin =-0.033,ymax = -0.033)+
  coord_cartesian( clip="off") #this keeps it from clipping off the stuff outside the plot
seed.num.pc1.pgls.plot

#ggsave("Germination.Fitness/Results/seed.num.pc1.pgls.pdf", height = 7, width = 7)

all.data.year.final$pred.year1fit.pc1 = predict(fit.pc1.pgls)

#define annotations to axis labels for PCA axes
text_x_low = text_grob("Hot & dry", size = 12,col="gray35")
text_x_high = text_grob("Cool & wet", size = 12,col="gray35")

year1fit.pc1.pgls.plot= ggplot(all.data.year.final, aes(y = year1fit, x = PC1, label = Pop.2))+
  geom_point(size=4)+
  geom_label_repel(size=5)+
  theme_classic(base_size=22)+
  labs(y = "Year 1 fitness\n (p(flower)*seed number)")+
  annotation_custom(text_x_high,xmin=3.8,xmax= 4,ymin =-0.03,ymax = -0.03) +
  annotation_custom(text_x_low,xmin=-2.8,xmax=-3,ymin =-0.03,ymax = -0.03)+
  coord_cartesian( clip="off") #this keeps it from clipping off the stuff outside the plot
year1fit.pc1.pgls.plot

#ggsave("Germination.Fitness/Results/yr1fit.pc1.pgls.pdf", height = 7, width = 7)

Figure6 = plot_grid(days2bud.pc1.pgls.plot,pflwr.pc1.pgls.plot,seed.num.pc1.pgls.plot,
                    year1fit.pc1.pgls.plot, nrow = 2, ncol = 2, labels = c("A.","B.","C.","D."))



#### Figure 6 REVISED ####

all.data.year.final$Pop.2 = c("CAIN","CACO","CAAN","STGL","STIN","STPO","STDI",
                              "STTO","STBR","STDR")
all.data.year.final$pred.days2bud.pc1 = predict(days2bud.pc1.pgls.final)

#define annotations to axis labels for PCA axes
text_x_low = text_grob("Hot & dry", size = 12,col="gray35")
text_x_high = text_grob("Cool & wet", size = 12,col="gray35")

days2bud.pc1.pgls.plot= ggplot(all.data.year.final, aes(y = days.2.bud.final, x = PC1, label = Pop.2))+
  geom_line(aes(y = pred.days2bud.pc1))+
  geom_point(size = 4)+
  geom_label_repel(size = 5)+
  theme_classic(base_size=22)+
  labs(y = "Time to First Bud (days)")+
  annotation_custom(text_x_high,xmin=3.8,xmax= 4,ymin =-1.07,ymax = -1.07) +
  annotation_custom(text_x_low,xmin=-2.8,xmax=-3,ymin =-1.07,ymax = -1.07)+
  coord_cartesian( clip="off") #this keeps it from clipping off the stuff outside the plot
days2bud.pc1.pgls.plot

ggsave("Germination.Fitness/Results/days2bud.pc1.pgls.final.pdf", height = 7, width = 7)

all.data.year.final$pred.flwprob.pc1 = predict(flwprob.pc1.pgls.final)

#define annotations to axis labels for PCA axes
text_x_low = text_grob("Hot & dry", size = 12,col="gray35")
text_x_high = text_grob("Cool & wet", size = 12,col="gray35")

# no line because not significant
pflwr.pc1.pgls.plot= ggplot(all.data.year.final, aes(y = pflwr.final, x = PC1, label = Pop.2))+
  geom_point(size = 4)+
  geom_label_repel(size = 5)+
  theme_classic(base_size=22)+
  labs(y = "Probability of Flowering")+
  annotation_custom(text_x_high,xmin=3.8,xmax= 4,ymin =-0.04,ymax = -0.04) +
  annotation_custom(text_x_low,xmin=-2.8,xmax=-3,ymin =-0.04,ymax = -0.04)+
  coord_cartesian( clip="off") #this keeps it from clipping off the stuff outside the plot
pflwr.pc1.pgls.plot

#ggsave("Germination.Fitness/Results/pflwr.pc1.pgls.pdf", height = 7, width = 7)

all.data.year.final$pred.seed.num.pc1 = predict(seed.num.pc1.pgls.final)

#define annotations to axis labels for PCA axes
text_x_low = text_grob("Hot & dry", size = 12,col="gray35")
text_x_high = text_grob("Cool & wet", size = 12,col="gray35")

# no line because not significant
seed.num.pc1.pgls.plot= ggplot(all.data.year.final, aes(y = seed_num_nb.final, x = PC1, label = Pop.2))+
  geom_line(aes(y = pred.seed.num.pc1))+
  geom_point(size=4)+
  geom_label_repel(size=5)+
  theme_classic(base_size=22)+
  labs(y = "Number of Seeds")+
  annotation_custom(text_x_high,xmin=3.8,xmax= 4,ymin =-0.033,ymax = -0.033) +
  annotation_custom(text_x_low,xmin=-2.8,xmax=-3,ymin =-0.033,ymax = -0.033)+
  coord_cartesian( clip="off") #this keeps it from clipping off the stuff outside the plot
seed.num.pc1.pgls.plot

#ggsave("Germination.Fitness/Results/seed.num.pc1.pgls.pdf", height = 7, width = 7)

all.data.year.final$pred.year1fit.pc1 = predict(fit.pc1.pgls.final)

#define annotations to axis labels for PCA axes
text_x_low = text_grob("Hot & dry", size = 12,col="gray35")
text_x_high = text_grob("Cool & wet", size = 12,col="gray35")

year1fit.pc1.pgls.plot= ggplot(all.data.year.final, aes(y = year1fit.final, x = PC1, label = Pop.2))+
  geom_point(size=4)+
  geom_label_repel(size=5)+
  theme_classic(base_size=22)+
  labs(y = "Year 1 fitness\n (p(flower)*seed number)")+
  annotation_custom(text_x_high,xmin=3.8,xmax= 4,ymin =-0.03,ymax = -0.03) +
  annotation_custom(text_x_low,xmin=-2.8,xmax=-3,ymin =-0.03,ymax = -0.03)+
  coord_cartesian( clip="off") #this keeps it from clipping off the stuff outside the plot
year1fit.pc1.pgls.plot

#ggsave("Germination.Fitness/Results/yr1fit.pc1.pgls.pdf", height = 7, width = 7)

Figure6 = plot_grid(days2bud.pc1.pgls.plot,pflwr.pc1.pgls.plot,seed.num.pc1.pgls.plot,
                    year1fit.pc1.pgls.plot, nrow = 2, ncol = 2, labels = c("A","B","C","D"))
Figure6

#ggsave("Germination.Fitness/Results/Figure6.final.pdf", height = 10, width = 11)

#### Testing for differences in slopes between field and screenhouse species ####

# read in slopes dataframe, also includes standard error of slopes
all.slopes=read.csv("./all.slopes.csv")

# add column of designation
all.slopes$local = c("field","SC","both","field","SC","field","SC",
                       "field","SC","SC")
# remove CAAN since one population is from the field and the other from SC

all.slopes.2 = all.slopes[c(1:2,4:10),]


# t.test to see if mean slope differs between groups

t.test(days.2.bud.final~local, data = all.slopes.2) # p = 0.6479
t.test(sept.1.bud.final~local, data = all.slopes.2) # p = 0.6479

t.test(pflwr.final~local, data = all.slopes.2)  # p = 0.3418
t.test(sizebud.final~local, data = all.slopes.2) # p = 0.2277
t.test(seed_num_nb.final~local, data = all.slopes.2) # p = 0.1302
t.test(year1fit.final~local, data = all.slopes.2) # p = 0.5513
t.test(seed_mass.final~local, data = all.slopes.2) # p = 0.08

# levenne's test to see if variance differs between groups

leveneTest(days.2.bud.final~local, data = all.slopes.2) # p = 0.7036
leveneTest(sept.1.bud.final~local, data = all.slopes.2) # p = 0.7036

leveneTest(pflwr.final~local, data = all.slopes.2)  # p = 0.9089
leveneTest(sizebud.final~local, data = all.slopes.2) # p = 0.7878
leveneTest(seed_num_nb.final~local, data = all.slopes.2) # p = 0.3099
leveneTest(year1fit.final~local, data = all.slopes.2) # p = 0.8160
leveneTest(seed_mass.final~local, data = all.slopes.2) # p = 0.9379






