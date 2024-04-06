# This script analyses data to address the question: 
# How have responses of flowering time and fitness to germination timing evolved across the clade?
# Do species with greater compensatory flowering plasticity have greater fitness stability in the face of variable germination timing?
# Is there evidence of climate adaptation?

# load libraries
library(ape)
library(phytools)
library(RColorBrewer)
library(ggnewscale)
library(caper)

remotes::install_github("clauswilke/colorblindr")
library(colorblindr)

BiocManager::install("ggtree", force = TRUE)
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
model.slopes=read.csv("Germination.Fitness/Results/pheno.slopes.csv")

# days.2.bud and sept.1.bud are direct opposites of each other, difference in p is due to randomization

phylosig(phylo, model.slopes$days.2.bud, method = "K", test = TRUE, se = model.slopes$days.2.bud.se, nsim=1000)
# K = 1.18966, p = 0.136
phylosig(phylo, model.slopes$sept.1.bud, method = "K", test = TRUE, se = model.slopes$sept.1.bud.se,nsim=1000)
# K = 1.18966, p = 0.151

# fitness slopes
fitness.slopes=read.csv("Germination.Fitness/Results/fitness.slopes.2.csv")

phylosig(phylo, fitness.slopes$pflwr, method = "K", test = TRUE, se = fitness.slopes$pflwr.se,nsim=1000)
# K = 0.973273, p = 0.908
phylosig(phylo, fitness.slopes$seed_num_nb, method = "K", test = TRUE, se = fitness.slopes$seed_num_nb.se,nsim=1000)
# K = 1.06775, p = 0.842
phylosig(phylo, fitness.slopes$year1fit, method = "K", test = TRUE, se = fitness.slopes$year1fit.se,nsim=1000)
# K = 0.682429, p = 0.9997
phylosig(phylo, fitness.slopes$seed_mass, method = "K", se = fitness.slopes$seed_mass.se,test = TRUE, nsim=1000)
# K = 0.771932, p = 0.442

#### Figure 5 #####

all.phylo <- ggtree::read.tree("./Germination.Fitness/Raw.Data/tree_pruned.new")

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

#### Figure S8 ####

all.phylo <- ggtree::read.tree("./Germination.Fitness/Raw.Data/tree_pruned.new")

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



#### PGLS to relate phenology and fitness slopes ####
# read in tree
all.phylo <- read.tree("./Germination.Fitness/Raw.Data/tree_pruned.new")

sp.list=c("Caulanthus_anceps","Caulanthus_coulteri","Caulanthus_inflatus",
          "Streptanthus_breweri","Streptanthus_diversifolius","Streptanthus_drepanoides",
          "Streptanthus_glandulosus","Streptanthus_insignis","Streptanthus_polygaloides",
          "Streptanthus_tortuosus")

# prune phylogeny to only includes our species
phylo=keep.tip(all.phylo, sp.list)

# pgls using caper package

comp.data.1<-comparative.data(phylo, all.slopes, names.col="sp", vcv.dim=2, warn.dropped=TRUE, vcv = TRUE)

d2bud.seed.num.pgls = pgls(days.2.bud~seed_num_nb, data=comp.data.1)
summary(d2bud.seed.num.pgls) # not-significant, p = 0.07

d2bud.pflwr.pgls = pgls(days.2.bud~pflwr, data=comp.data.1)
summary(d2bud.pflwr.pgls) # not significant p = 0.78

d2bud.yr1fit.pgls = pgls(days.2.bud~year1fit, data=comp.data.1)
summary(d2bud.yr1fit.pgls) # not significant p = 0.49

#### Figure S9 ####

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

#### PGLS to relate germination slopes to fitness slopes ####
# Germination timing slopes come from Worthy et al. 2023 bioRxiv

all.slopes=read.csv("Germination.Fitness/Results/all.slopes.csv")

# add slopes from germ pheno: germination_fraction~cohort
# had to average two CAIN and two CAAN slopes
all.slopes$pheno.slopes = c(-0.45898,-0.1238,-0.021865,-0.59355,-0.33908,-0.25232,-0.40204,-0.2401,
                            -0.36986,-0.36455)

all.phylo <- read.tree("./Germination.Fitness/Raw.Data/tree_pruned.new")

sp.list=c("Caulanthus_anceps","Caulanthus_coulteri","Caulanthus_inflatus",
          "Streptanthus_breweri","Streptanthus_diversifolius","Streptanthus_drepanoides",
          "Streptanthus_glandulosus","Streptanthus_insignis","Streptanthus_polygaloides",
          "Streptanthus_tortuosus")

# prune phylogeny to only includes our species
phylo=keep.tip(all.phylo, sp.list)

comp.data<-comparative.data(phylo, all.slopes, names.col="sp", vcv.dim=2, warn.dropped=TRUE, vcv = TRUE)

pheno.slopes.seed.num.pgls = pgls(seed_num_nb~pheno.slopes, data=comp.data)
summary(pheno.slopes.seed.num.pgls) # not-significant, p = 0.50

pheno.slopes.year1fit.pgls = pgls(year1fit~pheno.slopes, data=comp.data)
summary(pheno.slopes.year1fit.pgls) # not-significant, p = 0.71

#### Figure S10 ####

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


