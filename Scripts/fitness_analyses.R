# This script analyses data to address the question: 
# How does germination timing affect fitness, either directly or indirectly through flowering phenology?

# This script analyses the fitness metrics: flowering probability, number of seeds, total seed mass, year 1 fitness
# Script produces figures 2,3,S4,S6,S7

# libraries
library(tidyverse)
library(lubridate)
library(MASS) #for glm.nb
library(emmeans)
library(car) #for logit function

# inverse logit function
inv.logit <- function(f,a) {
  a <- (1-2*a)
  (a*(1+exp(f))+(exp(f)-1))/(2*a*(1+exp(f)))
}

# load data
fitdat = read.csv("Germination.Fitness/Formatted.Data/final.data.csv")

fitdat = fitdat %>%
  mutate(long.flower.bud.count = na_if(long.flower.bud.count, "N/A"), #need to remove any N/A in data, in lots of columns
         long.mature.fruit.count = na_if(long.mature.fruit.count, "N/A"),
         long.pedicel.count = na_if(long.pedicel.count, "N/A"),
         long.STTO.BH.height.cm = na_if(long.STTO.BH.height.cm, "N/A")) %>%
  mutate(long.flower.bud.count = if_else(long.flower.bud.count == "v", "NA", long.flower.bud.count)) %>%
  separate(Cohort, " ", into = c("text", "cohort")) %>%
  mutate(long.repro.census.date = as.Date(long.repro.census.date, "%m/%d/%y"), 
         short.repro.census.date = as.Date( short.repro.census.date, "%m/%d/%y"),
         First.Bud.Date = as.Date(First.Bud.Date, "%m/%d/%y"),
         First.Flower.Date= as.Date(First.Flower.Date, "%m/%d/%y"),
         First.Fruit.Date = as.Date(First.Fruit.Date, "%m/%d/%y"),
         Death.Date = as.Date(Death.Date, "%m/%d/%y"),
         plant.date = as.Date(plant.date, "%m/%d/%y"),
         transplant.date = as.Date(transplant.date, "%m/%d/%y"),
         rando.date = as.Date(rando.date, "%m/%d/%y")) %>% 
  mutate(cohort= as.numeric(cohort),
         long.flower.bud.count= as.numeric(long.flower.bud.count),
         long.mature.fruit.count = as.numeric(long.mature.fruit.count),
         long.pedicel.count = as.numeric(long.pedicel.count),
         first.flower.height.cm = as.numeric(first.flower.height.cm),
         first.fruit.height.cm = as.numeric(first.fruit.height.cm)) %>%
  mutate(short.flower.bud.count = as.numeric(short.flower.bud.count),
         short.mature.fruit.count = as.numeric(short.mature.fruit.count),
         short.pedicel.count = as.numeric(short.pedicel.count),
         short.STTO.BH.height.cm= as.numeric(short.STTO.BH.height.cm),
         long.STTO.BH.height.cm = as.numeric(long.STTO.BH.height.cm)) %>%
  mutate(cohortfact = as.factor(cohort),
         Mid.Season.Height.Pheno.Stage = as.factor(Mid.Season.Height.Pheno.Stage)) %>%
  mutate(uniqueID = as.factor(paste(Pop, Bench, Block,cohort,Position, sep="_"))) %>% #in case we need to keep track of plants
  mutate(Bench = as.factor(Bench)) %>%
  mutate(Pop = if_else(Pop == "STIN", "STIN1", Pop)) %>%                   
  mutate(Pop = if_else(Pop == "CAAM", "CAAM1", Pop)) %>%       
  mutate(Pop = if_else(Pop == "STDI", "STDI1", Pop)) %>%       
  mutate(Pop = if_else(Pop == "STPO", "STPO1", Pop)) %>%       
  #change Species to be just 4 letter abbrev, 
  mutate(Species = as.factor(substr(Species,1,4))) %>%
  mutate(Pop = as.factor(Pop), Position = as.factor(Position),
         Species = as.factor(Species)) %>%
  #code transplant date as day since Sept 1 
  mutate(transplant.daySept1 = difftime(transplant.date, ymd("2021-09-01"))) %>%
  mutate(transplant.daySept1 = as.numeric(transplant.daySept1))

fitdat = fitdat %>%
  filter(Pop != "CAAM1") %>% # not correct species
  filter(Pop != "STTO_BH") %>% # None of those flowered, population is more biennial than we thought
  #create pflower variable
  mutate(planted = 1,
         flowered = if_else(is.na(First.Flower.Date), 0,1),
         bud_yn = if_else(is.na(First.Bud.Date), 0,1)) %>%
  filter(Season == "LONG") %>%
  filter(is.na(transplant.height.cm) == FALSE)

# add list of transplant dates
transplantdates = c("7-Oct", "28-Oct", "18-Nov", "09-Dec", "30-Dec", "27-Jan", "24-Feb", "24-Mar")

#### Model of probability of flowering ####
# generalized linear model with binomial error and logit link
# predictors are exactly the same as phenology models

pflwr_sep1 = glm(flowered ~ transplant.daySept1*Species + transplant.height.cm + Bench + Pop, 
                 family = binomial(link = "logit"), data=fitdat)
summary(pflwr_sep1)

# test for significant of fixed effects
anova(pflwr_sep1, test = "Chisq") # bench and pop not significant

# estimates species-specific slopes
pflwr_emtrends = emtrends(pflwr_sep1, pairwise ~ Species, var = "transplant.daySept1")
pflwr_emtrends

#write.csv(pflwr_emtrends$emtrends, file = "./Germination.Fitness/Results/prob.flower.emtrends.csv")

#### Figure 2 ####
# new data to predict 
mylist <- list(transplant.daySept1= seq(min(fitdat$transplant.daySept1, na.rm=T),
                                        max(fitdat$transplant.daySept1, na.rm=T)+2, by=10), 
               Species=c("CAAN","CACO","CAIN","STBR", "STDI","STDR","STGL","STIN","STPO","STTO"))

# make predictions from model
pflwr.dat=emmip(pflwr_sep1, Species~transplant.daySept1, at=mylist, CIs=TRUE, plotit = FALSE, type = "response", nesting = F) %>%
  mutate(sigreg = if_else(Species %in% c("CAAN", "CAIN", "STBR", 
                                         "STDR", "STPO", "STTO"), 
                          "Significant", "Nonsignificant"))%>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN")) 

flwrprop = fitdat %>%
  group_by(cohort, Species, Bench, transplant.daySept1) %>%
  summarize(totflwr = sum(flowered),
            totplanted = sum(planted),
            transplant.height.cm = mean(transplant.height.cm, na.rm=T)) %>%
  mutate(yvar = totflwr/totplanted) %>%
  mutate(logitpropflwr = car::logit(yvar, adjust =0.025))%>%
  mutate(sigreg = if_else(Species %in% c("CAAN", "CAIN", "STBR", 
                                         "STDR", "STPO", "STTO"), 
                          "Significant", "Nonsignificant")) %>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))

# plot removing lines for nonsignificant values
pflwr.plot = ggplot(flwrprop, aes(x=transplant.daySept1,y=yvar)) +
  geom_point(show.legend = FALSE, size = 2, alpha = 0.5)+
  theme_classic(base_size=22) + theme(legend.position = "none")  +
  labs(y="Flowering proportion",x="Cohort")+
  scale_x_continuous(breaks = c(unique(fitdat$transplant.daySept1)),labels = transplantdates)+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~phy_order)
pflwr.plot

# only plot line for significant individuals
pflwr.dat.2 = subset(pflwr.dat, pflwr.dat$sigreg == "Significant")

pflwr.plot.2 = pflwr.plot +
  geom_line(data = pflwr.dat.2, aes(x = transplant.daySept1, y = yvar),
            show.legend = FALSE, color = "gray40", lwd = 1)+
  facet_wrap(~phy_order)
pflwr.plot.2

#ggsave("./Germination.Fitness/Results/Prob.flower_byspecies_glm.all.pts.pdf", height = 10, width = 12)

#### Model of size at first bud ####
# filter for individuals that did bud
buddat = fitdat %>%
  filter(bud_yn == 1)
summary(buddat)

sizebud_glob_sep1 = lm(first.bud.height.cm ~transplant.daySept1*Species + transplant.height.cm + Bench + Pop, data=buddat)
summary(sizebud_glob_sep1) 
anova(sizebud_glob_sep1, test = "Chisq")

# estimating species-specific slopes
sizebud_emtrends = emtrends(sizebud_glob_sep1, pairwise ~ Species, var = "transplant.daySept1") 
sizebud_emtrends

# write.csv(sizebud_emtrends$emtrends, file = "./Germination.Fitness/Results/bud.size.emtrends.csv")

#### Figure S4 ####
mylist <- list(transplant.daySept1= seq(min(fitdat$transplant.daySept1, na.rm=T),
                                        max(fitdat$transplant.daySept1, na.rm=T)+2, by=10), 
               Species=c("CAAN","CACO","CAIN","STBR", "STDI","STDR","STGL","STIN","STPO","STTO"))

sizebud.dat= emmip(sizebud_glob_sep1,Species~transplant.daySept1,at=mylist, CIs=TRUE, plotit = FALSE, type = "response", nesting = F) %>%
  mutate(sigreg = if_else(Species %in% c("STBR", "STDI", "STIN", "STPO"), 
                          "Significant", "Nonsignificant")) %>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN")) 

# subset predictions so don't plot line beyond observed data
sizebud.dat = sizebud.dat %>%
  #truncate lines to only cohorts that did have a bud
  mutate(yvar = if_else(Species == "STTO" & transplant.daySept1 >176 , 0, yvar)) %>%
  mutate(yvar = if_else(Species == "STGL" & transplant.daySept1 >186 , 0, yvar)) %>%
  mutate(yvar = if_else(Species == "STPO" & transplant.daySept1 >186 , 0, yvar)) %>%
  mutate(yvar = na_if(yvar,0))

size.all = buddat %>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))

# plot removing lines for nonsignificant values
sizebud.plot = ggplot(size.all, aes(x=transplant.daySept1,y=first.bud.height.cm)) +
  geom_point(show.legend = FALSE, size = 2, alpha = 0.5)+
  theme_classic(base_size=22) + theme(legend.position = "none")  +
  labs(y="Size at budding\n height (cm)",x="Cohort")+
  scale_x_continuous(breaks = c(unique(fitdat$transplant.daySept1)),labels = transplantdates)+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~phy_order)
sizebud.plot

# only plot line for significant individuals
sizebud.dat.2 = subset(sizebud.dat, sizebud.dat$sigreg == "Significant")

sizebud.plot.2 = sizebud.plot +
  geom_line(data = sizebud.dat.2, aes(x = transplant.daySept1, y = yvar),
            show.legend = FALSE, color = "gray40", lwd = 1)+
  facet_wrap(~phy_order)
sizebud.plot.2

#ggsave("./Germination.Fitness/Results/bud.size.all.pts.pdf", height = 10, width = 12)

#### Model of number of seeds ####

flwrdat = fitdat %>%
  filter(flowered == 1)
summary(flwrdat)

#remove STTO outlier
flwrdat2 = flwrdat %>%
  filter(long.seed.counts < 100)

# negative binomial model
seed_nb_glob_sep1 = glm.nb(long.seed.counts~transplant.daySept1*Species+ transplant.height.cm+ Bench + Pop, 
                           data= droplevels(flwrdat2))

summary(seed_nb_glob_sep1)

# test significance of fixed effects
anova(seed_nb_glob_sep1, test = "Chisq")

# estimating speices-specific slopes
seednb_emtrends= emtrends(seed_nb_glob_sep1,pairwise ~ Species, var = "transplant.daySept1") 
seednb_emtrends

# write.csv(seednb_emtrends$emtrends, file = "./Germination.Fitness/Results/seed_number.emtrends.csv")

#### Figure 3 ####

flwrdat3 = flwrdat2 %>%
  mutate(sigreg = if_else(Species %in% c("STDR", "STPO", "STGL", "STTO"), "Nonignificant", "Significant"))%>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))  
# new data to predict from
mylist <- list(transplant.daySept1= seq(min(fitdat$transplant.daySept1, na.rm=T),
                                        max(fitdat$transplant.daySept1, na.rm=T)+2, by=10), 
               Species=c("CAAN","CACO","CAIN","STBR", "STDI","STDR","STGL","STIN","STPO","STTO"))

# predictions from model
seedcount.dat=emmip(seed_nb_glob_sep1,Species~transplant.daySept1,at=mylist, CIs=TRUE, plotit = FALSE, type = "response") %>% # I think response will back transform from logit
  mutate(sigreg = if_else(Species %in% c("STDR", "STPO", "STGL", "STTO"), "Nonignificant", "Significant"))%>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN")) 

# removing predictions so line doesn't go beyond observed data
seedcount.dat = seedcount.dat%>%
  mutate(yvar = if_else(Species == "STTO" & transplant.daySept1 >120 , 0, yvar)) %>%
  mutate(yvar = if_else(Species == "CAIN" & transplant.daySept1 >148 , 0, yvar)) %>%
  mutate(yvar = if_else(Species == "STPO" & transplant.daySept1 >148 , 0, yvar)) %>%
  mutate(yvar = na_if(yvar,0))

# plot removing lines for nonsignificant values
seedcount.plot = ggplot(flwrdat3, aes(x=transplant.daySept1,y=long.seed.counts)) +
  geom_point(show.legend = FALSE, size = 2, alpha = 0.5)+
  theme_classic(base_size=22) + theme(legend.position = "none")  +
  labs(y="Number of seeds",x="Cohort")+
  scale_x_continuous(breaks = c(unique(fitdat$transplant.daySept1)),labels = transplantdates)+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~phy_order)
seedcount.plot

# only plot line for significant individuals
seedcount.dat.2 = subset(seedcount.dat, seedcount.dat$sigreg == "Significant")

# remove last STBR line since no data
seedcount.dat.3 = seedcount.dat.2[c(1:99,101:105,107,108),]

seedcount.plot.2 = seedcount.plot +
  geom_line(data = seedcount.dat.3, aes(x = transplant.daySept1, y = yvar),
            show.legend = FALSE, color = "gray40", lwd = 1)+
  facet_wrap(~phy_order)
seedcount.plot.2

#ggsave("./Germination.Fitness/Results/seed.number.all.pts.pdf", height = 10, width = 12)

#### Model of total seed mass ####

seedmass_glob_sep1 = lm(long.seed.weight.g~transplant.daySept1*Species + transplant.height.cm+ Bench + Pop, data= droplevels(flwrdat2))
summary(seedmass_glob_sep1)

# testing significance of fixed effects
anova(seedmass_glob_sep1, test = "Chisq")

# estimating species-specific slopes
seedmass_emtrends = emtrends(seedmass_glob_sep1,pairwise ~ Species, var = "transplant.daySept1") 
seedmass_emtrends

# write.csv(seedmass_emtrends$emtrends, file = "./Germination.Fitness/Results/seed.mass.emtrends.csv")

#### Figure S6 ####

flwrdat3 = flwrdat2 %>%
  mutate(sigreg = if_else(Species %in% c("CAAN", "CACO", "STDI", "STIN"), "Significant", "Nonsignificant"))%>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))  

# new data to predict from
mylist <- list(transplant.daySept1= seq(min(fitdat$transplant.daySept1, na.rm=T),
                                        max(fitdat$transplant.daySept1, na.rm=T)+2, by=10), 
               Species=c("CAAN","CACO","CAIN","STBR", "STDI","STDR","STGL","STIN","STPO","STTO"))

# make predictions from model
seedmass.dat=emmip(seedmass_glob_sep1,Species~transplant.daySept1,at=mylist, CIs=TRUE, plotit = FALSE, type = "response") %>% # I think response will back transform from logit
  mutate(sigreg = if_else(Species %in% c("CAAN","CACO", "STDI", "STIN"), "Significant", "Nonsignificant")) 

# remove predictions that go beyond observed data
seedmass.dat = seedmass.dat%>%
  mutate(yvar = if_else(Species == "STTO" & transplant.daySept1 >120 , 0, yvar)) %>%
  mutate(yvar = if_else(Species == "CAIN" & transplant.daySept1 >148 , 0, yvar)) %>%
  mutate(yvar = if_else(Species == "STPO" & transplant.daySept1 >148 , 0, yvar)) %>%
  mutate(yvar = na_if(yvar,0))%>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN")) 

# plot removing lines for nonsignificant values
seedmass.plot = ggplot(flwrdat3, aes(x=transplant.daySept1,y=long.seed.weight.g)) +
  geom_point(show.legend = FALSE, size = 2, alpha = 0.5)+
  theme_classic(base_size=22) + theme(legend.position = "none")  +
  labs(y="Seed mass (g)",x="Cohort")+
  scale_x_continuous(breaks = c(unique(fitdat$transplant.daySept1)),labels = transplantdates)+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~phy_order)
seedmass.plot

# only plot line for significant individuals
seedmass.dat.2 = subset(seedmass.dat, seedmass.dat$sigreg == "Significant")

seedmass.plot.2 = seedmass.plot +
  geom_line(data = seedmass.dat.2, aes(x = transplant.daySept1, y = yvar),
            show.legend = FALSE, color = "gray40", lwd = 1)+
  facet_wrap(~phy_order)
seedmass.plot.2

#ggsave("./Germination.Fitness/Results/seed.mass.all.pts.pdf", height = 10, width = 12)

#### Model of first year fitness ####
# calculate first year fitness
fitdat = fitdat %>%
  mutate(year1fit = flowered * long.seed.counts)

# removed STTO outlier
fitdat2 = fitdat %>%
  filter(long.seed.counts < 100)

year1_nb_glob_sep1 = glm.nb(year1fit~transplant.daySept1*Species + transplant.height.cm + Bench + Pop, data= droplevels(fitdat2))
summary(year1_nb_glob_sep1)

# testing significance of fixed effects
anova(year1_nb_glob_sep1, test = "Chisq")

# estimating species-specific slopes
year1_emtrends = emtrends(year1_nb_glob_sep1,pairwise ~ Species, var = "transplant.daySept1") 
year1_emtrends

# write.csv(year1_emtrends$emtrends, file = "./Germination.Fitness/Results/year1fit.emtrends.csv")

#### Figure S7 ####

fitdat3 = fitdat2 %>%
  mutate(sigreg = if_else(Species %in% c("STDR", "STGL", "STPO"), "Nonsignificant", "Significant"))%>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))  

# make predictions
year1.dat=emmip(year1_nb_glob_sep1,Species~transplant.daySept1,at=mylist, CIs=TRUE, plotit = FALSE, type = "response") %>% # I think response will back transform from logit
  mutate(sigreg = if_else(Species %in% c("STDR", "STGL", "STPO"), "Nonsignificant", "Significant"))%>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN")) 

# plot removing lines for nonsignificant values
year1.plot = ggplot(fitdat3, aes(x=transplant.daySept1,y=year1fit)) +
  geom_point(show.legend = FALSE, size = 2, alpha = 0.5)+
  theme_classic(base_size=22) + theme(legend.position = "none")  +
  labs(y="Year 1 fitness\n (p(flower)*seed number)",x="Cohort")+
  scale_x_continuous(breaks = c(unique(fitdat$transplant.daySept1)),labels = transplantdates)+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~phy_order)
year1.plot

# only plot line for significant individuals
year1.dat.2 = subset(year1.dat, year1.dat$sigreg == "Significant")

year1.plot.2 = year1.plot +
  geom_line(data = year1.dat.2, aes(x = transplant.daySept1, y = yvar),
            show.legend = FALSE, color = "gray40", lwd = 1)+
  facet_wrap(~phy_order)
year1.plot.2

#ggsave("./Germination.Fitness/Results/year1fitness.allpoints.pdf", height = 10, width = 12)

