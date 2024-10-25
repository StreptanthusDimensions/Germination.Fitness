# This script analyses data to address the question: 
# How does germination timing affect fitness, either directly or indirectly through flowering phenology?

# This script analyses the fitness metrics: flowering probability, number of seeds, total seed mass, year 1 fitness
# Script produces figures 2,3,S4,S6,S7

# Users, please specify your own path to load data in this script.

# libraries
library(tidyverse) # version 2.0.0
library(lubridate) # version 1.9.2 
library(MASS) # for glm.nb version 7.3.58.2
library(emmeans) # version 1.8.7
library(car) # for logit function version 3.1.2

# inverse logit function
inv.logit <- function(f,a) {
  a <- (1-2*a)
  (a*(1+exp(f))+(exp(f)-1))/(2*a*(1+exp(f)))
}

# load data
fitdat = read.csv("./final.data.csv") # final.data.csv generated using data_cleaning.R script

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

fitdat$seed.source = if_else(fitdat$Pop %in% c("CAAN1","CACO1","STBR3","STDI","STDR2","STIN"), "SH", "Field")
fitdat$seed.source = as.factor(fitdat$seed.source)

# add list of transplant dates
transplantdates = c("7-Oct", "28-Oct", "18-Nov", "09-Dec", "30-Dec", "27-Jan", "24-Feb", "24-Mar")

#### Model of probability of flowering ####
# generalized linear model with binomial error and logit link
# predictors are exactly the same as phenology models

pflwr_sep1 = glm(flowered ~ transplant.daySept1*Species + transplant.height.cm + Bench + Pop, 
                 family = binomial(link = "logit"), data=droplevels(fitdat))
summary(pflwr_sep1)

# test for significant of fixed effects
anova(pflwr_sep1, test = "Chisq") # bench and pop not significant

# estimates species-specific slopes
pflwr_emtrends = emtrends(pflwr_sep1, pairwise ~ Species, var = "transplant.daySept1")
pflwr_emtrends

# write.csv(pflwr_emtrends$emtrends, file = "./Germination.Fitness/Results/prob.flower.emtrends.csv")

#### Model of probability of flowering with weights ####
# create data frame that has successes (flowered) and total numbers of individuals on each bench, in each cohort, for each pop

n.cohort = fitdat %>%
  group_by(Pop,Bench,transplant.daySept1) %>%
  count()

n.flowered = fitdat %>%
  filter(flowered == 1) %>%
  group_by(Pop,Bench,transplant.daySept1,flowered) %>%
  count()

colnames(n.flowered)[5] = "Yes.Flower"

no.flowered = fitdat %>%
  filter(flowered == 0) %>%
  group_by(Pop,Bench,transplant.daySept1,flowered) %>%
  count()

colnames(no.flowered)[4] = "no.flowered"
colnames(no.flowered)[5] = "No.Flower"

all.data = full_join(n.flowered,no.flowered)
all.data$total.n = rowSums(all.data[,c("Yes.Flower", "No.Flower")], na.rm=TRUE)

# change NA is Yes.Flower to 0

all.data$Yes.Flower[is.na(all.data$Yes.Flower)] = 0

all.data$seed.source = if_else(all.data$Pop %in% c("CAAN1","CACO1","STBR3","STDI","STDR2","STIN"), "SH", "Field")

all.data = all.data %>%
  mutate(Species = dplyr::recode(Pop, "CAAN1" = "CAAN", "CAAN2" = "CAAN","CAIN3" = "CAIN", "CAIN4" = "CAIN","CACO1" = "CACO", 
                                 "STTO_TM2" = "STTO", "STDR2" = "STDR", "STPO1" = "STPO", "STBR3" = "STBR", 
                                 "STDI1" = "STDI", "STGL1" = "STGL", "STIN1" = "STIN")) %>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))

# generalized linear model with binomial error and logit link
# can't include transplant height as a predictor since now data is successes/failures

y=cbind(all.data$Yes.Flower, (all.data$total.n-all.data$Yes.Flower))

pflwr_sep1_weighted = glm(y ~ transplant.daySept1*Species + Bench + Pop, 
                 family = binomial(link = "logit"), data=droplevels(all.data))
summary(pflwr_sep1_weighted)

# test for significant of fixed effects
anova(pflwr_sep1_weighted, test = "Chisq") # bench and pop not significant

# estimates species-specific slopes
pflwr_weighted_emtrends = emtrends(pflwr_sep1_weighted, pairwise ~ Species, var = "transplant.daySept1")
pflwr_weighted_emtrends

# write.csv(pflwr_weighted_emtrends$emtrends, file = "./Germination.Fitness/Results/prob.flower.weighted.emtrends.csv")

#### Model of probability of flowering with weights FINAL ####

# generalized linear model with binomial error and logit link
# can't include transplant height as a predictor since now data is successes/failures

y=cbind(all.data$Yes.Flower, (all.data$total.n-all.data$Yes.Flower))

pflwr_sep1_weighted.final = glm(y ~ transplant.daySept1*Species + Bench, 
                          family = binomial(link = "logit"), data=all.data)
summary(pflwr_sep1_weighted.final)

# test for significant of fixed effects
anova(pflwr_sep1_weighted.final, test = "Chisq") # bench not significant

# estimates species-specific slopes
pflwr_weighted_emtrends.final = emtrends(pflwr_sep1_weighted.final, pairwise ~ Species, var = "transplant.daySept1")
pflwr_weighted_emtrends.final

# write.csv(pflwr_weighted_emtrends.final$emtrends, file = "./Germination.Fitness/Results/prob.flower.weighted.emtrends.final.csv")

#### Model of probability of flowering with weights with CAAN1 and CAIN3 ####

all.data.sub = all.data %>%
  filter(!Pop %in% c("CAAN2","CAIN4"))

# generalized linear model with binomial error and logit link
# can't include transplant height as a predictor since now data is successes/failures

y=cbind(all.data.sub$Yes.Flower, (all.data.sub$total.n-all.data.sub$Yes.Flower))

pflwr_sep1_weighted.nopop.ss = glm(y ~ transplant.daySept1*Species + Bench, 
                             family = binomial(link = "logit"), data=all.data.sub)
summary(pflwr_sep1_weighted.nopop.ss)

# test for significant of fixed effects
anova(pflwr_sep1_weighted.nopop.ss, test = "Chisq") # bench and pop not significant

# estimates species-specific slopes
pflwr_weighted_emtrends.nopop.ss = emtrends(pflwr_sep1_weighted.nopop.ss, pairwise ~ Species, var = "transplant.daySept1")
pflwr_weighted_emtrends.nopop.ss

#write.csv(pflwr_weighted_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/prob.flower.weighted.emtrends.CAAN1.CAIN3.csv")

#### Model of probability of flowering with weights with CAAN1 and CAIN4 ####

all.data.sub = all.data %>%
  filter(!Pop %in% c("CAAN2","CAIN3"))

# generalized linear model with binomial error and logit link
# can't include transplant height as a predictor since now data is successes/failures

y=cbind(all.data.sub$Yes.Flower, (all.data.sub$total.n-all.data.sub$Yes.Flower))

pflwr_sep1_weighted.nopop.ss = glm(y ~ transplant.daySept1*Species + Bench, 
                                   family = binomial(link = "logit"), data=all.data.sub)
summary(pflwr_sep1_weighted.nopop.ss)

# test for significant of fixed effects
anova(pflwr_sep1_weighted.nopop.ss, test = "Chisq") # bench and pop not significant

# estimates species-specific slopes
pflwr_weighted_emtrends.nopop.ss = emtrends(pflwr_sep1_weighted.nopop.ss, pairwise ~ Species, var = "transplant.daySept1")
pflwr_weighted_emtrends.nopop.ss

# write.csv(pflwr_weighted_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/prob.flower.weighted.emtrends.CAAN1.CAIN4.csv")

#### Model of probability of flowering with weights with CAAN2 and CAIN3 ####

all.data.sub = all.data %>%
  filter(!Pop %in% c("CAAN1","CAIN4"))

# generalized linear model with binomial error and logit link
# can't include transplant height as a predictor since now data is successes/failures

y=cbind(all.data.sub$Yes.Flower, (all.data.sub$total.n-all.data.sub$Yes.Flower))

pflwr_sep1_weighted.nopop.ss = glm(y ~ transplant.daySept1*Species + Bench, 
                                   family = binomial(link = "logit"), data=all.data.sub)
summary(pflwr_sep1_weighted.nopop.ss)

# test for significant of fixed effects
anova(pflwr_sep1_weighted.nopop.ss, test = "Chisq") # bench and pop not significant

# estimates species-specific slopes
pflwr_weighted_emtrends.nopop.ss = emtrends(pflwr_sep1_weighted.nopop.ss, pairwise ~ Species, var = "transplant.daySept1")
pflwr_weighted_emtrends.nopop.ss

write.csv(pflwr_weighted_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/prob.flower.weighted.emtrends.CAAN2.CAIN3.csv")

#### Model of probability of flowering with weights with CAAN2 and CAIN4 ####

all.data.sub = all.data %>%
  filter(!Pop %in% c("CAAN1","CAIN3"))

# generalized linear model with binomial error and logit link
# can't include transplant height as a predictor since now data is successes/failures

y=cbind(all.data.sub$Yes.Flower, (all.data.sub$total.n-all.data.sub$Yes.Flower))

pflwr_sep1_weighted.nopop.ss = glm(y ~ transplant.daySept1*Species + Bench, 
                                   family = binomial(link = "logit"), data=all.data.sub)
summary(pflwr_sep1_weighted.nopop.ss)

# test for significant of fixed effects
anova(pflwr_sep1_weighted.nopop.ss, test = "Chisq") # bench and pop not significant

# estimates species-specific slopes
pflwr_weighted_emtrends.nopop.ss = emtrends(pflwr_sep1_weighted.nopop.ss, pairwise ~ Species, var = "transplant.daySept1")
pflwr_weighted_emtrends.nopop.ss

#write.csv(pflwr_weighted_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/prob.flower.weighted.emtrends.CAAN2.CAIN4.csv")

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

#### Figure 2 REVISED ####
# new data to predict 
mylist <- list(transplant.daySept1= seq(min(fitdat$transplant.daySept1, na.rm=T),
                                        max(fitdat$transplant.daySept1, na.rm=T)+2, by=10), 
               Species=c("CAAN","CACO","CAIN","STBR", "STDI","STDR","STGL","STIN","STPO","STTO"))

# make predictions from model
pflwr.dat=emmip(pflwr_sep1_weighted.final, Species~transplant.daySept1, at=mylist, CIs=TRUE, plotit = FALSE, type = "response", nesting = F) %>%
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
  mutate(logitpropflwr = car::logit(yvar, adjust = 0.025))%>%
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

# Slight differences in lines on each plot

# ggsave("./Germination.Fitness/Results/Prob.flower_byspecies_glm.all.pts.final.pdf", height = 10, width = 12)

#### Model of size at first bud ####
# filter for individuals that did bud
buddat = fitdat %>%
  filter(bud_yn == 1)
summary(buddat)

sizebud_glob_sep1 = lm(first.bud.height.cm ~transplant.daySept1*Species + transplant.height.cm + Bench + Pop, data=buddat)
summary(sizebud_glob_sep1) 
anova(sizebud_glob_sep1)

# estimating species-specific slopes
sizebud_emtrends = emtrends(sizebud_glob_sep1, pairwise ~ Species, var = "transplant.daySept1") 
sizebud_emtrends

# write.csv(sizebud_emtrends$emtrends, file = "./Germination.Fitness/Results/bud.size.emtrends.csv")

#### Model of size at first bud FINAL ####
# filter for individuals that did bud
buddat = fitdat %>%
  filter(bud_yn == 1)
summary(buddat)

sizebud_glob_sep1.final = lm(first.bud.height.cm ~transplant.daySept1*Species + transplant.height.cm + Bench, data=buddat)
summary(sizebud_glob_sep1.final) 
anova(sizebud_glob_sep1.final)

# estimating species-specific slopes
sizebud_emtrends.final = emtrends(sizebud_glob_sep1.final, pairwise ~ Species, var = "transplant.daySept1") 
sizebud_emtrends.final

#write.csv(sizebud_emtrends.final$emtrends, file = "./Germination.Fitness/Results/bud.size.emtrends.final.csv")

#### Model of size at first bud with CAAN1 and CAIN3 ####
buddat.sub = fitdat %>%
  filter(bud_yn == 1) %>%
  filter(!Pop %in% c("CAAN2","CAIN4"))
summary(buddat.sub)

sizebud_glob_sep1.nopop.ss = lm(first.bud.height.cm ~ transplant.daySept1*Species + transplant.height.cm + Bench, data=buddat.sub)
summary(sizebud_glob_sep1.nopop.ss) 
anova(sizebud_glob_sep1.nopop.ss)

# estimating species-specific slopes
sizebud_emtrends.nopop.ss = emtrends(sizebud_glob_sep1.nopop.ss, pairwise ~ Species, var = "transplant.daySept1") 
sizebud_emtrends.nopop.ss

write.csv(sizebud_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/bud.size.emtrends.CAAN1.CAIN3.csv")

#### Model of size at first bud with CAAN1 and CAIN4 ####
buddat.sub = fitdat %>%
  filter(bud_yn == 1) %>%
  filter(!Pop %in% c("CAAN2","CAIN3"))
summary(buddat.sub)

sizebud_glob_sep1.nopop.ss = lm(first.bud.height.cm ~ transplant.daySept1*Species + transplant.height.cm + Bench, data=buddat.sub)
summary(sizebud_glob_sep1.nopop.ss) 
anova(sizebud_glob_sep1.nopop.ss)

# estimating species-specific slopes
sizebud_emtrends.nopop.ss = emtrends(sizebud_glob_sep1.nopop.ss, pairwise ~ Species, var = "transplant.daySept1") 
sizebud_emtrends.nopop.ss

write.csv(sizebud_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/bud.size.emtrends.CAAN1.CAIN4.csv")

#### Model of size at first bud with CAAN2 and CAIN3 ####
buddat.sub = fitdat %>%
  filter(bud_yn == 1) %>%
  filter(!Pop %in% c("CAAN1","CAIN4"))
summary(buddat.sub)

sizebud_glob_sep1.nopop.ss = lm(first.bud.height.cm ~ transplant.daySept1*Species + transplant.height.cm + Bench, data=buddat.sub)
summary(sizebud_glob_sep1.nopop.ss) 
anova(sizebud_glob_sep1.nopop.ss)

# estimating species-specific slopes
sizebud_emtrends.nopop.ss = emtrends(sizebud_glob_sep1.nopop.ss, pairwise ~ Species, var = "transplant.daySept1") 
sizebud_emtrends.nopop.ss

#write.csv(sizebud_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/bud.size.emtrends.CAAN2.CAIN3.csv")

#### Model of size at first bud with CAAN2 and CAIN4 ####
buddat.sub = fitdat %>%
  filter(bud_yn == 1) %>%
  filter(!Pop %in% c("CAAN1","CAIN3"))
summary(buddat.sub)

sizebud_glob_sep1.nopop.ss = lm(first.bud.height.cm ~ transplant.daySept1*Species + transplant.height.cm + Bench, data=buddat.sub)
summary(sizebud_glob_sep1.nopop.ss) 
anova(sizebud_glob_sep1.nopop.ss)

# estimating species-specific slopes
sizebud_emtrends.nopop.ss = emtrends(sizebud_glob_sep1.nopop.ss, pairwise ~ Species, var = "transplant.daySept1") 
sizebud_emtrends.nopop.ss

write.csv(sizebud_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/bud.size.emtrends.CAAN2.CAIN4.csv")

#### Model of size at first flower ####
flwrdat = fitdat %>%
  filter(flowered == 1)
summary(flwrdat)

sizeflwr_glob_sep1 = lm(first.flower.height.cm ~transplant.daySept1*Species + transplant.height.cm + Bench + Pop, data=flwrdat)
summary(sizeflwr_glob_sep1) 
anova(sizeflwr_glob_sep1)

# estimating species-specific slopes
sizeflwr_emtrends = emtrends(sizeflwr_glob_sep1, pairwise ~ Species, var = "transplant.daySept1") 
sizeflwr_emtrends

#write.csv(sizeflwr_emtrends$emtrends, file = "./Germination.Fitness/Results/Flower/flower.size.emtrends.csv")

#### Model of size at first flower with CAAN1 and CAIN3 ####
flwrdat.sub = fitdat %>%
  filter(flowered == 1) %>%
  filter(!Pop %in% c("CAAN2","CAIN4"))
summary(flwrdat.sub)

sizeflwr_glob_sep1.nopop.ss = lm(first.flower.height.cm ~ transplant.daySept1*Species + transplant.height.cm + Bench, data=flwrdat.sub)
summary(sizeflwr_glob_sep1.nopop.ss) 
anova(sizeflwr_glob_sep1.nopop.ss)

# estimating species-specific slopes
sizeflwr_emtrends.nopop.ss = emtrends(sizeflwr_glob_sep1.nopop.ss, pairwise ~ Species, var = "transplant.daySept1") 
sizeflwr_emtrends.nopop.ss

write.csv(sizeflwr_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Flower/flower.size.emtrends.CAAN1.CAIN3.csv")

#### Model of size at first flower with CAAN1 and CAIN4 ####
flwrdat.sub = fitdat %>%
  filter(flowered == 1) %>%
  filter(!Pop %in% c("CAAN2","CAIN3"))
summary(flwrdat.sub)

sizeflwr_glob_sep1.nopop.ss = lm(first.flower.height.cm ~ transplant.daySept1*Species + transplant.height.cm + Bench, data=flwrdat.sub)
summary(sizeflwr_glob_sep1.nopop.ss) 
anova(sizeflwr_glob_sep1.nopop.ss)

# estimating species-specific slopes
sizeflwr_emtrends.nopop.ss = emtrends(sizeflwr_glob_sep1.nopop.ss, pairwise ~ Species, var = "transplant.daySept1") 
sizeflwr_emtrends.nopop.ss

# write.csv(sizeflwr_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Flower/flower.size.emtrends.CAAN1.CAIN4.csv")

#### Model of size at first flower with CAAN2 and CAIN3 ####
flwrdat.sub = fitdat %>%
  filter(flowered == 1) %>%
  filter(!Pop %in% c("CAAN1","CAIN4"))
summary(flwrdat.sub)

sizeflwr_glob_sep1.nopop.ss = lm(first.flower.height.cm ~ transplant.daySept1*Species + transplant.height.cm + Bench, data=flwrdat.sub)
summary(sizeflwr_glob_sep1.nopop.ss) 
anova(sizeflwr_glob_sep1.nopop.ss)

# estimating species-specific slopes
sizeflwr_emtrends.nopop.ss = emtrends(sizeflwr_glob_sep1.nopop.ss, pairwise ~ Species, var = "transplant.daySept1") 
sizeflwr_emtrends.nopop.ss

# write.csv(sizeflwr_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Flower/flower.size.emtrends.CAAN2.CAIN3.csv")

#### Model of size at first flower with CAAN2 and CAIN4 ####
flwrdat.sub = fitdat %>%
  filter(flowered == 1) %>%
  filter(!Pop %in% c("CAAN1","CAIN3"))
summary(flwrdat.sub)

sizeflwr_glob_sep1.nopop.ss = lm(first.flower.height.cm ~ transplant.daySept1*Species + transplant.height.cm + Bench, data=flwrdat.sub)
summary(sizeflwr_glob_sep1.nopop.ss) 
anova(sizeflwr_glob_sep1.nopop.ss)

# estimating species-specific slopes
sizeflwr_emtrends.nopop.ss = emtrends(sizeflwr_glob_sep1.nopop.ss, pairwise ~ Species, var = "transplant.daySept1") 
sizeflwr_emtrends.nopop.ss

#write.csv(sizeflwr_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Flower/flower.size.emtrends.CAAN2.CAIN4.csv")


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

#### Figure S4 REVISED ####
mylist <- list(transplant.daySept1= seq(min(fitdat$transplant.daySept1, na.rm=T),
                                        max(fitdat$transplant.daySept1, na.rm=T)+2, by=10), 
               Species=c("CAAN","CACO","CAIN","STBR", "STDI","STDR","STGL","STIN","STPO","STTO"))

sizebud.dat= emmip(sizebud_glob_sep1.final,Species~transplant.daySept1,at=mylist, CIs=TRUE, plotit = FALSE, type = "response", nesting = F) %>%
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

# Same as original plot

#ggsave("./Germination.Fitness/Results/bud.size.all.pts.final.pdf", height = 10, width = 12)

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

#### Model of number of seeds FINAL ####
flwrdat = fitdat %>%
  filter(flowered == 1)
summary(flwrdat)

#remove STTO outlier
flwrdat2 = flwrdat %>%
  filter(long.seed.counts < 100)

# negative binomial model
seed_nb_glob_sep1.final = glm.nb(long.seed.counts~transplant.daySept1*Species+ transplant.height.cm + Bench, 
                           data= droplevels(flwrdat2))

summary(seed_nb_glob_sep1.final)

# test significance of fixed effects
anova(seed_nb_glob_sep1.final, test = "Chisq")

# estimating speices-specific slopes
seednb_emtrends.final = emtrends(seed_nb_glob_sep1.final,pairwise ~ Species, var = "transplant.daySept1") 
seednb_emtrends.final

#write.csv(seednb_emtrends.final$emtrends, file = "./Germination.Fitness/Results/seed_number.emtrends.final.csv")

#### Model of number of seeds with CAAN1 and CAIN3 ####
flwrdat.sub = fitdat %>%
  filter(flowered == 1) %>%
  filter(!Pop %in% c("CAAN2","CAIN4"))
summary(flwrdat.sub)

#remove STTO outlier
flwrdat2 = flwrdat.sub %>%
  filter(long.seed.counts < 100)

# negative binomial model
seed_nb_glob_sep1.nopop.ss = glm.nb(long.seed.counts~transplant.daySept1*Species+ transplant.height.cm+ Bench, 
                           data= droplevels(flwrdat2))

summary(seed_nb_glob_sep1.nopop.ss)

# test significance of fixed effects
anova(seed_nb_glob_sep1.nopop.ss, test = "Chisq")

# estimating speices-specific slopes
seednb_emtrends.nopop.ss= emtrends(seed_nb_glob_sep1.nopop.ss,pairwise ~ Species, var = "transplant.daySept1") 
seednb_emtrends.nopop.ss

#write.csv(seednb_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/seed_number.emtrends.CAAN1.CAIN3.csv")

#### Model of number of seeds with CAAN1 and CAIN4 ####
flwrdat.sub = fitdat %>%
  filter(flowered == 1) %>%
  filter(!Pop %in% c("CAAN2","CAIN3"))
summary(flwrdat.sub)

#remove STTO outlier
flwrdat2 = flwrdat.sub %>%
  filter(long.seed.counts < 100)

# negative binomial model
seed_nb_glob_sep1.nopop.ss = glm.nb(long.seed.counts~transplant.daySept1*Species+ transplant.height.cm+ Bench, 
                                    data= droplevels(flwrdat2))

summary(seed_nb_glob_sep1.nopop.ss)

# test significance of fixed effects
anova(seed_nb_glob_sep1.nopop.ss, test = "Chisq")

# estimating speices-specific slopes
seednb_emtrends.nopop.ss= emtrends(seed_nb_glob_sep1.nopop.ss,pairwise ~ Species, var = "transplant.daySept1") 
seednb_emtrends.nopop.ss

#write.csv(seednb_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/seed_number.emtrends.CAAN1.CAIN4.csv")

#### Model of number of seeds with CAAN2 and CAIN3 ####
flwrdat.sub = fitdat %>%
  filter(flowered == 1) %>%
  filter(!Pop %in% c("CAAN1","CAIN4"))
summary(flwrdat.sub)

#remove STTO outlier
flwrdat2 = flwrdat.sub %>%
  filter(long.seed.counts < 100)

# negative binomial model
seed_nb_glob_sep1.nopop.ss = glm.nb(long.seed.counts~transplant.daySept1*Species+ transplant.height.cm+ Bench, 
                                    data= droplevels(flwrdat2))

summary(seed_nb_glob_sep1.nopop.ss)

# test significance of fixed effects
anova(seed_nb_glob_sep1.nopop.ss, test = "Chisq")

# estimating speices-specific slopes
seednb_emtrends.nopop.ss= emtrends(seed_nb_glob_sep1.nopop.ss,pairwise ~ Species, var = "transplant.daySept1") 
seednb_emtrends.nopop.ss

write.csv(seednb_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/seed_number.emtrends.CAAN2.CAIN3.csv")

#### Model of number of seeds with CAAN2 and CAIN4 ####
flwrdat.sub = fitdat %>%
  filter(flowered == 1) %>%
  filter(!Pop %in% c("CAAN1","CAIN3"))
summary(flwrdat.sub)

#remove STTO outlier
flwrdat2 = flwrdat.sub %>%
  filter(long.seed.counts < 100)

# negative binomial model
seed_nb_glob_sep1.nopop.ss = glm.nb(long.seed.counts~transplant.daySept1*Species+ transplant.height.cm+ Bench, 
                                    data= droplevels(flwrdat2))

summary(seed_nb_glob_sep1.nopop.ss)

# test significance of fixed effects
anova(seed_nb_glob_sep1.nopop.ss, test = "Chisq")

# estimating speices-specific slopes
seednb_emtrends.nopop.ss= emtrends(seed_nb_glob_sep1.nopop.ss,pairwise ~ Species, var = "transplant.daySept1") 
seednb_emtrends.nopop.ss

#write.csv(seednb_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/seed_number.emtrends.CAAN2.CAIN4.csv")

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

#### Figure 3 REVISED ####

flwrdat3 = flwrdat2 %>%
  mutate(sigreg = if_else(Species %in% c("STDR", "STPO", "STGL", "STTO"), "Nonignificant", "Significant"))%>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))  
# new data to predict from
mylist <- list(transplant.daySept1= seq(min(fitdat$transplant.daySept1, na.rm=T),
                                        max(fitdat$transplant.daySept1, na.rm=T)+2, by=10), 
               Species=c("CAAN","CACO","CAIN","STBR", "STDI","STDR","STGL","STIN","STPO","STTO"))

# predictions from model
seedcount.dat=emmip(seed_nb_glob_sep1.final,Species~transplant.daySept1,at=mylist, CIs=TRUE, plotit = FALSE, type = "response") %>% # I think response will back transform from logit
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

# Same plot as original

#ggsave("./Germination.Fitness/Results/seed.number.all.pts.final.pdf", height = 10, width = 12)

#### Model of total seed mass ####
flwrdat = fitdat %>%
  filter(flowered == 1)
summary(flwrdat)

#remove STTO outlier
flwrdat2 = flwrdat %>%
  filter(long.seed.counts < 100)

seedmass_glob_sep1 = lm(long.seed.weight.g~transplant.daySept1*Species + transplant.height.cm+ Bench + Pop, data= droplevels(flwrdat2))
summary(seedmass_glob_sep1)

# testing significance of fixed effects
anova(seedmass_glob_sep1, test = "Chisq")

# estimating species-specific slopes
seedmass_emtrends = emtrends(seedmass_glob_sep1,pairwise ~ Species, var = "transplant.daySept1") 
seedmass_emtrends

# write.csv(seedmass_emtrends$emtrends, file = "./Germination.Fitness/Results/seed.mass.emtrends.csv")

#### Model of total seed mass FINAL ####

flwrdat = fitdat %>%
  filter(flowered == 1)
summary(flwrdat)

#remove STTO outlier
flwrdat2 = flwrdat %>%
  filter(long.seed.counts < 100)

seedmass_glob_sep1.final = lm(long.seed.weight.g~transplant.daySept1*Species + transplant.height.cm+ Bench, data= droplevels(flwrdat2))
summary(seedmass_glob_sep1.final)

# testing significance of fixed effects
anova(seedmass_glob_sep1.final)

# estimating species-specific slopes
seedmass_emtrends.final = emtrends(seedmass_glob_sep1.final,pairwise ~ Species, var = "transplant.daySept1") 
seedmass_emtrends.final

#write.csv(seedmass_emtrends.final$emtrends, file = "./Germination.Fitness/Results/seed.mass.emtrends.final.csv")

#### Model of total seed mass with CAAN1 and CAIN3 ####

flwrdat.sub = fitdat %>%
  filter(flowered == 1) %>%
  filter(!Pop %in% c("CAAN2","CAIN4"))
summary(flwrdat.sub)

#remove STTO outlier
flwrdat2 = flwrdat.sub %>%
  filter(long.seed.counts < 100)

seedmass_glob_sep1.nopop.ss = lm(long.seed.weight.g~transplant.daySept1*Species + transplant.height.cm+ Bench, data= droplevels(flwrdat2))
summary(seedmass_glob_sep1.nopop.ss)

# testing significance of fixed effects
anova(seedmass_glob_sep1.nopop.ss)

# estimating species-specific slopes
seedmass_emtrends.nopop.ss = emtrends(seedmass_glob_sep1.nopop.ss,pairwise ~ Species, var = "transplant.daySept1") 
seedmass_emtrends.nopop.ss

# write.csv(seedmass_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/seed.mass.emtrends.CAAN1.CAIN3.csv")

#### Model of total seed mass with CAAN1 and CAIN4 ####
flwrdat.sub = fitdat %>%
  filter(flowered == 1) %>%
  filter(!Pop %in% c("CAAN2","CAIN3"))
summary(flwrdat.sub)

#remove STTO outlier
flwrdat2 = flwrdat.sub %>%
  filter(long.seed.counts < 100)

seedmass_glob_sep1.nopop.ss = lm(long.seed.weight.g~transplant.daySept1*Species + transplant.height.cm+ Bench, data= droplevels(flwrdat2))
summary(seedmass_glob_sep1.nopop.ss)

# testing significance of fixed effects
anova(seedmass_glob_sep1.nopop.ss)

# estimating species-specific slopes
seedmass_emtrends.nopop.ss = emtrends(seedmass_glob_sep1.nopop.ss,pairwise ~ Species, var = "transplant.daySept1") 
seedmass_emtrends.nopop.ss

write.csv(seedmass_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/seed.mass.emtrends.CAAN1.CAIN4.csv")

#### Model of total seed mass with CAAN2 and CAIN3 ####
flwrdat.sub = fitdat %>%
  filter(flowered == 1) %>%
  filter(!Pop %in% c("CAAN1","CAIN4"))
summary(flwrdat.sub)

#remove STTO outlier
flwrdat2 = flwrdat.sub %>%
  filter(long.seed.counts < 100)

seedmass_glob_sep1.nopop.ss = lm(long.seed.weight.g~transplant.daySept1*Species + transplant.height.cm+ Bench, data= droplevels(flwrdat2))
summary(seedmass_glob_sep1.nopop.ss)

# testing significance of fixed effects
anova(seedmass_glob_sep1.nopop.ss)

# estimating species-specific slopes
seedmass_emtrends.nopop.ss = emtrends(seedmass_glob_sep1.nopop.ss,pairwise ~ Species, var = "transplant.daySept1") 
seedmass_emtrends.nopop.ss

# write.csv(seedmass_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/seed.mass.emtrends.CAAN2.CAIN3.csv")

#### Model of total seed mass with CAAN2 and CAIN4 ####
flwrdat.sub = fitdat %>%
  filter(flowered == 1) %>%
  filter(!Pop %in% c("CAAN1","CAIN3"))
summary(flwrdat.sub)

#remove STTO outlier
flwrdat2 = flwrdat.sub %>%
  filter(long.seed.counts < 100)

seedmass_glob_sep1.nopop.ss = lm(long.seed.weight.g~transplant.daySept1*Species + transplant.height.cm+ Bench, data= droplevels(flwrdat2))
summary(seedmass_glob_sep1.nopop.ss)

# testing significance of fixed effects
anova(seedmass_glob_sep1.nopop.ss, test = "Chisq")

# estimating species-specific slopes
seedmass_emtrends.nopop.ss = emtrends(seedmass_glob_sep1.nopop.ss,pairwise ~ Species, var = "transplant.daySept1") 
seedmass_emtrends.nopop.ss

#write.csv(seedmass_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/seed.mass.emtrends.CAAN2.CAIN4.csv")

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

#### Figure S6 REVISED ####

flwrdat3 = flwrdat2 %>%
  mutate(sigreg = if_else(Species %in% c("CAAN", "CACO", "STDI", "STIN"), "Significant", "Nonsignificant"))%>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))  

# new data to predict from
mylist <- list(transplant.daySept1= seq(min(fitdat$transplant.daySept1, na.rm=T),
                                        max(fitdat$transplant.daySept1, na.rm=T)+2, by=10), 
               Species=c("CAAN","CACO","CAIN","STBR", "STDI","STDR","STGL","STIN","STPO","STTO"))

# make predictions from model
seedmass.dat=emmip(seedmass_glob_sep1.final,Species~transplant.daySept1,at=mylist, CIs=TRUE, plotit = FALSE, type = "response") %>% # I think response will back transform from logit
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

# Same plot as original

# ggsave("./Germination.Fitness/Results/seed.mass.all.pts.final.pdf", height = 10, width = 12)

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

#### Model of first year fitness FINAL ####
# calculate first year fitness
fitdat = fitdat %>%
  mutate(year1fit = flowered * long.seed.counts)

# removed STTO outlier
fitdat2 = fitdat %>%
  filter(long.seed.counts < 100)

year1_nb_glob_sep1.final = glm.nb(year1fit~transplant.daySept1*Species + transplant.height.cm + Bench, data= droplevels(fitdat2))
summary(year1_nb_glob_sep1.final)

# testing significance of fixed effects
anova(year1_nb_glob_sep1.final, test = "Chisq")

# estimating species-specific slopes
year1_emtrends.final = emtrends(year1_nb_glob_sep1.final,pairwise ~ Species, var = "transplant.daySept1") 
year1_emtrends.final

#write.csv(year1_emtrends.final$emtrends, file = "./Germination.Fitness/Results/year1fit.emtrends.final.csv")

#### Model of first year fitness with seed source with CAAN1 and CAIN3 ####
# calculate first year fitness
fitdat = fitdat %>%
  mutate(year1fit = flowered * long.seed.counts)

# removed STTO outlier
fitdat2 = fitdat %>%
  filter(long.seed.counts < 100) %>%
  filter(!Pop %in% c("CAAN2","CAIN4"))

year1_nb_glob_sep1.nopop.ss = glm.nb(year1fit~transplant.daySept1*Species + transplant.height.cm + Bench, data= droplevels(fitdat2))
summary(year1_nb_glob_sep1.nopop.ss)

# testing significance of fixed effects
anova(year1_nb_glob_sep1.nopop.ss, test = "Chisq")

# estimating species-specific slopes
year1_emtrends.nopop.ss = emtrends(year1_nb_glob_sep1.nopop.ss,pairwise ~ Species, var = "transplant.daySept1") 
year1_emtrends.nopop.ss

#write.csv(year1_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/year1fit.emtrends.CAAN1.CAIN3.csv")

#### Model of first year fitness with CAAN1 and CAIN4 ####
# calculate first year fitness
fitdat = fitdat %>%
  mutate(year1fit = flowered * long.seed.counts)

# removed STTO outlier
fitdat2 = fitdat %>%
  filter(long.seed.counts < 100) %>%
  filter(!Pop %in% c("CAAN2","CAIN3"))

year1_nb_glob_sep1.nopop.ss = glm.nb(year1fit~transplant.daySept1*Species + transplant.height.cm + Bench, data= droplevels(fitdat2))
summary(year1_nb_glob_sep1.nopop.ss)

# testing significance of fixed effects
anova(year1_nb_glob_sep1.nopop.ss, test = "Chisq")

# estimating species-specific slopes
year1_emtrends.nopop.ss = emtrends(year1_nb_glob_sep1.nopop.ss,pairwise ~ Species, var = "transplant.daySept1") 
year1_emtrends.nopop.ss

#write.csv(year1_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/year1fit.emtrends.CAAN1.CAIN4.csv")

#### Model of first year fitness with CAAN2 and CAIN3 ####
# calculate first year fitness
fitdat = fitdat %>%
  mutate(year1fit = flowered * long.seed.counts)

# removed STTO outlier
fitdat2 = fitdat %>%
  filter(long.seed.counts < 100) %>%
  filter(!Pop %in% c("CAAN1","CAIN4"))

year1_nb_glob_sep1.nopop.ss = glm.nb(year1fit~transplant.daySept1*Species + transplant.height.cm + Bench, data= droplevels(fitdat2))
summary(year1_nb_glob_sep1.nopop.ss)

# testing significance of fixed effects
anova(year1_nb_glob_sep1.nopop.ss, test = "Chisq")

# estimating species-specific slopes
year1_emtrends.nopop.ss = emtrends(year1_nb_glob_sep1.nopop.ss,pairwise ~ Species, var = "transplant.daySept1") 
year1_emtrends.nopop.ss

write.csv(year1_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/year1fit.emtrends.CAAN2.CAIN3.csv")

#### Model of first year fitness with CAAN2 and CAIN4 ####
# calculate first year fitness
fitdat = fitdat %>%
  mutate(year1fit = flowered * long.seed.counts)

# removed STTO outlier
fitdat2 = fitdat %>%
  filter(long.seed.counts < 100) %>%
  filter(!Pop %in% c("CAAN1","CAIN3"))

year1_nb_glob_sep1.nopop.ss = glm.nb(year1fit~transplant.daySept1*Species + transplant.height.cm + Bench, data= droplevels(fitdat2))
summary(year1_nb_glob_sep1.nopop.ss)

# testing significance of fixed effects
anova(year1_nb_glob_sep1.nopop.ss, test = "Chisq")

# estimating species-specific slopes
year1_emtrends.nopop.ss = emtrends(year1_nb_glob_sep1.nopop.ss,pairwise ~ Species, var = "transplant.daySept1") 
year1_emtrends.nopop.ss

#write.csv(year1_emtrends.nopop.ss$emtrends, file = "./Germination.Fitness/Results/Pop.Sensitivity/year1fit.emtrends.CAAN2.CAIN4.csv")

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


#### Figure S7 REVISED ####

fitdat3 = fitdat2 %>%
  mutate(sigreg = if_else(Species %in% c("STDR", "STGL", "STPO"), "Nonsignificant", "Significant"))%>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))  

# make predictions
year1.dat=emmip(year1_nb_glob_sep1.final,Species~transplant.daySept1,at=mylist, CIs=TRUE, plotit = FALSE, type = "response") %>% # I think response will back transform from logit
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

# Same as original plot

#ggsave("./Germination.Fitness/Results/year1fitness.allpoints.final.pdf", height = 10, width = 12)

