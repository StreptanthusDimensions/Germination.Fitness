# This script analyses data to address the question: 
 # How does flowering phenology respond to germination timing?

# This script analyses the phenology metrics: days to first bud, first bud date
# script calculate dead, flowering, non-flowering proportions of individuals
# Script produces Figures 1,S3,S5

# Users, please specify your own path to load data in this script.

# load libraries
library(tidyverse) # version 2.0.0
library(lubridate) # version 1.9.2
library(emmeans) # 1.8.7

#### Calculating days to first bud, flower, fruit ####
# code includes calculation of days to first flower and fruit, but these are correlated with days to first bud

final.data=read.csv("./final.data.csv") # final.data.csv generated using data_cleaning.R script

# remove CAAM because don't actually have it in the study
# remove STTO_BH because didn't have any phenology/perennial

final.data.2=subset(final.data, !final.data$Pop %in% "CAAM")
final.data.2=subset(final.data.2, !final.data.2$Pop %in% "STTO_BH")

final.data.3 = final.data.2 %>%
  mutate(bud.date = mdy(First.Bud.Date)) %>%  # convert to date
  mutate(flower.date = mdy(First.Flower.Date)) %>%
  mutate(fruit.date = mdy(First.Fruit.Date)) %>%
  mutate(transplant.date.2 = mdy(transplant.date)) %>%
  mutate(budjul = yday(bud.date)) %>% # covert to Julian day
  mutate(flowerjul = yday(flower.date)) %>%
  mutate(fruitjul = yday(fruit.date)) %>%
  mutate(transplant.jul = yday(transplant.date.2)) %>%
  mutate(transplantjul.std = difftime(transplant.date.2, ymd("2021-09-01"))) %>% # cohort as transplants days since Sept 1
  mutate(transplantjul.std = as.numeric(transplantjul.std))

# days to event
final.data.3 = final.data.3 %>%
  mutate(day2bud = difftime(bud.date, transplant.date.2, units = c("days"))) %>%
  mutate(day2flower = difftime(flower.date, transplant.date.2, units = c("days"))) %>%
  mutate(day2fruit = difftime(fruit.date, transplant.date.2, units = c("days")))

# change to numeric
final.data.3$day2bud=as.numeric(final.data.3$day2bud)
final.data.3$day2flower=as.numeric(final.data.3$day2flower)
final.data.3$day2fruit=as.numeric(final.data.3$day2fruit)

# make names nice and put into phylogenetic order
final.data.4 = final.data.3 %>%
  mutate(Species = dplyr::recode(Species, "CAAN" = "CAAN", "CAIN" = "CAIN", "CACO1" = "CACO", 
                                 "STTO_TM2" = "STTO", "STDR2" = "STDR", "STPO1" = "STPO", "STBR3" = "STBR", 
                                 "STDI" = "STDI", "STGL1" = "STGL", "STIN" = "STIN")) %>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN")) 

# make bench and pop factors
final.data.4$bench.factor = as.factor(final.data.4$Bench)
final.data.4$pop_factor = as.factor(final.data.4$Pop)

# add seed source column to data frame and make a factor

final.data.4$seed.source = if_else(final.data.4$Pop %in% c("CAAN1","CACO1","STBR3","STDI","STDR2","STIN"), "SH", "Field")
final.data.4$seed.source = as.factor(final.data.4$seed.source)

# get columns we need for models
final.data.5 = final.data.4[,c(1:5,18,50:57)]

#### Correlation among phenophases ####

cor.test(final.data.5$day2bud,final.data.5$day2flower)
# r = 0.98, p < 2.2e-16
cor.test(final.data.5$day2bud,final.data.5$day2fruit)
# r = 0.96, p < 2.2e-16
cor.test(final.data.5$day2flower,final.data.5$day2fruit)
# r = 0.97, p < 2.2e-16

# sample sizes
# days2bud = 1008
# days2flower = 922
# days2fruit = 778

#### Model of days to first bud ####
# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench and pop included as fixed effects because only 4 benches and 2 pops, less than required for random effects

days.2.bud.model.pop = lm(day2bud~transplantjul.std*phy_order + transplant.height.cm + bench.factor + pop_factor, data=droplevels(final.data.5))

# estimate slope of relationship for each species
days.2.bud.model.pop.emm=emtrends(days.2.bud.model.pop, pairwise~phy_order, var="transplantjul.std")
# All slopes significant

# visualize slopes
emmip(days.2.bud.model.pop, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(days.2.bud.model.pop)
# bench and pop not significant

#write.csv(days.2.bud.model.pop.emm$emtrends, file="./Germination.Fitness/Results/days.2.bud.pheno.emtrends.csv")
#write.csv(days.2.bud.model.pop.emm$contrasts, file="./Germination.Fitness/Results/days.2.bud.pheno.contrasts.csv")

#### Model of days to first bud with seed source ####
# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects
# seed source included as fixed effect during revision, 2 levels Screenhouse (SH) and Field

days.2.bud.model.ss = lm(day2bud~transplantjul.std*phy_order + transplant.height.cm + bench.factor + seed.source, data=droplevels(final.data.5))
# seed.source confounded for all species except CAAN

# estimate slope of relationship for each species
days.2.bud.model.ss.emm=emtrends(days.2.bud.model.ss, pairwise~phy_order, var="transplantjul.std")

# visualize slopes
emmip(days.2.bud.model.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(days.2.bud.model.ss)

#### Model of days to first bud FINAL ####
# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

days.2.bud.model.final = lm(day2bud~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.5))

# estimate slope of relationship for each species
days.2.bud.model.emm.final=emtrends(days.2.bud.model.final, pairwise~phy_order, var="transplantjul.std")

# visualize slopes
emmip(days.2.bud.model.final, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(days.2.bud.model.final)
# bench NS, all others are

#write.csv(days.2.bud.model.emm.final$emtrends, file="./Germination.Fitness/Results/days.2.bud.pheno.emtrends.final.csv")

#### Model of days to first bud with with CAAN1 and CAIN3 ####

# data subset to only include CAAN1 and CAIN3 populations

final.data.sub = final.data.5 %>%
  filter(!Pop %in% c("CAAN2","CAIN4"))

# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

days.2.bud.model.nopop.ss = lm(day2bud~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.sub))

# estimate slope of relationship for each species
days.2.bud.model.nopop.ss.emm=emtrends(days.2.bud.model.nopop.ss, pairwise~phy_order, var="transplantjul.std")

# visualize slopes
emmip(days.2.bud.model.nopop.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(days.2.bud.model.nopop.ss)

# write.csv(days.2.bud.model.nopop.ss.emm$emtrends, file="./Germination.Fitness/Results/Pop.Sensitivity/days.2.bud.pheno.emtrends.CAAN1.CAIN3.csv")

#### Model of days to first bud with with CAAN1 and CAIN4 ####

# data subset to only include CAAN1 and CAIN4 populations

final.data.sub = final.data.5 %>%
  filter(!Pop %in% c("CAAN2","CAIN3"))

# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

days.2.bud.model.nopop.ss = lm(day2bud~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.sub))

# estimate slope of relationship for each species
days.2.bud.model.nopop.ss.emm=emtrends(days.2.bud.model.nopop.ss, pairwise~phy_order, var="transplantjul.std")

# visualize slopes
emmip(days.2.bud.model.nopop.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(days.2.bud.model.nopop.ss)

#write.csv(days.2.bud.model.nopop.ss.emm$emtrends, file="./Germination.Fitness/Results/Pop.Sensitivity/days.2.bud.pheno.emtrends.CAAN1.CAIN4.csv")

#### Model of days to first bud with with CAAN2 and CAIN3 ####

# data subset to only include CAAN2 and CAIN3 populations

final.data.sub = final.data.5 %>%
  filter(!Pop %in% c("CAAN1","CAIN4"))

# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

days.2.bud.model.nopop.ss = lm(day2bud~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.sub))

# estimate slope of relationship for each species
days.2.bud.model.nopop.ss.emm=emtrends(days.2.bud.model.nopop.ss, pairwise~phy_order, var="transplantjul.std")

# visualize slopes
emmip(days.2.bud.model.nopop.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(days.2.bud.model.nopop.ss)

#write.csv(days.2.bud.model.nopop.ss.emm$emtrends, file="./Germination.Fitness/Results/Pop.Sensitivity/days.2.bud.pheno.emtrends.CAAN2.CAIN3.csv")

#### Model of days to first bud with with CAAN2 and CAIN4 ####

# data subset to only include CAAN2 and CAIN4 populations

final.data.sub = final.data.5 %>%
  filter(!Pop %in% c("CAAN1","CAIN3"))

# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

days.2.bud.model.nopop.ss = lm(day2bud~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.sub))

# estimate slope of relationship for each species
days.2.bud.model.nopop.ss.emm=emtrends(days.2.bud.model.nopop.ss, pairwise~phy_order, var="transplantjul.std")

# visualize slopes
emmip(days.2.bud.model.nopop.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(days.2.bud.model.nopop.ss)

#write.csv(days.2.bud.model.nopop.ss.emm$emtrends, file="./Germination.Fitness/Results/Pop.Sensitivity/days.2.bud.pheno.emtrends.CAAN2.CAIN4.csv")

#### Figure 1 ####

# new data for prediction
mylist <- list(transplantjul.std=unique(final.data.5$transplantjul.std), 
               phy_order=c("CAAN","CACO","CAIN","STBR","STDI","STDR","STGL","STIN","STPO","STTO"))

# make predictions from model
bud.dat=emmip(days.2.bud.model.pop,phy_order~transplantjul.std,at=mylist, CIs=TRUE, plotit = FALSE)

bud.dat.2 = bud.dat %>%
  mutate(phy_order = fct_relevel(phy_order, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN")) 

# remove predicted values where no data actually occurs

bud.dat.2[c(60,70,80), c(3,6,7)] = NA # Remove STTO
bud.dat.2[c(69,79), c(3,6,7)] = NA # Remove STPO
bud.dat.2[c(67,77), c(3,6,7)] = NA # Remove STGL

days.2.bud.plot = ggplot(bud.dat.2, aes(x = transplantjul.std, y = yvar)) +
  geom_line(show.legend = FALSE,  lwd = 1, color = "gray40")+
  facet_wrap(.~phy_order)+
  theme_classic(base_size=22)+
  labs(y="Time to First Bud (days)",x="Cohort")+
  scale_x_continuous(breaks = c(36,57,78,99,120,148,176,204),labels = c("07-Oct","28-Oct","18-Nov","09-Dec","30-Dec","27-Jan","24-Feb","24-Mar"))+
  theme(axis.text.x = element_text(angle = 90))

days.2.bud.plot.2 = days.2.bud.plot +
  geom_point(data=final.data.5, aes(x=transplantjul.std,y=day2bud), alpha = 0.5,
             show.legend = FALSE, size = 2)+
  facet_wrap(.~phy_order)
days.2.bud.plot.2

#ggsave("Germination.Fitness/Results/bud.day2pheno.plot.all.pts.pdf", height = 10, width = 12)


#### Figure 1 REVISED ####

# new data for prediction
mylist <- list(transplantjul.std=unique(final.data.5$transplantjul.std), 
               phy_order=c("CAAN","CACO","CAIN","STBR","STDI","STDR","STGL","STIN","STPO","STTO"))

# make predictions from model
bud.dat=emmip(days.2.bud.model.final,phy_order~transplantjul.std,at=mylist, CIs=TRUE, plotit = FALSE)

bud.dat.2 = bud.dat %>%
  mutate(phy_order = fct_relevel(phy_order, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN")) 

# remove predicted values where no data actually occurs

bud.dat.2[c(60,70,80), c(3,6,7)] = NA # Remove STTO
bud.dat.2[c(69,79), c(3,6,7)] = NA # Remove STPO
bud.dat.2[c(67,77), c(3,6,7)] = NA # Remove STGL

days.2.bud.plot = ggplot(bud.dat.2, aes(x = transplantjul.std, y = yvar)) +
  geom_line(show.legend = FALSE,  lwd = 1, color = "gray40")+
  facet_wrap(.~phy_order)+
  theme_classic(base_size=22)+
  labs(y="Time to First Bud (days)",x="Cohort")+
  scale_x_continuous(breaks = c(36,57,78,99,120,148,176,204),labels = c("07-Oct","28-Oct","18-Nov","09-Dec","30-Dec","27-Jan","24-Feb","24-Mar"))+
  theme(axis.text.x = element_text(angle = 90))

days.2.bud.plot.2 = days.2.bud.plot +
  geom_point(data=final.data.5, aes(x=transplantjul.std,y=day2bud), alpha = 0.5,
             show.legend = FALSE, size = 2)+
  facet_wrap(.~phy_order)
days.2.bud.plot.2

# produces the exact same plot as the original Figure 1

#ggsave("Germination.Fitness/Results/bud.day2pheno.plot.all.pts.revised.pdf", height = 10, width = 12)

#### Model of days to first flower ####
# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench and pop included as fixed effects because only 4 benches and 2 pops, less than required for random effects

days.2.flower.model.pop = lm(day2flower~transplantjul.std*phy_order + transplant.height.cm + bench.factor + pop_factor, data=droplevels(final.data.5))

# estimate slope of relationship for each species
days.2.flower.model.pop.emm=emtrends(days.2.flower.model.pop, pairwise~phy_order, var="transplantjul.std")
# All slopes significant

# visualize slopes
emmip(days.2.flower.model.pop, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(days.2.flower.model.pop)
# bench not significant
# pop significant

# write.csv(days.2.flower.model.pop.emm$emtrends, file="./Germination.Fitness/Results/Flower/days.2.flower.pheno.emtrends.csv")

#### Model of days to first flower with CAAN1 and CAIN3 ####
# data subset to only include CAAN1 and CAIN3 populations

final.data.sub = final.data.5 %>%
  filter(!Pop %in% c("CAAN2","CAIN4"))

# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

days.2.flower.model.nopop.ss = lm(day2flower~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.sub))

# estimate slope of relationship for each species
days.2.flower.model.nopop.ss.emm=emtrends(days.2.flower.model.nopop.ss, pairwise~phy_order, var="transplantjul.std")

# visualize slopes
emmip(days.2.flower.model.nopop.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(days.2.flower.model.nopop.ss)

#write.csv(days.2.flower.model.nopop.ss.emm$emtrends, file="./Germination.Fitness/Results/Flower/days.2.flower.pheno.emtrends.CAAN1.CAIN3.csv")

#### Model of days to first flower with CAAN1 and CAIN4 ####
# data subset to only include CAAN1 and CAIN4 populations

final.data.sub = final.data.5 %>%
  filter(!Pop %in% c("CAAN2","CAIN3"))

# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

days.2.flower.model.nopop.ss = lm(day2flower~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.sub))

# estimate slope of relationship for each species
days.2.flower.model.nopop.ss.emm=emtrends(days.2.flower.model.nopop.ss, pairwise~phy_order, var="transplantjul.std")

# visualize slopes
emmip(days.2.flower.model.nopop.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(days.2.flower.model.nopop.ss)

#write.csv(days.2.flower.model.nopop.ss.emm$emtrends, file="./Germination.Fitness/Results/Flower/days.2.flower.pheno.emtrends.CAAN1.CAIN4.csv")

#### Model of days to first flower with CAAN2 and CAIN3 ####

# data subset to only include CAAN2 and CAIN3 populations

final.data.sub = final.data.5 %>%
  filter(!Pop %in% c("CAAN1","CAIN4"))

# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

days.2.flower.model.nopop.ss = lm(day2flower~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.sub))

# estimate slope of relationship for each species
days.2.flower.model.nopop.ss.emm=emtrends(days.2.flower.model.nopop.ss, pairwise~phy_order, var="transplantjul.std")

# visualize slopes
emmip(days.2.flower.model.nopop.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(days.2.flower.model.nopop.ss)


#write.csv(days.2.flower.model.nopop.ss.emm$emtrends, file="./Germination.Fitness/Results/Flower/days.2.flower.pheno.emtrends.CAAN2.CAIN3.csv")

#### Model of days to first flower with CAAN2 and CAIN4 ####

# data subset to only include CAAN2 and CAIN4 populations

final.data.sub = final.data.5 %>%
  filter(!Pop %in% c("CAAN1","CAIN3"))

# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

days.2.flower.model.nopop.ss = lm(day2flower~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.sub))

# estimate slope of relationship for each species
days.2.flower.model.nopop.ss.emm=emtrends(days.2.flower.model.nopop.ss, pairwise~phy_order, var="transplantjul.std")

# visualize slopes
emmip(days.2.flower.model.nopop.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(days.2.flower.model.nopop.ss)

#write.csv(days.2.flower.model.nopop.ss.emm$emtrends, file="./Germination.Fitness/Results/Flower/days.2.flower.pheno.emtrends.CAAN2.CAIN4.csv")

#### Calculating first bud, flower, fruit date ####
# days since Sept. 1 to first bud to standardize among cohorts
# when in the season the first bud occurs

final.data=read.csv("./final.data.csv")

# remove CAAM because don't actually have it in the study
# remove STTO_BH because didn't have any phenology/perennial

final.data.2=subset(final.data, !final.data$Pop %in% "CAAM")
final.data.2=subset(final.data.2, !final.data.2$Pop %in% "STTO_BH")

# days to phenology stage from Sept 1

final.data.3 = final.data.2 %>%
  mutate(bud.date = mdy(First.Bud.Date)) %>%  # convert to date
  mutate(flower.date = mdy(First.Flower.Date)) %>%
  mutate(fruit.date = mdy(First.Fruit.Date)) %>%
  mutate(transplant.date.2 = mdy(transplant.date)) %>%
  mutate(budjul = yday(bud.date)) %>% # covert to Julian day
  mutate(flowerjul = yday(flower.date)) %>%
  mutate(fruitjul = yday(fruit.date)) %>%
  mutate(budjul.std = difftime(bud.date, ymd("2021-09-01"))) %>% # calculates julian days since Sept. 1
  mutate(flowerjul.std = difftime(flower.date, ymd("2021-09-01"))) %>%
  mutate(fruitjul.std = difftime(fruit.date, ymd("2021-09-01"))) %>%
  mutate(transplant.jul = yday(transplant.date.2)) %>%
  mutate(transplantjul.std = difftime(transplant.date.2, ymd("2021-09-01"))) %>% # cohort as transplants days since Sept 1
  mutate(transplantjul.std = as.numeric(transplantjul.std)) %>% # change to numeric
  mutate(budjul.std = as.numeric(budjul.std)) %>% 
  mutate(flowerjul.std = as.numeric(flowerjul.std)) %>%
  mutate(fruitjul.std = as.numeric(fruitjul.std))

final.data.4 = final.data.3 %>%
  mutate(Species = dplyr::recode(Species, "CAAN" = "CAAN", "CAIN" = "CAIN", "CACO1" = "CACO", 
                                 "STTO_TM2" = "STTO", "STDR2" = "STDR", "STPO1" = "STPO", "STBR3" = "STBR", 
                                 "STDI" = "STDI", "STGL1" = "STGL", "STIN" = "STIN")) %>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN")) 

# make bench and pop factors
final.data.4$bench.factor = as.factor(final.data.4$Bench)
final.data.4$pop_factor = as.factor(final.data.4$Pop)

# add seed source column to data frame and make a factor

final.data.4$seed.source = if_else(final.data.4$Pop %in% c("CAAN1","CACO1","STBR3","STDI","STDR2","STIN"), "SH", "Field")
final.data.4$seed.source = as.factor(final.data.4$seed.source)

# get columns we need
final.data.5 = final.data.4[,c(1:5,18,49:51,53:57)]

#### Correlation among phenophases ####

cor.test(final.data.5$budjul.std,final.data.5$flowerjul.std)
# r = 0.97, p < 2.2e-16
cor.test(final.data.5$budjul.std,final.data.5$fruitjul.std)
# r = 0.93, p < 2.2e-16
cor.test(final.data.5$flowerjul.std,final.data.5$fruitjul.std)
# r = 0.95, p < 2.2e-16

# sample sizes
# budjul.std = 1008
# flowerjul.std = 922
# fruitjul.std = 778

#### Model of first bud date ####
# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench and pop included as fixed effects because only 4 benches and 2 pops, less than required for random effects

bud.sept1.model.pop = lm(budjul.std~transplantjul.std*phy_order + transplant.height.cm + bench.factor + pop_factor, data=droplevels(final.data.5))

# estimate slope of relationship for each species
bud.sept1.model.pop.emm=emtrends(bud.sept1.model.pop, pairwise~phy_order, var="transplantjul.std")
# All slopes significant, except for STTO

# visualize slopes
emmip(bud.sept1.model.pop, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(bud.sept1.model.pop)
# bench and pop not sig

#write.csv(bud.sept1.model.pop.emm$emtrends, file="./Germination.Fitness/Results/bud.sept1.pheno.emtrends.csv")
#write.csv(bud.sept1.model.pop.emm$contrasts, file="./Germination.Fitness/Results/bud.sept1.pheno.contrasts.csv")

#### Model of first bud date FINAL ####
# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

bud.sept1.model.pop.final = lm(budjul.std~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.5))

# estimate slope of relationship for each species
bud.sept1.model.pop.emm.final=emtrends(bud.sept1.model.pop.final, pairwise~phy_order, var="transplantjul.std")
# All slopes significant, except for STTO

# visualize slopes
emmip(bud.sept1.model.pop.final, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(bud.sept1.model.pop.final)
# bench not sig

#write.csv(bud.sept1.model.pop.emm.final$emtrends, file="./Germination.Fitness/Results/bud.sept1.pheno.emtrends.final.csv")

#### Model of first bud date with CAAN1 and CAIN3 ####

# data subset to only include CAAN1 and CAIN3 populations

final.data.sub = final.data.5 %>%
  filter(!Pop %in% c("CAAN2","CAIN4"))

# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

bud.sept1.model.nopop.ss = lm(budjul.std~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.sub))

# estimate slope of relationship for each species
bud.sept1.model.nopop.ss.emm=emtrends(bud.sept1.model.nopop.ss, pairwise~phy_order, var="transplantjul.std")
# All slopes significant, except for STTO

# visualize slopes
emmip(bud.sept1.model.nopop.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(bud.sept1.model.nopop.ss)
# bench and pop not sig

#write.csv(bud.sept1.model.nopop.ss.emm$emtrends, file="./Germination.Fitness/Results/Pop.Sensitivity/bud.sept1.pheno.emtrends.CAAN1.CAIN3.csv")

#### Model of first bud date with CAAN1 and CAIN4 ####

# data subset to only include CAAN1 and CAIN4 populations

final.data.sub = final.data.5 %>%
  filter(!Pop %in% c("CAAN2","CAIN3"))

# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

bud.sept1.model.nopop.ss = lm(budjul.std~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.sub))

# estimate slope of relationship for each species
bud.sept1.model.nopop.ss.emm=emtrends(bud.sept1.model.nopop.ss, pairwise~phy_order, var="transplantjul.std")
# All slopes significant, except for STTO

# visualize slopes
emmip(bud.sept1.model.nopop.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(bud.sept1.model.nopop.ss)
# bench and pop not sig

#write.csv(bud.sept1.model.nopop.ss.emm$emtrends, file="./Germination.Fitness/Results/Pop.Sensitivity/bud.sept1.pheno.emtrends.CAAN1.CAIN4.csv")

#### Model of first bud date with CAAN2 and CAIN3 ####

# data subset to only include CAAN2 and CAIN3 populations

final.data.sub = final.data.5 %>%
  filter(!Pop %in% c("CAAN1","CAIN4"))

# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

bud.sept1.model.nopop.ss = lm(budjul.std~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.sub))

# estimate slope of relationship for each species
bud.sept1.model.nopop.ss.emm=emtrends(bud.sept1.model.nopop.ss, pairwise~phy_order, var="transplantjul.std")
# All slopes significant, except for STTO

# visualize slopes
emmip(bud.sept1.model.nopop.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(bud.sept1.model.nopop.ss)
# bench and pop not sig

#write.csv(bud.sept1.model.nopop.ss.emm$emtrends, file="./Germination.Fitness/Results/Pop.Sensitivity/bud.sept1.pheno.emtrends.CAAN2.CAIN3.csv")


#### Model of first bud date with CAAN2 and CAIN4 ####

# data subset to only include CAAN2 and CAIN4 populations

final.data.sub = final.data.5 %>%
  filter(!Pop %in% c("CAAN1","CAIN3"))

# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

bud.sept1.model.nopop.ss = lm(budjul.std~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.sub))

# estimate slope of relationship for each species
bud.sept1.model.nopop.ss.emm=emtrends(bud.sept1.model.nopop.ss, pairwise~phy_order, var="transplantjul.std")
# All slopes significant, except for STTO

# visualize slopes
emmip(bud.sept1.model.nopop.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(bud.sept1.model.nopop.ss)
# bench and pop not sig

#write.csv(bud.sept1.model.nopop.ss.emm$emtrends, file="./Germination.Fitness/Results/Pop.Sensitivity/bud.sept1.pheno.emtrends.CAAN2.CAIN4.csv")

#### Figure S4 ####

# new data for prediction
mylist <- list(transplantjul.std=unique(final.data.5$transplantjul.std), 
               phy_order=c("CAAN","CACO","CAIN","STBR","STDI","STDR","STGL","STIN","STPO","STTO"))

# predict from model
bud.dat=emmip(bud.sept1.model.pop,phy_order~transplantjul.std,at=mylist, CIs=TRUE, plotit = FALSE)

bud.dat.2 = bud.dat %>%
  mutate(phy_order = fct_relevel(phy_order, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN")) 

# remove predicted values where no data actually occurs
bud.dat.2[c(60,70,80),c(3,6,7)] = NA # Remove STTO, remove all STTO b/c not significant
bud.dat.2[c(69,79), c(3,6,7)] = NA # Remove STPO
bud.dat.2[c(67,77), c(3,6,7)] = NA # Remove STGL

bud_daySep1.plot = ggplot(bud.dat.2, aes(x = transplantjul.std, y = yvar)) +
  geom_line(show.legend = FALSE,  lwd = 1, color = "gray40")+
  facet_wrap(.~phy_order)+
  geom_line(data = subset(bud.dat.2, bud.dat.2$phy_order == "STTO"),show.legend = FALSE,  lwd = 1, color = "white")+
  geom_ribbon(data = subset(bud.dat.2, bud.dat.2$phy_order == "STTO"),aes(ymin = LCL, ymax =UCL),fill = "white", colour = "white")+
  theme_classic(base_size=22)+
  labs(y="Days since Sept. 1 to First Bud",x="Cohort")+
  scale_x_continuous(breaks = c(36,57,78,99,120,148,176,204),labels = c("07-Oct","28-Oct","18-Nov","09-Dec","30-Dec","27-Jan","24-Feb","24-Mar"))+
  theme(axis.text.x = element_text(angle = 90))

bud_daySep1.plot.2 = bud_daySep1.plot +
  geom_point(data=final.data.5, aes(x=transplantjul.std,y=budjul.std), alpha = 0.5,
             show.legend = FALSE, size = 2)+
  facet_wrap(.~phy_order)
bud_daySep1.plot.2

#ggsave("Germination.Fitness/Results/bud.date.sept1.pheno.plot.all.pts.pdf", height = 10, width = 12)


#### Figure S4 REVISED ####

# new data for prediction
mylist <- list(transplantjul.std=unique(final.data.5$transplantjul.std), 
               phy_order=c("CAAN","CACO","CAIN","STBR","STDI","STDR","STGL","STIN","STPO","STTO"))

# predict from model
bud.dat=emmip(bud.sept1.model.pop.final,phy_order~transplantjul.std,at=mylist, CIs=TRUE, plotit = FALSE)

bud.dat.2 = bud.dat %>%
  mutate(phy_order = fct_relevel(phy_order, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN")) 

# remove predicted values where no data actually occurs
bud.dat.2[c(60,70,80),c(3,6,7)] = NA # Remove STTO, remove all STTO b/c not significant
bud.dat.2[c(69,79), c(3,6,7)] = NA # Remove STPO
bud.dat.2[c(67,77), c(3,6,7)] = NA # Remove STGL

bud_daySep1.plot = ggplot(bud.dat.2, aes(x = transplantjul.std, y = yvar)) +
  geom_line(show.legend = FALSE,  lwd = 1, color = "gray40")+
  facet_wrap(.~phy_order)+
  geom_line(data = subset(bud.dat.2, bud.dat.2$phy_order == "STTO"),show.legend = FALSE,  lwd = 1, color = "white")+
  geom_ribbon(data = subset(bud.dat.2, bud.dat.2$phy_order == "STTO"),aes(ymin = LCL, ymax =UCL),fill = "white", colour = "white")+
  theme_classic(base_size=22)+
  labs(y="Days since Sept. 1 to First Bud",x="Cohort")+
  scale_x_continuous(breaks = c(36,57,78,99,120,148,176,204),labels = c("07-Oct","28-Oct","18-Nov","09-Dec","30-Dec","27-Jan","24-Feb","24-Mar"))+
  theme(axis.text.x = element_text(angle = 90))

bud_daySep1.plot.2 = bud_daySep1.plot +
  geom_point(data=final.data.5, aes(x=transplantjul.std,y=budjul.std), alpha = 0.5,
             show.legend = FALSE, size = 2)+
  facet_wrap(.~phy_order)
bud_daySep1.plot.2

# plot is identical to original

#ggsave("Germination.Fitness/Results/bud.date.sept1.pheno.plot.all.pts.final.pdf", height = 10, width = 12)


#### Model of first flower date ####
# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench and pop included as fixed effects because only 4 benches and 2 pops, less than required for random effects

flower.sept1.model.pop = lm(flowerjul.std~transplantjul.std*phy_order + transplant.height.cm + bench.factor + pop_factor, data=droplevels(final.data.5))

# estimate slope of relationship for each species
flower.sept1.model.pop.emm=emtrends(flower.sept1.model.pop, pairwise~phy_order, var="transplantjul.std")
# All slopes significant, except for STTO

# visualize slopes
emmip(flower.sept1.model.pop, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(flower.sept1.model.pop)
# bench and pop not sig

# write.csv(flower.sept1.model.pop.emm$emtrends, file="./Germination.Fitness/Results/Flower/flower.sept1.pheno.emtrends.csv")

#### Model of first flower date with CAAN1 and CAIN3 ####
# data subset to only include CAAN1 and CAIN3 populations

final.data.sub = final.data.5 %>%
  filter(!Pop %in% c("CAAN2","CAIN4"))

# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

flower.sept1.model.nopop.ss = lm(flowerjul.std~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.sub))

# estimate slope of relationship for each species
flower.sept1.model.nopop.ss.emm=emtrends(flower.sept1.model.nopop.ss, pairwise~phy_order, var="transplantjul.std")
# All slopes significant, except for STTO

# visualize slopes
emmip(flower.sept1.model.nopop.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(flower.sept1.model.nopop.ss)
# bench not sig

# write.csv(flower.sept1.model.nopop.ss.emm$emtrends, file="./Germination.Fitness/Results/Flower/flower.sept1.pheno.emtrends.CAAN1.CAIN3.csv")

#### Model of first flower date with CAAN1 and CAIN4 ####
# data subset to only include CAAN1 and CAIN4 populations

final.data.sub = final.data.5 %>%
  filter(!Pop %in% c("CAAN2","CAIN3"))

# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

flower.sept1.model.nopop.ss = lm(flowerjul.std~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.sub))

# estimate slope of relationship for each species
flower.sept1.model.nopop.ss.emm=emtrends(flower.sept1.model.nopop.ss, pairwise~phy_order, var="transplantjul.std")
# All slopes significant, except for STTO

# visualize slopes
emmip(flower.sept1.model.nopop.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(flower.sept1.model.nopop.ss)
# bench not sig

# write.csv(flower.sept1.model.nopop.ss.emm$emtrends, file="./Germination.Fitness/Results/Flower/flower.sept1.pheno.emtrends.CAAN1.CAIN4.csv")

#### Model of first flower date with CAAN2 and CAIN3 ####
# data subset to only include CAAN2 and CAIN3 populations

final.data.sub = final.data.5 %>%
  filter(!Pop %in% c("CAAN1","CAIN4"))

# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

flower.sept1.model.nopop.ss = lm(flowerjul.std~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.sub))

# estimate slope of relationship for each species
flower.sept1.model.nopop.ss.emm=emtrends(flower.sept1.model.nopop.ss, pairwise~phy_order, var="transplantjul.std")
# All slopes significant, except for STTO

# visualize slopes
emmip(flower.sept1.model.nopop.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(flower.sept1.model.nopop.ss)
# bench not sig

# write.csv(flower.sept1.model.nopop.ss.emm$emtrends, file="./Germination.Fitness/Results/Flower/flower.sept1.pheno.emtrends.CAAN2.CAIN3.csv")

#### Model of first flower date with CAAN2 and CAIN4 ####
# data subset to only include CAAN2 and CAIN4 populations

final.data.sub = final.data.5 %>%
  filter(!Pop %in% c("CAAN1","CAIN3"))

# linear model
# standardized transplant date (cohort date as continuous from Sept. 1)
# phy_order (species as factor)
# bench included as fixed effects because only 4 benches, less than required for random effects

flower.sept1.model.nopop.ss = lm(flowerjul.std~transplantjul.std*phy_order + transplant.height.cm + bench.factor, data=droplevels(final.data.sub))

# estimate slope of relationship for each species
flower.sept1.model.nopop.ss.emm=emtrends(flower.sept1.model.nopop.ss, pairwise~phy_order, var="transplantjul.std")
# All slopes significant, except for STTO

# visualize slopes
emmip(flower.sept1.model.nopop.ss, phy_order ~ transplantjul.std, cov.reduce = range)

# evaluate significance of fixed effects in model
anova(flower.sept1.model.nopop.ss)
# bench not sig

#write.csv(flower.sept1.model.nopop.ss.emm$emtrends, file="./Germination.Fitness/Results/Flower/flower.sept1.pheno.emtrends.CAAN2.CAIN4.csv")

#### Calculate dead, Flowering, non-flowering proportions ####
# Proportions of individuals of each species in each cohort that fall into 
# three categories: flowered during the experiment (black), lived until the end 
# of the experiment but never flowered (dark gray), died before reaching the 
# first reproductive stage, budding (light gray)

final.data=read.csv("./final.data.csv")

# correct date for one
final.data[294,29] = "12/30/21"

# remove CAAM because don't actually have it in the study
# remove STTO_BH because didn't have any phenology/perennial

final.data.2=subset(final.data, !final.data$Pop %in% "CAAM")
final.data.2=subset(final.data.2, !final.data.2$Pop %in% "STTO_BH")

final.data.3 = final.data.2 %>%
  mutate(bud.date = mdy(First.Bud.Date)) %>%  # convert to date
  mutate(flower.date = mdy(First.Flower.Date)) %>%
  mutate(fruit.date = mdy(First.Fruit.Date)) %>%
  mutate(transplant.date.2 = mdy(transplant.date)) %>%
  mutate(death.date.2 = mdy(Death.Date)) %>%
  mutate(budjul = yday(bud.date)) %>% # covert to Julian day
  mutate(flowerjul = yday(flower.date)) %>%
  mutate(fruitjul = yday(fruit.date)) %>%
  mutate(deathjul = yday(death.date.2)) %>%
  mutate(transplant.jul = yday(transplant.date.2)) %>%
  mutate(transplantjul.std = difftime(transplant.date.2, ymd("2021-09-01"))) %>% # cohort as transplants days since Sept 1
  mutate(transplantjul.std = as.numeric(transplantjul.std))

final.data.4 = final.data.3 %>%
  mutate(Species = dplyr::recode(Species, "CAAN" = "CAAN", "CAIN" = "CAIN", "CACO1" = "CACO", 
                                 "STTO_TM2" = "STTO", "STDR2" = "STDR", "STPO1" = "STPO", "STBR3" = "STBR", 
                                 "STDI" = "STDI", "STGL1" = "STGL", "STIN" = "STIN")) %>%
  mutate(phy_order = fct_relevel(Species, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN")) 

# get total number of each species in each cohort
total.numb = final.data.4 %>%
  group_by(phy_order,Cohort) %>%
  tally()

# subset for individuals with a death date and no bud date
final.data.5 = subset(final.data.4, final.data.4$deathjul > 0)
final.data.6 = subset(final.data.5, is.na(final.data.5$budjul))

# get totals of dead individuals of each species/cohort
death.num = final.data.6 %>%
  group_by(phy_order, Cohort) %>%
  dplyr::summarize(dead.number = sum(is.na(budjul)))

# merge total with death
total.numb.2 = left_join(total.numb,death.num)

# individuals alive til end but didn't flower
final.data.7 = subset(final.data.4, is.na(final.data.4$flowerjul)) # no flower
final.data.8 = subset(final.data.7, is.na(final.data.7$deathjul)) # no death
# 291 individuals never flowered but were alive 

flower.alive = final.data.8 %>%
  group_by(phy_order, Cohort) %>%
  tally()

colnames(flower.alive)[3] = "no.flower.alive"

# merge with other data
total.numb.3 = left_join(total.numb.2,flower.alive)

# check it all worked
test = final.data.4 %>%
  filter(phy_order == "STDR",
         Cohort == "cohort 1")

total.numb.3[is.na(total.numb.3)]=0
total.numb.3$flowered = (total.numb.3$n-(total.numb.3$dead.number+total.numb.3$no.flower.alive))

# turn numbers into proportions
total.numb.3$dead.proportion = total.numb.3$dead.number/total.numb.3$n
total.numb.3$no.flower.proportion = total.numb.3$no.flower.alive/total.numb.3$n
total.numb.3$flower.proportion = total.numb.3$flowered/total.numb.3$n

# write.csv(total.numb.3, file = "./Germination.Fitness/Results/proportions.csv")

# stacked bar chart with flowering proportion on bottom
# made it long format so read back in (proportions.long.csv)

#### Figure S6 ####
# read in data
proportions = read.csv("./proportions.long.csv") # found in the results folder and generated in the above chuck of code

proportions = proportions %>%
  mutate(phy_order = fct_relevel(phy_order, "STTO", "STDI", "STPO", "STDR", "STBR", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"),
         Cohort = fct_relevel(Cohort, "7-Oct","28-Oct","18-Nov","9-Dec","30-Dec","27-Jan","24-Feb","24-Mar"))


proportions.fig=ggplot(proportions, aes(fill=factor(Status,levels = c("Flowered","No.Flower","Dead")), y=Proportions, x=factor(Cohort))) + 
  geom_bar(position=position_stack(reverse = TRUE), stat="identity")+
  theme_classic(base_size = 15)+scale_fill_manual(values = c("black","gray48","gray90"))+
  facet_wrap(~phy_order)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")+
  labs(x="Cohort", y="Proportions")
proportions.fig

#ggsave("Germination.Fitness/Results/proportions.legend.pdf", height = 10, width = 12)
#ggsave("Germination.Fitness/Results/proportions.pdf", height = 10, width = 12)


#### Table S1 Sample Sizes ####

table(final.data$Pop, final.data$Cohort)

