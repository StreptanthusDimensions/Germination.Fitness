# This script analyses data to address the question: 
# How does germination timing affect fitness, either directly or indirectly through flowering phenology?

# This script fits and plots structural equation models for time to first bud and number of seeds

# https://jslefche.github.io/sem_book/coefficients.html
# https://cran.r-project.org/web/packages/DiagrammeR/DiagrammeR.pdf

# libraries
library(lubridate)
library(tidyverse)
library(piecewiseSEM)
library(MASS)
library(semEff)
library(DiagrammeR)
library(rsvg)
library(DiagrammeRsvg)

#### Prepare the Data ####

pheno.conditions = read.csv("./Germination.Fitness/Formatted.Data/final.data.pheno.conditions.csv", row.names = 1)
cohort.conditions = read.csv("./Germination.Fitness/Formatted.Data/cohort.growing.conditions.csv")

# calculate difference in time since September-01-2021

pheno.conditions.2 = pheno.conditions[c(1:5,18,23,25,35,38,40)]

pheno.conditions.3 = pheno.conditions.2 %>%
  mutate(First.Bud.Date = mdy(First.Bud.Date)) %>%
  mutate(bud.sept1.days = difftime(First.Bud.Date, ymd("2021-09-01"))) %>%
  mutate(sept1.bud = as.numeric(bud.sept1.days)) %>%
  mutate(transplant.date = mdy(transplant.date)) %>%
  mutate(transplant.sept1.days = difftime(transplant.date, ymd("2021-09-01"))) %>%
  mutate(transplantjul.std = as.numeric(transplant.sept1.days)) %>%
  mutate(day2bud = difftime(First.Bud.Date, transplant.date, units = c("days"))) %>%
  mutate(day2bud = as.numeric(day2bud)) %>%
  mutate(First.Flower.Date = mdy(First.Flower.Date)) %>%
  mutate(day2flower = difftime(First.Flower.Date, transplant.date, units = c("days"))) %>%
  mutate(day2flower = as.numeric(day2flower))
  

# subset for long season conditions
cohort.conditions.long = subset(cohort.conditions, cohort.conditions$Season.Length == "long")

# subset pheno data for only long season
pheno.conditions.4 = subset(pheno.conditions.3, pheno.conditions.3$Season == "LONG")

# merge with the rest of pheno data
all.data = merge(pheno.conditions.4, cohort.conditions.long, by = "Cohort")

#### Subset data for each species ####

cohort.mod_new = all.data[,c(1,2,6,10,15,16)]

# removing NAs from data
table(is.na(cohort.mod_new$day2bud)) # 425
cohort.mod.new.2 = subset(cohort.mod_new, (!is.na(cohort.mod_new$day2bud)))
table(is.na(cohort.mod.new.2$transplant.height.cm)) # 3
cohort.mod.new.3 = subset(cohort.mod.new.2, (!is.na(cohort.mod.new.2$transplant.height.cm)))

caan.cohort.new = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% c("CAAN1","CAAN2"))
caco.cohort.new  = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% "CACO1")
cain.cohort.new  = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% c("CAIN3","CAIN4"))
stbr.cohort.new  = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% "STBR3")
stdi.cohort.new = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% "STDI")
stdr.cohort.new  = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% "STDR2")
stgl.cohort.new  = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% "STGL1")
stin.cohort.new  = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% "STIN")
stpo.cohort.new  = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% "STPO")
stto.cohort.new  = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% "STTO_TM2")

#### CAAN ####
caan.cohort.full.mediation = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = caan.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplant.height.cm, 
         data = caan.cohort.new))
# says we are missing the significant direct path 

summary(caan.cohort.full.mediation)
basisSet(caan.cohort.full.mediation)
fisherC(caan.cohort.full.mediation)
LLchisq(caan.cohort.full.mediation)

# add significant direct path
caan.cohort.sem.new = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = caan.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplantjul.std + transplant.height.cm, 
         data = caan.cohort.new))

# compare the models
anova(caan.cohort.full.mediation,caan.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(caan.cohort.sem.new)
fisherC(caan.cohort.sem.new)
basisSet(caan.cohort.sem.new)
rsquared(caan.cohort.sem.new)
coefs(caan.cohort.sem.new)
plot(caan.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
caan.seed.counts.glm.nb.new <- caan.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(caan.cohort.new$long.seed.counts, predict(caan.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(caan.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(caan.seed.counts.glm.nb.new)[2]*sd(caan.cohort.new$day2bud) / cohort.sd.yhat
coef(caan.seed.counts.glm.nb.new)[3]*sd(caan.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(caan.seed.counts.glm.nb.new)[4]*sd(caan.cohort.new$transplant.height.cm) / cohort.sd.yhat

# diagrammeR code test
caan.day2bud =
  "digraph SEM { 
  graph [label = 'CAAN', labelloc = 't', fontsize = 20]
  
  node [shape = rectangle]
  rec1 [label = 'Cohort']
  rec2 [label = 'Time to First bud']
  rec3 [label = 'Seed Count']
  rec4 [label = 'Transplant Height']
  
  rec1 -> rec2 [label = '-0.61'] [style = 'dashed'] [penwidth = 6]
  rec2 -> rec3 [label = '  -0.69'] [style = 'dashed'] [penwidth = 6]
  rec1 -> rec3 [label = '-0.69'] [style = 'dashed'] [penwidth = 6]
  rec4 -> rec2 [label = '-0.08'] [style = 'dashed'] [color = 'gray']
  rec4 -> rec3 [label = '  0.15']  [penwidth = 3]
  }"

grViz(caan.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("CAAN.sem.day2bud.pdf")
grViz(caan.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_png("CAAN.sem.day2bud.png")

#### CAAN1 ####

caan1.cohort.new = subset(caan.cohort.new, caan.cohort.new$Pop == "CAAN1")

caan1.cohort.full.mediation = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = caan1.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplant.height.cm, 
         data = caan1.cohort.new))
# says we are missing the significant direct path 

summary(caan1.cohort.full.mediation)
basisSet(caan1.cohort.full.mediation)
fisherC(caan1.cohort.full.mediation)
LLchisq(caan1.cohort.full.mediation)

# add significant direct path
caan1.cohort.sem.new = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = caan1.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplantjul.std + transplant.height.cm, 
         data = caan1.cohort.new))

# compare the models
anova(caan1.cohort.full.mediation,caan1.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(caan1.cohort.sem.new)
fisherC(caan1.cohort.sem.new)
basisSet(caan1.cohort.sem.new)
rsquared(caan1.cohort.sem.new)
coefs(caan1.cohort.sem.new)
plot(caan1.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
caan1.seed.counts.glm.nb.new <- caan1.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(caan1.cohort.new$long.seed.counts, predict(caan1.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(caan1.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(caan1.seed.counts.glm.nb.new)[2]*sd(caan1.cohort.new$day2bud) / cohort.sd.yhat
coef(caan1.seed.counts.glm.nb.new)[3]*sd(caan1.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(caan1.seed.counts.glm.nb.new)[4]*sd(caan1.cohort.new$transplant.height.cm) / cohort.sd.yhat



#### CAAN2 ####

caan2.cohort.new = subset(caan.cohort.new, caan.cohort.new$Pop == "CAAN2")

caan2.cohort.full.mediation = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = caan2.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplant.height.cm, 
         data = caan2.cohort.new))
# says we are missing the significant direct path 

summary(caan2.cohort.full.mediation)
basisSet(caan2.cohort.full.mediation)
fisherC(caan2.cohort.full.mediation)
LLchisq(caan2.cohort.full.mediation)

# add significant direct path
caan2.cohort.sem.new = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = caan2.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplantjul.std + transplant.height.cm, 
         data = caan2.cohort.new))

# compare the models
anova(caan2.cohort.full.mediation,caan2.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(caan2.cohort.sem.new)
fisherC(caan2.cohort.sem.new)
basisSet(caan2.cohort.sem.new)
rsquared(caan2.cohort.sem.new)
coefs(caan2.cohort.sem.new)
plot(caan2.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
caan2.seed.counts.glm.nb.new <- caan2.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(caan2.cohort.new$long.seed.counts, predict(caan2.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(caan2.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(caan2.seed.counts.glm.nb.new)[2]*sd(caan2.cohort.new$day2bud) / cohort.sd.yhat
coef(caan2.seed.counts.glm.nb.new)[3]*sd(caan2.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(caan2.seed.counts.glm.nb.new)[4]*sd(caan2.cohort.new$transplant.height.cm) / cohort.sd.yhat
#### CACO ####
caco.cohort.full.mediation = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = caco.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplant.height.cm, 
         data = caco.cohort.new))
# says we are missing the significant direct path 

summary(caco.cohort.full.mediation)
basisSet(caco.cohort.full.mediation)
fisherC(caco.cohort.full.mediation)
LLchisq(caco.cohort.full.mediation)

# add significant direct path
caco.cohort.sem.new = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = caco.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplantjul.std + transplant.height.cm, 
         data = caco.cohort.new))

# compare the models
anova(caco.cohort.full.mediation,caco.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(caco.cohort.sem.new)
fisherC(caco.cohort.sem.new)
basisSet(caco.cohort.sem.new)
rsquared(caco.cohort.sem.new)
coefs(caco.cohort.sem.new)
plot(caco.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
caco.seed.counts.glm.nb.new <- caco.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(caco.cohort.new$long.seed.counts, predict(caco.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(caco.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(caco.seed.counts.glm.nb.new)[2]*sd(caco.cohort.new$day2bud) / cohort.sd.yhat
coef(caco.seed.counts.glm.nb.new)[3]*sd(caco.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(caco.seed.counts.glm.nb.new)[4]*sd(caco.cohort.new$transplant.height.cm) / cohort.sd.yhat

# diagrammeR code test
caco.day2bud =
  "digraph SEM { 
  graph [label = 'CACO', labelloc = 't', fontsize = 20]
  
  node [shape = rectangle]
  rec1 [label = 'Cohort']
  rec2 [label = 'Time to First bud']
  rec3 [label = 'Seed Count']
  rec4 [label = 'Transplant Height']
  
  rec1 -> rec2 [label = '-0.81'] [style = 'dashed'] [penwidth = 6]
  rec2 -> rec3 [label = '  -0.22'] [style = 'dashed'] [penwidth = 6]
  rec1 -> rec3 [label = '-0.26'] [style = 'dashed'] [penwidth = 6]
  rec4 -> rec2 [label = '0.03'] [color = 'gray']
  rec4 -> rec3 [label = '  0.004']  [color = 'gray']
  }"

grViz(caco.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("CACO.sem.day2bud.pdf")
grViz(caco.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_png("CACO.sem.day2bud.png")

#### CAIN ####
cain.cohort.full.mediation = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = cain.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplant.height.cm, 
         data = cain.cohort.new))
# says we are missing the significant direct path 

summary(cain.cohort.full.mediation)
basisSet(cain.cohort.full.mediation)
fisherC(cain.cohort.full.mediation)
LLchisq(cain.cohort.full.mediation)

# add significant direct path
cain.cohort.sem.new = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = cain.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplantjul.std + transplant.height.cm, 
         data = cain.cohort.new))

# compare the models
anova(cain.cohort.full.mediation,cain.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(cain.cohort.sem.new)
fisherC(cain.cohort.sem.new)
basisSet(cain.cohort.sem.new)
rsquared(cain.cohort.sem.new)
coefs(cain.cohort.sem.new)
plot(cain.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
cain.seed.counts.glm.nb.new <- cain.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(cain.cohort.new$long.seed.counts, predict(cain.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(cain.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(cain.seed.counts.glm.nb.new)[2]*sd(cain.cohort.new$day2bud) / cohort.sd.yhat
coef(cain.seed.counts.glm.nb.new)[3]*sd(cain.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(cain.seed.counts.glm.nb.new)[4]*sd(cain.cohort.new$transplant.height.cm) / cohort.sd.yhat

# diagrammeR code test
cain.day2bud =
  "digraph SEM { 
  graph [label = 'CAIN', labelloc = 't', fontsize = 20]
  
  node [shape = rectangle]
  rec1 [label = 'Cohort']
  rec2 [label = 'Time to First bud']
  rec3 [label = 'Seed Count']
  rec4 [label = 'Transplant Height']
  
  rec1 -> rec2 [label = '-0.75'] [style = 'dashed'] [penwidth = 6]
  rec2 -> rec3 [label = '  -0.37'] [style = 'dashed'] [penwidth = 6]
  rec1 -> rec3 [label = '-0.43'] [style = 'dashed'] [penwidth = 6]
  rec4 -> rec2 [label = '-0.01'] [style = 'dashed'] [color = 'gray']
  rec4 -> rec3 [label = '  0.08'] [color = 'gray']
  }"

grViz(cain.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("CAIN.sem.day2bud.pdf")
grViz(cain.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_png("CAIN.sem.day2bud.png")

#### CAIN3 ####

cain3.cohort.new = subset(cain.cohort.new, cain.cohort.new$Pop == "CAIN3")

cain3.cohort.full.mediation = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = cain3.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplant.height.cm, 
         data = cain3.cohort.new))
# says we are missing the significant direct path 

summary(cain3.cohort.full.mediation)
basisSet(cain3.cohort.full.mediation)
fisherC(cain3.cohort.full.mediation)
LLchisq(cain3.cohort.full.mediation)

# add significant direct path
cain3.cohort.sem.new = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = cain3.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplantjul.std + transplant.height.cm, 
         data = cain3.cohort.new))

# compare the models
anova(cain3.cohort.full.mediation,cain3.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(cain3.cohort.sem.new)
fisherC(cain3.cohort.sem.new)
basisSet(cain3.cohort.sem.new)
rsquared(cain3.cohort.sem.new)
coefs(cain3.cohort.sem.new)
plot(cain3.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
cain3.seed.counts.glm.nb.new <- cain3.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(cain3.cohort.new$long.seed.counts, predict(cain3.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(cain3.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(cain3.seed.counts.glm.nb.new)[2]*sd(cain3.cohort.new$day2bud) / cohort.sd.yhat
coef(cain3.seed.counts.glm.nb.new)[3]*sd(cain3.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(cain3.seed.counts.glm.nb.new)[4]*sd(cain3.cohort.new$transplant.height.cm) / cohort.sd.yhat

#### CAIN4 ####

cain4.cohort.new = subset(cain.cohort.new, cain.cohort.new$Pop == "CAIN4")

cain4.cohort.full.mediation = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = cain4.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplant.height.cm, 
         data = cain4.cohort.new))
# says we are missing the significant direct path 

summary(cain4.cohort.full.mediation)
basisSet(cain4.cohort.full.mediation)
fisherC(cain4.cohort.full.mediation)
LLchisq(cain4.cohort.full.mediation)

# add significant direct path
cain4.cohort.sem.new = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = cain4.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplantjul.std + transplant.height.cm, 
         data = cain4.cohort.new))

# compare the models
anova(cain4.cohort.full.mediation,cain4.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(cain4.cohort.sem.new)
fisherC(cain4.cohort.sem.new)
basisSet(cain4.cohort.sem.new)
rsquared(cain4.cohort.sem.new)
coefs(cain4.cohort.sem.new)
plot(cain4.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
cain4.seed.counts.glm.nb.new <- cain4.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(cain4.cohort.new$long.seed.counts, predict(cain4.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(cain4.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(cain4.seed.counts.glm.nb.new)[2]*sd(cain4.cohort.new$day2bud) / cohort.sd.yhat
coef(cain4.seed.counts.glm.nb.new)[3]*sd(cain4.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(cain4.seed.counts.glm.nb.new)[4]*sd(cain4.cohort.new$transplant.height.cm) / cohort.sd.yhat

#### STBR ####
stbr.cohort.full.mediation = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = stbr.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplant.height.cm, 
         data = stbr.cohort.new))
# says we are missing the significant direct path 

summary(stbr.cohort.full.mediation)
basisSet(stbr.cohort.full.mediation)
fisherC(stbr.cohort.full.mediation)
LLchisq(stbr.cohort.full.mediation)

# add significant direct path
stbr.cohort.sem.new = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = stbr.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplantjul.std + transplant.height.cm, 
         data = stbr.cohort.new))

# compare the models
anova(stbr.cohort.full.mediation,stbr.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(stbr.cohort.sem.new)
fisherC(stbr.cohort.sem.new)
basisSet(stbr.cohort.sem.new)
rsquared(stbr.cohort.sem.new)
coefs(stbr.cohort.sem.new)
plot(stbr.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
stbr.seed.counts.glm.nb.new <- stbr.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(stbr.cohort.new$long.seed.counts, predict(stbr.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(stbr.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(stbr.seed.counts.glm.nb.new)[2]*sd(stbr.cohort.new$day2bud) / cohort.sd.yhat
coef(stbr.seed.counts.glm.nb.new)[3]*sd(stbr.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(stbr.seed.counts.glm.nb.new)[4]*sd(stbr.cohort.new$transplant.height.cm) / cohort.sd.yhat

# diagrammeR code test
stbr.day2bud =
  "digraph SEM { 
  graph [label = 'STBR', labelloc = 't', fontsize = 20]
  
  node [shape = rectangle]
  rec1 [label = 'Cohort']
  rec2 [label = 'Time to First bud']
  rec3 [label = 'Seed Count']
  rec4 [label = 'Transplant Height']
  
  rec1 -> rec2 [label = '-0.89'] [style = 'dashed'] [penwidth = 6]
  rec2 -> rec3 [label = '  -0.15'] [style = 'dashed'] [penwidth = 6]
  rec1 -> rec3 [label = '-0.16'] [style = 'dashed'] [penwidth = 6]
  rec4 -> rec2 [label = '-0.10'] [style = 'dashed'] [penwidth = 1]
  rec4 -> rec3 [label = '  -0.01']  [style = 'dashed'] [color = 'gray']
  }"

grViz(stbr.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("STBR.sem.day2bud.pdf")
grViz(stbr.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_png("STBR.sem.day2bud.png")

#### STDI ####
stdi.cohort.full.mediation = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = stdi.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplant.height.cm, 
         data = stdi.cohort.new))
# says we are missing the significant direct path 

summary(stdi.cohort.full.mediation)
basisSet(stdi.cohort.full.mediation)
fisherC(stdi.cohort.full.mediation)
LLchisq(stdi.cohort.full.mediation)

# add significant direct path
stdi.cohort.sem.new = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = stdi.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplantjul.std + transplant.height.cm, 
         data = stdi.cohort.new))

# compare the models
anova(stdi.cohort.full.mediation,stdi.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(stdi.cohort.sem.new)
fisherC(stdi.cohort.sem.new)
basisSet(stdi.cohort.sem.new)
rsquared(stdi.cohort.sem.new)
coefs(stdi.cohort.sem.new)
plot(stdi.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
stdi.seed.counts.glm.nb.new <- stdi.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(stdi.cohort.new$long.seed.counts, predict(stdi.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(stdi.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(stdi.seed.counts.glm.nb.new)[2]*sd(stdi.cohort.new$day2bud) / cohort.sd.yhat
coef(stdi.seed.counts.glm.nb.new)[3]*sd(stdi.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(stdi.seed.counts.glm.nb.new)[4]*sd(stdi.cohort.new$transplant.height.cm) / cohort.sd.yhat

# diagrammeR code test
stdi.day2bud =
  "digraph SEM { 
  graph [label = 'STDI', labelloc = 't', fontsize = 20]
  
  node [shape = rectangle]
  rec1 [label = 'Cohort']
  rec2 [label = 'Time of First bud']
  rec3 [label = 'Seed Count']
  rec4 [label = 'Transplant Height']
  
  rec1 -> rec2 [label = '-0.96'] [style = 'dashed'] [penwidth = 6]
  rec2 -> rec3 [label = '  -0.54'] [style = 'dashed'] [penwidth = 6]
  rec1 -> rec3 [label = '-0.88'] [style = 'dashed'] [penwidth = 6]
  rec4 -> rec2 [label = '-0.16'] [style = 'dashed'] [penwidth = 3]
  rec4 -> rec3 [label = '  -0.03'] [style = 'dashed'] [color = 'gray']
  }"

grViz(stdi.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("STDI.sem.day2bud.pdf")
grViz(stdi.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_png("STDI.sem.day2bud.png")

#### STDR ####
stdr.cohort.full.mediation = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = stdr.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplant.height.cm, 
         data = stdr.cohort.new))
# says we are missing the significant direct path 

summary(stdr.cohort.full.mediation)
basisSet(stdr.cohort.full.mediation)
fisherC(stdr.cohort.full.mediation)
LLchisq(stdr.cohort.full.mediation)

# add significant direct path
stdr.cohort.sem.new = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = stdr.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplantjul.std + transplant.height.cm, 
         data = stdr.cohort.new))

# compare the models
anova(stdr.cohort.full.mediation,stdr.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(stdr.cohort.sem.new)
fisherC(stdr.cohort.sem.new)
basisSet(stdr.cohort.sem.new)
rsquared(stdr.cohort.sem.new)
coefs(stdr.cohort.sem.new)
plot(stdr.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
stdr.seed.counts.glm.nb.new <- stdr.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(stdr.cohort.new$long.seed.counts, predict(stdr.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(stdr.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(stdr.seed.counts.glm.nb.new)[2]*sd(stdr.cohort.new$day2bud) / cohort.sd.yhat
coef(stdr.seed.counts.glm.nb.new)[3]*sd(stdr.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(stdr.seed.counts.glm.nb.new)[4]*sd(stdr.cohort.new$transplant.height.cm) / cohort.sd.yhat

# diagrammeR code test
stdr.day2bud =
  "digraph SEM { 
  graph [label = 'STDR', labelloc = 't', fontsize = 20]
  
  node [shape = rectangle]
  rec1 [label = 'Cohort']
  rec2 [label = 'Time to First bud']
  rec3 [label = 'Seed Count']
  rec4 [label = 'Transplant Height']
  
  rec1 -> rec2 [label = '-0.87'] [style = 'dashed'] [penwidth = 6]
  rec2 -> rec3 [label = '  -0.75'] [style = 'dashed'] [penwidth = 6]
  rec1 -> rec3 [label = '-0.63'] [style = 'dashed'] [penwidth = 6]
  rec4 -> rec2 [label = '-0.01'] [style = 'dashed'] [color = 'gray']
  rec4 -> rec3 [label = '  0.17']  [style = 'dashed'] [penwidth = 1]
  }"

grViz(stdr.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("STDR.sem.day2bud.pdf")
grViz(stdr.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_png("STDR.sem.day2bud.png")

#### STGL ####
stgl.cohort.full.mediation = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = stgl.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplant.height.cm, 
         data = stgl.cohort.new))

summary(stgl.cohort.full.mediation)
basisSet(stgl.cohort.full.mediation)
fisherC(stgl.cohort.full.mediation)
LLchisq(stgl.cohort.full.mediation)

# add significant direct path
stgl.cohort.sem.new = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = stgl.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplantjul.std + transplant.height.cm, 
         data = stgl.cohort.new))

# compare the models
anova(stgl.cohort.full.mediation,stgl.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(stgl.cohort.sem.new)
fisherC(stgl.cohort.sem.new)
basisSet(stgl.cohort.sem.new)
rsquared(stgl.cohort.sem.new)
coefs(stgl.cohort.sem.new)
plot(stgl.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
stgl.seed.counts.glm.nb.new <- stgl.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(stgl.cohort.new$long.seed.counts, predict(stgl.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(stgl.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(stgl.seed.counts.glm.nb.new)[2]*sd(stgl.cohort.new$day2bud) / cohort.sd.yhat
coef(stgl.seed.counts.glm.nb.new)[3]*sd(stgl.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(stgl.seed.counts.glm.nb.new)[4]*sd(stgl.cohort.new$transplant.height.cm) / cohort.sd.yhat

# diagrammeR code test
stgl.day2bud =
  "digraph SEM { 
  graph [label = 'STGL', labelloc = 't', fontsize = 20]
  
  node [shape = rectangle]
  rec1 [label = 'Cohort']
  rec2 [label = 'Time to First bud']
  rec3 [label = 'Seed Count']
  rec4 [label = 'Transplant Height']
  
  rec1 -> rec2 [label = '-0.78'] [style = 'dashed'] [penwidth = 6]
  rec2 -> rec3 [label = '  -0.08'] [style = 'dashed'] [penwidth = 6]
  rec1 -> rec3 [label = '-0.06'] [style = 'dashed'] [penwidth = 6]
  rec4 -> rec2 [label = '0.01'] [color = 'gray']
  rec4 -> rec3 [label = '  -0.02'] [style = 'dashed'] [color = 'gray']
  }"

grViz(stgl.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("STGL.sem.day2bud.pdf")
grViz(stgl.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_png("STGL.sem.day2bud.png")

#### STIN ####
stin.cohort.full.mediation = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = stin.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplant.height.cm, 
         data = stin.cohort.new))
# direct path is not significant

summary(stin.cohort.full.mediation)
basisSet(stin.cohort.full.mediation)
fisherC(stin.cohort.full.mediation)
LLchisq(stin.cohort.full.mediation)

# add significant direct path
stin.cohort.sem.new = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = stin.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplantjul.std + transplant.height.cm, 
         data = stin.cohort.new))

# compare the models
anova(stin.cohort.full.mediation,stin.cohort.sem.new)
# model without direct path is better fit

summary(stin.cohort.sem.new)
fisherC(stin.cohort.sem.new)
basisSet(stin.cohort.sem.new)
rsquared(stin.cohort.sem.new)
coefs(stin.cohort.sem.new)
plot(stin.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
stin.seed.counts.glm.nb.new <- stin.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(stin.cohort.new$long.seed.counts, predict(stin.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(stin.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(stin.seed.counts.glm.nb.new)[2]*sd(stin.cohort.new$day2bud) / cohort.sd.yhat
coef(stin.seed.counts.glm.nb.new)[3]*sd(stin.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(stin.seed.counts.glm.nb.new)[4]*sd(stin.cohort.new$transplant.height.cm) / cohort.sd.yhat

# diagrammeR code test
stin.day2bud =
  "digraph SEM { 
  graph [label = 'STIN', labelloc = 't', fontsize = 20]
  
  node [shape = rectangle]
  rec1 [label = 'Cohort']
  rec2 [label = 'Time to First bud']
  rec3 [label = 'Seed Count']
  rec4 [label = 'Transplant Height']
  
  rec1 -> rec2 [label = '-0.76'] [style = 'dashed'] [penwidth = 6]
  rec2 -> rec3 [label = '  -0.02'] [style = 'dashed'] [color = 'gray']
  rec1 -> rec3 [label = '-0.42'] [style = 'dashed'] [penwidth = 6]
  rec4 -> rec2 [label = '-0.19'] [style = 'dashed'] [penwidth = 6]
  rec4 -> rec3 [label = '  -0.03'] [style = 'dashed'] [color = 'gray']
  }"

grViz(stin.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("STIN.sem.day2bud.pdf")
grViz(stin.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_png("STIN.sem.day2bud.png")

#### STPO ####
stpo.cohort.full.mediation = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = stpo.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplant.height.cm, 
         data = stpo.cohort.new))
# direct path is not significant

summary(stpo.cohort.full.mediation)
basisSet(stpo.cohort.full.mediation)
fisherC(stpo.cohort.full.mediation)
LLchisq(stpo.cohort.full.mediation)

# add significant direct path
stpo.cohort.sem.new = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = stpo.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplantjul.std + transplant.height.cm, 
         data = stpo.cohort.new))

# compare the models
anova(stpo.cohort.full.mediation,stpo.cohort.sem.new)
# model without direct path is better fit

summary(stpo.cohort.sem.new)
fisherC(stpo.cohort.sem.new)
basisSet(stpo.cohort.sem.new)
rsquared(stpo.cohort.sem.new)
coefs(stpo.cohort.sem.new)
plot(stpo.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
stpo.seed.counts.glm.nb.new <- stpo.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(stpo.cohort.new$long.seed.counts, predict(stpo.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(stpo.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(stpo.seed.counts.glm.nb.new)[2]*sd(stpo.cohort.new$day2bud) / cohort.sd.yhat
coef(stpo.seed.counts.glm.nb.new)[3]*sd(stpo.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(stpo.seed.counts.glm.nb.new)[4]*sd(stpo.cohort.new$transplant.height.cm) / cohort.sd.yhat

# diagrammeR code test
stpo.day2bud =
  "digraph SEM { 
  graph [label = 'STPO', labelloc = 't', fontsize = 20]
  
  node [shape = rectangle]
  rec1 [label = 'Cohort']
  rec2 [label = 'Time to First bud']
  rec3 [label = 'Seed Count']
  rec4 [label = 'Transplant Height']
  
  rec1 -> rec2 [label = '-0.69'] [style = 'dashed'] [penwidth = 6]
  rec2 -> rec3 [label = '  -0.01'] [style = 'dashed'] [color = 'gray']
  rec1 -> rec3 [label = '-0.01'] [style = 'dashed'] [color = 'gray']
  rec4 -> rec2 [label = '-0.22'] [style = 'dashed'] [penwidth = 1]
  rec4 -> rec3 [label = '  0.01'] [color = 'gray']
  }"

grViz(stpo.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("STPO.sem.day2bud.pdf")
grViz(stpo.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_png("STPO.sem.day2bud.png")

#### STTO ####
stto.cohort.new = stto.cohort.new %>%
filter(long.seed.counts < 100)

stto.cohort.full.mediation = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = stto.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplant.height.cm, 
         data = stto.cohort.new))
# says we are missing the significant direct path 

summary(stto.cohort.full.mediation)
basisSet(stto.cohort.full.mediation)
fisherC(stto.cohort.full.mediation)
LLchisq(stto.cohort.full.mediation)

# add significant direct path
stto.cohort.sem.new = psem(
  lm(day2bud ~ transplantjul.std + transplant.height.cm, data = stto.cohort.new),
  glm.nb(long.seed.counts ~ day2bud + transplantjul.std + transplant.height.cm, 
         data = stto.cohort.new))

# compare the models
anova(stto.cohort.full.mediation,stto.cohort.sem.new)
# no difference in model fit

summary(stto.cohort.sem.new)
fisherC(stto.cohort.sem.new)
basisSet(stto.cohort.sem.new)
rsquared(stto.cohort.sem.new)
coefs(stto.cohort.sem.new)
plot(stto.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
stto.seed.counts.glm.nb.new <- stto.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(stto.cohort.new$long.seed.counts, predict(stto.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(stto.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(stto.seed.counts.glm.nb.new)[2]*sd(stto.cohort.new$day2bud) / cohort.sd.yhat
coef(stto.seed.counts.glm.nb.new)[3]*sd(stto.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(stto.seed.counts.glm.nb.new)[4]*sd(stto.cohort.new$transplant.height.cm) / cohort.sd.yhat

# diagrammeR code test
stto.day2bud =
  "digraph SEM { 
  graph [label = 'STTO', labelloc = 't', fontsize = 20]
  
  node [shape = rectangle]
  rec1 [label = 'Cohort']
  rec2 [label = 'Time to First bud']
  rec3 [label = 'Seed Count']
  rec4 [label = 'Transplant Height']
  
  rec1 -> rec2 [label = '-1.01'] [style = 'dashed'] [penwidth = 6]
  rec2 -> rec3 [label = '  -0.85'] [style = 'dashed'] [penwidth = 1]
  rec1 -> rec3 [label = '-0.67'] [style = 'dashed'] [color = 'gray']
  rec4 -> rec2 [label = '1.46'] [color = 'gray']
  rec4 -> rec3 [label = ' 0.17'] [color = 'gray']
  }"

grViz(stto.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("STTO.sem.day2bud.revised.pdf")
grViz(stto.day2bud) %>%
  export_svg %>% charToRaw %>% rsvg_png("STTO.sem.day2bud.revised.png")

#### Subset data for each species for days 2 flower ####

cohort.mod_new = all.data[,c(1,2,6,10,15,17)]

# removing NAs from data
table(is.na(cohort.mod_new$day2flower)) # 504
cohort.mod.new.2 = subset(cohort.mod_new, (!is.na(cohort.mod_new$day2flower)))
table(is.na(cohort.mod.new.2$transplant.height.cm)) # 3
cohort.mod.new.3 = subset(cohort.mod.new.2, (!is.na(cohort.mod.new.2$transplant.height.cm)))

caan.cohort.new = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% c("CAAN1","CAAN2"))
caco.cohort.new  = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% "CACO1")
cain.cohort.new  = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% c("CAIN3","CAIN4"))
stbr.cohort.new  = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% "STBR3")
stdi.cohort.new = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% "STDI")
stdr.cohort.new  = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% "STDR2")
stgl.cohort.new  = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% "STGL1")
stin.cohort.new  = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% "STIN")
stpo.cohort.new  = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% "STPO")
stto.cohort.new  = subset(cohort.mod.new.3, cohort.mod.new.3$Pop %in% "STTO_TM2")

#### CAAN Flower ####
caan.cohort.full.mediation = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = caan.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplant.height.cm, 
         data = caan.cohort.new))
# says we are missing the significant direct path 

summary(caan.cohort.full.mediation)
basisSet(caan.cohort.full.mediation)
fisherC(caan.cohort.full.mediation)
LLchisq(caan.cohort.full.mediation)

# add significant direct path
caan.cohort.sem.new = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = caan.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplantjul.std + transplant.height.cm, 
         data = caan.cohort.new))

# compare the models
anova(caan.cohort.full.mediation,caan.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(caan.cohort.sem.new)
fisherC(caan.cohort.sem.new)
basisSet(caan.cohort.sem.new)
rsquared(caan.cohort.sem.new)
coefs(caan.cohort.sem.new)
plot(caan.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
caan.seed.counts.glm.nb.new <- caan.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(caan.cohort.new$long.seed.counts, predict(caan.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(caan.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(caan.seed.counts.glm.nb.new)[2]*sd(caan.cohort.new$day2flower) / cohort.sd.yhat
coef(caan.seed.counts.glm.nb.new)[3]*sd(caan.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(caan.seed.counts.glm.nb.new)[4]*sd(caan.cohort.new$transplant.height.cm) / cohort.sd.yhat

#### CAAN1 Flower ####

caan1.cohort.new = subset(caan.cohort.new, caan.cohort.new$Pop == "CAAN1")

caan1.cohort.full.mediation = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = caan1.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplant.height.cm, 
         data = caan1.cohort.new))
# says we are missing the significant direct path 

summary(caan1.cohort.full.mediation)
basisSet(caan1.cohort.full.mediation)
fisherC(caan1.cohort.full.mediation)
LLchisq(caan1.cohort.full.mediation)

# add significant direct path
caan1.cohort.sem.new = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = caan1.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplantjul.std + transplant.height.cm, 
         data = caan1.cohort.new))

# compare the models
anova(caan1.cohort.full.mediation,caan1.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(caan1.cohort.sem.new)
fisherC(caan1.cohort.sem.new)
basisSet(caan1.cohort.sem.new)
rsquared(caan1.cohort.sem.new)
coefs(caan1.cohort.sem.new)
plot(caan1.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
caan1.seed.counts.glm.nb.new <- caan1.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(caan1.cohort.new$long.seed.counts, predict(caan1.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(caan1.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(caan1.seed.counts.glm.nb.new)[2]*sd(caan1.cohort.new$day2flower) / cohort.sd.yhat
coef(caan1.seed.counts.glm.nb.new)[3]*sd(caan1.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(caan1.seed.counts.glm.nb.new)[4]*sd(caan1.cohort.new$transplant.height.cm) / cohort.sd.yhat

#### CAAN2 Flower ####

caan2.cohort.new = subset(caan.cohort.new, caan.cohort.new$Pop == "CAAN2")

caan2.cohort.full.mediation = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = caan2.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplant.height.cm, 
         data = caan2.cohort.new))
# says we are missing the significant direct path 

summary(caan2.cohort.full.mediation)
basisSet(caan2.cohort.full.mediation)
fisherC(caan2.cohort.full.mediation)
LLchisq(caan2.cohort.full.mediation)

# add significant direct path
caan2.cohort.sem.new = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = caan2.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplantjul.std + transplant.height.cm, 
         data = caan2.cohort.new))

# compare the models
anova(caan2.cohort.full.mediation,caan2.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(caan2.cohort.sem.new)
fisherC(caan2.cohort.sem.new)
basisSet(caan2.cohort.sem.new)
rsquared(caan2.cohort.sem.new)
coefs(caan2.cohort.sem.new)
plot(caan2.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
caan2.seed.counts.glm.nb.new <- caan2.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(caan2.cohort.new$long.seed.counts, predict(caan2.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(caan2.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(caan2.seed.counts.glm.nb.new)[2]*sd(caan2.cohort.new$day2flower) / cohort.sd.yhat
coef(caan2.seed.counts.glm.nb.new)[3]*sd(caan2.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(caan2.seed.counts.glm.nb.new)[4]*sd(caan2.cohort.new$transplant.height.cm) / cohort.sd.yhat
#### CACO Flower ####
caco.cohort.full.mediation = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = caco.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplant.height.cm, 
         data = caco.cohort.new))
# says we are missing the significant direct path 

summary(caco.cohort.full.mediation)
basisSet(caco.cohort.full.mediation)
fisherC(caco.cohort.full.mediation)
LLchisq(caco.cohort.full.mediation)

# add significant direct path
caco.cohort.sem.new = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = caco.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplantjul.std + transplant.height.cm, 
         data = caco.cohort.new))

# compare the models
anova(caco.cohort.full.mediation,caco.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(caco.cohort.sem.new)
fisherC(caco.cohort.sem.new)
basisSet(caco.cohort.sem.new)
rsquared(caco.cohort.sem.new)
coefs(caco.cohort.sem.new)
plot(caco.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
caco.seed.counts.glm.nb.new <- caco.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(caco.cohort.new$long.seed.counts, predict(caco.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(caco.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(caco.seed.counts.glm.nb.new)[2]*sd(caco.cohort.new$day2flower) / cohort.sd.yhat
coef(caco.seed.counts.glm.nb.new)[3]*sd(caco.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(caco.seed.counts.glm.nb.new)[4]*sd(caco.cohort.new$transplant.height.cm) / cohort.sd.yhat

#### CAIN Flower ####
cain.cohort.full.mediation = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = cain.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplant.height.cm, 
         data = cain.cohort.new))
# says we are missing the significant direct path 

summary(cain.cohort.full.mediation)
basisSet(cain.cohort.full.mediation)
fisherC(cain.cohort.full.mediation)
LLchisq(cain.cohort.full.mediation)

# add significant direct path
cain.cohort.sem.new = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = cain.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplantjul.std + transplant.height.cm, 
         data = cain.cohort.new))

# compare the models
anova(cain.cohort.full.mediation,cain.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(cain.cohort.sem.new)
fisherC(cain.cohort.sem.new)
basisSet(cain.cohort.sem.new)
rsquared(cain.cohort.sem.new)
coefs(cain.cohort.sem.new)
plot(cain.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
cain.seed.counts.glm.nb.new <- cain.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(cain.cohort.new$long.seed.counts, predict(cain.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(cain.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(cain.seed.counts.glm.nb.new)[2]*sd(cain.cohort.new$day2flower) / cohort.sd.yhat
coef(cain.seed.counts.glm.nb.new)[3]*sd(cain.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(cain.seed.counts.glm.nb.new)[4]*sd(cain.cohort.new$transplant.height.cm) / cohort.sd.yhat

#### CAIN3 Flower ####

cain3.cohort.new = subset(cain.cohort.new, cain.cohort.new$Pop == "CAIN3")

cain3.cohort.full.mediation = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = cain3.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplant.height.cm, 
         data = cain3.cohort.new))
# says we are missing the significant direct path 

summary(cain3.cohort.full.mediation)
basisSet(cain3.cohort.full.mediation)
fisherC(cain3.cohort.full.mediation)
LLchisq(cain3.cohort.full.mediation)

# add significant direct path
cain3.cohort.sem.new = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = cain3.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplantjul.std + transplant.height.cm, 
         data = cain3.cohort.new))

# compare the models
anova(cain3.cohort.full.mediation,cain3.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(cain3.cohort.sem.new)
fisherC(cain3.cohort.sem.new)
basisSet(cain3.cohort.sem.new)
rsquared(cain3.cohort.sem.new)
coefs(cain3.cohort.sem.new)
plot(cain3.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
cain3.seed.counts.glm.nb.new <- cain3.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(cain3.cohort.new$long.seed.counts, predict(cain3.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(cain3.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(cain3.seed.counts.glm.nb.new)[2]*sd(cain3.cohort.new$day2flower) / cohort.sd.yhat
coef(cain3.seed.counts.glm.nb.new)[3]*sd(cain3.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(cain3.seed.counts.glm.nb.new)[4]*sd(cain3.cohort.new$transplant.height.cm) / cohort.sd.yhat

#### CAIN4 Flower ####

cain4.cohort.new = subset(cain.cohort.new, cain.cohort.new$Pop == "CAIN4")

cain4.cohort.full.mediation = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = cain4.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplant.height.cm, 
         data = cain4.cohort.new))
# says we are missing the significant direct path 

summary(cain4.cohort.full.mediation)
basisSet(cain4.cohort.full.mediation)
fisherC(cain4.cohort.full.mediation)
LLchisq(cain4.cohort.full.mediation)

# add significant direct path
cain4.cohort.sem.new = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = cain4.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplantjul.std + transplant.height.cm, 
         data = cain4.cohort.new))

# compare the models
anova(cain4.cohort.full.mediation,cain4.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(cain4.cohort.sem.new)
fisherC(cain4.cohort.sem.new)
basisSet(cain4.cohort.sem.new)
rsquared(cain4.cohort.sem.new)
coefs(cain4.cohort.sem.new)
plot(cain4.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
cain4.seed.counts.glm.nb.new <- cain4.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(cain4.cohort.new$long.seed.counts, predict(cain4.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(cain4.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(cain4.seed.counts.glm.nb.new)[2]*sd(cain4.cohort.new$day2flower) / cohort.sd.yhat
coef(cain4.seed.counts.glm.nb.new)[3]*sd(cain4.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(cain4.seed.counts.glm.nb.new)[4]*sd(cain4.cohort.new$transplant.height.cm) / cohort.sd.yhat

#### STBR Flower ####
stbr.cohort.full.mediation = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = stbr.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplant.height.cm, 
         data = stbr.cohort.new))
# says we are missing the significant direct path 

summary(stbr.cohort.full.mediation)
basisSet(stbr.cohort.full.mediation)
fisherC(stbr.cohort.full.mediation)
LLchisq(stbr.cohort.full.mediation)

# add significant direct path
stbr.cohort.sem.new = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = stbr.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplantjul.std + transplant.height.cm, 
         data = stbr.cohort.new))

# compare the models
anova(stbr.cohort.full.mediation,stbr.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(stbr.cohort.sem.new)
fisherC(stbr.cohort.sem.new)
basisSet(stbr.cohort.sem.new)
rsquared(stbr.cohort.sem.new)
coefs(stbr.cohort.sem.new)
plot(stbr.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
stbr.seed.counts.glm.nb.new <- stbr.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(stbr.cohort.new$long.seed.counts, predict(stbr.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(stbr.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(stbr.seed.counts.glm.nb.new)[2]*sd(stbr.cohort.new$day2flower) / cohort.sd.yhat
coef(stbr.seed.counts.glm.nb.new)[3]*sd(stbr.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(stbr.seed.counts.glm.nb.new)[4]*sd(stbr.cohort.new$transplant.height.cm) / cohort.sd.yhat

#### STDI Flower ####
stdi.cohort.full.mediation = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = stdi.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplant.height.cm, 
         data = stdi.cohort.new))
# says we are missing the significant direct path 

summary(stdi.cohort.full.mediation)
basisSet(stdi.cohort.full.mediation)
fisherC(stdi.cohort.full.mediation)
LLchisq(stdi.cohort.full.mediation)

# add significant direct path
stdi.cohort.sem.new = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = stdi.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplantjul.std + transplant.height.cm, 
         data = stdi.cohort.new))

# compare the models
anova(stdi.cohort.full.mediation,stdi.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(stdi.cohort.sem.new)
fisherC(stdi.cohort.sem.new)
basisSet(stdi.cohort.sem.new)
rsquared(stdi.cohort.sem.new)
coefs(stdi.cohort.sem.new)
plot(stdi.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
stdi.seed.counts.glm.nb.new <- stdi.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(stdi.cohort.new$long.seed.counts, predict(stdi.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(stdi.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(stdi.seed.counts.glm.nb.new)[2]*sd(stdi.cohort.new$day2flower) / cohort.sd.yhat
coef(stdi.seed.counts.glm.nb.new)[3]*sd(stdi.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(stdi.seed.counts.glm.nb.new)[4]*sd(stdi.cohort.new$transplant.height.cm) / cohort.sd.yhat

#### STDR Flower ####
stdr.cohort.full.mediation = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = stdr.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplant.height.cm, 
         data = stdr.cohort.new))
# says we are missing the significant direct path 

summary(stdr.cohort.full.mediation)
basisSet(stdr.cohort.full.mediation)
fisherC(stdr.cohort.full.mediation)
LLchisq(stdr.cohort.full.mediation)

# add significant direct path
stdr.cohort.sem.new = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = stdr.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplantjul.std + transplant.height.cm, 
         data = stdr.cohort.new))

# compare the models
anova(stdr.cohort.full.mediation,stdr.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(stdr.cohort.sem.new)
fisherC(stdr.cohort.sem.new)
basisSet(stdr.cohort.sem.new)
rsquared(stdr.cohort.sem.new)
coefs(stdr.cohort.sem.new)
plot(stdr.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
stdr.seed.counts.glm.nb.new <- stdr.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(stdr.cohort.new$long.seed.counts, predict(stdr.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(stdr.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(stdr.seed.counts.glm.nb.new)[2]*sd(stdr.cohort.new$day2flower) / cohort.sd.yhat
coef(stdr.seed.counts.glm.nb.new)[3]*sd(stdr.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(stdr.seed.counts.glm.nb.new)[4]*sd(stdr.cohort.new$transplant.height.cm) / cohort.sd.yhat

#### STGL Flower ####
stgl.cohort.full.mediation = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = stgl.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplant.height.cm, 
         data = stgl.cohort.new))

summary(stgl.cohort.full.mediation)
basisSet(stgl.cohort.full.mediation)
fisherC(stgl.cohort.full.mediation)
LLchisq(stgl.cohort.full.mediation)

# add significant direct path
stgl.cohort.sem.new = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = stgl.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplantjul.std + transplant.height.cm, 
         data = stgl.cohort.new))

# compare the models
anova(stgl.cohort.full.mediation,stgl.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(stgl.cohort.sem.new)
fisherC(stgl.cohort.sem.new)
basisSet(stgl.cohort.sem.new)
rsquared(stgl.cohort.sem.new)
coefs(stgl.cohort.sem.new)
plot(stgl.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
stgl.seed.counts.glm.nb.new <- stgl.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(stgl.cohort.new$long.seed.counts, predict(stgl.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(stgl.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(stgl.seed.counts.glm.nb.new)[2]*sd(stgl.cohort.new$day2flower) / cohort.sd.yhat
coef(stgl.seed.counts.glm.nb.new)[3]*sd(stgl.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(stgl.seed.counts.glm.nb.new)[4]*sd(stgl.cohort.new$transplant.height.cm) / cohort.sd.yhat

#### STIN Flower ####
stin.cohort.full.mediation = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = stin.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplant.height.cm, 
         data = stin.cohort.new))

summary(stin.cohort.full.mediation)
basisSet(stin.cohort.full.mediation)
fisherC(stin.cohort.full.mediation)
LLchisq(stin.cohort.full.mediation)

# add significant direct path
stin.cohort.sem.new = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = stin.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplantjul.std + transplant.height.cm, 
         data = stin.cohort.new))

# compare the models
anova(stin.cohort.full.mediation,stin.cohort.sem.new)
# model without direct path is better fit

summary(stin.cohort.sem.new)
fisherC(stin.cohort.sem.new)
basisSet(stin.cohort.sem.new)
rsquared(stin.cohort.sem.new)
coefs(stin.cohort.sem.new)
plot(stin.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
stin.seed.counts.glm.nb.new <- stin.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(stin.cohort.new$long.seed.counts, predict(stin.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(stin.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(stin.seed.counts.glm.nb.new)[2]*sd(stin.cohort.new$day2flower) / cohort.sd.yhat
coef(stin.seed.counts.glm.nb.new)[3]*sd(stin.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(stin.seed.counts.glm.nb.new)[4]*sd(stin.cohort.new$transplant.height.cm) / cohort.sd.yhat

#### STPO Flower ####
stpo.cohort.full.mediation = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = stpo.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplant.height.cm, 
         data = stpo.cohort.new))
# direct path is not significant

summary(stpo.cohort.full.mediation)
basisSet(stpo.cohort.full.mediation)
fisherC(stpo.cohort.full.mediation)
LLchisq(stpo.cohort.full.mediation)

# add significant direct path
stpo.cohort.sem.new = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = stpo.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplantjul.std + transplant.height.cm, 
         data = stpo.cohort.new))

# compare the models
anova(stpo.cohort.full.mediation,stpo.cohort.sem.new)
# model without direct path is better fit

summary(stpo.cohort.sem.new)
fisherC(stpo.cohort.sem.new)
basisSet(stpo.cohort.sem.new)
rsquared(stpo.cohort.sem.new)
coefs(stpo.cohort.sem.new)
plot(stpo.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
stpo.seed.counts.glm.nb.new <- stpo.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(stpo.cohort.new$long.seed.counts, predict(stpo.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(stpo.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(stpo.seed.counts.glm.nb.new)[2]*sd(stpo.cohort.new$day2flower) / cohort.sd.yhat
coef(stpo.seed.counts.glm.nb.new)[3]*sd(stpo.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(stpo.seed.counts.glm.nb.new)[4]*sd(stpo.cohort.new$transplant.height.cm) / cohort.sd.yhat

#### STTO Flower ####
stto.cohort.full.mediation = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = stto.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplant.height.cm, 
         data = stto.cohort.new))
# says we are missing the significant direct path 
# very little data to fit this, n = 19

summary(stto.cohort.full.mediation)
basisSet(stto.cohort.full.mediation)
fisherC(stto.cohort.full.mediation)
LLchisq(stto.cohort.full.mediation)

# add significant direct path
stto.cohort.sem.new = psem(
  lm(day2flower ~ transplantjul.std + transplant.height.cm, data = stto.cohort.new),
  glm.nb(long.seed.counts ~ day2flower + transplantjul.std + transplant.height.cm, 
         data = stto.cohort.new))

# compare the models
anova(stto.cohort.full.mediation,stto.cohort.sem.new)
# model with direct path is better fit
# partial mediation

summary(stto.cohort.sem.new)
fisherC(stto.cohort.sem.new)
basisSet(stto.cohort.sem.new)
rsquared(stto.cohort.sem.new)
coefs(stto.cohort.sem.new)
plot(stto.cohort.sem.new, show = "unstd", ns_dashed = TRUE)

# Calculate standardized coefficients for poisson distribution

# isolate the glm.nb
stto.seed.counts.glm.nb.new <- stto.cohort.sem.new[[2]]

# Observation Empirical 
R.cohort <- cor(stto.cohort.new$long.seed.counts, predict(stto.seed.counts.glm.nb.new, type = "response"))^2
cohort.sd.yhat <- sqrt(var(predict(stto.seed.counts.glm.nb.new, type = "link")) / R.cohort)

# Get standardized coefficient
coef(stto.seed.counts.glm.nb.new)[2]*sd(stto.cohort.new$day2flower) / cohort.sd.yhat
coef(stto.seed.counts.glm.nb.new)[3]*sd(stto.cohort.new$transplantjul.std) / cohort.sd.yhat
coef(stto.seed.counts.glm.nb.new)[4]*sd(stto.cohort.new$transplant.height.cm) / cohort.sd.yhat
