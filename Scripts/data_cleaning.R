# Script to clean Germ Fitness data

# Users, please specify your own path to load data in this script. All files are in the Raw data folder.

library(tidyverse) # version 2.0.0
library(lubridate) # version 1.9.2

# files to merge 
# 1. short season reproductive data
# 2. long season reproductive data
# 3. height at transplant data
# 4. mid-season height survey data
# 5. phenology and height survey data
# 6. biomass
# 7. fruit counts, seed counts and weights

# Short season reproductive survey
repro.survey.short = read.csv('./Germ Fitness Short Season Reproductive Success Survey Datasheets.final.csv')
colnames(repro.survey.short)[6]="short.repro.census.date"
colnames(repro.survey.short)[7]="short.flower.bud.count"
colnames(repro.survey.short)[8]="short.mature.fruit.count"
colnames(repro.survey.short)[9]="short.pedicel.count"
colnames(repro.survey.short)[10]="short.STTO.BH.height.cm"
colnames(repro.survey.short)[11]="short.repro.notes"

# Long season reproductive survey
repro.survey.long = read.csv('./Germination Fitness Long Season Reproductive Success Survey Datasheet.final.csv')
colnames(repro.survey.long)[6]="long.repro.census.date"
colnames(repro.survey.long)[7]="long.flower.bud.count"
colnames(repro.survey.long)[8]="long.mature.fruit.count"
colnames(repro.survey.long)[9]="long.pedicel.count"
colnames(repro.survey.long)[10]="long.STTO.BH.height.cm"
colnames(repro.survey.long)[11]="long.repro.notes"

repro.survey.merge = left_join(repro.survey.long, repro.survey.short, by=c("Pop","Bench",
                                                                           "Block","Position","Cohort"))
# Transplant height

transplant.height=read.csv('./Transplant.height.final.csv')

# covert height to cm to match units of other height
transplant.height$transplant.height.cm=transplant.height$Height..mm./10

colnames(transplant.height)[3]="Pop"
colnames(transplant.height)[5]="Position"
colnames(transplant.height)[7]="transplant.height.notes"

transplant.height.2=transplant.height[,c(1:3,5,8,7)]

# change bench from chr to int for merge. Makes Extra into NA
transplant.height.2$Bench=as.integer(transplant.height.2$Bench)

merge.2 = left_join(repro.survey.merge, transplant.height.2, by=c("Pop","Bench",
                                                                  "Position","Cohort"))

# mid-season height data #
mid.height=read.csv('./Midseason.height.survey.final.csv')
colnames(mid.height)[6]="Mid.Season.Height.Pheno.Stage"
colnames(mid.height)[7]="Mid.Season.Height.cm"
colnames(mid.height)[8]="Mid.Season.Height.Notes"

merge.3 = left_join(merge.2, mid.height, by = c("Pop","Bench","Position",
                                                "Cohort","Block"))

#  phenology data

phenology=read.csv('./Germination.Fitness/Raw.Data/Final.Germ Fitness 2.0 Phenology and Height Survey Datasheet.Sam.csv')
colnames(phenology)[7]="first.bud.height.cm"
colnames(phenology)[9]="first.flower.height.cm"
colnames(phenology)[11]="first.fruit.height.cm"
colnames(phenology)[13]="phenology.notes"

merge.4=left_join(merge.3, phenology, by=c("Pop","Bench","Position","Cohort","Block"))

# biomass data #

biomass=read.csv('./Germ_Fitness_round_2_biomass_data_entry_final.csv')
biomass$Envelope.1.Biomass..g.= as.numeric(biomass$Envelope.1.Biomass..g.)
colnames(biomass)[5]="biomass.AGB.envelope1.g"
colnames(biomass)[6]="biomass.AGB.envelope2.g"
colnames(biomass)[8]="biomass.BGB.STTO.BH.g"
colnames(biomass)[9]="biomass.notes"

# add biomass from 2 envelopes
biomass$total.ABG.biomass = rowSums(biomass[,c(5:6)], na.rm=TRUE)

biomass.2=biomass[,c(1:4,8,10,9)]

merge.5=left_join(merge.4, biomass.2, by=c("Pop","Bench","Position","Cohort"))

# write.csv(merge.5, file="test.csv")

# Adding cohort planting date to file

# load lookup table with dates for each cohort
# 3 dates per cohort: planting date, transplant data, randomization date
# This data is in the raw data file

planted = read.csv("./Germ Fitness 2.0 Timeline.csv")
planted = mutate(planted, plant.date = mdy(Planting.Deployment.Date))
planted = mutate(planted, transplant.date = mdy(Cone.Transplant.Date))
planted = mutate(planted, rando.date = mdy(Randomization.Date))

planted.2=planted[c(1,5:7)]

# merge germ data file and planted date file
# change upper case to lower case Cohort to match b/t dataframes
planted.2$Cohort=c("cohort 1", "cohort 2", "cohort 3", "cohort 4", "cohort 5", "cohort 6",
                   "cohort 7", "cohort 8")
merge.6 = left_join(merge.5, planted.2, by = c("Cohort"))

# make dates into same format

merge.6 = mutate(merge.6, long.repro.census.date = mdy(long.repro.census.date))
merge.6 = mutate(merge.6, short.repro.census.date = mdy(short.repro.census.date))
merge.6 = mutate(merge.6, First.Bud.Date = mdy(First.Bud.Date))
merge.6 = mutate(merge.6, First.Flower.Date = mdy(First.Flower.Date))
merge.6 = mutate(merge.6, First.Fruit.Date = mdy(First.Fruit.Date))
merge.6 = mutate(merge.6, Death.Date = mdy(Death.Date))

# seeds and fruit
seeds=read.csv('./Germination.Fitness/Raw.Data/Germ_Fitness_seeds.2.csv')
seeds$long.fruit.count.seeds=rowSums(seeds[,c("Envelope.1.Fruits", "Envelope.2.Fruits","Envelope.3.Fruits","Envelope.4.Fruits")], na.rm=TRUE)
seeds$long.seed.counts=rowSums(seeds[,c("Envelope.1.Seed.Count", "Envelope.2.Seed.Count","Envelope.3.Seed.Count","Envelope.4.Seed.Count")], na.rm=TRUE)
seeds$long.seed.weight.g=rowSums(seeds[,c("Envelope.1.Seed.Weight..g.", "Envelope.2.Seed.Weight..g.","Envelope.3.Seed.Weight..g.","Envelope.4.Seed.Weight..g.")], na.rm=TRUE)
seeds.2=seeds[,c(1:4,21:23)]

merge.7=left_join(merge.6, seeds.2, by=c("Pop","Bench","Position","Cohort"))

write.csv(merge.7, file="./final.data.temp.csv")

# Issue manually removed from final.data.temp to make final.data.csv
# remove STGL1	4	7	D08, it was pulled b/c it was a weed
# some issues with missing data, but verified with data sheets

#write.csv(merge.8, file="all.data.csv")






