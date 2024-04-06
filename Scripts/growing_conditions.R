# Script to evaluate growing conditions of the different cohorts
# Script produces Figure S2

# libraries
library(chillR) # to get chilling units and growing degree hours
library(lubridate)
library(tidyverse)

# approximate lat long of orchard park = 38.543267, -121.763047
# lat in degrees = 38.5

##### Getting hourly temps organized ####
# temps come from iButtons in the screenhouse

ibuttonall=read.csv("./Germination.Fitness/Formatted.Data/ibutton_all.csv")

#### Get hourly temps for missing days ####

# missing temperatures at 1/22/22 through 1/25/22
# downloaded climate data for missing dates: http://atm.ucdavis.edu/weather/uc-davis-weather-climate-station/
# use 2 m air temperature max and min = closest to matching iButton data

# also missing some of the hours for 
# 6/22/22, 10/28/2021, 1/26/22, 1/21/22, 3/13/2022, 5/5/2022

missing.temps=read.csv("./Germination.Fitness/Raw.Data/missing.temps.csv")
missing.temps$Date=mdy(missing.temps$Date)

colnames(missing.temps)[4] = "Tmax"
colnames(missing.temps)[5] = "Tmin"

missing.temps.2 = missing.temps %>%
  mutate(Year = year(Date)) %>% 
  mutate(Month = month(Date)) %>%
  mutate(Day = day(Date))

missing.temps.3 = missing.temps.2[c(1,4:8)]
hourtemps.missing=stack_hourly_temps(missing.temps.3, latitude = 38.5)$hourtemps
hourtemps.missing$DATE<-ISOdate(hourtemps.missing$Year,hourtemps.missing$Month,hourtemps.missing$Day,hourtemps.missing$Hour)

hourtemps.missing.2 = hourtemps.missing[,c(10,2:9)]

# write.csv(hourtemps.missing.2, file="./Germination.Fitness/Formatted.Data/hourtemps.missing.csv")

#### Merge ibutton data ####

# average ibutton data for each hour of each day
ibutton.2 = ibuttonall %>%
  group_by(Date,hour) %>%
  summarize(Temp = mean(temp))

# remove dates because only partial data so full data in missing.temps
# 10/28/2021, 6/22/22, 1/26/22, 1/21/22, 3/13/2022, 5/5/2022

ibutton.3 = subset(ibutton.2, !ibutton.2$Date %in% c("2021-10-28","2022-06-22","2022-01-26","2022-01-21","2022-03-13","2022-05-05"))

# get ibutton data into same format

ibutton.4 = ibutton.3 %>%
  mutate(Year = year(Date)) %>% 
  mutate(Month = month(Date)) %>%
  mutate(Day = day(Date))

ibutton.5 = ibutton.4 %>%
  group_by(Year, Month, Day) %>%
  mutate(Tmax = max(Temp),
         Tmin = min(Temp))

ibutton.5$DATE = ISOdate(ibutton.5$Year,ibutton.5$Month,ibutton.5$Day, ibutton.5$hour)
ibutton.5$JDay = yday(ibutton.5$Date)

ibutton.6=ibutton.5[,c(9,7,8,4,5,6,10,2,3)]
colnames(ibutton.6)[8]="Hour"

all.hourtemps = rbind(hourtemps.missing.2,ibutton.6)

# write.csv(all.hourtemps, file = "./Germination.Fitness/Formatted.Data/all.hourtemps.csv")

# add temperatures that go back to transplant date
# data come from iButtons in screenhouse deployed for Worthy et al. 2023 bioRxiv

ibutton.hourly.early = read.csv("./Germination.Fitness/Formatted.Data/ibutton_hourlytemps_round_2.csv")

# subset for the dates we need 10/7/2021 - 10/27/2021
ibutton.early.2 = ibutton.hourly.early[c(1312:1815),]

# average ibutton data for each hour of each day
ibutton.2 = ibutton.early.2 %>%
  group_by(Date,hour) %>%
  summarize(Temp = mean(temp))

ibutton.2$Date = mdy(ibutton.2$Date)

ibutton.3 = ibutton.2 %>%
  mutate(Year = year(Date)) %>% 
  mutate(Month = month(Date)) %>%
  mutate(Day = day(Date))

ibutton.4 = ibutton.3 %>%
  group_by(Year, Month, Day) %>%
  mutate(Tmax = max(Temp),
         Tmin = min(Temp))

ibutton.4$DATE = ISOdate(ibutton.4$Year,ibutton.4$Month,ibutton.4$Day, ibutton.4$hour)
ibutton.4$JDay = yday(ibutton.4$Date)

ibutton.5=ibutton.4[,c(9,7,8,4,5,6,10,2,3)]
colnames(ibutton.5)[8]="Hour"

ibutton.6 = as.data.frame(ibutton.5)

# merge new temps with hourtemps
hour.temps=read.csv("./Germination.Fitness/Formatted.Data/all.hourtemps.csv", row.names = 1)

new.all.hourtemps = rbind(ibutton.6,hour.temps)

# write.csv(new.all.hourtemps, file="./Germination.Fitness/Formatted.Data/transplant.all.hourtemps.csv")

#### Calculate chilling hours, Utah Chill units, Chill portions, Growing degree hours ####

# ends on 6/22/22
hour.temps=read.csv("./Germination.Fitness/Formatted.Data/transplant.all.hourtemps.csv", row.names = 1)

cohort.1.chill=chilling(hourtemps=hour.temps, Start_JDay = 280, End_JDay = 173)
cohort.2.chill=chilling(hourtemps=hour.temps, Start_JDay = 301, End_JDay = 173)
cohort.3.chill=chilling(hourtemps=hour.temps, Start_JDay = 322, End_JDay = 173)
cohort.4.chill=chilling(hourtemps=hour.temps, Start_JDay = 343, End_JDay = 173)
cohort.5.chill=chilling(hourtemps=hour.temps, Start_JDay = 364, End_JDay = 173)
cohort.6.chill=chilling(hourtemps=hour.temps, Start_JDay = 27, End_JDay = 173)
cohort.7.chill=chilling(hourtemps=hour.temps, Start_JDay = 55, End_JDay = 173)
cohort.8.chill=chilling(hourtemps=hour.temps, Start_JDay = 83, End_JDay = 173)

growing.conditions.long=rbind(cohort.1.chill,cohort.2.chill,cohort.3.chill,cohort.4.chill,
                              cohort.5.chill,cohort.6.chill,cohort.7.chill,cohort.8.chill)

# write.csv(growing.conditions.long, file="./Germination.Fitness/Formatted.Data/long.cohort.growing.conditions.csv")

# calculations for short season length
# ends on 6/1/22
cohort.1.chill.short=chilling(hourtemps=hour.temps, Start_JDay = 280, End_JDay = 152)
cohort.2.chill.short=chilling(hourtemps=hour.temps, Start_JDay = 301, End_JDay = 152)
cohort.3.chill.short=chilling(hourtemps=hour.temps, Start_JDay = 322, End_JDay = 152)
cohort.4.chill.short=chilling(hourtemps=hour.temps, Start_JDay = 343, End_JDay = 152)
cohort.5.chill.short=chilling(hourtemps=hour.temps, Start_JDay = 364, End_JDay = 152)
cohort.6.chill.short=chilling(hourtemps=hour.temps, Start_JDay = 27, End_JDay = 152)
cohort.7.chill.short=chilling(hourtemps=hour.temps, Start_JDay = 55, End_JDay = 152)
cohort.8.chill.short=chilling(hourtemps=hour.temps, Start_JDay = 83, End_JDay = 152)

growing.conditions.short=rbind(cohort.1.chill.short,cohort.2.chill.short,cohort.3.chill.short,cohort.4.chill.short,
                               cohort.5.chill.short,cohort.6.chill.short,cohort.7.chill.short,cohort.8.chill.short)

# write.csv(growing.conditions.short, file="./Germination.Fitness/Formatted.Data/short.cohort.growing.conditions.csv")

# long.cohort.growing.conditions.csv and short.cohort.growing.conditions.csv
# combined manually into cohort.growing.conditions.csv

#### Calculate heat accumulation per hour ####

hour.temps=read.csv("./Germination.Fitness/Formatted.Data/transplant.all.hourtemps.csv", row.names = 1)

all_daylengths=cbind(JDay=1:365, sapply(daylength(latitude = 38.5, JDay = 1:365),
                                        cbind))
all_daylengths.2=as.data.frame(all_daylengths)

all.data=merge(hour.temps,all_daylengths.2, by="JDay")

# need to make all values between Sunrise and Sunset = 1
# all values < sunrise = 0
# all values > sunset = 0

# add new blank column
all.data$d.hours=NA
all.data$new.hour=all.data$Hour-all.data$Sunrise
all.data$new.sunset=all.data$Hour-all.data$Sunset

heat.hour=all.data %>%
  mutate(hours.2 = ifelse(all.data$Hour > all.data$Sunset, 0, all.data$new.sunset)) %>%
  mutate(hours.3 = ifelse(all.data$Hour < all.data$Sunrise-1, 0, all.data$new.hour)) %>%
  mutate(hours.4 = ifelse(all.data$Hour < all.data$Sunset-1 & all.data$Hour > all.data$Sunrise, 1, 0))

heat.hour.2=heat.hour %>%
  mutate(new.hours.2 = ifelse(heat.hour$hours.2 < -1, "0", abs(heat.hour$hours.2))) %>%
  mutate(new.hours.3 = ifelse(heat.hour$hours.3 > 0.000005, "0", abs(heat.hour$hours.3)))

heat.hour.2$new.hours.2=as.numeric(heat.hour.2$new.hours.2)
heat.hour.2$new.hours.3=as.numeric(heat.hour.2$new.hours.3)

# subtract sunrise contribution from 1 because get the rest of the time
heat.hour.2$new.hours.4=ifelse(heat.hour.2$new.hours.3 == 0, "0", 1-heat.hour.2$new.hours.3)
heat.hour.2$new.hours.4=as.numeric(heat.hour.2$new.hours.4)

heat.hour.2$d.hours=heat.hour.2$new.hours.2+heat.hour.2$new.hours.4+heat.hour.2$hours.4

heat.hours.output=heat.hour.2[,c(1,3:13)]

# write.csv(heat.hours.output, file="./Germination.Fitness/Formatted.Data/temp.daylight.hours.csv")

#### Calculate sum of PTUs for each cohort ####
# formula: (Temperature - base) *d_range when T > base otherwise 0
# called thermal progress toward flowering in Appendix A, Burghardt et al. 2015

heat.hours.output=read.csv("./Germination.Fitness/Formatted.Data/temp.daylight.hours.csv", row.names = 1)

cohort.1.PTU = sum(ifelse(heat.hours.output[,8] - 4 > 0, (heat.hours.output[,8] - 4)*heat.hours.output$d.hours, 0))
# 50734.07, Dates:10/07/21-6/22/22
cohort.2.PTU = sum(ifelse(heat.hours.output[c(1:4152,4657:6216),8] - 4 > 0, (heat.hours.output[c(1:4152,4657:6216),8] - 4)*heat.hours.output[c(1:4152,4657:6216),12], 0)) 
# 47564.94 # Dates:10/28/21-6/22/22
cohort.3.PTU = sum(ifelse(heat.hours.output[c(1:4152,5161:6216),8] - 4 > 0, (heat.hours.output[c(1:4152,5161:6216),8] - 4)*heat.hours.output[c(1:4152,5161:6216),12], 0)) 
# 44472.18 # Dates:11/18/21-6/22/22
cohort.4.PTU = sum(ifelse(heat.hours.output[c(1:4152,5665:6216),8] - 4 > 0, (heat.hours.output[c(1:4152,5665:6216),8] - 4)*heat.hours.output[c(1:4152,5665:6216),12], 0)) 
# 42014.57 # Dates:12/09/21-6/22/22
cohort.5.PTU = sum(ifelse(heat.hours.output[c(1:4152,6169:6216),8] - 4 > 0, (heat.hours.output[c(1:4152,6169:6216),8] - 4)*heat.hours.output[c(1:4152,6169:6216),12], 0)) 
# 40708.16 # Dates:12/30/22-6/22/22
cohort.6.PTU = sum(ifelse(heat.hours.output[c(625:4152),8] - 4 > 0, (heat.hours.output[c(625:4152),8] - 4)*heat.hours.output[c(625:4152),12], 0)) 
# 37989.99 # Dates:01/27/22-6/22/22
cohort.7.PTU = sum(ifelse(heat.hours.output[c(1297:4152),8] - 4 > 0, (heat.hours.output[c(1297:4152),8] - 4)*heat.hours.output[c(1297:4152),12], 0)) 
# 33922.3 # Dates:02/24/22-6/22/22
cohort.8.PTU = sum(ifelse(heat.hours.output[c(1969:4152),8] - 4 > 0, (heat.hours.output[c(1969:4152),8] - 4)*heat.hours.output[c(1969:4152),12], 0)) 
# 28723.5 # Dates:03/24/22-6/22/22

# short season
cohort.1.PTU = sum(ifelse(heat.hours.output[c(1:3625,4153:6216),8] - 4 > 0, (heat.hours.output[c(1:3625,4153:6216),8] - 4)*heat.hours.output[c(1:3625,4153:6216),12], 0))
# 41438.95, Dates:10/07/21-6/1/22
cohort.2.PTU = sum(ifelse(heat.hours.output[c(1:3625,4657:6216),8] - 4 > 0, (heat.hours.output[c(1:3625,4657:6216),8] - 4)*heat.hours.output[c(1:3625,4657:6216),12], 0)) 
# 38269.83 # Dates:10/28/21-6/1/22
cohort.3.PTU = sum(ifelse(heat.hours.output[c(1:3625,5161:6216),8] - 4 > 0, (heat.hours.output[c(1:3625,5161:6216),8] - 4)*heat.hours.output[c(1:3625,5161:6216),12], 0)) 
# 35177.06 # Dates:11/18/21-6/1/22
cohort.4.PTU = sum(ifelse(heat.hours.output[c(1:3625,5665:6216),8] - 4 > 0, (heat.hours.output[c(1:3625,5665:6216),8] - 4)*heat.hours.output[c(1:3625,5665:6216),12], 0)) 
# 32719.46 # Dates:12/09/21-6/1/22
cohort.5.PTU = sum(ifelse(heat.hours.output[c(1:3625,6169:6216),8] - 4 > 0, (heat.hours.output[c(1:3625,6169:6216),8] - 4)*heat.hours.output[c(1:3625,6169:6216),12], 0)) 
# 31413.05 # Dates:12/30/22-6/1/22
cohort.6.PTU = sum(ifelse(heat.hours.output[c(625:3625),8] - 4 > 0, (heat.hours.output[c(625:3625),8] - 4)*heat.hours.output[c(625:3625),12], 0)) 
# 28694.88 # Dates:01/27/22-6/1/22
cohort.7.PTU = sum(ifelse(heat.hours.output[c(1297:3625),8] - 4 > 0, (heat.hours.output[c(1297:3625),8] - 4)*heat.hours.output[c(1297:3625),12], 0)) 
# 24627.18 # Dates:02/24/22-6/1/22
cohort.8.PTU = sum(ifelse(heat.hours.output[c(1969:3625),8] - 4 > 0, (heat.hours.output[c(1969:3625),8] - 4)*heat.hours.output[c(1969:3625),12], 0)) 
# 19428.39 # Dates:03/24/22-6/1/22

# values manually added to cohort.growing.conditions.csv

#### Figure S2 ####

cohort.conditions=read.csv("./Germination.Fitness/Formatted.Data/cohort.growing.conditions.csv")
long.conditions = subset(cohort.conditions, cohort.conditions$Season.Length == "long")
long.conditions$cohort.date = c("07-Oct","28-Oct","18-Nov","09-Dec",
                                "30-Dec","27-Jan","24-Feb","24-Mar")

# creating 1 plot with 2 y-axes
Figure_S2A = ggplot(long.conditions, aes(x = fct_inorder(cohort.date))) +
  geom_point(aes(y=Chill_portions, color = "Chill Portions"), size = 3) +
  geom_line(aes(y=Chill_portions, group = 1, color = "Chill Portions"), linewidth = 1.5)+
  geom_point(aes(y = PTU/1000, color = "Photothermal Units"), size = 3) + # divide to get the same range as chill data
  geom_line(aes(y= PTU/1000, group = 1,color = "Photothermal Units"), linewidth = 1.5)+
  scale_y_continuous(
    name = "Chill Portion",
    sec.axis = sec_axis(~.*1000, name = "Photothermal Units"))+
  theme_classic(base_size = 22)+
  labs(x = "Cohort", color = "")+
  scale_color_manual(values = c("black","gray70"))+
  theme(legend.position = "top")
Figure_S2A

#ggsave("./Germination.Fitness/Results/Figure_S2A.pdf", height = 6, width = 12)

all_daylengths=cbind(JDay=1:365, sapply(daylength(latitude = 38.5, JDay = 1:365),
                                        cbind))
all_daylengths.2=as.data.frame(all_daylengths)

# get rid of day not in experiment

all_daylengths.3 = all_daylengths.2[c(1:173,280:365),]

# write.csv(all_daylengths.3, file = "./Germination.Fitness/Formatted.Data/daylengths.csv")

daylengths = read.csv("./Germination.Fitness/Formatted.Data/daylengths.csv", row.names = 1)
daylengths$JDay.fact = as.factor(daylengths$JDay)

Figure_S2B = ggplot(daylengths, aes(x = fct_inorder(JDay.fact))) +
  geom_point(aes(y = Daylength))+
  theme_classic(base_size = 22)+
  scale_x_discrete(breaks = c(280,301,322,343,364,27,55,83,164),labels = c("07-Oct","28-Oct","18-Nov","09-Dec","30-Dec","27-Jan","24-Feb","24-Mar","13-Jun"))+
  labs(x = "Cohort", y = "Daylength (hours)")
Figure_S2B

#ggsave("Germination.Fitness/Results/Figure_S2B.pdf", height = 6, width = 12)





