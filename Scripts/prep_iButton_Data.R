# Script to prep iButton data from Germination Fitness experiment

# Users, please specify your own path to load data in this script. All files are in the Raw data folder/iButton Data

# load libraries
library("tidyverse") # version 2.0.0
library("lubridate") # version 1.9.2

# have to adjust mdy_hms vs mdy_ms depending if ibutton records seconds or not

##### reading in ibutton data for each bench and block and adding bench and block columns ####
#Bench 1
#Block 5
B1B5.1 = read.csv('./Bench 1 Block 5 (1).csv', skip = 14)
head(B1B5.1)
tail(B1B5.1)
B1B5.1$bench <- 1
B1B5.1$block <- 5
B1B5.1$Date.Time <- mdy_hms(B1B5.1$Date.Time) #change this because original format didn't have seconds and the rest of the files do,
#even though format says hm, check the values, it gives seconds (as 00)

#Bench 1
#Block 5 Part 2
B1B5.2 = read.csv("./Bench 1 Block 5 (2).csv", skip = 14)
B1B5.2$bench <- 1
B1B5.2$block <- 5
B1B5.2$Date.Time <- mdy_hm(B1B5.2$Date.Time)

#Bench 1
#Block 5 Part 3
B1B5.3 = read.csv("./Bench 1 Block 5 (3).csv", skip = 14)
B1B5.3$bench <- 1
B1B5.3$block <- 5
B1B5.3$Date.Time <- mdy_hm(B1B5.3$Date.Time)

#Bench 1
#Block 5 Part 4
B1B5.4 = read.csv("./Bench 1 Block 5 (4).csv", skip = 14)
B1B5.4$bench <- 1
B1B5.4$block <- 5
B1B5.4$Date.Time <- mdy_hm(B1B5.4$Date.Time)

#Bench 1
#Block 9
B1B9.1 = read.csv('./Bench 1 Block 9 (1).csv', skip = 14)
B1B9.1$bench <- 1
B1B9.1$block <- 9
B1B9.1$Date.Time <- mdy_hms(B1B9.1$Date.Time)

#Bench 1
#Block 9 Part 2
B1B9.2 = read.csv("./Bench 1 Block 9 (2).csv", skip = 14)
B1B9.2$bench <- 1
B1B9.2$block <- 9
B1B9.2$Date.Time <- mdy_hm(B1B9.2$Date.Time)

#Bench 1
#Block 9 Part 3
B1B9.3 = read.csv("./Bench 1 Block 9 (3).csv", skip = 14)
B1B9.3$bench <- 1
B1B9.3$block <- 9
B1B9.3$Date.Time <- mdy_hm(B1B9.3$Date.Time)

#Bench 1
#Block 13
B1B13.1 = read.csv('./Bench 1 Block 13 (1).csv', skip = 14)
B1B13.1$bench <- 1
B1B13.1$block <- 13
B1B13.1$Date.Time <- mdy_hms(B1B13.1$Date.Time)

#Bench 1
#Block 13 Part 2
B1B13.2 = read.csv("./Bench 1 Block 13 (2).csv", skip = 14)
B1B13.2$bench <- 1
B1B13.2$block <- 13
B1B13.2$Date.Time <- mdy_hm(B1B13.2$Date.Time)

#Bench 1
#Block 13 Part 3
B1B13.3 = read.csv("./Bench 1 Block 13 (3).csv", skip = 14)
B1B13.3$bench <- 1
B1B13.3$block <- 13
B1B13.3$Date.Time <- mdy_hm(B1B13.3$Date.Time)

#Bench 2
#Block 2
B2B2.1 = read.csv('./Bench 2 Block 2 (1).csv', skip = 14)
B2B2.1$bench <- 2
B2B2.1$block <- 2
B2B2.1$Date.Time <- mdy_hms(B2B2.1$Date.Time) 

#Bench 2
#Block 2 Part 2
B2B2.2 = read.csv("./Bench 2 Block 2 (2).csv", skip = 14)
B2B2.2$bench <- 2
B2B2.2$block <- 2
B2B2.2$Date.Time <- mdy_hm(B2B2.2$Date.Time)

#Bench 2
#Block 2 Part 3
B2B2.3 = read.csv("./Bench 2 Block 2 (3).csv", skip = 14)
B2B2.3$bench <- 2
B2B2.3$block <- 2
B2B2.3$Date.Time <- mdy_hm(B2B2.3$Date.Time)

#Bench 2
#Block 5
B2B5.1 = read.csv('./Bench 2 Block 5 (1).csv', skip = 14)
B2B5.1$bench <- 2
B2B5.1$block <- 5
B2B5.1$Date.Time <- mdy_hms(B2B5.1$Date.Time) 

#Bench 2
#Block 5 Part 2
B2B5.2 = read.csv("./Bench 2 Block 5 (2).csv", skip = 14)
B2B5.2$bench <- 2
B2B5.2$block <- 5
B2B5.2$Date.Time <- mdy_hm(B2B5.2$Date.Time)

#Bench 2
#Block 5 Part 3
B2B5.3 = read.csv("./Bench 2 Block 5 (3).csv", skip = 14)
B2B5.3$bench <- 2
B2B5.3$block <- 5
B2B5.3$Date.Time <- mdy_hm(B2B5.3$Date.Time)

#Bench 2
#Block 12
B2B12.1 = read.csv('./Bench 2 Block 12 (1).csv', skip = 14)
B2B12.1$bench <- 2
B2B12.1$block <- 12
B2B12.1$Date.Time <- mdy_hms(B2B12.1$Date.Time) 

#Bench 2
#Block 12 Part 2
B2B12.2 = read.csv("./Bench 2 Block 12 (2).csv", skip = 14)
B2B12.2$bench <- 2
B2B12.2$block <- 12
B2B12.2$Date.Time <- mdy_hm(B2B12.2$Date.Time)

#Bench 2
#Block 12 Part 3
B2B12.3 = read.csv("./Bench 2 Block 12 (3).csv", skip = 14)
B2B12.3$bench <- 2
B2B12.3$block <- 12
B2B12.3$Date.Time <- mdy_hm(B2B12.3$Date.Time)

#Bench 2
#Block 12 Part 4
B2B12.4 = read.csv("./Bench 2 Block 12 (4).csv", skip = 14)
B2B12.4$bench <- 2
B2B12.4$block <- 12
B2B12.4$Date.Time <- mdy_hm(B2B12.4$Date.Time)

#Bench 3
#Block 4
B3B4.1 = read.csv('./Bench 3 Block 4 (1).csv', skip = 14)
B3B4.1$bench <- 3
B3B4.1$block <- 4
B3B4.1$Date.Time <- mdy_hms(B3B4.1$Date.Time) 

#Bench 3
#Block 4 Part 2
B3B4.2 = read.csv("./Bench 3 Block 4 (2).csv", skip = 14)
B3B4.2$bench <- 3
B3B4.2$block <- 4
B3B4.2$Date.Time <- mdy_hm(B3B4.2$Date.Time)

#Bench 3
#Block 8
B3B8.1 = read.csv('./Bench 3 Block 8 (1).csv', skip = 14)
B3B8.1$bench <- 3
B3B8.1$block <- 8
B3B8.1$Date.Time <- mdy_hms(B3B8.1$Date.Time) 

#Bench 3
#Block 8 Part 2
B3B8.2 = read.csv("./Bench 3 Block 8 (2).csv", skip = 14)
B3B8.2$bench <- 3
B3B8.2$block <- 8
B3B8.2$Date.Time <- mdy_hm(B3B8.2$Date.Time)

#Bench 3
#Block 8 Part 3
B3B8.3 = read.csv("./Bench 3 Block 8 (3).csv", skip = 14)
B3B8.3$bench <- 3
B3B8.3$block <- 8
B3B8.3$Date.Time <- mdy_hm(B3B8.3$Date.Time)

#Bench 3
#Block 14
B3B14.1 = read.csv('./Bench 3 Block 14 (1).csv', skip = 14)
B3B14.1$bench <- 3
B3B14.1$block <- 14
B3B14.1$Date.Time <- mdy_hms(B3B14.1$Date.Time) 

#Bench 3
#Block 14 Part 2
B3B14.2 = read.csv("./Bench 3 Block 14 (2).csv", skip = 14)
B3B14.2$bench <- 3
B3B14.2$block <- 14
B3B14.2$Date.Time <- mdy_hm(B3B14.2$Date.Time)

#Bench 3
#Block 14 Part 3
B3B14.3 = read.csv("./Bench 3 Block 14 (3).csv", skip = 14)
B3B14.3$bench <- 3
B3B14.3$block <- 14
B3B14.3$Date.Time <- mdy_hm(B3B14.3$Date.Time)

#Bench 4
#Block 2
B4B2.1 = read.csv('./Bench 4 Block 2 (1).csv', skip = 14)
B4B2.1$bench <- 4
B4B2.1$block <- 2
B4B2.1$Date.Time <- mdy_hms(B4B2.1$Date.Time) 

#Bench 4
#Block 2 Part 2
B4B2.2 = read.csv("./Bench 4 Block 2 (2).csv", skip = 14)
B4B2.2$bench <- 4
B4B2.2$block <- 2
B4B2.2$Date.Time <- mdy_hm(B4B2.2$Date.Time)

#Bench 4
#Block 2 Part 3
B4B2.3 = read.csv("./Bench 4 Block 2 (3).csv", skip = 14)
B4B2.3$bench <- 4
B4B2.3$block <- 2
B4B2.3$Date.Time <- mdy_hm(B4B2.3$Date.Time)

#Bench 4
#Block 3
B4B3.1 = read.csv('./Bench 4 Block 3 (1).csv', skip = 14)
B4B3.1$bench <- 4
B4B3.1$block <- 3
B4B3.1$Date.Time <- mdy_hms(B4B3.1$Date.Time) 

#Bench 4
#Block 3 Part 2
B4B3.2 = read.csv("./Bench 4 Block 3 (2).csv", skip = 14)
B4B3.2$bench <- 4
B4B3.2$block <- 3
B4B3.2$Date.Time <- mdy_hm(B4B3.2$Date.Time)

#Bench 4
#Block 3 Part 3
B4B3.3 = read.csv(".//Bench 4 Block 3 (3).csv", skip = 14)
B4B3.3$bench <- 4
B4B3.3$block <- 3
B4B3.3$Date.Time <- mdy_hm(B4B3.3$Date.Time)

#Bench 4
#Block 10
B4B10.1 = read.csv('./Bench 4 Block 10 (1).csv', skip = 14)
B4B10.1$bench <- 4
B4B10.1$block <- 10
B4B10.1$Date.Time <- mdy_hms(B4B10.1$Date.Time) 

#Bench 4
#Block 10 Part 2
B4B10.2 = read.csv("./Bench 4 Block 10 (2).csv", skip = 14)
B4B10.2$bench <- 4
B4B10.2$block <- 10
B4B10.2$Date.Time <- mdy_hm(B4B10.2$Date.Time)

#Bench 4
#Block 10 Part 3
B4B10.3 = read.csv("./Bench 4 Block 10 (3).csv", skip = 14)
B4B10.3$bench <- 4
B4B10.3$block <- 10
B4B10.3$Date.Time <- mdy_hm(B4B10.3$Date.Time)

#Bench 4
#Block 10 Part 4
B4B10.4 = read.csv("./Bench 4 Block 10 (4).csv", skip = 14)
B4B10.4$bench <- 4
B4B10.4$block <- 10
B4B10.4$Date.Time <- mdy_hm(B4B10.4$Date.Time)

##### joining ibutton block datasets ####

ibuttonall=list(B1B13.1,B1B13.2,B1B13.3,B1B5.1,B1B5.2,B1B5.3,B1B5.4,B1B9.1,B1B9.2,B1B9.3,
                B2B12.1,B2B12.2,B2B12.3,B2B12.4,B2B2.1,B2B2.2,B2B2.3,B2B5.1,B2B5.2,B2B5.3,
                B3B14.1,B3B14.2,B3B14.3,B3B4.1,B3B4.2,B3B8.1,B3B8.2,B3B8.3,B4B10.1,B4B10.2,
                B4B10.3,B4B10.4,B4B2.1,B2B2.2,B4B2.3,B4B3.1,B4B3.2,B4B3.3) %>% 
  reduce(full_join)

dim(ibuttonall)
summary(ibuttonall)

#Jenny adding dplyr/tidyverse way to calculate/reformat variables

ibuttonall = ibuttonall %>%
  rename( temp = 'Value')   %>% 
  mutate(hour = hour(Date.Time)) %>% #this gives just the hour
  mutate(minute = minute(Date.Time)) %>%
  mutate(Date = as.Date(Date.Time))
summary(ibuttonall)

head(ibuttonall)
str(ibuttonall)
dim(ibuttonall)
summary(ibuttonall)

# dates covered 10-28-2021 to 01-21-2022, 01-26-2022 to 06-22-2022
# days skipped: 01-22-2022 to 01-25-2022

#### write csv files #####
#writing full dataset to a CSV
write.csv(ibuttonall, file = "Germination.Fitness/Formatted.Data/ibutton_all.csv", row.names = FALSE)

