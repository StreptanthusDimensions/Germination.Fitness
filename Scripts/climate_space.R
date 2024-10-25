# script to plot populations included in the study in the climate space of each species

# Users, please specify your own path to load data in this script.


# load libraries
library(tidyverse) # version 2.0.0
library(cowplot) # version 1.1.1


#### Population in Climate Space ####

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

##### Get climate variables of populations from Flint ####

flint.data = read.csv("./HTG_climate_data.csv")

flint.data.2 = flint.data %>%
  filter(id %in% c("CAAN1","CAAN2","CACO1","CAIN3","CAIN4","STBR3", "STDI","STDR2",
                   "STGL1","STIN","STPO1","STTO-TM2")) %>%
  filter(clim_year > 1990) %>% 
  filter(clim_year < 2016) %>%
  group_by(id) %>%
  dplyr::summarize(Cwd = mean(cwd),PPT_CV = (sd(ppt_mm)/mean(ppt_mm)),PPT = mean(ppt_mm), 
                   Tmin_SD=sd(tmin), Tmax_SD=sd(tmax),Tmin = mean(tmin), Tmax = mean(tmax)) %>%
  mutate(type= "seedpop")

##### temp,precip,variance months #####
climsummaries.2 = subset(climsummaries, climsummaries$folder %in% c("c_anceps","c_coulteri","c_inflatus",
                                                                    "s_breweri","s_diversifolius",
                                                                    "s_drepanoides","s_glandulosus",
                                                                    "s_insignis","s_polygaloides",
                                                                    "s_tortuosus"))

all.data = bind_rows(climsummaries.2, flint.data.2)

all.data.4pc = all.data %>%
  ungroup() %>%
  dplyr:: select(Cwd,PPT,PPT_CV,Tmin,Tmax,Tmin_SD,Tmax_SD)


all.data.pc = prcomp(all.data.4pc, scale  = TRUE, center = TRUE)
all.pc.dat = data.frame(all.data.pc$x)
all_locs_pc = cbind(all.data, all.pc.dat)
all_loadings = data.frame(varnames=rownames(all.data.pc$rotation), all.data.pc$rotation)

tibble(var_explained = (all.data.pc$sdev^2) / (sum(all.data.pc$sdev^2))) %>%
  mutate(pc = seq(1,length(var_explained), 1)) %>%
  mutate(cumvarexplained = cumsum(var_explained))
# PC 1: 52%, neg. assoc. with Cwd,PPT_CV,Tmax,Tmin, pos. assoc. with PPT
# PC 2: 22%, pos. assoc. Tmin_SD, Tmax_SD
# cumulative 74%

ggplot() +
  geom_point(data = filter(all_locs_pc, type == "herbarium"), aes(x = PC1, y = PC2, color = folder), alpha = 0.3, size = 1) +
  geom_point(data = filter(all_locs_pc, type == "seedpop"), aes(x = PC1, y = PC2, color = folder)) +
  geom_text_repel(data = filter(all_locs_pc, type == "seedpop"), aes(x = PC1, y = PC2, label = id))+
  theme_classic()+
  labs(x = paste0("Standardized PC1\n (52% explained var.)"), 
       y = paste0("Standardized PC2\n (22% explained var.)"))+
  ggtitle("cwd,ppt,temp,variance")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

all_locs_pc[c(1657:1668),10] = c("c_anceps","c_anceps",
                                 "c_coulteri","c_inflatus","c_inflatus",
                                 "s_breweri","s_diversifolius",
                                 "s_drepanoides","s_glandulosus",
                                 "s_insignis","s_polygaloides",
                                 "s_tortuosus")

range(all_locs_pc$PC1) # -4.1 to 6.9
range(all_locs_pc$PC2) # -6.7 to 3.1

caan = subset(all_locs_pc, all_locs_pc$folder == "c_anceps")
caco = subset(all_locs_pc, all_locs_pc$folder == "c_coulteri")
cain = subset(all_locs_pc, all_locs_pc$folder == "c_inflatus")
stbr = subset(all_locs_pc, all_locs_pc$folder == "s_breweri")
stdi = subset(all_locs_pc, all_locs_pc$folder == "s_diversifolius")
stdr = subset(all_locs_pc, all_locs_pc$folder == "s_drepanoides")
stgl = subset(all_locs_pc, all_locs_pc$folder == "s_glandulosus")
stin = subset(all_locs_pc, all_locs_pc$folder == "s_insignis")
stpo = subset(all_locs_pc, all_locs_pc$folder == "s_polygaloides")
stto = subset(all_locs_pc, all_locs_pc$folder == "s_tortuosus")

##### plotting #####

caan.climate = ggplot() +
  geom_point(data = filter(caan, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(caan, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4,color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("CAAN")+
  xlim(-4.1,6.9)+
  ylim(-6.7,3.1)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

caco.climate = ggplot() +
  geom_point(data = filter(caco, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(caco, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("CACO")+
  xlim(-4.1,6.9)+
  ylim(-6.7,3.1)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

cain.climate = ggplot() +
  geom_point(data = filter(cain, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(cain, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("CAIN")+
  xlim(-4.1,6.9)+
  ylim(-6.7,3.1)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

stbr.climate = ggplot() +
  geom_point(data = filter(stbr, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(stbr, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("STBR")+
  xlim(-4.1,6.9)+
  ylim(-6.7,3.1)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

stdi.climate = ggplot() +
  geom_point(data = filter(stdi, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(stdi, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("STDI")+
  xlim(-4.1,6.9)+
  ylim(-6.7,3.1)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

stdr.climate = ggplot() +
  geom_point(data = filter(stdr, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(stdr, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("STDR")+
  xlim(-4.1,6.9)+
  ylim(-6.7,3.1)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

stgl.climate = ggplot() +
  geom_point(data = filter(stgl, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(stgl, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("STGL")+
  xlim(-4.1,6.9)+
  ylim(-6.7,3.1)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

stin.climate = ggplot() +
  geom_point(data = filter(stin, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(stin, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("STIN")+
  xlim(-4.1,6.9)+
  ylim(-6.7,3.1)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

stpo.climate = ggplot() +
  geom_point(data = filter(stpo, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(stpo, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("STPO")+
  xlim(-4.1,6.9)+
  ylim(-6.7,3.1)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

stto.climate = ggplot() +
  geom_point(data = filter(stto, type == "herbarium"), aes(x = PC1, y = PC2, color = minimumElevationInMeters), size = 1) +
  scale_color_gradient(low = "black", high = "gray", "Elevation (m)")+
  geom_point(data = filter(stto, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("STTO")+
  xlim(-4.1,6.9)+
  ylim(-6.7,3.1)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")

stto.climate.legend = ggplot() +
  geom_point(data = filter(stto, type == "herbarium"), aes(x = PC1, y = PC2, color = minimumElevationInMeters), size = 1) +
  scale_color_gradient(low = "black", high = "gray", "Elevation (m)")+
  geom_point(data = filter(stto, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("STTO")+
  xlim(-4.1,6.9)+
  ylim(-6.7,3.1)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

climate.plots = plot_grid(stdr.climate,stbr.climate,stto.climate,stdi.climate,
                          stpo.climate,stin.climate,stgl.climate,caan.climate,
                          caco.climate,cain.climate)
climate.plots.legend = plot_grid(stdr.climate,stbr.climate,stto.climate.legend,stdi.climate,
                                 stpo.climate,stin.climate,stgl.climate,caan.climate,
                                 caco.climate,cain.climate)

climate.plots = plot_grid(stto.climate,stdi.climate,stpo.climate,stdr.climate,
                          stbr.climate,stin.climate,stgl.climate,caan.climate,
                          caco.climate,cain.climate)
climate.plots.legend = plot_grid(stto.climate.legend,stdi.climate,stpo.climate,stdr.climate,
                                 stbr.climate,stin.climate,stgl.climate,caan.climate,
                                 caco.climate,cain.climate)

#ggsave("Germination.Fitness/Results/climate.plots.no.legend.pdf", height = 10, width = 12)
#ggsave("Germination.Fitness/Results/climate.plots.legend.pdf", height = 10, width = 12)


