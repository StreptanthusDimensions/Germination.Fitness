## Repository for the paper "*Accelerated phenology fails to buffer fitness loss from delayed rain onset in a clade of wildflowers*"

## Abstract
The timing of early life cycle events has cascading effects on phenology and fitness. These effects may be critical for climate resilience of plant populations, especially in Mediterranean environments, where delayed rainfall onset causes delayed germination. To examine impacts of germination timing on ten species of the Streptanthus/Caulanthus clade, we induced germination across a range of dates in ambient seasonal conditions and recorded phenological and fitness traits. Later germination cohorts accelerated flowering, partially stabilizing flowering date, but the degree of this compensatory plasticity differed across species. Fitness declined with later germination; the magnitude of this decline depended on the balance between direct negative effects of later germination and compensatory positive effects of accelerated flowering. The resulting species’ differences in fitness responses suggest differential vulnerability to climate change. Species from wetter, cooler, less variable habitats exhibited greater phenological plasticity, accelerating flowering more and declining less in seed set with later germination than desert species. However, other fitness responses to germination timing, like first year fitness, were evolutionarily labile across the clade and unrelated to climate. Although compensatory phenological plasticity may buffer the impacts of delayed germination, it cannot prevent long term declines in population fitness as fall rains come later with climate change.

<a href="https://doi.org/10.5281/zenodo.14362672"><img src="https://zenodo.org/badge/782690980.svg" alt="DOI"></a>

## Authors

-   Samantha J. Worthy
-   Sarah R. Ashlock
-   Arquel Miller
-   Julin N. Maloof
-   Sharon Y. Strauss
-   Jennifer R. Gremer
-   Johanna Schmitt

### Corresponding authors contact

-   sjworthy@ucdavis.edu or sworthy2@unl.edu

### Software version and package information

This code has been verified to run on R version 4.2.3. To run the code, you will need packages: tidyverse version 2.0.0, lubridate version 1.9.2, emmeans version 1.8.7, car version 3.1.2, MASS version 7.3.58.2, chillR version 0.75, weathermetrics version 1.2.2, piecewiseSEM version 2.3.0, DiagrammeRsvg version 0.1, rsvg version 2.6.0, ape version 5.7.1, phytools version 1.9.16, ggtree version 3.6.2, colorblindr version 0.1.0, cowplot version 1.1.1, ggnewscale version 0.5.0, caper version 1.0.3, ggpubr version 0.6.0, ggrepel version 0.9.3.

### The files

Folders in this repository:

-   Raw.Data
-   Results
-   Scripts

#### Raw.Data

This folder includes alldata collected and used for analyses in this study. The iButton Data subfolder includes raw files of temperature data recorded during the experiment using Thermochron DS1921G iButtons buried in soil-filled cones.

Format is **Column Name**:Explanation

**MetaData for Germ Fitness Short Season Reproductive Success Survey Datasheets.final.csv**

-   **Pop**:which population the plant belong to
-   **Bench**:which of the 4 benches the plant for on
-   **Block**:which of the 14 blocks the plant was in
-   **Position**:which position in a block the plant was in
-   **Cohort**:which seed germination cohort the plant was part of
-   **Census.Date**:the date this data was collected
-   **Flower.Bud.Count**:total number of viable or active flowers and buds
-   **Mature.Fruit.Count**:total number of mature fruits on the plant
-   **Pedicel.Count**:total number of pedicels on the plant
-   **STTO.BH.Height..cm.**:height of STTO.BH plants in centimeters
-   **Notes..AC...Already.Collected.Plant.**:notes column were AC mean the plant was already collected.

**MetaData for Germination Fitness Long Season Reproductive Success Survey Datasheet.final.csv**

-   **Bench**:which of the 4 benches the plant for on
-   **Block**:which of the 14 blocks the plant was in
-   **Pop**:which population the plant belong to
-   **Position**:which position in a block the plant was in
-   **Cohort**:which seed germination cohort the plant was part of
-   **Census.Date**:the date this data was collected
-   **Flower.Bud.Count**:total number of viable or active flowers and buds
-   **Mature.Fruit.Count**:total number of mature fruits on the plant
-   **Pedicel.Count**:total number of pedicels on the plant
-   **STTO.BH.Height..cm.**:height of STTO.BH plants in centimeters
-   **Notes..AC...Already.Collected.Plant.**:notes column were AC mean the plant was already collected.

**MetaData for Transplant.height.final.csv**

-   **Cohort**:which seed germination cohort the plant was part of
-   **Bench**:which of the 4 benches the plant for on
-   **Population**:which population the plant belong to
-   **Block**:which of the 14 blocks the plant was in
-   **Positions**:which position in a block the plant was in
-   **Height..mm.**:height of the individual at transplant in millimeters
-   **Notes**:notes 

**MetaData for Midseason.height.survey.final.csv**

-   **Bench**:which of the 4 benches the plant for on
-   **Block**:which of the 14 blocks the plant was in
-   **Pop**:which population the plant belong to
-   **Position**:which position in a block the plant was in
-   **Cohort**:which seed germination cohort the plant was part of
-   **Pheno.Stage**:which phenological stage the individual was in at the time of the measurement
-   **Height..cm.**:height of the individual at transplant in centimeters
-   **Notes**:notes 

**MetaData for Final.Germ Fitness 2.0 Phenology and Height Survey Datasheet.Sam.csv**

-   **Bench**:which of the 4 benches the plant for on
-   **Block**:which of the 14 blocks the plant was in
-   **Pop**:which population the plant belong to
-   **Position**:which position in a block the plant was in
-   **Cohort**:which seed germination cohort the plant was part of
-   **First.Bud.Date**:date first bud appeared
-   **Height..cm....First.Bud**:height on first bud date
-   **First.Flower.Date**:date first flower appeared
-   **Height..cm....First.Flower**:height on first flower date
-   **First.Fruit.Date**:date first fruit appeared
-   **Height..cm....First.Fruit**:height on first fruit date
-   **Death.Date**:date of plant's death
-   **Notes**:notes 

**MetaData for Germ_Fitness_round_2_biomass_data_entry_final.csv**

-   **Pop**:which population the plant belong to
-   **Cohort**:which seed germination cohort the plant was part of
-   **Bench**:which of the 4 benches the plant for on
-   **Position**:which position in a block the plant was in
-   **Envelope.1.Biomass..g.**:weight of biomass in first envelope in grams
-   **Envelope.2.Biomass..g.**:weight of biomass in second envelope in grams
-   **Envelope.3.Biomass..g.**:weight of biomass in third envelope in grams
-   **Below.Ground.Biomass..g..STTO.BH**:weight of belowground biomass in grams for STTO.BH
-   **Notes**:notes 

**MetaData for Raw.Data/Germ Fitness 2.0 Timeline.csv**

-   **Cohort**:seed germination cohort
-   **Planting.Deployment.Date**:date seeds were sown for each cohort
-   **Cone.Transplant.Date**:date seedlings were transplanted into cones for each cohort
-   **Randomization.Date**:date cones were randomized on the benches for each cohort

**MetaData for Germ_Fitness_seeds.2.csv**

-   **Pop**:which population the plant belong to
-   **Cohort**:which seed germination cohort the plant was part of
-   **Bench**:which of the 4 benches the plant for on
-   **Position**:which position in a block the plant was in
-   **Envelope.1.Fruits**:number of fruits in the envelope
-   **Envelope.1.Seed.Count**:number of seeds in the envelope
-   **Envelope.1.Seed.Weight..g.**:weight of seeds in the envelope in grams
-   **Envelope.1.date**:date the fruits were collected from the plant
-   **Envelope.2.Fruits**:number of fruits in the envelope
-   **Envelope.2.Seed.Count**:number of seeds in the envelope
-   **Envelope.2.Seed.Weight..g.**:weight of seeds in the envelope in grams
-   **Envelope.2.date**:date the fruits were collected from the plant
-   **Envelope.3.Fruits**:number of fruits in the envelope
-   **Envelope.3.Seed.Count**:number of seeds in the envelope
-   **Envelope.3.Seed.Weight..g.**:weight of seeds in the envelope in grams
-   **Envelope.3.date**:date the fruits were collected from the plant
-   **Envelope.4.Fruits**:number of fruits in the envelope
-   **Envelope.4.Seed.Count**:number of seeds in the envelope
-   **Envelope.4.Seed.Weight..g.**:weight of seeds in the envelope in grams
-   **Envelope.4.date**:date the fruits were collected from the plant

**MetaData for all .csv files inside the iButton Data folder**

-   **Date.Time**:date and time of the temperature recording
-   **Unit**:unit of temperature measurement, C for Celsius 
-   **Value**:temperature value

**MetaData for missing.temps.csv**

-   **Date**:Date of temperature value
-   **temp.max.F**:maximum temperature in Fahrenheit 
-   **temp.min.F**:minimum temperature in Fahrenheit 
-   **temp.max.C**:maximum temperature in Celsius 
-   **temp.min.C**:minimum temperature in Celsius 

**MetaData for ibutton_hourlytemps_round_2.csv**

-   **Date.Time**:date and time of the temperature recording
-   **Unit**:unit of temperature measurement, C for Celsius 
-   **temp**:temperature value
-   **block**:block the ibutton was located in 
-   **hour**:hour the temperature was recorded
-   **minute**:minute the temperature was recorded
-   **Date**:date the temperature was recorded

**MetaData for cohort.growing.conditions.csv**

-   **Cohort**:Seed germination cohorts
-   **Season**:Years of dates included in calculations
-   **End_year**:the year when the period ended
-   **Season_days**:duration of the period
-   **Data_days**:days in the period with data
-   **Perc_complete**:percent of days with data
-   **Chilling_Hours**:chilling hours
-   **Utah_Model**:utah model
-   **Chill_portions**:chill portions, the dynamic model
-   **GDH**:growing degree hours
-   **PTU**:photothermal units
-   **Season.Length**:long or short season

**MetaData for HTG_climate_data.csv**

-   **id**:population id
-   **clim_year**:year of climate variable
-   **clim_month**:month of climate variable
-   **clim_date**:year, month combined 
-   **cwd**:total climate water deficit for the month and year
-   **pck**:total snow pack for the month and year
-   **ppt_mm**:total precipitation measured in millimeters for the month and year
-   **tmin**:average minimum temperature for the month and year in Celsius
-   **tmax**:average maximum temperature for the month and year in Celsius

**MetaData for weather.station.temp.csv**
##### weather.station.temp.csv is compiled from data in .txt files found in Raw.Data/Davis.Climate

-   **Month**:Month
-   **Day**:Day of the month
-   **Max.2m**:maximum air temperature 2 meters above the ground
-   **Min.2m**:minimum air temperature 2 meters above the ground
-   **Max.15m**:maximum air temperature 15 meters above the ground
-   **Min.15m**:minimum air temperature 15 meters above the ground

**MetaData for prism.data.1991_2020.csv**

-   **Name**:population name
-   **Date**:Date
-   **ppt.(inches)**:amount of precipitation measured in inches

**MetaData for all .csv files inside the Germ Fitness 2.0 Pollination Datasheets folder**

-   **Bench**:bench where the plant is located
-   **Block**:block that the plant is in 
-   **Pop**:population of the plant
-   **Cohort**:seed germination cohort of the plant
-   **Position**:position in the block the plant is located
-   **Date**:Date columns of when pollination was done or attempted. The actual dates differ within each file. Each file is a different population.

**MetaData for tree_pruned.new**

# newick formatted phylogeny

**MetaData for pheno.slopes.csv**

-   **sp**:species
-   **Pop**:population
-   **days.2.bud**:slopes of relationship between days to first bud and cohort
-   **days.2.flower**:slopes of relationship between days to first flower and cohort
-   **sept.1.bud**:slopes of relationship between days since Sept. 1 to first bud and cohort
-   **sept.1.flower**:slopes of relationship between days since Sept. 1 to first flower and cohort
-   **days.2.bud.se**:standard error of slopes of relationship between days to first bud and cohort
-   **days.2.flower.se**:standard error of slopes of relationship between days to first flower and cohort
-   **sept.1.bud.se**:standard error of slopes of relationship between days since Sept. 1 to first bud and cohort
-   **sept.1.flower.se**:standard error of slopes of relationship between days since Sept. 1 to first flower and cohort
-   **days.2.bud.CAAN1.CAIN3**:slopes of relationship between days to first bud and cohort with specific populations
-   **days.2.bud.CAAN1.CAIN3.se**:standard error of slopes of relationship between days to first bud and cohort with specific populations
-   **days.2.bud.CAAN1.CAIN4**:slopes of relationship between days to first bud and cohort with specific populations
-   **days.2.bud.CAAN1.CAIN4.se**:standard error of slopes of relationship between days to first bud and cohort with specific populations
-   **days.2.bud.CAAN2.CAIN3**:slopes of relationship between days to first bud and cohort with specific populations
-   **days.2.bud.CAAN2.CAIN3.se**:standard error of slopes of relationship between days to first bud and cohort with specific populations
-   **days.2.bud.CAAN2.CAIN4**:slopes of relationship between days to first bud and cohort with specific populations
-   **days.2.bud.CAAN2.CAIN4.se**:standard error of slopes of relationship between days to first bud and cohort with specific populations
-   **days.2.flower.CAAN1.CAIN3**:slopes of relationship between days to first flower and cohort with specific populations
-   **days.2.flower.CAAN1.CAIN3.se**:standard error of slopes of relationship between days to first flower and cohort with specific populations
-   **days.2.flower.CAAN1.CAIN4**:slopes of relationship between days to first flower and cohort with specific populations
-   **days.2.flower.CAAN1.CAIN4.se**:standard error of slopes of relationship between days to first flower and cohort with specific populations
-   **days.2.flower.CAAN2.CAIN3**:slopes of relationship between days to first flower and cohort with specific populations
-   **days.2.flower.CAAN2.CAIN3.se**:standard error of slopes of relationship between days to first flower and cohort with specific populations
-   **days.2.flower.CAAN2.CAIN4**:slopes of relationship between days to first flower and cohort with specific populations
-   **days.2.flower.CAAN2.CAIN4.se**:standard error of slopes of relationship between days to first flower and cohort with specific populations
-   **sept.1.bud.CAAN1.CAIN3**:slopes of relationship between days since Sept. 1 to first bud and cohort with specific populations
-   **sept.1.bud.CAAN1.CAIN3.se**:standard error of slopes of relationship between days since Sept. 1 to first bud and cohort with specific populations
-   **sept.1.bud.CAAN1.CAIN4**:slopes of relationship between days since Sept. 1 to first bud and cohort with specific populations
-   **sept.1.bud.CAAN1.CAIN4.se**:standard error of slopes of relationship between days since Sept. 1 to first bud and cohort with specific populations
-   **sept.1.bud.CAAN2.CAIN3**:slopes of relationship between days since Sept. 1 to first bud and cohort with specific populations
-   **sept.1.bud.CAAN2.CAIN3.se**:standard error of slopes of relationship between days since Sept. 1 to first bud and cohort with specific populations
-   **sept.1.bud.CAAN2.CAIN4**:slopes of relationship between days since Sept. 1 to first bud and cohort with specific populations
-   **sept.1.bud.CAAN2.CAIN4.se**:standard error of standard error of slopes of relationship between days to first bud and cohort with specific populations
-   **sept.1.flower.CAAN1.CAIN3**:slopes of relationship between days since Sept. 1 to first flower and cohort with specific populations
-   **sept.1.flower.CAAN1.CAIN3.se**:standard error of slopes of relationship between days since Sept. 1 to first flower and cohort with specific populations
-   **sept.1.flower.CAAN1.CAIN4**:slopes of relationship between days since Sept. 1 to first flower and cohort with specific populations
-   **sept.1.flower.CAAN1.CAIN4.se**:standard error of slopes of relationship between days since Sept. 1 to first flower and cohort with specific populations
-   **sept.1.flower.CAAN2.CAIN3**:slopes of relationship between days since Sept. 1 to first flower and cohort with specific populations
-   **sept.1.flower.CAAN2.CAIN3.se**:standard error of slopes of relationship between days since Sept. 1 to first flower and cohort with specific populations
-   **sept.1.flower.CAAN2.CAIN4**:slopes of relationship between days since Sept. 1 to first flower and cohort with specific populations
-   **sept.1.flower.CAAN2.CAIN4.se**:standard error of slopes of relationship between days since Sept. 1 to first flower and cohort with specific populations populations
-   **days.2.bud.final**:final slopes of relationship between days to first bud and cohort
-   **days.2.bud.final.se**:final standard error of slopes of relationship between days to first bud and cohort
-   **sept.1.bud.final**:final slopes of relationship between days since Sept. 1 to first bud and cohort
-   **sept.1.bud.final.se**:final standard error of slopes of relationship between days since Sept. 1 to first bud and cohort

**MetaData for fitness.slopes.2.csv**

-   **sp**:species
-   **Pop**:population
-   **pflwr**:slopes of relationship between probability of flowering and cohort
-   **sizeflw**:slopes of relationship between size at first flower and cohort
-   **seed_num_nb**:slopes of relationship between seed number and cohort
-   **year1fit**:slopes of relationship between first year fitness and cohort
-   **seed_mass**:slopes of relationship between seed mass and cohort
-   **pflwr.se**:standard error of slopes of relationship between probability of flowering and cohort
-   **sizeflw.se**:standard error of slopes of relationship between size at first flower and cohort
-   **seed_num_nb.se**:standard error of slopes of relationship between seed number and cohort
-   **year1fit.se**:standard error of slopes of relationship between first year fitness and cohort
-   **seed_mass.se**:standard error of slopes of relationship between seed mass and cohort
-   **sizeflw.CAAN1.CAIN3**:slopes of relationship between size at first flower and cohort with specific populations
-   **sizeflw.CAAN1.CAIN3.se**:standard error of slopes of relationship between size at first flower and cohort with specific populations
-   **sizeflw.CAAN1.CAIN4**:slopes of relationship between size at first flower and cohort with specific populations
-   **sizeflw.CAAN1.CAIN4.se**:standard error of slopes of relationship between size at first flower and cohort with specific populations
-   **sizeflw.CAAN2.CAIN3**:slopes of relationship between size at first flower and cohort with specific populations
-   **sizeflw.CAAN2.CAIN3.se**:standard error of slopes of relationship between size at first flower and cohort with specific populations
-   **sizeflw.CAAN2.CAIN4**:slopes of relationship between size at first flower and cohort with specific populations
-   **sizeflw.CAAN2.CAIN4.se**:standard error of slopes of relationship between size at first flower and cohort with specific populations
-   **seed_num_nb.CAAN1.CAIN3**:slopes of relationship between number of seeds and cohort with specific populations
-   **seed_num_nb.CAAN1.CAIN3.se**:standard error of slopes of relationship between number of seeds and cohort with specific populations
-   **seed_num_nb.CAAN1.CAIN4**:slopes of relationship between number of seeds and cohort with specific populations
-   **seed_num_nb.CAAN1.CAIN4.se**:standard error of slopes of relationship between number of seeds and cohort with specific populations
-   **seed_num_nb.CAAN2.CAIN3**:slopes of relationship between number of seeds and cohort with specific populations
-   **seed_num_nb.CAAN2.CAIN3.se**:standard error of slopes of relationship between number of seeds and cohort with specific populations
-   **seed_num_nb.CAAN2.CAIN4**:slopes of relationship between number of seeds and cohort with specific populations
-   **seed_num_nb.CAAN2.CAIN4.se**:standard error of slopes of relationship between number of seeds and cohort with specific populations
-   **year1fit.CAAN1.CAIN3**:slopes of relationship between first year fitness and cohort with specific populations
-   **year1fit.CAAN1.CAIN3.se**:standard error of slopes of relationship between first year fitness and cohort with specific populations
-   **year1fit.CAAN1.CAIN4**:slopes of relationship between first year fitness and cohort with specific populations
-   **year1fit.CAAN1.CAIN4.se**:standard error of slopes of relationship between first year fitness and cohort with specific populations
-   **year1fit.CAAN2.CAIN3**:slopes of relationship between first year fitness and cohort with specific populations
-   **year1fit.CAAN2.CAIN3.se**:standard error of slopes of relationship between first year fitness and cohort with specific populations
-   **year1fit.CAAN2.CAIN4**:slopes of relationship between first year fitness and cohort with specific populations
-   **year1fit.CAAN2.CAIN4.se**:standard error of slopes of relationship between first year fitness and cohort with specific populations
-   **seed_mass.CAAN1.CAIN3**:slopes of relationship between seed mass and cohort with specific populations
-   **seed_mass.CAAN1.CAIN3.se**:standard error of slopes of relationship between seed mass and cohort with specific populations
-   **seed_mass.CAAN1.CAIN4**:slopes of relationship between seed mass and cohort with specific populations
-   **seed_mass.CAAN1.CAIN4.se**:standard error of slopes of relationship between seed mass and cohort with specific populations
-   **seed_mass.CAAN2.CAIN3**:slopes of relationship between seed mass and cohort with specific populations
-   **seed_mass.CAAN2.CAIN3.se**:standard error of slopes of relationship between seed mass and cohort with specific populations
-   **seed_mass.CAAN2.CAIN4**:slopes of relationship between seed mass and cohort with specific populations
-   **seed_mass.CAAN2.CAIN4.se**:standard error of slopes of relationship between seed mass and cohort with specific populations
-   **pflwr.weighted**:slopes of relationship between probability of flowering and cohort using weighted model
-   **pflwr.CAAN1.CAIN3**:slopes of relationship between probability of flowering and cohort with specific populations
-   **pflwr.CAAN1.CAIN3.se**:standard error of slopes of relationship between probability of flowering and cohort with specific populations
-   **pflwr.CAAN1.CAIN4**:slopes of relationship between probability of flowering and cohort with specific populations
-   **pflwr.CAAN1.CAIN4.se**:standard error of slopes of relationship between probability of flowering and cohort with specific populations
-   **pflwr.CAAN2.CAIN3**:slopes of relationship between probability of flowering and cohort with specific populations
-   **pflwr.CAAN2.CAIN3.se**:standard error of slopes of relationship between probability of flowering and cohort with specific populations
-   **pflwr.CAAN2.CAIN4**:slopes of relationship between probability of flowering and cohort with specific populations
-   **pflwr.CAAN2.CAIN4.se**:standard error of slopes of relationship between probability of flowering and cohort with specific populations
-   **pflwr.final**:final slopes of relationship between probability of flowering and cohort using weighted model
-   **pflwr.final.se**:final standard error of slopes of relationship between probability of flowering and cohort using weighted model
-   **sizebud.final**:final slopes of relationship between size at first bud and cohort
-   **sizebud.final.se**:final standard error of slopes of relationship between size at first bud  and cohort
-   **seed_num_nb.final**:final slopes of relationship between number of seeds and cohort
-   **seed_num_nb.final.se**:final standard error of slopes of relationship between number of seeds and cohort
-   **year1fit.final**:final slopes of relationship between first year fitness and cohort
-   **year1fit.final.se**:final standard error of slopes of relationship between first year fitness and cohort
-   **seed_mass.final**:final slopes of relationship between total seed mass and cohort
-   **seed_mass.final.se**:final standard error of slopes of relationship between total seed mass and cohort 

**MetaData for all.slopes.csv**

-   **sp**:species
-   **Pop**:population
-   **days.2.bud**:slopes of relationship between days to first bud and cohort
-   **days.2.flower**:slopes of relationship between days to first flower and cohort
-   **sept.1.bud**:slopes of relationship between days since Sept. 1 to first bud and cohort
-   **sept.1.flower**:slopes of relationship between days since Sept. 1 to first flower and cohort
-   **pflwr**:slopes of relationship between probability of flowering and cohort
-   **sizeflw**:slopes of relationship between size at first flower and cohort
-   **seed_num_nb**:slopes of relationship between seed number and cohort
-   **year1fit**:slopes of relationship between first year fitness and cohort
-   **seed_mass**:slopes of relationship between seed mass and cohort
-   **seed_num_nb.CAAN1.CAIN3**:slopes of relationship between number of seeds and cohort with specific populations
-   **seed_num_nb.CAAN1.CAIN4**:slopes of relationship between number of seeds and cohort with specific populations
-   **seed_num_nb.CAAN2.CAIN3**:slopes of relationship between number of seeds and cohort with specific populations
-   **seed_num_nb.CAAN2.CAIN4**:slopes of relationship between number of seeds and cohort with specific populations
-   **days.2.bud.CAAN1.CAIN3**:slopes of relationship between days to first bud and cohort with specific populations
-   **days.2.bud.CAAN1.CAIN4**:slopes of relationship between days to first bud and cohort with specific populations
-   **days.2.bud.CAAN2.CAIN3**:slopes of relationship between days to first bud and cohort with specific populations
-   **days.2.bud.CAAN2.CAIN4**:slopes of relationship between days to first bud and cohort with specific populations
-   **days.2.flower.CAAN1.CAIN3**:slopes of relationship between days to first flower and cohort with specific populations
-   **days.2.flower.CAAN1.CAIN4**:slopes of relationship between days to first flower and cohort with specific populations
-   **days.2.flower.CAAN2.CAIN3**:slopes of relationship between days to first flower and cohort with specific populations
-   **days.2.flower.CAAN2.CAIN4**:slopes of relationship between days to first flower and cohort with specific populations
-   **year1fit.CAAN1.CAIN3**:slopes of relationship between first year fitness and cohort with specific populations
-   **year1fit.CAAN1.CAIN4**:slopes of relationship between first year fitness and cohort with specific populations
-   **year1fit.CAAN2.CAIN3**:slopes of relationship between first year fitness and cohort with specific populations
-   **year1fit.CAAN2.CAIN4**:slopes of relationship between first year fitness and cohort with specific populations
-   **pflwr.weighted**:slopes of relationship between probability of flowering and cohort using weighted model
-   **pflwr.CAAN1.CAIN3**:slopes of relationship between probability of flowering and cohort with specific populations
-   **pflwr.CAAN1.CAIN4**:slopes of relationship between probability of flowering and cohort with specific populations
-   **pflwr.CAAN2.CAIN3**:slopes of relationship between probability of flowering and cohort with specific populations
-   **pflwr.CAAN2.CAIN4**:slopes of relationship between probability of flowering and cohort with specific populations
-   **days.2.bud.final**:final slopes of relationship between days to first bud and cohort
-   **days.2.bud.final.se**:final standard error of slopes of relationship between days to first bud and cohort
-   **sept.1.bud.final**:final slopes of relationship between days since Sept. 1 to first bud and cohort
-   **sept.1.bud.final.se**:final standard error of slopes of relationship between days since Sept. 1 to first bud and cohort
-   **pflwr.final**:final slopes of relationship between probability of flowering and cohort using weighted model
-   **pflwr.final.se**:final standard error of slopes of relationship between probability of flowering and cohort using weighted model
-   **sizebud.final**:final slopes of relationship between size at first bud and cohort
-   **sizebud.final.se**:final standard error of slopes of relationship between size at first bud  and cohort
-   **seed_num_nb.final**:final slopes of relationship between number of seeds and cohort
-   **seed_num_nb.final.se**:final standard error of slopes of relationship between number of seeds and cohort
-   **year1fit.final**:final slopes of relationship between first year fitness and cohort
-   **year1fit.final.se**:final standard error of slopes of relationship between first year fitness and cohort
-   **seed_mass.final**:final slopes of relationship between total seed mass and cohort
-   **seed_mass.final.se**:final standard error of slopes of relationship between total seed mass and cohort 

**MetaData for georeferencing_clean.csv**

-   **folder**:species
-   **specimen**:herbarium specimen number
-   **day**:day of collection
-   **month**:month of collection
-   **year**:year of collection
-   **decimalLatitude**:latitude of specimen location
-   **decimalLongitude**:longitude of specimen location
-   **coordinateUncertaintyInMeters**:uncertainity in specimen location in meters
-   **minimumElevationInMeters**:minimum elevation in meters
-   **maximumElevationInMeters**:maximum elevation in meters

**MetaData for all_herbarium_climate.csv**

-   **specimen**:herbarium specimen number
-   **clim_year**:year of climate variable
-   **clim_month**:month of climate variable
-   **clim_date**:year, month combined 
-   **cwd**:total climate water deficit for the month and year
-   **pck**:total snow pack for the month and year
-   **ppt_mm**:total precipitation measured in millimeters for the month and year
-   **snw**:total snow for the month and year
-   **str**:total soil water storage for the month and year
-   **tmin**:average minimum temperature for the month and year in Celsius
-   **tmax**:average maximum temperature for the month and year in Celsius

#### Results

This folder includes results of analyses presented in the manuscript along with graphs that were imported into word to make figures found in the manuscript and supplemental material. The subfolder SEM.plots includes species-specific graphs of the SEM diagrams.

#### Scripts

This folder includes all scripts used for analyses in the manuscript. Within each script file, there is a brief description of what the script was used to do.

**Order to Use Scripts**

1. data_cleaning.R
2. phenology_analyses.R
3. fitness_analyses.R
4. prep_iButton_Data.R
5. growing_conditions.R
6. SEM_analyses.R
7. phylogenetic_analyses.R
8. climate_space.R

