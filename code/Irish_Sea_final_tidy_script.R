# Script to reproduce the analyses and figures presented in the manuscript:
# "Functional diversity buffers multi-decadal shifts in fish-mediated bioturbation within a heavily exploited sea"
# Written by Mara Fischer (m.fischer2@exeter.ac.uk)
# April 2026

# Libraries ----
library(tidyverse)
library(maps)
library(sf)
library(DHARMa)
library(mgcv)
library(visreg)
library(visibly)
library(vegan)
library(glmmTMB)
library(emmeans)
library(ggrepel)
library(purrr)
library(broom)
library(ggpubr)
library(ggeffects)
library(cowplot)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)

# Load raw data ----

# ICES Trawl Survey data:

## Data downloaded from DATRAS (https://datras.ices.dk/)
## Specifications: Survey: BTS, Quarter: 3, Years: all, countries: GB, species: all

## Two raw datasets:
## HH: contains metadata, incl. haul duration
## HL: contains biological info

hl <- read.csv("raw_data/IBTS_data/DATRAS_HL_all_years.csv", header = T)
hh <- read.csv("raw_data/IBTS_data/DATRAS_HH_all_years.csv", header = T)

# Bioturbation data

## UK fish species bioturbation classification (described in Fischer et al., 2025 and 
## downloaded from: doi: 10.17632/t8kg43f8kn.1 
biot <- read.csv("raw_data/bioturbation.csv", header = T)

# EMODnet substrate data
# based on a hierarchy of seven Folk classes from EMODnet Geology (Kaskela et al., 2019; EMODnet, 2025). 
substrate <- st_read("raw_data/seabed_substrate_250k.shp")

## ICES stat recs shapefile
stat_recs <- read_sf("raw_data/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp")

# Load species code file
spp <- read.csv("raw_data/landings_species.csv")

# Custom plot theme ----

theme_convex <- function(){
  theme_bw() +
    theme(text = element_text(family = "Arial"),
          axis.text = element_text(size = 9), 
          axis.title = element_text(size = 11),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 16, vjust = 1, hjust = 0),
          legend.text = element_text(size = 11),          
          #legend.title = element_blank(),                              
          legend.position = c(0.95, 0.15), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           linewidth = 2, linetype = "blank"))
}

# Data manipulation ----

## Trawl data tidying ----

# Select BT4AI gear for Irish Sea only

hl <- hl %>% 
  filter(Gear == "BT4AI")

hh <- hh %>% 
  filter(Gear == "BT4AI")

# Create unique haul ID, made up of year, station, haul

hl <- hl %>% 
  mutate(ID = paste(Year, StNo, HaulNo, sep = "_"))

length(unique(hl$ID)) # 4068, should match the length of rows in hh

hh <- hh %>% 
  mutate(ID = paste(Year, StNo, HaulNo, sep = "_"),
         Date = paste(Day, Month, Year, sep = "/")) # add date column

# Calculate Julian day within the year
hh$DayInYear = strptime(hh$Date, "%d/%m/%Y")$yday + 1

length(unique(hh$ID)) # 4068

sum(hh$SurSa > 0)

# Only keep relevant columns from hh
hh_tidy <- hh %>% 
  dplyr::select(ID, Quarter, Country, Ship, Gear, Year, Month, Day,
                Date, DayInYear, TimeShot,
                StNo, HaulNo, HaulDur, Depth, StatRec, 
                HaulLat, HaulLong)

# Tidy HL
hl_tidy <- hl %>% 
  dplyr::select(ID, Year, StNo, HaulNo,
                TotalNo, CatCatchWgt, LngtCode, LngtClass, HLNoAtLngt, ScientificName_WoRMS) %>% 
  rename(species = ScientificName_WoRMS,
         weight = CatCatchWgt,
         no.at.lngt = HLNoAtLngt) %>% 
  filter(no.at.lngt != "-9")  # remove observations with missing catch data

# Join together
data <- hl_tidy %>% 
  left_join(hh_tidy, by = c("ID", "Year", "StNo", "HaulNo"))

unique(data$StatRec)

# Subset to Irish Sea statistical rectangles
# calculate CPUE
data <- data %>% 
  filter(StatRec %in% c("38E4","38E5","38E6",
                        "37E3","37E4","37E5","37E6","37E7",
                        "36E3","36E4","36E5","36E6","36E7",
                        "35E3","35E4","35E5","35E6",
                        "34E3","34E4","34E5", 
                        "33E2","33E3","33E4","33E5"),                    # subset to Irish Sea
         HaulDur >= 15) %>%                                              # remove haul duration of < 15 mins
  mutate(CPUE = no.at.lngt / (HaulDur / 60),                             # calculate cpue per hour for each length class
         Length_mm = ifelse(LngtCode == "1", LngtClass * 10, LngtClass)) # add length class in mm

write.csv(data, "processed_data/TIDY_IBTS_DATA.csv")

## Add bioturbation to trawl data ----

## Tidy bioturbation data

biot <- biot %>% 
  dplyr::select(- X, - Suborder, - Incertae.sedis, - weight.g, - weight.kg, 
                - log10.size.class, - size.category, - size.category.score,
                - bioturbation.x.size.score, - Description, - bioturbation.score, 
                - commercial, - bioturbation.mode)

# Subset trawl dataset to just bioturbating species

biot_data <- data %>%
  inner_join(biot, by = c("species" = "Species"))

length(unique(biot_data$species)) 
unique(biot_data$species)

## Calculate BPi ----

## Calculate BPi using individual square-root transformed weights in each length class.
## This results in a per-capita BPi value for a individual in a given length class.
## The square root transformation follows the approach by Solan et al (2004) and is
## to ensure that we linearize the relationship between size and impact.
## If we were using the untransformed weights, we would essentially assume that a
## fish that weighs 100g has 100x higher impact than a fish that weighs 1g
## but this is not biologically realistic. So instead with square-root transformed weight,
## the larger fish has 10x higher impact (sqrt(100) = 10). 

head(biot_data)

biot_final <- biot_data %>%
  mutate(
    weight.g = a * (Length_mm / 10)^b,                               # Calculate individual weight in g
    weight.cpue = weight.g * CPUE,                                   # weight in grams per hour
    weight.cpue.kg = weight.cpue/1000,                               # weight in kg per hour
    BPi = bioturbation.mode.score * frequency.score *sqrt(weight.g), # calculate BPi using indiv. weights
    BPclass = BPi * CPUE                                             # get BP per length class by multiplying BPi by CPUE
  )

head(biot_final)

# Reorder bioturbation modes
biot_final$Bioturbation.mode <- factor(biot_final$Bioturbation.mode,
                                       levels = c("Burrower", "Vertical excavator",
                                                  "Nest-builder", "Lateral excavator",
                                                  "Sediment sifter"))
# add unique site column
biot_final <- biot_final %>% 
  mutate(Site = paste(StatRec, StNo, sep = "_"))

length(unique(biot_final$Site))

## Calculate BPp ----

# Calculate population-level BP per species × haul by summing across length classes
biot_haul <- biot_final %>%
  group_by(ID, Year, StNo, HaulNo, species) %>%
  summarise(
    Family = unique(Family),
    bioturbation.mode = unique(Bioturbation.mode),
    Depth = first(Depth),
    Date = first(Date),
    DayInYear = first(DayInYear),
    TimeIn = unique(TimeShot),
    Site = unique(Site),
    HaulLat = unique(HaulLat),
    HaulLong = unique(HaulLong),
    StatRec = unique(StatRec),
    BPp = sum(BPclass, na.rm = TRUE),
    total_abun_cpue = sum(CPUE, na.rm = TRUE),
    total_biomass_cpue = sum(weight.cpue, na.rm = TRUE)
  )

length(unique(biot_haul$Site)) # 188

# Extract site coordinates
sites <- biot_haul %>% 
  group_by(Site) %>% 
  summarise(Lat = mean(HaulLat),
            Long = mean(HaulLong),
            StatRec = unique(StatRec))

# Add ecological groups
biot_haul <- biot_haul %>% 
  mutate(Group = case_when(
    species %in% c("Raja clavata", "Raja montagui", "Raja brachyura", 
                   "Raja microocellata", "Raja undulata") ~ "Skates",
    species %in% c("Buglossidium luteum", "Glyptocephalus cynoglossus",
                   "Hippoglossoides platessoides", "Limanda limanda", 
                   "Microchirus variegatus", "Microstomus kitt", 
                   "Pegusa lascaris", "Platichthys flesus", 
                   "Pleuronectes platessa", "Lepidorhombus whiffiagonis",
                   "Scophthalmus maximus", "Scophthalmus rhombus", 
                   "Zeugopterus punctatus", "Zeugopterus regius", 
                   "Solea solea", "Arnoglossus laterna", 
                   "Arnoglossus imperialis") ~ "Flatfish",
    species %in% c("Gadus morhua", "Merlangius merlangus", 
                   "Melanogrammus aeglefinus", "Pollachius pollachius", 
                   "Trisopterus minutus", "Trisopterus luscus", 
                   "Enchelyopus cimbrius", "Raniceps raninus", 
                   "Molva molva") ~ "Gadoids",
    species %in% c("Scyliorhinus canicula", "Scyliorhinus stellaris", 
                   "Squalus acanthias", "Mustelus asterias", 
                   "Mustelus mustelus", "Galeorhinus galeus") ~ "Sharks",
    species %in% c("Callionymus lyra", "Callionymus maculatus", "Callionymus reticulatus",
                   "Ammodytes marinus", "Ammodytes tobianus", 
                   "Hyperoplus immaculatus", "Hyperoplus lanceolatus",
                   "Gobius niger", "Gobius paganellus", "Lesueurigobius friesii",
                   "Parablennius gattorugine", "Blennius ocellaris",
                   "Pholis gunnellus", "Chirolophis ascanii", 
                   "Liparis liparis", "Liparis montagui",
                   "Taurulus bubalis", "Artediellus atlanticus",
                   "Micrenophrys lilljeborgii", "Myoxocephalus scorpius",
                   "Myoxocephalus scorpioides", "Entelurus aequoreus",
                   "Syngnathus acus", "Syngnathus typhle",
                   "Diplecogaster bimaculata", "Spinachia spinachia",
                   "Ciliata mustela", "Ciliata septentrionalis",
                   "Gaidropsarus mediterraneus", "Gaidropsarus vulgaris") ~ "Small benthic",
    TRUE ~ "Other"
  ))

## Add substrate data ----

# Benthic substrate data from EMODnet
# Relevant scale: Folk-7

# After quick visual inspection, there are two coordinate outliers:in the biot_haul data
# the one in the Atlantic the Latitude entered twice... delete that one
# the one in the North Sea is just missing the minus sign in front of the Long - add that

# Remove 
biot_haul2 <- biot_haul %>% 
  filter(ID != "2006_408_80")

# Flip coordinate on North Sea point
biot_haul2 <- biot_haul2 %>% 
  mutate(HaulLong = ifelse(HaulLong == 3.8567, -3.8567, HaulLong))

# Convert hauls to sf object

hauls_sf <- st_as_sf(
  biot_haul2,
  coords = c("HaulLong", "HaulLat"),
  crs = 4326  # WGS84
)

st_crs(hauls_sf)
st_crs(substrate)

# Select relevant columns
substrate_simple <- substrate %>% 
  dplyr::select(folk_7cl, folk_7cl_t, geometry)

# Extract substrate class to each haul
hauls_sf <- st_join(hauls_sf, substrate_simple)

# This will attach the substrate polygon’s attributes to every haul point.

hauls_sf <- hauls_sf %>% 
  rename(substrate_class = folk_7cl,
         substrate_name = folk_7cl_t)

unique(hauls_sf$substrate_name)

# count NAs
sum(is.na(hauls_sf$substrate_name)) # 52 NAs
sum(is.na(hauls_sf$substrate_class))

# Subset rows with NA in substrate
na_substrate <- hauls_sf %>% 
  filter(is.na(substrate_class))

# How many unique hauls?
length(unique(na_substrate$ID)) # 4 hauls

# Get UK/Ireland map data
uk_map <- map_data("world") %>% 
  filter(region %in% c("UK", "Ireland"))

# Plot
ggplot() +
  geom_polygon(data = uk_map, aes(x = long, y = lat, group = group),
               fill = "lightgrey", color = "black") +
  #geom_point(data = na_substrate, aes(x = HaulLong, y = HaulLat),
  #          color = "red", size = 2, alpha = 0.6) +
  geom_sf(data = na_substrate) +
  #coord_fixed(xlim = c(-10, 5), ylim = c(49, 61)) +
  coord_sf(xlim = c(-8, -3), ylim = c(51.5, 55)) +
  theme_minimal() +
  labs(title = "Trawl Survey Haul Locations",
       x = "Longitude", y = "Latitude")

# one point on land...
# remove those four hauls

# Get the unique haul IDs with NA
hauls_to_remove <- unique(hauls_sf$ID[is.na(hauls_sf$substrate_name)])

hauls_clean <- hauls_sf[!(hauls_sf$ID %in% hauls_to_remove), ]

length(unique(hauls_clean$ID))

# Tidy (remove numbers from front of substrate categories)
hauls_clean <- hauls_clean %>% 
  separate(substrate_name, into = c(NA, "substrate_name"), sep = "\\.\\s+")

unique(hauls_clean$substrate_name)

# Turn sf object back into normal dataframe for modelling
hauls_clean <- st_drop_geometry(hauls_clean)

str(hauls_clean)

# Add Long and Lat again

# site coords by haul
coords <- biot_haul2 %>% 
  group_by(ID) %>% 
  summarise(lat = unique(HaulLat),
            lon = unique(HaulLong))

# Join coords to haul level data
hauls_clean <- hauls_clean %>% 
  left_join(coords, by = "ID")

# Save updated haul-level dataset with substrate data
# tidy data to use in temporal analysis
write.csv(hauls_clean, "processed_data/haul_level_data_with_BPp_and_substrate.csv")

## Calculate BPc ----

# Community-level BP per haul
biot_comm_haul <- hauls_clean %>%
  group_by(ID, Year) %>%
  summarise(Depth = first(Depth),
            Date = first(Date),
            DayInYear = first(DayInYear),
            TimeIn = unique(TimeIn),
            Site = unique(Site),
            lat = mean(lat),
            lon = mean(lon),
            StatRec = unique(StatRec),
            BPc = sum(BPp, na.rm = TRUE),
            Log10_BPc = log10(BPc + 1),
            abun_cpue = sum(total_abun_cpue, na.rm = TRUE),
            biomass_cpue = sum(total_biomass_cpue, na.rm = TRUE))

hist(biot_comm_haul$BPc) # right-tailed, lots of low values
hist(biot_comm_haul$Log10_BPc)  # more normal

# Annual mean (trend)
biot_year <- biot_comm_haul %>%
  group_by(Year) %>%
  summarise(mean_BPc = mean(BPc, na.rm = TRUE),
            se_BPc = sd(BPc)/sqrt(n()),
            n = n(),
            mean_abun_cpue = mean(abun_cpue, na.rm = TRUE),
            mean_biomass_cpue = mean(biomass_cpue, na.rm = TRUE))

# Temporal analysis ----

# Analysis on subset of sites surveyed consistently over time

## Subset sites ----

# Subset to sites surveyed in >= 5 years and spanning >= 10 years

head(hauls_clean)

# Aggregate Community-level per haul
hauls_com <- hauls_clean %>%
  group_by(ID, Year, Site) %>%
  summarise(Depth = first(Depth),
            Date = first(Date),
            DayInYear = first(DayInYear),
            TimeIn = unique(TimeIn),
            StatRec = unique(StatRec),
            substrate = unique(substrate_name),
            lat = mean(lat),
            lon = mean(lon),
            BPc = sum(BPp, na.rm = TRUE),
            abun_cpue = sum(total_abun_cpue, na.rm = TRUE),
            biomass_cpue = sum(total_biomass_cpue, na.rm = TRUE))

# this is the same as biot_comm_haul

length(unique(hauls_com$Site)) # 188 sites
length(unique(hauls_com$ID)) # 2529 hauls, so each row is a haul

# Are there multiple hauls per site per year?
hauls_per_site_year <- hauls_com  %>% 
  count(Site, Year, name = "n_hauls")

summary(hauls_per_site_year$n_hauls)

# only 1 haul per sites per year so we are good

# Check site coverage over time
site_coverage <- hauls_com  %>% 
  group_by(Site)  %>% 
  summarise(
    n_years = n(),
    first_year = min(Year),
    last_year  = max(Year),
    span_years = last_year - first_year + 1,
    .groups = "drop"
  )

# Exclude sites surveyed in fewer than 5 years and spanning less than 10 years
keep_sites <- site_coverage %>% 
  filter(n_years >= 5, span_years >= 10)

length(unique(keep_sites$Site)) # 74 sites left

hauls_subset <- hauls_com  %>% 
  filter(Site %in% keep_sites$Site)

unique(hauls_subset$Site) # 74 sites

# log10 transform BPc and scale depth
hauls_subset <- hauls_subset %>% 
  ungroup() %>% 
  mutate(Log10_BPc = log10(BPc),
         Depth_s = as.numeric(scale(Depth)))

## Check outliers ----

# IQR: look at global outliers

# Using haul-level dataset (hauls_subset)
# Calculate the global threshold
Q1 <- quantile(hauls_subset$BPc, 0.25)
Q3 <- quantile(hauls_subset$BPc, 0.75)
IQR_val <- Q3 - Q1

# Define the upper cutoff
upper_limit <- Q3 + (3 * IQR_val) 

# Identify the outliers
outlier_hauls <- hauls_subset %>%
  filter(BPc > upper_limit)

# Winsorization (cap the data at upper threshold)
# outliers identified are replaced by the upper limit rather than deleted
hauls_subset <- hauls_subset %>%
  mutate(BPc_capped = ifelse(BPc > upper_limit, upper_limit, BPc))

# compare
p1 <- ggplot(hauls_subset, aes(x = Year, y = BPc)) + 
  geom_point(alpha = 0.3) + geom_smooth() + labs(title = "Raw Data")

p2 <- ggplot(hauls_subset, aes(x = Year, y = BPc_capped)) + 
  geom_point(alpha = 0.3) + geom_smooth() + labs(title = "Winsorized (Capped) Data")

p1 + p2

# use capped data in analysis

## GAMMs ----

# Ensure covariates are factors
hauls_subset$Site <- as.factor(hauls_subset$Site)
hauls_subset$substrate <- factor(hauls_subset$substrate)

# OG model structure: allow year trend to vary by substrate

# Determine best distribution:

# run with ML to compare
GAM1.1 <- gam(BPc_capped ~ s(Year, by = substrate, bs = "tp") + 
                s(Year, bs = "tp") + 
                substrate,
              data = hauls_subset, method = "ML")

GAM1.2 <- gam(BPc_capped ~ s(Year, by = substrate, bs = "tp") + 
                s(Year, bs = "tp") + 
                substrate,
              family = Gamma(link="log"),
              data = hauls_subset, method = "ML")

AIC(GAM1.1, GAM1.2)
54116.47 - 52603.37 # delta AIC = 1513.1, so gamma is better model

# Run different combinations and do backwards stepwise selection to find best model
# using maximum restricted likelihood (ML)

# The Interaction Model (Most complex) - and different combinations of covariates
# Does the trend vary by habitat?
m1 <- gam(BPc_capped ~ s(Year, by = substrate, bs = "tp") + substrate + Depth_s + s(Site, bs = "re"), 
          family = Gamma(link="log"), data = hauls_subset, method = "ML")

m1_sum <- summary(m1)

m1.2 <- gam(BPc_capped ~ s(Year, by = substrate, bs = "tp") + s(Year, bs = "tp") + substrate + Depth_s + s(Site, bs = "re"), 
            family = Gamma(link="log"), data = hauls_subset, method = "ML")

m1.2_sum <- summary(m1.2)

m1.3 <- gam(BPc_capped ~ s(Year, by = substrate, bs = "tp") + + s(Year, bs = "tp") + substrate + s(Site, bs = "re"), 
            family = Gamma(link="log"), data = hauls_subset, method = "ML")

m1.3_sum <- summary(m1.3)

m1.4 <- gam(BPc_capped ~ s(Year, by = substrate, bs = "tp") + s(Year, bs = "tp") + Depth_s + s(Site, bs = "re"), 
            family = Gamma(link="log"), data = hauls_subset, method = "ML")

m1.4_sum <- summary(m1.4)

m1.5 <- gam(BPc_capped ~ s(Year, by = substrate, bs = "tp") + s(Year, bs = "tp") + s(Site, bs = "re"), 
            family = Gamma(link="log"), data = hauls_subset, method = "ML")

m1.5_sum <- summary(m1.5)

m1.6 <- gam(BPc_capped ~ s(Year, by = substrate, bs = "tp") + s(Site, bs = "re"), 
            family = Gamma(link="log"), data = hauls_subset, method = "ML")

m1.6_sum <- summary(m1.6)

m1.7 <- gam(BPc_capped ~ s(Year, by = substrate, bs = "tp") + substrate + s(Site, bs = "re"), 
            family = Gamma(link="log"), data = hauls_subset, method = "ML")

m1.7_sum <- summary(m1.7)

m1.8 <- gam(BPc_capped ~ s(Year, by = substrate, bs = "tp") + Depth_s + s(Site, bs = "re"), 
            family = Gamma(link="log"), data = hauls_subset, method = "ML")

m1.8_sum <- summary(m1.8)

# The Additive Model
# Is there one shared trend across all habitats?
m2 <- gam(BPc_capped ~ s(Year, bs = "tp") + substrate + Depth_s + s(Site, bs = "re"), 
          family = Gamma(link="log"), data = hauls_subset, method = "ML")

m2_sum <- summary(m2)

# No Depth
m3 <- gam(BPc_capped ~ s(Year, bs = "tp") + substrate + s(Site, bs = "re"), 
          family = Gamma(link="log"), data = hauls_subset, method = "ML")

m3_sum <- summary(m3)

# No Substrate
m4 <- gam(BPc_capped ~ s(Year, bs = "tp") + Depth_s + s(Site, bs = "re"), 
          family = Gamma(link="log"), data = hauls_subset, method = "ML")

m4_sum <- summary(m4)

# Year only
m5 <- gam(BPc_capped ~ s(Year, bs = "tp") + s(Site, bs = "re"), 
          family = Gamma(link="log"), data = hauls_subset, method = "ML")

m5_sum <- summary(m5)

# No Year
m6 <- gam(BPc_capped ~ substrate + Depth_s + s(Site, bs = "re"), 
          family = Gamma(link="log"), data = hauls_subset, method = "ML")


m6_sum <- summary(m6)

# 5. Null model
m7 <- gam(BPc_capped ~ 1 + s(Site, bs = "re"), 
          family = Gamma(link="log"), data = hauls_subset, method = "ML")

m7_sum <- summary(m7)

# Compare AIC
AIC_table <- AIC(m1, m1.2, m1.3, m1.4, m1.5, m1.6, m1.7, m1.8, m2, m3, m4, m5, m6, m7) %>%
  arrange(AIC) %>%
  mutate(delta_AIC = AIC - min(AIC))

head(AIC_table, 10) 

# m1.7 is best (but only 1.2 lower than 1.6, but makes sense to include substrate as a parametric term)

model_summaries <- list(m1_sum, m1.2_sum, m1.3_sum, m1.4_sum, m1.5_sum, m1.6_sum, m1.7_sum, m1.8_sum,
                        m2_sum, m3_sum, m4_sum, m5_sum, m6_sum, m7_sum)

AIC_scores <- AIC(m1, m1.2, m1.3, m1.4, m1.5, m1.6, m1.7, m1.8, m2, m3, m4, m5, m6, m7)
AIC_scores 

## interestingly, the null model including only site as random effect, still has an R2 of 0.6 (compared to 0.64 for the best model)
## this suggests that variability between sites actually explains most of the change in BPc?

### Best model ----

# refit with REML:

m_final <- gam(BPc_capped ~ s(Year, by = substrate, bs = "tp") + substrate + s(Site, bs = "re"), 
               family = Gamma(link="log"), 
               data = hauls_subset, 
               method = "REML")

summary(m_final)

# Plot overall year effect
# Figure S1

# Predict BPc across the years, averaging over substrate and site
net_trend <- ggpredict(m_final, terms = "Year [all]")

# Get predicted mean starting and end values
start_val <- net_trend$predicted[1]
end_val <- net_trend$predicted[nrow(net_trend)]

# Calculate percentage change
pct_change <- ((end_val - start_val) / start_val) * 100 # 20.6%

# Plot the "Net" Irish Sea Trend
(GAM_net_year_plot <-
    ggplot(net_trend, aes(x = x, y = predicted)) +
    geom_line(size = 1, color = "darkblue") +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "darkblue") +
    labs(#title = "Net Bioturbation Potential: Irish Sea (1988-2024)",
      x = "Year",
      y = expression(BP[c])) +
    scale_y_continuous(expand = c(0, 0.1), limits = c(0,70000)) +
    # lims(y = c(0,70000)) +
    #scale_y_log10() +
    theme_convex()
)

ggsave(GAM_net_year_plot, filename = "outputs/Figure_S1.png", width = 5, height = 4, dpi = 600)

### Diagnostics ----

gam.check(m_final, k.rep = 1000) # k too low?

# refit with higher k
m_final2 <- gam(BPc_capped ~ s(Year, by = substrate, bs = "tp", k = 20) + substrate + s(Site, bs = "re"), 
                family = Gamma(link="log"), 
                data = hauls_subset, 
                method = "REML")

summary(m_final2)

gam.check(m_final2, k.rep = 1000)

# increasing k only increases the wigglyness of sandy mud, all others remain unchanged

plot_gam_check(m_final) # fairly good fit

concurvity(m_final, full = TRUE)

plot(m_final2, pages = 1, all.terms = TRUE, residuals = TRUE, seWithMean = TRUE, 
     shift = coef(m1.7)[1], trans = exp)

# look at spatial autocorrelation
res <- simulateResiduals(m_final)

# need a vector of X and Y coordinates for hauls

# Recalculate residuals by location
# This groups the errors by their Lat/Lon so there is only 1 point per coordinate
res_spatial <- recalculateResiduals(res, group = interaction(hauls_subset$lon, hauls_subset$lat))

# Extract the unique coordinates that match the grouped residuals
coords_unique <- hauls_subset %>%
  group_by(lon, lat) %>%
  summarize(.groups = "drop")

# Run the test
testSpatialAutocorrelation(res_spatial, x = coords_unique$lon, y = coords_unique$lat)

# no spatial autocorrelation after accounting for habitat and site-level random effects
# Moran's I test: p = 0.62

# table of model outputs:

# Get parametric terms
parametric_table <- tidy(m_final, parametric = TRUE)

# Get smooth terms
smooth_table <- tidy(m_final, parametric = FALSE)

# export to csv
write.csv(parametric_table, "outputs/Table_GAMM_Parametric.csv")
write.csv(smooth_table, "outputs/Table_GAMM_Smooths.csv")

### Visualise final GAMM ----

# Create a data frame for prediction
# create a sequence of years for every substrate type
plot_data <- expand.grid(
  Year = seq(min(hauls_subset$Year), max(hauls_subset$Year), length.out = 200),
  substrate = unique(hauls_subset$substrate)
)

# Add a dummy Site (the 'exclude' argument in predict handles the random effect)
plot_data$Site <- hauls_subset$Site[1] 

# Predict values and Standard Errors
# type = "link" gives us the log scale, which we then transform
preds <- predict(m_final, newdata = plot_data, se.fit = TRUE, exclude = "s(Site)")

plot_data$fit <- preds$fit
plot_data$se  <- preds$se.fit

# Transform back to the original scale (since it's a log link)
# use exp() because Gamma(link="log") uses natural logs
plot_data <- plot_data %>%
  mutate(
    BPc_pred = exp(fit),
    lower    = exp(fit - 1.96 * se),
    upper    = exp(fit + 1.96 * se)
  )

# re-order substrate for plotting
unique(plot_data$substrate)
plot_data$substrate <- factor(plot_data$substrate, 
                              levels = c("Mud", "sandy Mud", "muddy Sand", "Sand", 
                                         "Coarse-grained sediment", "Mixed sediment",
                                         "Rock & boulders"),
                              labels = c("Mud", "Sandy mud", "Muddy sand", "Sand", 
                                         "Coarse-grained sediment", "Mixed sediment",
                                         "Rock & boulders"))

# Grab 9 colors from the palette but skip the middle (white) one
my_colors <- brewer.pal(11, "BrBG")[c(1, 2, 3, 4, 8, 9, 10, 11)]

# Plot GAMM temporal trend by substrate
(GAM_substrate_plot <-
    ggplot(plot_data, aes(x = Year, y = BPc_pred, 
                          color = substrate, 
                          fill = substrate)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
    geom_line(size = 1) +
    facet_wrap(~ substrate, scales = "free_y", ncol = 2) + 
    scale_color_manual(values = my_colors) +
    scale_fill_manual(values = my_colors) +
    scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k")) +
    theme_bw() +
    labs(
      x = "Year",
      y = expression(BP[c])
    ) +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) 
)

# Save as Figure 2
ggsave(GAM_substrate_plot, filename = "outputs/Figure_2.png", width = 6, height = 7, dpi = 600)

# Site-level analysis ----

# fit one regression per site to look at spatial differences in temporal trends in BPc

# use gamma glm with log-link
# better suited for right-skewed positive data

site_slopes <- hauls_subset %>%
  group_by(Site) %>%
  # Using GLM with log link for exponential rates
  do(tidy(glm(BPc_capped ~ Year, 
              family = Gamma(link = "log"), 
              data = .))) %>%
  filter(term == "Year") %>%
  select(Site, estimate, std.error, p.value) %>%
  mutate(sig = ifelse(p.value < 0.05, "significant", "not significant"),
         abs_estimate = abs(estimate),
         trajectory = case_when(
           p.value < 0.05 & estimate > 0 ~ "Increasing",
           p.value < 0.05 & estimate < 0 ~ "Decreasing",
           TRUE ~ "Stable"
         )) %>%
  left_join(sites, by = "Site") %>%
  filter(Site != "36E4_408")

summary(site_slopes)

table(site_slopes$p.value < 0.05) # 36 sig

# back-transform slopes to get annual rate of change
site_slopes <- site_slopes %>%
  mutate(annual_change = (exp(estimate) - 1) * 100)

# Map annual change across Irish Sea

# Get coastline

coast <- ne_countries(scale = "medium", returnclass = "sf")

sites_sf <- st_as_sf(site_slopes, coords = c("Long", "Lat"), crs = 4326)

# remove one rogue point on land
sites_sf <- sites_sf %>% 
  filter(Site != "36E4_408")

(slope_map_final <- ggplot() +
    geom_sf(data = coast, color = "black", fill = "lightgray") + 
    geom_sf(data = sites_sf %>% arrange(sig), 
            aes(color = annual_change,  # Color = Direction & Strength
                shape = sig),           # Shape = Significance
            alpha = 0.9,
            size = 5) +             
    scale_color_gradient2(low = "#c1121f", mid = "white", high = "#669bbc", 
                          midpoint = 0, name = "Annual % change \nin BPc") +
    scale_shape_manual(values = c(1, 16), name = "Significance (p < 0.05)") +
    coord_sf(xlim = c(-8, -3), ylim = c(51.5, 55)) +
    theme_bw() +
    labs(x = "Longitude", y = "Latitude") +
    theme(legend.position = "right",
          legend.box = "vertical") + 
    guides(size = guide_legend(order = 1), 
           color = guide_colorbar(order = 2), 
           shape = guide_legend(order = 3))
)

# Save as figure 3
ggsave(slope_map_final, filename = "outputs/Figure_3.png",
       width = 170, height = 140, units = "mm", dpi = 600)

# Multivariate analysis ----

## Prep data ----

# Subset raw data to subset sites

hauls_spp_sub <- hauls_clean %>% 
  filter(Site %in% keep_sites$Site)

length(unique(hauls_spp_sub$Site)) # 74 sites
unique(hauls_spp_sub$species) # 93, so lose one species by subsetting

# How many species per mode?

hauls_spp_sub %>% group_by(bioturbation.mode) %>% summarise(length(unique(species)))

# Burrower                                   8
# Lateral excavator                         18
# Nest-builder                               9
# Sediment sifter                           29
# Vertical excavator                        29

# Add site slopes info to raw data so I know which ones are sig.

hauls_spp_sub <- hauls_spp_sub %>% 
  left_join(site_slopes, by = "Site")

# look at composition of increasing vs decreasing sites

sig_sites_spp <- hauls_spp_sub %>% 
  filter(sig == "significant") %>% 
  mutate(direction = if_else(estimate > 0, "increase", "decrease"))

head(sig_sites_spp)

# this is the key dataframe of sites showing significant increase or decrease over time

# aggregate species data across years per site
# since direction groups are already a temporal signal

site_comm <- sig_sites_spp %>%
  group_by(Site, species, bioturbation.mode, direction, trajectory) %>%
  summarise(Depth = mean(Depth),
            Substrate = first(substrate_name),
            abun = mean(total_abun_cpue), .groups = "drop")

# convert to wide format
site_wide <- site_comm %>% 
  select(Site, species, abun) %>%
  pivot_wider(names_from = species,
              values_from = abun,
              values_fill = 0)

# Metadata
meta_site <- site_comm %>%
  group_by(Site) %>%
  summarise(
    direction = first(direction),
    trajectory = first(trajectory),
    Depth = mean(Depth),
    Substrate = first(Substrate)
  )

meta_site %>% group_by(trajectory) %>% tally()

## PERMANOVA ----

# Tests whether sites with increasing vs decreasing bioturbation differ in 
# species composition, after accounting for environmental variation.

# Extract species matrix
species_matrix <- site_wide %>%
  select(-Site) %>%
  as.data.frame()

# Square-root transform
species_sqrt <- sqrt(species_matrix)

set.seed(8)

# PERMANOVA
perm <- adonis2(
  species_sqrt ~ trajectory + Depth + Substrate,
  data = meta_site,
  method = "bray",
  by = "margin",
  permutations = 999
)

perm

### dbRDA ----

# Visualise increasing/decreasing sites with env variables

species_bray <- vegdist(species_sqrt, method = "bray")

dbRDA <- capscale(
  species_sqrt ~ trajectory + Depth + Substrate,
  data = meta_site
)

anova(dbRDA, by = "margin")

# Extract site scores (samples)
site_scores <- scores(dbRDA, display = "sites")
site_df <- as.data.frame(site_scores)
site_df$trajectory <- meta_site$trajectory 

# Extract environmental variables (biplot arrows)
env_vectors <- scores(dbRDA, display = "bp")
rownames(env_vectors)

# Clean biplot scores: keep only continuous variables
env_vectors <- env_vectors[!rownames(env_vectors) %in% c("trajectoryIncreasing"), ]
# Rename labels
rownames(env_vectors)[rownames(env_vectors) == "SubstrateMixed sediment"] <- "Mixed sediment"
rownames(env_vectors)[rownames(env_vectors) == "Substratemuddy Sand"] <- "Muddy sand"
rownames(env_vectors)[rownames(env_vectors) == "SubstrateRock & boulders"] <- "Rock & boulders"
rownames(env_vectors)[rownames(env_vectors) == "SubstrateSand"] <- "Sand"
rownames(env_vectors)[rownames(env_vectors) ==  "Substratesandy Mud"] <- "Sandy mud"

env_df <- as.data.frame(env_vectors)
env_df$var <- rownames(env_df)

# Scale arrows
env_df_scaled <- env_df
scaling_factor <- 3  # adjust as needed
env_df_scaled$CAP1 <- env_df$CAP1 * scaling_factor
env_df_scaled$CAP2 <- env_df$CAP2 * scaling_factor

# Extract species scores
species_scores <- scores(dbRDA, display = "species")
species_df <- as.data.frame(species_scores)
species_df$species <- rownames(species_df)

# Percent variance for axis labels
ve <- summary(dbRDA)$cont$importance
cap1_pct <- round(ve[2,1]*100, 1)
cap2_pct <- round(ve[2,2]*100, 1)

(dbRDA_plot_site <- ggplot(data = site_df, 
                           aes(CAP1, CAP2)) +
    geom_point(aes(CAP1, CAP2, color = trajectory), 
               size = 4, alpha = 0.8) +
    #scale_shape_manual(values = c(16, 17, 18)) +
    geom_segment(data = env_df_scaled, aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
                 arrow = arrow(length = unit(0.3, "cm")), 
                 color = "black") +
    geom_text_repel(data = env_df_scaled, aes(x = CAP1, y = CAP2, label = var), 
                    size = 5, segment.color = "black", fontface = "italic",
                    point.padding = 0.2 ) +
    #add_phylopic(x = -1.1, y = -1, img = fish_img, height = 0.3) +
    geom_hline(yintercept= 0, linetype = "dashed", colour = "grey") +
    geom_vline(xintercept= 0, linetype = "dashed", colour = "grey") +
    scale_fill_manual(values = c("#c1121f",  "#669bbc")) +
    scale_colour_manual(values = c("#c1121f",  "#669bbc")) +
    labs(x = paste0("CAP1 (", cap1_pct, "%)"),
         y = paste0("CAP2 (", cap2_pct, "%)"),
         colour = "Site trajectory") +
    theme_bw() +
    theme(legend.position = c(0.15, 0.9),
          legend.frame = element_rect(),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 20),
          legend.text = element_text(size = 16),          
          legend.title = element_text(size = 16, face = "bold"),
          legend.key.size = unit(1.2,"lines"),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank())
)

# Save as figure 4
ggsave(dbRDA_plot_site, filename = "outputs/Figure_4.png", 
       width = 8, height = 5.5, units = "in", dpi = 600)

# Functional contributions----

# plot relative contribution of each bioturbation mode to total BPc at both site trajectories

# Annual
mode_contribution_annual <- hauls_spp_sub %>%
  # Sum BPp by Year and Mode
  group_by(Year, bioturbation.mode) %>%
  summarize(Mode_Total_BP = sum(BPp, na.rm = TRUE), .groups = "drop") %>%
  # Calculate % within each Year
  group_by(Year) %>%
  mutate(Percentage = (Mode_Total_BP / sum(Mode_Total_BP)) * 100) %>%
  ungroup()

# Quick check: Do the percentages in 1990 sum to 100?
mode_contribution_annual %>% filter(Year == 1990) %>% summarize(check = sum(Percentage))

mode_contribution_global <- hauls_spp_sub %>%
  # Sum everything by Mode across the whole dataset
  group_by(bioturbation.mode) %>%
  summarize(Total_Mode_BP = sum(BPp, na.rm = TRUE)) %>%
  # Calculate % of the grand total
  mutate(Percentage = (Total_Mode_BP / sum(Total_Mode_BP)) * 100) %>%
  arrange(desc(Percentage))

print(mode_contribution_global)

## Mode area plot ----

head(hauls_spp_sub)

# Prepare the data: Join trajectories and calculate % per mode per year
mode_plot_data <- hauls_spp_sub %>%
  #inner_join(site_slopes %>% dplyr::select(Site, trajectory), by = "Site") %>%
  filter(trajectory %in% c("Increasing", "Decreasing")) %>%
  group_by(Year, trajectory, bioturbation.mode) %>%
  summarize(Total_Mode_BP = sum(BPp, na.rm = TRUE), .groups = "drop") %>%
  # Fill in missing Year/Mode combinations with 0 
  complete(Year, trajectory, bioturbation.mode, fill = list(Total_Mode_BP = 0)) %>%
  group_by(Year, trajectory) %>%
  mutate(Percentage = (Total_Mode_BP / sum(Total_Mode_BP)) * 100) %>%
  # Replace any NaNs (from dividing by zero if a whole year is empty) with 0
  mutate(Percentage = ifelse(is.nan(Percentage), 0, Percentage)) %>% 
  ungroup()

# Reorder modes by their average contribution (biggest at the bottom)
mode_order <- mode_plot_data %>%
  group_by(bioturbation.mode) %>%
  summarize(avg = mean(Percentage)) %>%
  arrange(desc(avg)) %>%
  pull(bioturbation.mode)

mode_plot_data$bioturbation.mode <- factor(mode_plot_data$bioturbation.mode, 
                                           levels = c("Burrower", "Vertical excavator", "Nest-builder", 
                                                      "Lateral excavator", "Sediment sifter"))

# Create the Plot
(trajectory_mode_plot <-
    ggplot(mode_plot_data, aes(x = Year, y = Percentage, fill = bioturbation.mode)) +
    geom_area(alpha = 0.7, color = "white", linewidth = 0.1) +
    facet_wrap(~ trajectory) +
    scale_fill_viridis_d() +
    scale_y_continuous(expand = c(0,0.01)) +
    scale_x_continuous(expand = c(0,0.01)) +
    labs(
      x = "Year",
      y = "Contribution to \ntotal BPc (%)",
      fill = "Bioturbation mode"
    ) +
    theme_convex() +
    theme(
      #panel.spacing = unit(2, "lines"),
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(color="black", fill="white", alpha = 0.5, size=1),
      plot.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.box = "vertical", 
      legend.margin = margin(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.6),
      panel.background = element_rect(fill = "white", colour = NA)) 
)

## Species contributions ----

### Vertical excavators ----

target_mode <- "Vertical excavator"

# Join trajectories to main data first
spp_with_traj <- hauls_spp_sub %>%
  #inner_join(site_slopes %>% select(Site, trajectory), by = "Site") %>%
  filter(trajectory %in% c("Increasing", "Decreasing"))

# Get the top species for this mode per trajectory
top_species_by_traj <- spp_with_traj %>%
  filter(bioturbation.mode == target_mode) %>%
  group_by(trajectory, species) %>%
  summarize(total_contribution = sum(BPp, na.rm = TRUE), .groups = "drop") %>%
  group_by(trajectory) %>%
  slice_max(total_contribution, n = 5) %>% 
  mutate(is_top = TRUE)

# Prepare the plotting data
mode_zoom_traj <- spp_with_traj %>%
  filter(bioturbation.mode == target_mode) %>%
  left_join(top_species_by_traj %>% select(trajectory, species, is_top), 
            by = c("trajectory", "species")) %>%
  # Label anything not in the top 5 for THAT trajectory as "Other"
  mutate(Species_Label = ifelse(is.na(is_top), "Other", species)) %>%
  group_by(Year, trajectory, Species_Label) %>%
  summarize(Mode_BPc = sum(BPp, na.rm = TRUE), .groups = "drop") %>%
  complete(Year, trajectory, Species_Label, fill = list(Mode_BPc = 0)) %>%
  # Calculate Percentage
  group_by(Year, trajectory) %>%
  mutate(Total_Year_BP = sum(Mode_BPc, na.rm = TRUE)) %>%
  # Only calculate if there is actually data for that year/traj
  mutate(Rel_Contrib = ifelse(Total_Year_BP > 0, (Mode_BPc / Total_Year_BP) * 100, 0)) %>%
  ungroup()

# Add common names for labelling

# extract from biot dataframe
common <- biot %>% 
  select(Species, Common.name)

mode_zoom_traj <- mode_zoom_traj %>% 
  left_join(common, by = c("Species_Label" = "Species"))

mode_zoom_traj <- mode_zoom_traj %>%
  mutate(Common.name = replace_na(Common.name, "Other"), 
         Category = case_when(Common.name == "Lesser spotted dogfish/smallspotted catshark" ~ "Small-spotted catshark",
                              TRUE ~ Common.name))

# Calculate the total contribution of each species label to define the order
label_order <- mode_zoom_traj %>%
  group_by(Category) %>%
  summarize(total = sum(Mode_BPc, na.rm = TRUE)) %>%
  filter(Category != "Other") %>%
  arrange(desc(total)) %>%
  pull(Category)

# Reconstruct the factor levels: "Other" first (bottom), then species in order
mode_zoom_traj <- mode_zoom_traj %>%
  mutate(Category = factor(Category, levels = c(label_order, "Other")))

# Create custom colour palette:

# Get the levels exactly as they are in data
levels_in_order <- levels(mode_zoom_traj$Category)

# Define how many species need a color (Total minus other because other are grey)
num_actual_species <- length(levels_in_order) - 1

# Create the Blue Palette from RColourBrewer
blue_palette <- colorRampPalette(brewer.pal(7, "Blues"))(num_actual_species)

# Create the named vector
# "Other" gets a neutral grey, then the rest get the blues in order
spec_colors <- setNames(
  c(blue_palette, "#D3D3D3"), # The colors
  levels_in_order             # The names they attach to
)

# Vertical excavator plot
(ve_traj_plot <-
    ggplot(mode_zoom_traj, aes(x = Year, y = Rel_Contrib, fill = Category)) +
    geom_area(alpha = 0.8, 
              color = "white", 
              linewidth = 0.1, position = "stack") + 
    facet_wrap(~ trajectory) +
    scale_fill_manual(values = spec_colors) +
    scale_y_continuous(expand = c(0,0.01)) +
    scale_x_continuous(expand = c(0,0.01)) +
    labs(y = "Contribution to \ntotal BPc (%)",
         fill = "Vertical excavators") +
    theme_convex() +
    theme(
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(color="black", fill="white", alpha = 0.5, size=1),
      plot.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.box = "vertical", 
      legend.margin = margin(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.6
      ),
      panel.background = element_rect(fill = "white", colour = NA))
)

### Lateral excavators ----

target_mode <- "Lateral excavator"

# Join trajectories to main data first
spp_with_traj <- hauls_spp_sub %>%
  #inner_join(site_slopes %>% select(Site, trajectory), by = "Site") %>%
  filter(trajectory %in% c("Increasing", "Decreasing"))

# Get the top species for this mode per trajectory
top_species_by_traj <- spp_with_traj %>%
  filter(bioturbation.mode == target_mode) %>%
  group_by(trajectory, species) %>%
  summarize(total_contribution = sum(BPp, na.rm = TRUE), .groups = "drop") %>%
  group_by(trajectory) %>%
  slice_max(total_contribution, n = 5) %>% 
  mutate(is_top = TRUE)

# Prepare the plotting data
mode_zoom_traj <- spp_with_traj %>%
  filter(bioturbation.mode == target_mode) %>%
  left_join(top_species_by_traj %>% select(trajectory, species, is_top), 
            by = c("trajectory", "species")) %>%
  # Label anything not in the top 5 for THAT trajectory as "Other"
  mutate(Species_Label = ifelse(is.na(is_top), "Other", species)) %>%
  group_by(Year, trajectory, Species_Label) %>%
  summarize(Mode_BPc = sum(BPp, na.rm = TRUE), .groups = "drop") %>%
  complete(Year, trajectory, Species_Label, fill = list(Mode_BPc = 0)) %>%
  # Calculate Percentage
  group_by(Year, trajectory) %>%
  mutate(Total_Year_BP = sum(Mode_BPc, na.rm = TRUE)) %>%
  # Only calculate if there is actually data for that year/traj
  mutate(Rel_Contrib = ifelse(Total_Year_BP > 0, (Mode_BPc / Total_Year_BP) * 100, 0)) %>%
  ungroup()

# Add common names for labelling
mode_zoom_traj <- mode_zoom_traj %>% 
  left_join(common, by = c("Species_Label" = "Species"))

mode_zoom_traj <- mode_zoom_traj %>%
  mutate(Common.name = replace_na(Common.name, "Other"), 
         Category = case_when(#Common.name == "Lesser spotted dogfish/smallspotted catshark" ~ "Small-spotted catshark",
           TRUE ~ Common.name))

# Calculate the total contribution of each species label to define the order
label_order <- mode_zoom_traj %>%
  group_by(Category) %>%
  summarize(total = sum(Mode_BPc, na.rm = TRUE)) %>%
  filter(Category != "Other") %>%
  arrange(desc(total)) %>%
  pull(Category)

# Reconstruct the factor levels: "Other" first (bottom), then species in order
mode_zoom_traj <- mode_zoom_traj %>%
  mutate(Category = factor(Category, levels = c(label_order, "Other")))

# Create custom colour palette

# Get the levels exactly as they are in data
levels_in_order <- levels(mode_zoom_traj$Category)

# Define how many species need a green color (Total minus 'Other')
num_actual_species <- length(levels_in_order) - 1

# Create the Green Palette 
green_palette <- colorRampPalette(brewer.pal(9, "Greens"))(num_actual_species)

# Create the named vector
# "Other" gets a neutral grey, then the rest get the greens in order
spec_colors <- setNames(
  c(green_palette, "#D3D3D3"), # The colors
  levels_in_order             # The names they attach to
)

# Lateral excavator plot
(le_traj_plot <-
    ggplot(mode_zoom_traj, aes(x = Year, y = Rel_Contrib, fill = Category)) +
    geom_area(alpha = 0.8, 
              color = "white", 
              linewidth = 0.1, position = "stack") + 
    facet_wrap(~ trajectory) +
    scale_fill_manual(values = spec_colors) +
    scale_y_continuous(expand = c(0,0.01)) +
    scale_x_continuous(expand = c(0,0.01)) +
    labs(y = "Contribution to \ntotal BPc (%)",
         fill = "Lateral excavators") +
    theme_convex() +
    theme(
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(color="black", fill="white", alpha = 0.5, size=1),
      plot.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.box = "vertical", 
      legend.margin = margin(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.6),
      panel.background = element_rect(fill = "white", colour = NA))
)

### Multi-panel area plot ----

p1 <- trajectory_mode_plot + theme(plot.margin = unit(c(0.5, 0.5, 0.1, 0.5), units = , "cm"),
                                   legend.margin = margin(b = 20))
p1

p2 <- ve_traj_plot + theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.5), units = , "cm"),
                           legend.margin = margin(t = 20, b = 20))
p2

p3 <- le_traj_plot + theme(plot.margin = unit(c(0.5, 0.1, 0.5, 0.5), units = , "cm"),
                           legend.margin = margin(t = 20))
p3

area_panel <- p1 / p2 / p3 + plot_annotation(tag_levels = "A") + plot_layout(axes = "collect_x") + plot_layout(guides = "collect") 
area_panel

# Save as figure 5
ggsave(area_panel, filename = "outputs/Figure_5.png", width = 8, height = 8, dpi = 600)

# Supplementary materials ----

## Table S1 ----

## Species summary table

# for supplementary material:
# including mode, how many sites they are present,
# mean body size (across all sites), mean abundance across all sites,
# mean BPi, mean BPp

# first, subset raw biot data to only my final subset of sites
# subsetting is done lower down the script...

biot_spp <- biot_final %>% 
  filter(Site %in% keep_sites$Site)

length(unique(biot_spp$Site)) # 74
length(unique(biot_spp$species)) # 93
# all looks good

spp_table <- biot_spp %>% 
  group_by(species, Common.name, Family, Bioturbation.mode) %>% 
  summarise(n = length(unique(Site)),
            B = mean(weight.g),
            A = mean(CPUE),
            BPi = mean(BPi))

# Add mean BPp
BPp <- biot_haul %>% 
  filter(Site %in% keep_sites$Site)

length(unique(BPp$Site))

BPp2 <- BPp %>% 
  group_by(species) %>% 
  summarise(BPp = mean(BPp))

spp_table <- spp_table %>% 
  left_join(BPp2, by = "species")

write.csv(spp_table, "outputs/Table_S1.csv")

# Heatmaps (Figure S2) ----

# for supplementary materials

# Species abundance/biomass by trajectory over time

# this will tell me whether any species have increased/decreased at certain sites over time

spp_agg <- hauls_spp_sub %>% 
  filter(trajectory %in% c("Increasing", "Decreasing")) %>% 
  group_by(Year, trajectory, species) %>% 
  summarise(mean_abun = mean(total_abun_cpue, na.rm = TRUE),
            mean_biom = mean(total_biomass_cpue, na.rm = TRUE),
            .groups = "drop")

# Identify Top 22 species overall
top_species <- spp_agg %>%
  group_by(species) %>%
  summarize(total = sum(mean_abun)) %>%
  slice_max(total, n = 22) %>%
  pull(species)

# Filter and Prepare for Heatmap
heatmap_data <- spp_agg %>%
  filter(species %in% top_species) %>%
  # Log transform to see patterns across all species
  mutate(log_abun = log10(mean_abun)) %>% 
  left_join(common, by = c("species" = "Species")) %>% 
  mutate(Category = case_when(Common.name == "Common dab" ~ "Dab",
                              Common.name == "Lesser spotted dogfish/smallspotted catshark" ~ "Small-spotted catshark",
                              TRUE ~ Common.name))

# Reorder species so the most abundant are at the top
species_order <- heatmap_data %>%
  group_by(Category) %>%
  summarize(m = mean(log_abun)) %>%
  arrange(m) %>%
  pull(Category)

heatmap_data$species <- factor(heatmap_data$species, levels = species_order)
heatmap_data$Category <- factor(heatmap_data$Category, levels = species_order)

# abundance heatmap
(abun_heatmap <-
    ggplot(heatmap_data, aes(x = Year, y = Category, fill = log_abun)) +
    geom_tile(#color = "white", 
      size = 0.1) +
    facet_grid(. ~ trajectory) + # Side-by-side comparison
    scale_fill_viridis_c(option = "plasma", name = "log10 Abundance ") +
    labs(x = "Year", y = NULL) +
    scale_x_continuous(expand = c(0,0.01)) +
  theme_convex() +
    theme(
      panel.spacing = unit(1, "lines"),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank(),
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(color="black", fill="white", alpha = 0.5, size=1),
      plot.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.box = "vertical", 
      legend.margin = margin(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.6
      ),
      panel.background = element_rect(fill = "white", colour = NA))
)

## Biomass heatmap

# Identify Top 22 species overall
top_species <- spp_agg %>%
  group_by(species) %>%
  summarize(total = sum(mean_biom)) %>%
  slice_max(total, n = 22) %>%
  pull(species)

# Filter and Prepare for Heatmap
heatmap_data <- spp_agg %>%
  filter(species %in% top_species) %>%
  # Log transform to see patterns across all species
  mutate(log_biom = log10(mean_biom)) %>% 
  left_join(common, by = c("species" = "Species")) %>% 
  mutate(Category = case_when(Common.name == "Common dab" ~ "Dab",
                              Common.name == "Lesser spotted dogfish/smallspotted catshark" ~ "Small-spotted catshark",
                              TRUE ~ Common.name))

# Reorder species so the most abundant are at the top (better than alphabetical)
species_order <- heatmap_data %>%
  group_by(Category) %>%
  summarize(m = mean(log_biom)) %>%
  arrange(m) %>%
  pull(Category)

heatmap_data$species <- factor(heatmap_data$species, levels = species_order)
heatmap_data$Category <- factor(heatmap_data$Category, levels = species_order)

(biom_heatmap <-
    ggplot(heatmap_data, aes(x = Year, y = Category, fill = log_biom)) +
    geom_tile(size = 0.1) +
    facet_grid(. ~ trajectory) +
    scale_fill_viridis_c(option = "plasma", name = "log10 Biomass ") +
    labs(x = "Year", y = NULL) +
    scale_x_continuous(expand = c(0, 0.01)) +
    theme_convex() +
    theme(
      panel.spacing = unit(1, "lines"),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank(),
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(color="black", fill="white", alpha = 0.5, size=1),
      plot.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.box = "vertical", 
      legend.margin = margin(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.6
      ),
      panel.background = element_rect(fill = "white", colour = NA))
)

# Combine: 
abun_heatmap <- abun_heatmap + theme(plot.margin = unit(c(0.5, 0.5, 0.1, 0.5), units = , "cm"))

biom_heatmap <- biom_heatmap + theme(plot.margin = unit(c(0.5, 0.1, 0.5, 0.5), units = , "cm"))

heatmaps <- abun_heatmap/biom_heatmap + plot_annotation(tag_levels = "A")
#+ plot_layout(guides = "collect") 
heatmaps

# Save as Fig. S2
ggsave(heatmaps, filename = "outputs/Figure_S2.png", width = 7, height = 8, dpi = 600)

# Landings plots ----

# ICES provide catch data on their website:
## historical data 1903-1949
## historical records 1950-2010
## recent catch data 2006-2023

# Load raw landings data
landings1903_1949 <- read.csv("raw_data/Landings_data/1903-1949_Landings.csv")
landings1950_2010 <- read.csv("raw_data/Landings_data/ICES_1950-2010.csv")
landings2006_2023 <- read.csv("raw_data/Landings_data/ICESCatchDataset2006-2023.csv")

## Landings data tidying ----

### 1903-1949 ----

head(landings1903_1949, 10) # already in long format
str(landings1903_1949)

unique(landings1903_1949$FAO_Area)
# VIIa is the Irish Sea

# subset to Irish Sea (VIIa)
landings1903_1949 <- landings1903_1949 %>% 
  filter(FAO_Area == "VIIa")

length(unique(landings1903_1949$Species.scientific.name)) # 51 species

# Sort out missing values
## replace empty cells with 0s
## p means < half a ton, use half point of 0-0.5 t (same approach as for 1950-2006 data)

landings1903_1949 <- landings1903_1949 %>%
  mutate(Landings = case_when(
    is.na(Landings) | Landings == "" ~ "0",     # replace empty cells with 0
    Landings == "P"            ~ "0.25",        # midpoint of 0–0.5 t
    TRUE                  ~ Landings
  ),
  Landings = as.numeric(Landings))           

summary(landings1903_1949)

# Sum catches landed by different countries by year and species
landings1903_1949_sum <- landings1903_1949 %>% 
  group_by(Year, FAO_Species_Code, FAO_Species_Name, Species.scientific.name) %>% 
  summarise(Landings = sum(Landings, na.rm = T)) %>% 
  rename(Species = FAO_Species_Name)

summary(landings1903_1949_sum)

unique(landings1903_1949_sum$Species)

### 1950-2010 -----

# Sort out missing values
landings1950_2010_clean <- landings1950_2010 %>%
  mutate(across(4:64, ~ case_when(
    . == "-" ~ 0,                        # replace '-' with 0
    . == "." | . == "..." ~ NA_real_,    # replace '.' or '...' with NA
    . == "<0.5" ~ 0.25,                  # midpoint of 0–0.5 t
    TRUE ~ as.numeric(.))))              # convert numeric-looking cells to numbers

head(landings1950_2010_clean)

# Remove Xs
landings1950_2010_clean <- landings1950_2010_clean %>%
  rename_with(~ sub("^X", "", .x))

# Convert to long format
landings1950_2010_long <- landings1950_2010_clean %>%
  pivot_longer(
    cols = `1950`:`2010`,      # all the year columns
    names_to = "Year",         # new column name for the years
    values_to = "Landings"     # new column name for the landings
  ) %>%
  mutate(Year = as.integer(Year))  # convert year from text to number

# Remove missing catches
landings1950_2010_long <- landings1950_2010_long %>% 
  filter(!is.na(Landings))

# Subset to Irish Sea

unique(landings1950_2010_long$Division)

landings1950_2010_is <- landings1950_2010_long %>% 
  filter(Division == "VII a")

# Explore
length(unique(landings1950_2010_is$Species)) # 193 species
unique(landings1950_2010_is$Species)
length(unique(landings1950_2010_is$Country)) # 15 countries

## How do deal with these species - some are just groups and all common names?
## Need to find a way to match them with other landings data

# Sum landings across countries for each year and each species
landings1950_2010_sum <- landings1950_2010_is %>% 
  group_by(Year, Species) %>% 
  summarise(Landings = sum(Landings, na.rm = T))

unique(landings1950_2010_sum$Species)

### 2006-2023 ----

# Remove Xs
landings2006_2023 <- landings2006_2023 %>%
  rename_with(~ sub("^X", "", .x))

length(unique(landings2006_2023$Area)) # 131
unique(landings2006_2023$Area)

# Subset to Irish Sea
landings2006_2023_is <- landings2006_2023 %>% 
  filter(Area == "27.7.a")

# Convert 0 c to NAs
landings2006_2023_is <- landings2006_2023_is %>%
  mutate(across(5:22, ~ case_when(
    . == "0 c" ~ NA_real_,               # replace '0 c' with NA
    TRUE ~ as.numeric(.))))              # convert numeric-looking cells to numbers

# Convert to long format
landings2006_2023_long <- landings2006_2023_is %>%
  pivot_longer(
    cols = 5:22,           # all the year columns
    names_to = "Year",     # new column name for the years
    values_to = "Landings" # new column name for the landings
  ) %>%
  mutate(Year = as.integer(Year))  # convert year from text to number

# Check
unique(landings2006_2023_long$Units) # all TLW
unique(landings2006_2023_long$Area) # all 27.7.a

str(landings2006_2023_long)
summary(landings2006_2023_long)

# Sum across countries
landings2006_2023_sum <- landings2006_2023_long %>% 
  group_by(Year, Species) %>% 
  summarise(Landings = sum(Landings, na.rm = T))

summary(landings2006_2023_sum)

# Add species names
landings2006_2023_final <- landings2006_2023_sum %>% 
  left_join(spp, by = c("Species" = "FAO_code")) %>% 
  dplyr::select(- French, - Spanish) %>% 
  rename(FAO_Species_Code = Species,
         Species = English,
         Species.scientific.name = Latin.species.name)

length(unique(landings2006_2023_final$Species.scientific.name)) # 208 spp
unique(landings2006_2023_final$Species)

head(landings1950_2010_sum)
head(landings2006_2023_final)

# How to merge the three datasets with different species names?
unique(landings1903_1949_sum$Species)
unique(landings1950_2010_sum$Species)
unique(landings2006_2023_final$Species)

### Join landings data ----

# joining all landings data 1903-2023

unique(landings1903_1949_sum$Species)   # 51
unique(landings1950_2010_sum$Species)   # 193
unique(landings2006_2023_final$Species) # 211

# extract all unique species names across datasets
all_species <- unique(c(
  landings1903_1949_sum$Species,
  landings1950_2010_sum$Species,
  landings2006_2023_final$Species
))
# n = 253

# add raw species names to a master species key
species_key <- tibble(
  Species_raw = all_species,    # original name in dataset
  Species_std = NA_character_,  # standardized common name
  Latin_name  = NA_character_,  # Latin name if available
  FAO_group   = NA_character_,  # e.g. "Ray", "Flatfish", "Shark"
  Resolution  = NA_character_   # "species", "nei", "group"
)

# Classify resolution - species level or nei?
species_key <- species_key %>%
  mutate(
    Species_std = Species_raw,
    Resolution  = "species"
  )

species_key <- species_key %>%
  mutate(
    Species_std = if_else(
      str_detect(Species_raw, " nei"),
      Species_raw,
      Species_std
    ),
    Resolution = if_else(
      str_detect(Species_raw, " nei"),
      "nei",
      Resolution
    )
  )

table(species_key$Resolution)
# 92 nei
# 161 species

# Add Latin names where available

latin_lookup <- bind_rows(
  landings1903_1949_sum %>%
    ungroup() %>% 
    select(Species, Species.scientific.name),
  landings2006_2023_final %>%
    ungroup() %>% 
    select(Species, Species.scientific.name)
) %>%
  distinct() %>%
  filter(!is.na(Species.scientific.name))

# Join to species key
species_key <- species_key %>%
  left_join(
    latin_lookup,
    by = c("Species_raw" = "Species")
  ) %>% 
  select(- Latin_name)

species_key %>%
  summarise(
    total = n(),
    with_latin = sum(!is.na(Species.scientific.name))
  )

# 215 out of 253 have Latin names

head(species_key)

unique(species_key$Species_raw)
unique(species_key$Species.scientific.name)

# Fix spelling variants

species_key <- species_key %>%
  mutate(
    Species_std = case_when(
      Species_raw %in% c("Thickback sole", "Thickback soles") ~ "Thickback sole",
      Species_raw %in% c("Dogfish sharks nei", "Dogfish sharks, etc. nei") ~ "Dogfish sharks nei",
      Species_raw %in% c("Cuttlefish,bobtail squids nei", "Cuttlefish, bobtail squids nei") ~
        "Cuttlefish, bobtail squids nei",
      TRUE ~ Species_std
    )
  )

species_key <- species_key %>%
  mutate(
    Species_std = case_when(
      # Singular / plural
      Species_raw %in% c("Thickback sole", "Thickback soles") ~ "Thickback sole",
      Species_raw %in% c("Thickback sole nei", "Thickback soles nei") ~ "Thickback sole nei",
      Species_raw %in% c("Spiny lobster nei", "Spiny lobsters nei") ~ "Spiny lobsters nei",
      Species_raw %in% c("Flat oyster nei", "Flat oysters nei") ~ "Flat oysters nei",
      Species_raw %in% c("Sea mussel nei", "Sea mussels nei") ~ "Sea mussels nei",
      Species_raw %in% c("Cockle nei", "Cockles nei") ~ "Cockles nei",
      Species_raw %in% c("Sole nei", "Soles nei") ~ "Soles nei",
      # Punctuation
      Species_raw %in% c(
        "Cuttlefish,bobtail squids nei",
        "Cuttlefish, bobtail squids nei"
      ) ~ "Cuttlefish, bobtail squids nei",
      Species_raw %in% c(
        "Houndsharks,smoothhounds nei",
        "Houndsharks, smoothhounds nei"
      ) ~ "Houndsharks, smoothhounds nei",
      Species_raw %in% c(
        "Triggerfishes,durgons nei",
        "Triggerfishes, durgons nei"
      ) ~ "Triggerfishes, durgons nei",
      # Truncated wording
      Species_raw %in% c(
        "Dogfish sharks nei",
        "Dogfish sharks, etc. nei"
      ) ~ "Dogfish sharks nei",
      Species_raw %in% c(
        "Catsharks, etc. nei",
        "Catsharks, nursehounds nei"
      ) ~ "Catsharks, nursehounds nei",
      Species_raw %in% c(
        "Razor clams nei",
        "Razor clams, knife clams nei"
      ) ~ "Razor clams nei",
      Species_raw %in% c(
        "Squids nei",
        "Various squids nei"
      ) ~ "Squids nei",
      # Hyphen / spacing 
      Species_raw %in% c(
        "Thicklip grey mullet",
        "Thick-lip grey mullet"
      ) ~ "Thicklip grey mullet",
      Species_raw %in% c(
        "Blackspot (=red) seabream",
        "Blackspot(=red) seabream"
      ) ~ "Blackspot(=red) seabream",
      # Default 
      TRUE ~ Species_std
    )
  )

# How many standardised names?
length(unique(species_key$Species_std)) #247

species_key %>%
  count(Species_std) %>%
  filter(n > 1) %>%
  arrange(desc(n))

# Assign broader groups
species_key <- species_key %>%
  mutate(
    Group_wide = case_when(
      
      # ---- Flatfishes 
      str_detect(Species.scientific.name, 
                 "Pleuronectes|Solea|Limanda|Microstomus|Platichthys|Hippoglossus|Soleidae|Scophthalmus|Psetta|Pleuronectiformes") |
        str_detect(Species_std, "sole|flounder|plaice|turbot|brill|dab|flatfishes") ~
        "Flatfishes",
      
      # ---- Gadoids 
      str_detect(Species.scientific.name,
                 "Gadus|Melanogrammus|Merlangius|Pollachius|Trisopterus|Molva|Gadiforme|Merluccius") |
        str_detect(Species_std, "cod|haddock|whiting|pollack|saithe|ling|hake|hakes") ~
        "Gadoids",
      
      # ---- Pelagic fishes 
      str_detect(Species.scientific.name,
                 "Clupea|Sardina|Engraulis|Scomber|Thunnus|Scombridae|Scombroidei|Thunnini") |
        str_detect(Species_std, "herring|sardine|anchovy|mackerel|tuna|sprat") ~
        "Pelagic fishes",
      
      # ---- Skates and rays
      str_detect(Species.scientific.name,
                 "Raja|Rajidae|Rajiformes") |
        str_detect(Species_std, "skate|ray") ~
        "Skates and rays",
      
      # ---- Sharks 
      str_detect(Species.scientific.name,
                 "Squalus|Scyliorhinus|Carcharhinus|Lamna|Mustelus|Alopias ") |
        str_detect(Species_std, "shark|dogfish|Thresher|hammerhead") ~
        "Sharks",
      
      # ---- Crustaceans
      str_detect(Species.scientific.name,
                 "Nephrops|Homarus|Crangon|Penaeus|Cancer") |
        str_detect(Species_std,
                   "lobster|prawn|shrimp|crab|nephrops|Craylets|crustaceans|squat lobsters") ~
        "Crustaceans",
      
      # ---- Cephalopods 
      str_detect(Species.scientific.name,
                 "Loligo|Illex|Sepia|Octopus|Teuthida|Octopodidae") |
        str_detect(Species_std,
                   "squid|cuttlefish|octopus|Squids|octopuses") ~
        "Cephalopods",
      
      # ---- Bivalves 
      str_detect(Species.scientific.name,
                 "Mytilus|Ostrea|Pecten|Cerastoderma|Venus|Bivalvia") |
        str_detect(Species_std,
                   "mussel|oyster|scallop|cockle|clam|clams|scallops|shell|molluscs|Cockles") ~
        "Bivalves",
      
      # ---- Gastropods
      str_detect(Species.scientific.name,
                 "Buccinum|Littorina|Patella|Haliotis") |
        str_detect(Species_std,
                   "whelk|periwinkle|limpet") ~
        "Gastropods",
      
      # ---- Other invertebrates 
      str_detect(Species.scientific.name,
                 "Echinodermata|Holothuroidea|Invertebrata") |
        str_detect(Species_std,
                   "urchin|starfish|sea cucumber") ~
        "Other invertebrates",
      
      # ---- Seaweeds 
      str_detect(Species.scientific.name,
                 "Algae") |
        str_detect(Species_std,
                   "Seaweed|seaweeds") ~
        "Algae",
      
      # ---- Fallback 
      TRUE ~ "Other"
    )
  )

# Craylets, squat lobsters got wrongly classes as skates and rays, manually change:
species_key <- species_key %>% 
  mutate(Group_wide = if_else(Species_raw == "Craylets, squat lobsters", "Crustaceans", Group_wide))

table(species_key$Group_wide)

# Add higher level grouping

species_key <- species_key %>%
  mutate(
    Group_high = case_when(
      
      # ---- Elasmobranchs 
      Group_wide %in% c("Sharks", "Skates and rays") ~
        "Elasmobranchs",
      
      # ---- Shellfish 
      Group_wide %in% c("Crustaceans", "Bivalves", "Gastropods", "Cephalopods") ~
        "Shellfish",
      
      # ---- Pelagic fishes 
      Group_wide == "Pelagic fishes" ~
        "Pelagic fishes",
      
      # Demersal fishes 
      Group_wide %in% c("Flatfishes", "Gadoids") ~
        "Demersal fishes",
      
      # Fallback
      TRUE ~ "Other"
    )
  )

species_key %>%
  count(Group_high) %>%
  arrange(desc(n))

species_key <- species_key %>% 
  select(- FAO_group) %>% 
  rename(Latin.name = Species.scientific.name)

# Use the species key to join the three datasets together

# Join key to each dataset

## 1903-1949

landings1 <- landings1903_1949_sum %>%
  left_join(
    species_key,
    by = c("Species" = "Species_raw")
  )

# check all species are matched
landings1 %>%
  filter(is.na(Species_std)) %>%
  distinct(Species)

table(landings1$Species == landings1$Species_std) # 1 false

# Replace original species column with standardised one
landings1 <- landings1 %>% 
  mutate(Species = Species_std)

table(landings1$Species == landings1$Species_std) # all true now
table(landings1$Latin.name == landings1$Species.scientific.name) # Latin all match

# tidy
landings1 <- landings1 %>% 
  ungroup() %>% 
  select(- Species_std, - Species.scientific.name, - FAO_Species_Code)

## 1950-2010

landings2 <- landings1950_2010_sum %>% 
  left_join(
    species_key,
    by = c("Species" = "Species_raw")
  )

# check all species are matched
landings2 %>%
  filter(is.na(Species_std)) %>%
  distinct(Species)

table(landings2$Species == landings2$Species_std) #319 false

# Replace original species column with standardised one
landings2 <- landings2 %>% 
  mutate(Species = Species_std)

table(landings2$Species == landings2$Species_std) # now all true

# tidy
landings2 <- landings2 %>% 
  select(- Species_std)

## 2006-2023

landings3 <- landings2006_2023_final %>% 
  left_join(
    species_key,
    by = c("Species" = "Species_raw")
  )

# check all species are matched
landings3 %>%
  filter(is.na(Species_std)) %>%
  distinct(Species)

table(landings3$Species == landings3$Species_std) # 90 false

# Replace original species column with standardised one
landings3 <- landings3 %>% 
  mutate(Species = Species_std)

table(landings3$Species == landings3$Species_std) # all treu now
table(landings3$Latin.name == landings3$Species.scientific.name) # Latin all match

# tidy
landings3 <- landings3 %>% 
  select(- Species.scientific.name, - Species_std, - FAO_Species_Code)

# Check columns are the same

common_cols <- Reduce(
  intersect,
  list(
    names(landings1),
    names(landings2),
    names(landings3)
  )
)

# Combine all

landings_all <- bind_rows(
  landings1,
  landings2,
  landings3
)

summary(landings_all) 

unique(landings_all$Species) #246

# # Check for duplicate years
landings_all %>%
  count(Year) %>%
  filter(n > 1)

# Check out 2008
land_2008 <- landings_all %>% 
  filter(Year == 2008)
# good sign that duplicates have very similar if not exact same landings values!
# means that the two dataset match up

# remove the duplicates
landings_all <- landings_all %>%
  arrange(Species, Year, desc(Species)) %>%
  distinct(Year, Species, .keep_all = TRUE)

# Check out 2008 again
land_2008_v2 <- landings_all %>% 
  filter(Year == 2008)

# no more duplicates

## Save full landings data:
write.csv(landings_all, "processed_data/landings_1903_2023.csv")

### Landings plots ----

unique(landings_all$Species)

unique(landings_all$Group_wide)

## Filter out Nephrops, skates and rays, catcharks, key gadoid and flatfish species
## and scallops

landings_plot_df <- landings_all %>% 
  filter(Group_wide == "Skates and rays" | 
           Latin.name %in% c("Pleuronectes platessa", "Limanda limanda",
                             "Scyliorhinus canicula", "Solea solea", "Melanogrammus aeglefinus",
                             "Merlangius merlangus", "Nephrops norvegicus", "Gadus morhua") |
           Species %in% c("Catsharks, nursehounds nei", "Small-spotted catshark",  "Nursehound",
                          "Great Atlantic scallop", "Queen scallop", "Scallops nei"))

unique(landings_plot_df$Species)

landings_plot_df <- landings_plot_df %>% 
  mutate(Category = case_when(Group_wide == "Skates and rays" ~ "Skates and rays",
                              Species %in% c("Catsharks, nursehounds nei", "Small-spotted catshark",  "Nursehound") ~ "Catsharks",
                              Species %in% c("Great Atlantic scallop", "Queen scallop", "Scallops nei") ~ "Scallops",
                              TRUE ~ Species))

unique(landings_plot_df$Category)

# Reorder 
landings_plot_df$Category <- factor(landings_plot_df$Category, 
                                    levels = c("Atlantic cod", "Haddock", "Whiting",
                                               "Common dab", "Common sole", "European plaice", 
                                               "Catsharks", "Skates and rays", "Scallops", "Norway lobster"),
                                    labels = c("Cod", "Haddock", "Whiting", "Dab", "Sole", "Plaice", 
                                               "Catsharks", "Skates and rays", "Scallops", "Nephrops"))

# Aggregate by category for plotting
landings_plot_df_cat <- landings_plot_df %>% 
  group_by(Year, Category) %>%
  summarise(Landings = sum(Landings, na.rm = TRUE), .groups = "drop")

cols_9_cb <- c(
  # Blues 
  "#03045e",  # deep navy
  "#90e0ef",  # light blue
  "#0077b6",  # mid blue
  
  # Greens 
  #"#33A02C",  # deep green
  "#B2DF8A",  # light green
  "#132a13",  # deep forest
  "#90a955",  # light olive
  
  # Purples 
  "#CAB2D6",  # light purple
  "#6A3D9A",  # deep purple
  
  # Orange
  "#ffd670",
  #"#f4a261",
  "#e76f51"
)

# top commercial groups and top bioturbators
(landings_plot <- ggplot(landings_plot_df_cat, aes(x = Year, y = Landings, fill = Category)) +
    geom_area() +
    geom_vline(xintercept = 1988, colour = "lightgrey", 
               linetype = "dashed", size = 1) +
    labs(y = "Landings (tonnes)\n",
         fill = "") +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.1))) +
    scale_fill_manual(values = cols_9_cb) +
    theme_convex() +
    theme(#legend.position = "right",
      legend.position = c(0.15, 0.65),
      #legend.text = element_text(face = "italic")
    )
)


#### Landings by mode ----

# Subset landings dataset to just bioturbating species

biot_landings2 <- landings_all %>%
  inner_join(biot, by = c("Latin.name" = "Species"))

unique(biot_landings2$Species)

# Plot landings of all bioturbating species (split into bioturbation modes)

mode_landings <- biot_landings2 %>% 
  group_by(Year, Bioturbation.mode) %>% 
  summarise(Landings = sum(Landings))

# Reorder modes
mode_landings$Bioturbation.mode <- factor(mode_landings$Bioturbation.mode, levels = c("Burrower", "Vertical excavator",
                                                                                      "Nest-builder", "Lateral excavator",
                                                                                      "Sediment sifter"))

(landings_mode_plot <- ggplot(mode_landings, aes(x = Year, y = Landings, fill = Bioturbation.mode)) +
    geom_area() +
    geom_vline(xintercept = 1988, colour = "grey", 
               linetype = "dashed", size = 1) +
    labs(y = "Landings (tonnes)\n",
         fill = "") +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.1))) +
    # scale_fill_manual(values = cols_22) +
    scale_fill_viridis_d() +
    theme_convex() +
    theme(#legend.position = "right",
      legend.position = c(0.16, 0.85),
      #legend.text = element_text(face = "italic")
    )
)

# Combine 2 landings plots for Figure S3

landings_plots <- plot_grid(landings_plot, landings_mode_plot, 
                            labels = "AUTO",
                            nrow = 2,
                            align = "hv")
landings_plots

ggsave(landings_plots, filename = "outputs/Figure_S3.png", width = 9, height = 9)
