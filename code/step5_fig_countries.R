# clear existing workspace 
rm(list = ls())

# install necessary packages
if(!require(rcartocolor)){install.packages('rcartocolor'); library(rcartocolor)}
if(!require(rnaturalearth)){install.packages('rnaturalearth'); library(rnaturalearth)}

# get Africa shapefile
countries <- ne_countries(continent = 'Africa')

# add in the CV groups 
countries$CV_group <- NA
countries$CV_group[countries$adm0_a3_is %in% c('SSD', 'SLE', 'UGA')] <- 1
countries$CV_group[countries$adm0_a3_is %in% c('GNQ', 'TGO', 'BFA')] <- 2
countries$CV_group[countries$adm0_a3_is %in% c('ERI', 'GMB', 'NGA')] <- 3
countries$CV_group[countries$adm0_a3_is %in% c('SDN', 'GIN', 'TZA')] <- 4
countries$CV_group[countries$adm0_a3_is %in% c('GNB', 'KEN', 'CIV')] <- 5
countries$CV_group[countries$adm0_a3_is %in% c('NER', 'COD', 'TCD')] <- 6
countries$CV_group[countries$adm0_a3_is %in% c('GAB', 'MLI', 'ZMB', 'SOM')] <- 7
countries$CV_group[countries$adm0_a3_is %in% c('CMR', 'COG', 'GHA', 'AGO')] <- 8
countries$CV_group[countries$adm0_a3_is %in% c('RWA', 'ETH', 'SEN', 'CAF')] <- 9
countries$CV_group[countries$adm0_a3_is %in% c('BEN', 'MRT', 'LBR', 'BDI')] <- 10

# specify color palette
palette <- carto_pal(n = 10, name = 'Prism')

# plot 
jpeg(filename = '../figures/map_cross_validation.jpg', width = 6.5, height = 5, units = 'in', res = 500)
par(mar = c(0,0,0,0))
plot(countries, col = palette[countries$CV_group])
dev.off()