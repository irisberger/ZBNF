####### The code below estimates square-level yield and profit#####

###Download the data from 10.5281/zenodo.16687021

###load packages
library(dplyr)
library(plyr)
library(tidyr)

  #############   1. SQUARE LEVEL YIELD #########

##import yield data
yield.data <- read.csv("yield_v3.csv")

##import data needed to convert harvest in kg into energy
energy <- read.csv("energy.csv")
energy2 <- energy %>% dplyr::rename(Crop = food_item) 

##combine datasets
yield2 <- plyr::join_all(list(yield.data, energy2),
                         by = c("Crop"), type = 'full')

##calculate the food energy in KJ per acre for each harvest
yield2$food.enegry.kg_actual <- yield2$food.energy_KJ.1kg/100 ##get the food energy in KJ per kg per crop
yield3.1 <- yield2 %>% ##multiply proportion of harvest kept (excl husks, peels etc) by harvest in kg and per kg energy values
  dplyr::mutate(yield_kj_acre = (yield_kg.acre_annual*prop_kept)*food.enegry.kg_actual) %>% 
  ungroup()

## Calculate yield in GJ per ha
yield3.2 <- yield3.1 %>%
  dplyr::mutate(
    yield_gj_acre = yield_kj_acre / 1000000, ##convert to GJ per acre
    yield_gj_ha = yield_gj_acre / 0.404686 ##convert to GJ per ha
  ) %>%
  ungroup()

##get total number of crops grown per field per year
yield4 <- yield3.2 %>%
  dplyr::group_by(Field.ID) %>%
  dplyr::mutate(no_crops_year = n()) %>% 
  ungroup()

##get covariates:
###import data on the percentage of vegetation patch area in the 100m around each point
covar <- read.csv("covariates_combined_farm.csv")
covar2 <- covar%>% dplyr::rename(Point.ID = point_ID) 
###import data on the % vegetation patches in each square
wfh.square <- read.csv("wfh_squares_final.csv")
covar3 <- join_all(list(covar2, wfh.square),
                   by = c("site_ID"), type = 'full')
##import data on the travel time (min) to nearest city (>50k inhabitants)
travel <- read.csv("farms_traveltime.csv")
travel2 <- travel%>% dplyr::rename(Point.ID = point_ID)
####import data on farmers' total landholding sizes, ownership, and region (tribal vs plain area) data
f.size <- read.csv("field_details_results.csv")
f.size2 <- f.size%>% dplyr::select(Point.ID, total_acres_farmed, owner, plain_tribal)
####import data on irrigation
irrigation <- read.csv("irrigation.csv")

yield5 <- join_all(list(yield4, covar3, travel2, f.size2, irrigation),
                   by = c("Point.ID"), type = 'full')

# create tree and non-tree crop groups
yield6 <- yield5 %>%
  mutate(crop_type = case_when(
    Crop %in% c("rice", "sunflower", "maize", "tomato", "beans", "bottlegourd",
                "aubergine", "ragi", "groundnut", "blackgram", "greengram", "sesame", "sugarcane",
                "tapioca") ~ "ground",
    Crop %in% c("cashew", "coconut") ~ "tree",
    TRUE ~ NA_character_
  )) %>%
  ungroup()

#get total yield per point ID (field)
yield6.1 <- yield6.2 %>% 
  dplyr::group_by(Point.ID) %>%
  dplyr::mutate(yield_kj_point = sum(yield_kj)) %>%
  ungroup()
##keep only 1 row per point
yield6.3 <- yield6.1 %>% distinct(Point.ID, .keep_all = TRUE)

###calculate the weighted mean for the area under production per acre for each crop class in each square (weighted by field size)
yield7 <- yield6.3 %>% 
  dplyr::group_by(Square.ID, crop_type) %>%
  dplyr::mutate(weighted_yield_kjacre = weighted.mean(yield_kj_point, Field.size..acres., na.rm = TRUE)) %>%
  dplyr::mutate(non_weighted_yield_kjacre = mean(yield_kj_point)) %>%
  ungroup()

###keep 1 row per crop type per square
yield8 <- yield7 %>%
  distinct(Square.ID, crop_type, .keep_all = TRUE)

yield9 <- yield8 %>%
  tidyr::complete(Square.ID, crop_type = c("ground", "tree"), fill = list(weighted_yield_kjacre = 0))

ee2 <- yield9 %>% subset(select = c(Square.ID, Field.ID, Crop, crop_type, weighted_yield_kjacre, perc_ground_cropland, perc_plantation, perc_quality, yield_kg.acre_annual))

##calculate the total cropland area per square: total area of square minus the area 
##occupied by vegetation patches and minus the area of tree crops (plantation)
comp <- ee2 %>% 
  dplyr::mutate(perc_crop_new = 100 - (perc_quality + perc_plantation))
yield82 <- comp %>% 
  dplyr::mutate(crop_area = 61.7763 * perc_crop_new, ##each 500m x 500m square is 61.7763 acres
                tree_crop_area = 61.7763 * perc_plantation) %>%
  ungroup()

# Calculate total yield for each square
yield84 <- yield82 %>%
  dplyr::group_by(Square.ID) %>%
  dplyr::mutate(
    yield_ground_square = case_when(
      crop_type == "ground" ~ crop_area * weighted_yield_kjacre,
      TRUE ~ NA_real_
    ),
    yield_tree_square = case_when(
      crop_type == "tree" ~ tree_crop_area * weighted_yield_kjacre,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()

yield85 <- yield84 %>%
  dplyr::group_by(Square.ID) %>%
  dplyr::mutate(yield_ground_square = replace_na(yield_ground_square, 0),
                yield_tree_square = replace_na(yield_tree_square, 0),
                yield_square = sum(yield_ground_square, yield_tree_square)) %>%
  ungroup()
##keep one row per square
yield86 <- yield85 %>%
  distinct(Square.ID, .keep_all = TRUE)


#############   2. SQUARE LEVEL PROFIT #########
prof_dat <- read.csv("profit_data_points.csv")
wfh.square$Square.ID <- wfh.square$site_ID
profit1 <- join_all(list(prof_dat, wfh.square),
                    by = c("Square.ID"), type = 'full')

# create tree and non-tree groups (as above)
profit2 <- profit1 %>%
  dplyr::mutate(crop_type = case_when(
    Crop %in% c("rice", "sunflower", "maize", "tomato", "beans", "bottlegourd",
                "aubergine", "ragi", "groundnut", "blackgram", "greengram", "sesame", "sugarcane",
                "tapioca") ~ "ground",
    Crop %in% c("cashew", "coconut") ~ "tree",
    TRUE ~ NA_character_
  )) %>%
  ungroup()

###calculate the weighted mean per farmed bit acre for each crop class in each square
profit3 <- profit2 %>% 
  dplyr::group_by(Square.ID, crop_type) %>%
  dplyr::mutate(weighted_profit_acre = weighted.mean(profit_point_acre, Field.size..acres., na.rm = TRUE)) %>%
  dplyr::mutate(non_weighted_profit_acre = mean(profit_point_acre)) %>%
  ungroup()

profit5 <- profit3 %>%
  ungroup() %>%  # Ensure no grouping is interfering
  complete(Square.ID, crop_type = c("ground", "tree"), fill = list(weighted_profit_acre = 0))

ee2 <- profit5 %>% subset(select = c(Square.ID, Field.ID, Crop, crop_type, weighted_profit_acre, perc_ground_cropland, perc_plantation, perc_quality))

##calculate cropland area as total area (as above)
comp <- ee2 %>% 
  dplyr::group_by(Square.ID) %>%
  dplyr::mutate(perc_crop_new = 100 - (perc_quality + perc_plantation))
###keep unique Square ID and yield_ground_square values
comp2 <- comp %>%
  dplyr::distinct(Square.ID, crop_type, .keep_all = TRUE)
##get area under cropland in each square
yield82 <- comp2 %>% 
  dplyr::mutate(crop_area = 61.7763 * perc_crop_new, ##each 500m x 500m square is 61.7763 acres
                tree_crop_area = 61.7763 * perc_plantation) %>%
  ungroup()

# Calculate total profit for each square
yield84 <- yield82 %>%
  dplyr::group_by(Square.ID) %>%
  dplyr::mutate(
    profit_ground_square = case_when(
      crop_type == "ground" ~ crop_area * weighted_profit_acre,
      TRUE ~ NA_real_
    ),
    profit_tree_square = case_when(
      crop_type == "tree" ~ tree_crop_area * weighted_profit_acre,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()

yield85 <- yield84 %>%
  dplyr::group_by(Square.ID) %>%
  dplyr::mutate(profit_ground_square = replace_na(profit_ground_square, 0),
                profit_tree_square = replace_na(profit_tree_square, 0),
                profit_square = sum(profit_ground_square, profit_tree_square)) %>%
  ungroup()

##keep one row per square
yield86 <- yield85 %>%
  distinct(Square.ID, .keep_all = TRUE)
