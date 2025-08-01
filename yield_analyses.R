##COMPARISON OF HARVEST-LEVEL YIELD##
##The code below calculates harvest-level yield (in food energy), conducts matching (testing different
##algorithms and distance metrics and identifies the best matching run), and conducts g-computation 
##(for the main analysis of the paper)

###Download the data from 10.5281/zenodo.16687021

##Load all packages needed
library(dplyr)
library(plyr)
library(magrittr)
library(lattice)
library(car)
library(ggplot2)
library(lme4)
library(car)
library(DHARMa)
library(ape)
library(geosphere)
library(MatchIt)
library(cobalt)
library("marginaleffects")
library(lme4)
library("sandwich")
library("clubSandwich")
library("boot")

############# 1. PREPARE DATA#############

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

####check that ZBNF and chemical don't differ in the number of crops grown per year
test7 <- yield4 %>%
  dplyr::group_by(Farming.practice) %>%
  dplyr::summarise_at(vars(no_crops_year), mean, na.fm = T)
##all ok

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

##code farming practice as 0 and 1
yield10.1 <- yield5 %>%
  dplyr::mutate(farm_type = case_when(Farming.practice == "ZBNF" ~ 1,Farming.practice == "Chemical" ~0)) %>% 
  dplyr::group_by(farm_type) %>% 
  ungroup()

####create new cleaned harvest season column 
yield10 <- yield10.1 %>%
  mutate(season2 = ifelse(Season %in% c("kharif", "rabi"), Season, "multiyear"))

############# 2. EXPLORE DATA / TEST ASSUMPTIONS#############
###check for outliers
boxplot(yield10$yield_gj_ha,  ylab = "yield gj")
dotchart(yield10$yield_gj_ha, xlab = "yield gj", 
         ylab = "Order of the data")
bwplot(yield_gj_ha ~ Farming.practice, data = yield10,
       strip = strip.custom(bg = 'white'),
       cex = .5, layout = c(2, 1),
       xlab = "Farming practice", ylab = "yield per ha",
       par.settings = list(
         box.rectangle = list(col = 1),
         box.umbrella  = list(col = 1),
         plot.symbol   = list(cex = .5, col = 1)),
       scales = list(x = list(relation = "same"),
                     y = list(relation = "same")))
#remove outliers
yield13 <- yield10 %>% subset(Crop != "tapioca") ##remove tapioca
boxplot(yield13$yield_gj_ha,  ylab = "yield gj")
dotchart(yield13$yield_gj_ha, xlab = "yield gJ", ##ok now
         ylab = "Order of the data")
bwplot(yield_gj_ha~ Farming.practice, data = yield13,
       strip = strip.custom(bg = 'white'),
       cex = .5, layout = c(2, 1),
       xlab = "Farming practice", ylab = "yield per ha per harvest (gj)",
       par.settings = list(
         box.rectangle = list(col = 1),
         box.umbrella  = list(col = 1),
         plot.symbol   = list(cex = .5, col = 1)),
       scales = list(x = list(relation = "same"),
                     y = list(relation = "same")))
##check for normal distribution
hist(yield13$yield_gj_ha,
     xlab = "yield GJ/ha",
     main = "", ylab = "Frequency")
##data is right skewed --> transform by taking the sqrt
hist(sqrt(yield13$yield_gj_ha),
     xlab = "yield GJ/ha",
     main = "", ylab = "Frequency") #ok now
yield13$yield_gj_ha_sqrt <- sqrt(yield13$yield_gj_ha)

###check for homogeneity of variance (Levene's test)
var.test(yield_gj_ha_sqrt ~ Farming.practice, data = yield13) ##ok

###check for co-linearity
modt1 <- car::vif(lm(yield_gj_ha_sqrt ~ Farming.practice + travel50 + owner2+ irrigation2+ 
                       + perc_quality + suitability_index + Crop + prop_harvest_sold,
                     data = yield13))
modt1 ##ok

##check for overdispersion
DHARMa::testDispersion(modt1) ##ok

##check for spatial autocorrelation (Moran I distance)
##get coordinates for each point
coor <- read.csv("points_coordinates.csv") ###NB as this constitutes confidential information coordinates will only be supplied upon reasonable request 
coor2 <- coor %>% dplyr::rename(Point.ID = point_ID) 
##combine datasets
yield14 <- join_all(list(yield13, coor2),
                    by = c("Point.ID"), type = 'full')
##get mean per square
yield.mean <- yield14 %>% group_by(Point.ID) %>%
  dplyr::mutate(yield_kj2 = mean(yield_kj_sqrt)) %>%
  ungroup()
##to check for autocorrelation keep only one row per square
yield.sq <- yield.mean %>% distinct(Square.ID, .keep_all = TRUE) 
plot.dists <- as.matrix(distm(cbind(yield.sq$Longitude, yield.sq$Latitude), fun= distHaversine))
plot.dists.inv <- 1/plot.dists
plot.dists.inv[which(!is.finite(plot.dists.inv))] <- 0
plot.dists.inv
diag(plot.dists.inv) <- 0
Moran.I(yield.sq$yield_kj_sqrt, plot.dists.inv) ##no spatial autocorrelation


############# 3. Identify matching algorithm and distance that best suits the data#############

####check collinearity of covariates for matching (ie predicting farm type, not yield)
modn2 <- car::vif(lm(farm_type  ~ ##use farm type coded as 0 as 1 as can't use categorical values otherwise
                       + perc_quality + suitability_index + Crop + travel50 + owner2 + prop_harvest_sold + irrigation2,
                     data = yield13))
modn2 ##all ok - all GVIF^(1/(2*Df)) are < 2

###check covariate balance prior to matching
## Mahalanobis distance
mod1pf <- matchit(farm_type ~ perc_quality + suitability_index + Crop
                  + travel50 +owner2 + irrigation2+prop_harvest_sold, ##y needs to be in 0 and 1
                  data = yield13,
                  method = NULL, distance = "mahalanobis", estimand = "ATT")
mod1pf
summary(mod1pf)
###evaluate balanace
wfarm <- bal.tab(mod1pf, disp = c("means", "sds"), un = TRUE, 
                 stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                 thresholds = c(m = .25),
                 disp.means = TRUE, 
                 s.d.denom = "treated") 
wfarm
love.plot(wfarm, binary = "std", thresholds = c(m = .25)) 
## % WFH, suitability index, and travel time are imbalanced
######extract SDMs for comparison between matching methods
colSDM <- summary(mod1pf)[["sum.all"]][, "Std. Mean Diff."]
SDMun <- as.data.frame(colSDM)
SDMun$covariate <- rownames(SDMun)
SDMun$method <- "unadjusted"

###Nearest neighbour matching, mahalanobis distance, with repacement
mod_near_mahal <- matchit(farm_type ~ 
                            perc_quality + suitability_index + Crop + travel50 + owner2 + prop_harvest_sold +irrigation2,
                          data = yield13,
                          method = "nearest", distance = "mahalanobis", estimand = "ATT", replace = TRUE, pop.size = 1000)
mod_near_mahal
summary(mod_near_mahal)
w_near_mahal <- bal.tab(mod_near_mahal, disp = c("means", "sds"), un = TRUE, 
                        stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                        thresholds = c(m = .25),
                        disp.means = TRUE, 
                        s.d.denom = "treated")
w_near_mahal
love.plot(mod_near_mahal, binary = "std", thresholds = c(m = .25)) 
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_near_mahal)[["sum.matched"]][, "Std. Mean Diff."]
SDMnm <- as.data.frame(colSDM)
SDMnm$covariate <- rownames(SDMnm)
SDMnm$method <- "near_mahal_repl"

###Nearest neighbour matching, mahalanobis distance, without repacement
mod_near_mahal2 <- matchit(farm_type ~ 
                             perc_quality + suitability_index + Crop + travel50 + owner2 + prop_harvest_sold+irrigation2,
                           data = yield13,
                           method = "nearest", distance = "mahalanobis", estimand = "ATT", replace = FALSE, pop.size = 1000)
mod_near_mahal2
summary(mod_near_mahal2)
w_near_mahal2 <- bal.tab(mod_near_mahal2, disp = c("means", "sds"), un = TRUE, 
                         stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                         thresholds = c(m = .25),
                         disp.means = TRUE, 
                         s.d.denom = "treated")
w_near_mahal2
love.plot(mod_near_mahal2, binary = "std", thresholds = c(m = .25)) 
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_near_mahal2)[["sum.matched"]][, "Std. Mean Diff."]
SDMnm2 <- as.data.frame(colSDM)
SDMnm2$covariate <- rownames(SDMnm2)
SDMnm2$method <- "near_mahal_no_repl"

#####Nearest neighbour matching, propensity score distance, with replacement
mod_near_prop <- matchit(farm_type ~  perc_quality + suitability_index + Crop + travel50 + owner2 + prop_harvest_sold+irrigation2,
                         data = yield13,
                         method = "nearest", distance = "glm", estimand = "ATT", replace = TRUE, pop.size = 1000)
mod_near_prop
summary(mod_near_prop)
w_near_prop <- bal.tab(mod_near_prop, disp = c("means", "sds"), un = TRUE, 
                       stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                       thresholds = c(m = .25),
                       disp.means = TRUE, 
                       s.d.denom = "treated") 
w_near_prop
love.plot(mod_near_prop, binary = "std", thresholds = c(m = .25)) 
##both good/ok balance but ESS below 30 for both
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_near_prop)[["sum.matched"]][, "Std. Mean Diff."]
SDMnp <- as.data.frame(colSDM)
SDMnp$covariate <- rownames(SDMnp)
SDMnp$method <- "near_prop_repl"

#####Nearest neighbour matching, propensity score distance, without replacement
mod_near_prop2 <- matchit(farm_type ~  perc_quality + suitability_index + Crop + travel50 + owner2 + prop_harvest_sold+irrigation2,
                          data = yield13,
                          method = "nearest", distance = "glm", estimand = "ATT", replace = FALSE, pop.size = 1000)
mod_near_prop2
summary(mod_near_prop2)
w_near_prop2 <- bal.tab(mod_near_prop2, disp = c("means", "sds"), un = TRUE, 
                        stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                        thresholds = c(m = .25),
                        disp.means = TRUE, 
                        s.d.denom = "treated") 
w_near_prop2
love.plot(mod_near_prop2, binary = "std", thresholds = c(m = .25)) ##still as above, even worse than original
##both good/ok balance but ESS below 30 for both
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_near_prop2)[["sum.matched"]][, "Std. Mean Diff."]
SDMnp2 <- as.data.frame(colSDM)
SDMnp2$covariate <- rownames(SDMnp2)
SDMnp2$method <- "near_prop_no_repl"

#####genetic matching, with replacement
mod_genet <- matchit(farm_type ~ perc_quality + suitability_index + Crop + travel50 + owner2 + prop_harvest_sold+irrigation2,
                     data = yield13, pop.size=1000,
                     method = "genetic", estimand = "ATT", replace = TRUE)
mod_genet
summary(mod_genet)
w_genet <- bal.tab(mod_genet, disp = c("means", "sds"), un = TRUE, 
                   stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                   thresholds = c(m = .25),
                   disp.means = TRUE, 
                   s.d.denom = "treated") 
w_genet
love.plot(mod_genet, binary = "std", thresholds = c(m = .25)) 
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_genet)[["sum.matched"]][, "Std. Mean Diff."]
SDMg <- as.data.frame(colSDM)
SDMg$covariate <- rownames(SDMg)
SDMg$method <- "genet_repl"

#####genetic matching, without replacement
mod_genet2 <- matchit(farm_type ~ perc_quality + suitability_index + Crop + travel50 +
                        owner2 + prop_harvest_sold+irrigation2,
                      data = yield13, pop.size=1000,
                      method = "genetic", estimand = "ATT", replace = FALSE)
mod_genet2
summary(mod_genet2)
w_genet2 <- bal.tab(mod_genet2, disp = c("means", "sds"), un = TRUE, 
                    stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                    thresholds = c(m = .25),
                    disp.means = TRUE, 
                    s.d.denom = "treated") 
w_genet2
love.plot(mod_genet2, binary = "std", thresholds = c(m = .25)) 
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_genet2)[["sum.matched"]][, "Std. Mean Diff."]
SDMg2 <- as.data.frame(colSDM)
SDMg2$covariate <- rownames(SDMg2)
SDMg2$method <- "genet_no_repl"

#####genetic matching, with replacement, matched on propensity score first
ps_model <- glm(farm_type ~ perc_quality + suitability_index + Crop + travel50 + owner2 + prop_harvest_sold+irrigation2,
                data = yield13, family = "binomial")
propensity_scores <- predict(ps_model, type = "response")
mod_genet_prop <- matchit(farm_type ~ perc_quality + suitability_index + Crop + travel50 + owner2 
                          + prop_harvest_sold + propensity_scores,
                          data = yield13, pop.size=1000,
                          method = "genetic", estimand = "ATT", replace = TRUE)
mod_genet_prop
summary(mod_genet_prop)
w_genet_prop <- bal.tab(mod_genet_prop, disp = c("means", "sds"), un = TRUE, 
                        stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                        thresholds = c(m = .25),
                        disp.means = TRUE, 
                        s.d.denom = "treated") 
w_genet_prop
love.plot(mod_genet_prop, binary = "std", thresholds = c(m = .25)) 
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_genet_prop)[["sum.matched"]][, "Std. Mean Diff."]
SDMg3 <- as.data.frame(colSDM)
SDMg3$covariate <- rownames(SDMg3)
SDMg3$method <- "genet_repl_prop"

#####genetic matching, with replacement, matched on propensity score first but only use most important 
#covariates in the calculation of the distance and weights (mahavars) but optimise balance on all covarites see https://ngreifer.github.io/blog/genetic-matching/
mod_genet_mah <- matchit(farm_type ~ perc_quality + suitability_index + Crop + travel50 + owner2 + prop_harvest_sold+irrigation2,
                         data = yield13, pop.size=1000,
                         method = "genetic", mahvars = ~ perc_quality + suitability_index + Crop, estimand = "ATT", replace = TRUE)
mod_genet_mah
summary(mod_genet_mah)
w_genet_mah <- bal.tab(mod_genet_mah, disp = c("means", "sds"), un = TRUE, 
                       stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                       thresholds = c(m = .25),
                       disp.means = TRUE, 
                       s.d.denom = "treated") 
w_genet_mah
love.plot(mod_genet_mah, binary = "std", thresholds = c(m = .25))
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_genet)[["sum.matched"]][, "Std. Mean Diff."]
SDMg4 <- as.data.frame(colSDM)
SDMg4$covariate <- rownames(SDMg4)
SDMg4$method <- "genet_repl_mahvars"

#####full matching, mahalanobis distance
mod_full_mahal <- matchit(farm_type ~  perc_quality + suitability_index + Crop + travel50 + owner2 + prop_harvest_sold+irrigation2,
                          data = yield13,
                          method = "full", distance = "mahalanobis", estimand = "ATT", pop.size = 1000)
mod_full_mahal
summary(mod_full_mahal)
w_full_mahal <- bal.tab(mod_full_mahal, disp = c("means", "sds"), un = TRUE, 
                        stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                        thresholds = c(m = .25),
                        disp.means = TRUE, 
                        estimand = "ATT")
w_full_mahal
love.plot(mod_full_mahal, binary = "std", thresholds = c(m = .25))
##suitability and travel time remain slightly imbalanced (around 0.3), but ESS > 30
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_full_mahal)[["sum.matched"]][, "Std. Mean Diff."]
SDMfm <- as.data.frame(colSDM)
SDMfm$covariate <- rownames(SDMfm)
SDMfm$method <- "full_mahal"

#####full matching, propensity score distance
mod_full_prop <- matchit(farm_type ~ perc_quality+ suitability_index + Crop + travel50 + owner2 + prop_harvest_sold+irrigation2,
                         data = yield13,
                         method = "full", distance = "glm", estimand = "ATT")
mod_full_prop
summary(mod_full_prop)
w_full_prop <- bal.tab(mod_full_prop, disp = c("means", "sds"), un = TRUE, 
                       stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                       thresholds = c(m = .25),
                       disp.means = TRUE, 
                       s.d.denom = "treated") ##treated recommended for 
w_full_prop
love.plot(mod_full_prop, binary = "std", thresholds = c(m = .25)) 
###all matched now but ESS < 30
colSDM <- summary(mod_full_prop)[["sum.matched"]][, "Std. Mean Diff."]
SDMfp <- as.data.frame(colSDM)
SDMfp$covariate <- rownames(SDMfp)
SDMfp$method <- "full_prop"

##crate dataset containing SDMs from all matching methods
sdms.data <- purrr::reduce(
  list(SDMun, SDMnm, SDMnm2, SDMnp, SDMnp2, SDMg, SDMg2, SDMg3, SDMg4, SDMfm, SDMfp),
  function(left, right) {
    dplyr::full_join(left, right)
  }
)
####Create figure comparing SDMs #########
custom_colors <- c(
  full_mahal = "springgreen",
  full_prop = "chartreuse3",
  genet_no_repl = "lightblue",
  genet_repl = "darkblue",
  genet_repl_prop = "cyan",
  genet_repl_mahvars = "cornflowerblue", 
  near_mahal_no_repl = "gold1",
  near_mahal_repl = "goldenrod3",
  near_prop_no_repl = "orangered",
  near_prop_repl = "red4",
  unadjusted = "black"
)

# Create a custom y-axis label and order
y_axis_labels <- c(
  distance = "Propensity Score",
  irrigation2Y = "Irrigated",
  irrigation2N = "Rainfed",
  suitability_index = "Agroecological suitability",
  perc_quality = "Proportion of habitat patches",
  owner2N = "Land rights - owner",
  owner2Y = "Land rights - tenure",
  prop_harvest_sold = "Market ties",
  travel50 = "Access to cities",
  Croprice = "Crop - rice",
  Cropblackgram = "Crop - blackgram",
  Cropgreengram = "Crop - greengram",
  Croptomato = "Crop - tomato",
  Cropsunflower = "Crop - sunflower",
  Cropsugarcane = "Crop - sugarcane",
  Cropsesame = "Crop - sesame",
  Cropragi = "Crop - ragi",
  Cropmaize = "Crop - maize",
  Cropgroundnut = "Crop - groundnut",
  Cropcoconut = "Crop - coconut",
  Cropcashew = "Crop - cashew",
  Cropbottlegourd = "Crop - bottlegourd",
  Cropbeans = "Crop - beans",
  Cropaubergine = "Crop - aubergine"
)

sdms.data2 <- sdms.data %>% filter(covariate != "propensity_scores") ##already included under "distance" row
# Create the ggplot figure with custom colors, labels, and y-axis order
ggplot(data = sdms.data2, aes(x = colSDM, y = factor(covariate, levels = names(y_axis_labels)), color = method)) +
  geom_point() +
  geom_vline(xintercept = c(-0.25, 0, 0.25), linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +  # Add a black horizontal line at y = 0
  labs(
    x = "Standardised Difference in Means",
    y = NULL,
  ) +
  scale_color_manual(
    name = "Matching Method",
    values = custom_colors,
    labels = c(
      "full_mahal" = "Full matching, Mahalanobis (N=87/115, u=4)",
      "full_prop" = "Full matching, Propensity Score (N=87/115, u=6)",
      "genet_no_repl" = "Genetic matching, without replacement (N=87/87, u=5)",
      "genet_repl" = "Genetic matching, with replacement (N=28/115, u=1)",
      "genet_repl_prop" = "Genetic matching, with replacement and propensity score (N=31/115, u=1)",
      "genet_repl_mahvars" = "Genetic matching, with replacement, key covariates (N=46/115, u=4)",
      "near_mahal_no_repl" = "Nearest Neighbour, Mahalanobis, without replacement (N=87/87, u=8)",
      "near_mahal_repl" = "Nearest Neighbour, Mahalanobis, with replacement (N=33/115, u=3)",
      "near_prop_no_repl" = "Nearest Neighbour, Propensity Score, without replacement (N=87/87, u=5)",
      "near_prop_repl" = "Nearest Neighbour, Propensity Score, with replacement (N=24/115, u=5)",
      "unadjusted" = "Unadjusted (N=87/115, u=4)"
    )
  ) +
  scale_y_discrete(labels = y_axis_labels) +  # Customize y-axis labels and order
  theme_minimal() +  # Set a minimal theme
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.text.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14),# Increase y-axis label size
    legend.text = element_text(size = 14),    # Increase legend text size
    legend.title = element_text(size = 16)     # Increase legend title size
  )


############# 4. Conduct post-matching analyses#############
##only for matching analyses where n(control)>30 and a maximum of 3 unbalanced covariates (which must then be included as covariates in the model)

#####Main analysis: genetic matching with replacement, incl. matched on propensity score
boot_fun3 <- function(data, i) { ##bootstrapping
  boot_data3 <- data[i,]
  ps_model <- glm(farm_type ~ perc_quality + suitability_index + Crop + travel50 + owner2 + prop_harvest_sold + irrigation2,
                  data = boot_data3, family = "binomial")
  propensity_scores <- predict(ps_model, type = "response")
  mod_genet_prop <- matchit(farm_type ~ perc_quality + suitability_index + Crop + travel50 + owner2 
                            + prop_harvest_sold + propensity_scores,
                            data = boot_data3, pop.size=1000,
                            method = "genetic", estimand = "ATT", replace = TRUE) ###change for the different matching runs accordingly
  #Extract matched dataset
  md <- match.data(mod_genet_prop, data = boot_data3)
  #Fit outcome model
  fit <- lmer(yield_kj_sqrt ~ farm_type + perc_quality+
                (1|Square.ID), 
              data = md, weights = weights) 
  ## G-computation ##
  #Subset to treated units for ATT
  md1 <- subset(md, farm_type == 1)
  #Estimated potential outcomes under treatment
  p1 <- predict(fit, type = "response",
                newdata = transform(md1, farm_type = 1))
  Ep1 <- weighted.mean(p1, md1$weights)
  #Estimated potential outcomes under control
  p0 <- predict(fit, type = "response",
                newdata = transform(md1, farm_type = 0))
  Ep0 <- weighted.mean(p0, md1$weights)
  #Risk ratio (RR) - Ep1 / Ep0 - use to get % change in yield due to ZBNF programme  
  #change to Ep1 - Ep0 to get RD - risk difference, i.e. measure in outcome unit (yield)
  ##both are measures of marginal effect
  ##return Ep1 to get estimated yield for ZBNF, and Ep0 to get estimated yield for chemical farming
  return(Ep1/Ep0)
}
set.seed(54321)
boot_out3 <- boot(yield13, boot_fun3, R = 999) ##change to 999 minimum
boot_out3
boot.ci(boot_out3, type = "bca") #bias-corrected and accelerated (BCa) bootstrap confidence intervals


