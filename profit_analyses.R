###FIELD LEVEL PROFIT ANALYSES
##The code below calculates field-level profit, conducts matching (testing different
##algorithms and distance metrics and identifies the best matching run), and conducts g-computation 
##(for the main analysis of the paper)

###Download the data from 10.5281/zenodo.16687021

##load packages
library(dplyr)
library(plyr)
library(car)
library(ggplot2)
library(lme4)
library(car)
library(DHARMa)
library(ape)
library(geosphere)
library(MatchIt)
library(cobalt)
library("boot")

##########1. Prepare data (calculate field-level profit)##############

###import wholesale crop price data
price.data <- read.csv("wholesale_prices.csv")
head(price.data)
##remove rows/years for which there is no price data
price2 <- price.data %>% subset(!is.na(Prices.in.Rs.Quintal..january.))
###get mean price per crop
price3 <- price2 %>%
  dplyr::group_by(crop) %>%
  dplyr::summarize(price_kg = mean(price_kg_2023INR)) %>% ungroup()
price4 <- price3 %>%
  dplyr::rename("Crop" = "crop") %>%
  ungroup()

###import yield data
yield.data <- read.csv("yield_v3.csv")
# check number of fields
yield.data %>%
  dplyr::group_by(Farming.practice) %>%
  dplyr::summarize(unique_point_ID_count = n_distinct(Point.ID))
##get number of harvests per field
yield.data2 <- yield.data %>%
  dplyr::group_by(Field.ID) %>%
  dplyr::mutate(no_harvests = n())
##combine datasets
revenue <- join_all(list(yield.data2, price4),
                    by = c("Crop"), type = 'full')
##get revenue per harvest, which is sold/bought with inedible portion (husk, shell etc)
revenue3 <- revenue %>% 
  dplyr::mutate(revenue_harvest = harvest..kg. * price_kg) %>% ##use absolute (not per acre) harvest value in kg (need to subtract cost first before estimating per acre profit)
  ungroup()
##cost data is given per year (not per harvest), so need to summarise at the point (i.e. field) level
revenue4 <- revenue3 %>%
  dplyr::group_by(Point.ID) %>%
  dplyr::mutate(revenue_point = sum(revenue_harvest),
                num_harvests = n()) %>% 
  ungroup()
revenue5 <- revenue4 %>% 
  dplyr::distinct(Point.ID, .keep_all = TRUE)

##import cost data
energy <- read.csv("energy_costs.csv")
labour <- read.csv("labour_costs.csv")
rental <- read.csv("rental_costs.csv")
seeds <- read.csv("seed_costs.csv")
substances <- read.csv("substances.csv")
tools <- read.csv("tool_costs.csv")
energy2 <- energy %>%
  dplyr::group_by(Point.ID) %>%
  dplyr::summarize(energy_inr_point = sum(petrol_electr_cost_per_field)) %>% ungroup()
labour2 <- labour %>%
  dplyr::group_by(Point.ID) %>%
  dplyr::summarize(labour_inr_point = sum(Total..annual..amount.paid, Total.equivalent.family.labour.cost..self.repoted.data.)) %>% ungroup()
rental2 <- rental %>%
  dplyr::group_by(Point.ID) %>%
  dplyr::summarize(rental_inr_point = sum(cost)) %>% ungroup()
seeds2 <- seeds %>%
  dplyr::group_by(Point.ID) %>%
  dplyr::summarize(seeds_inr_point = sum(cost_per_harvest_per_field)) %>% ungroup()
substances2 <- substances %>%
  dplyr::group_by(Point.ID) %>%
  dplyr::summarize(substances_inr_point = sum(cost_per_field_per_substance)) %>% ungroup()
tools2 <- tools %>%
  dplyr::group_by(Point.ID) %>%
  dplyr::summarize(tools_inr_point = sum(yearly_cost_per_tool_per_acre)) %>% ungroup()

###combine all cost data (cost for whole year per point, ie field)
costs <- join_all(list(energy2, labour2, rental2, seeds2, substances2, tools2),
                  by = c("Point.ID"), type = 'full') %>%
  tidyr::replace_na(list(energy_inr_point = 0, 
                         labour_inr_point = 0,
                         rental_inr_point = 0,
                         seeds_inr_point = 0,
                         substances_inr_point = 0,
                         tools_inr_point = 0)) 
##calculate total cost
costs2 <- costs %>%
  dplyr::group_by(Point.ID) %>%
  dplyr::mutate(cost_all = sum(energy_inr_point, labour_inr_point, rental_inr_point, seeds_inr_point, substances_inr_point, tools_inr_point)) %>% ungroup()
###combine cost and yield data
profit <- join_all(list(revenue5, costs2),
                   by = c("Point.ID"), type = 'full')
profit2 <- dplyr::distinct(profit, Point.ID, .keep_all = TRUE)
###calculate profit per point
profit3 <- profit %>% 
  dplyr::group_by(Point.ID) %>% 
  dplyr::mutate(profit_point = revenue_point - cost_all)
##calculate per acre profit
profit4 <- profit3 %>%
  dplyr::group_by(Point.ID) %>%
  dplyr::mutate(revenue_point_acre = revenue_point / Field.size..acres.,
                profit_point_acre = profit_point / Field.size..acres.,
                cost_point_acre = cost_all / Field.size..acres.) %>% ungroup()

##import covariates
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
####import irrigation data
irrigation <- read.csv("irrigation.csv")
##combine datasets
yield5 <- join_all(list(profit4, covar3, travel2, f.size2, irrigation),
                   by = c("Point.ID"), type = 'full')
##create column whether or not at least one harvest per field was irrigated
yield6.1 <- yield5 %>%
  dplyr::group_by(Point.ID) %>%
  dplyr::mutate(irrigation2 = if_else(any(irrigation == "Y"), "Y", "N"))
###remove duplicates
yield9 <- dplyr::distinct(yield6.1, Point.ID, .keep_all = TRUE)
##code farming practice in 0 and 1
yield10 <- yield9 %>%
  dplyr::mutate(farm_type = case_when(Farming.practice == "ZBNF" ~ 1,Farming.practice == "Chemical" ~0)) %>%
  dplyr::group_by(farm_type) %>% ungroup()


##########2. Explore data and test assumptions##############
###check for outliers
boxplot(yield10$profit_point_acre,  ylab = "profit inr")
dotchart(yield10$profit_point_acre, xlab = "profit inr", 
         ylab = "Order of the data")
bwplot(profit_point_acre ~ Farming.practice, data = yield10,
       strip = strip.custom(bg = 'white'),
       cex = .5, layout = c(2, 1),
       xlab = "Farming practice", ylab = "Profit INR/acre",
       par.settings = list(
         box.rectangle = list(col = 1),
         box.umbrella  = list(col = 1),
         plot.symbol   = list(cex = .5, col = 1)),
       scales = list(x = list(relation = "same"),
                     y = list(relation = "same")))
#remove outliers
yield13 <- yield10 %>% subset(profit_point_acre < 200000)
boxplot(yield13$profit_point_acre,  ylab = "profit inr")
dotchart(yield13$profit_point_acre, xlab = "profit inr",
         ylab = "Order of the data")
bwplot(profit_point_acre ~ Farming.practice, data = yield13,
       strip = strip.custom(bg = 'white'),
       cex = .5, layout = c(2, 1),
       xlab = "Farming practice", ylab = "yield per acre (kg)",
       par.settings = list(
         box.rectangle = list(col = 1),
         box.umbrella  = list(col = 1),
         plot.symbol   = list(cex = .5, col = 1)),
       scales = list(x = list(relation = "same"),
                     y = list(relation = "same")))
###check if normally distributed
hist(yield13$profit_point_acre,
     xlab = "profit INR/acre ",
     main = "", ylab = "Frequency") ##ok
###check for homogeneity of variance (Levene's test)
var.test(profit_point_acre ~ Farming.practice, data = yield13) ##ok
###check for collinearity
modt <- car::vif(lm(profit_point_acre ~ Farming.practice + travel50 + owner2+ irrigation2+ 
                      + perc_quality + suitability_index + no_harvests,
                    data = yield13))
modt ##all GVIF^(1/(2*Df)) should be below 2; rainfall was greater and was thus excluded
###check for overdispersion
DHARMa::testOverdispersion(model1) ##no

##check for spatial autocorrelation
##get coordinates for each point
coor <- read.csv("points_coordinates.csv") ##NB due this data being sensitive we only provide it upon reasonable request
coor2 <- coor %>% dplyr::rename(Point.ID = point_ID)
##combine data
yield14 <- join_all(list(yield13, coor2),
                    by = c("Point.ID"), type = 'full')
##get mean per square
yield.mean <- yield14 %>% group_by(Point.ID) %>%
  dplyr::mutate(profit_acre2 = mean(profit_point_acre)) %>%
  ungroup()
##to check for autocorrelation keep only one row per square
yield.sq <- yield.mean %>% distinct(Square.ID, .keep_all = TRUE) 
plot.dists <- as.matrix(distm(cbind(yield.sq$Longitude, yield.sq$Latitude), fun= distHaversine))
plot.dists.inv <- 1/plot.dists
plot.dists.inv[which(!is.finite(plot.dists.inv))] <- 0
plot.dists.inv
diag(plot.dists.inv) <- 0
Moran.I(yield.sq$profit_acre2, plot.dists.inv) ##no spatial autocorrelation


############# 3. Identify matching algorithm and distance that best suits the data#############
###check covariate balance prior to matching (Mahalanobis distance)
mod1pf <- matchit(farm_type ~ suitability_index, ##y needs to be in 0 and 1
                  data = yield13,
                  method = NULL, distance = "mahalanobis", estimand = "ATT")
mod1pf
summary(mod1pf)
###evaluate balanace using cobalt
#install.packages("cobalt")
library(cobalt)
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

#######Nearest neighbour matching, mahalanobis distance, with replacement
mod_near_mahal <- matchit(farm_type ~ suitability_index,
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

#########Nearest neighbour matching, mahalanobis distance, without replacement
mod_near_mahal2 <- matchit(farm_type ~ suitability_index,
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

#######Nearest neighbour matching, propensity score distance, with replacement
mod_near_prop <- matchit(farm_type ~  suitability_index,
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

#######Nearest neighbour matching, propensity score distance, without replacement
mod_near_prop2 <- matchit(farm_type ~  suitability_index,
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
love.plot(mod_near_prop2, binary = "std", thresholds = c(m = .25)) 
##both good/ok balance but ESS below 30 for both
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_near_prop2)[["sum.matched"]][, "Std. Mean Diff."]
SDMnp2 <- as.data.frame(colSDM)
SDMnp2$covariate <- rownames(SDMnp2)
SDMnp2$method <- "near_prop_no_repl"

####Genetic matching, with replacement
mod_genet <- matchit(farm_type ~ suitability_index,
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

####Genetic matching, without replacement
mod_genet2 <- matchit(farm_type ~ suitability_index,
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

####Genetic matching, with replacement, matched on propensity score first
ps_model <- glm(farm_type ~ suitability_index,
                data = yield13, family = "binomial")
propensity_scores <- predict(ps_model, type = "response")

mod_genet_prop <- matchit(farm_type ~ suitability_index + 
                            propensity_scores,
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

####full matching, mahalanobis distance
mod_full_mahal <- matchit(farm_type ~  suitability_index,
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
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_full_mahal)[["sum.matched"]][, "Std. Mean Diff."]
SDMfm <- as.data.frame(colSDM)
SDMfm$covariate <- rownames(SDMfm)
SDMfm$method <- "full_mahal"

####full matching, propensity score distance
mod_full_prop <- matchit(farm_type ~ suitability_index,
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
  list(SDMun, SDMnm, SDMnm2, SDMnp, SDMnp2, SDMg, SDMg2, SDMg3, SDMfm, SDMfp),
  function(left, right) {
    dplyr::full_join(left, right)
  }
)
sdms.data$covariate <- ifelse(sdms.data$covariate == "1", "suitability_index", sdms.data$covariate)
####Create figure comapring SDMs #########
# Create a custom color palette
custom_colors <- c(
  full_mahal = "springgreen",
  full_prop = "chartreuse3",
  genet_no_repl = "lightblue",
  genet_repl = "darkblue",
  genet_repl_prop = "cyan",
  near_mahal_no_repl = "gold1",
  near_mahal_repl = "goldenrod3",
  near_prop_no_repl = "orangered",
  near_prop_repl = "red4",
  unadjusted = "black"
)

# Create a custom y-axis label and order
y_axis_labels <- c(
  distance = "Propensity Score",
  suitability_index = "Agroecological suitability"
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
      "full_mahal" = "Full matching, Mahalanobis (N=55/66, u=0)",
      "full_prop" = "Full matching, Propensity Score (N=55/66, u=0)",
      "genet_no_repl" = "Genetic matching, without replacement (N=55/55, u=1)",
      "genet_repl" = "Genetic matching, with replacement (N=30/66, u=0)",
      "genet_repl_prop" = "Genetic matching, with replacement and propensity score (N=29/66, u=0)",
      "near_mahal_no_repl" = "Nearest Neighbour, Mahalanobis, without replacement (N=55/55, u=1)",
      "near_mahal_repl" = "Nearest Neighbour, Mahalanobis, with replacement (N=11/66, u=0)",
      "near_prop_no_repl" = "Nearest Neighbour, Propensity Score, without replacement (N=55/55, u=1)",
      "near_prop_repl" = "Nearest Neighbour, Propensity Score, with replacement (N=12/66, u=0)",
      "unadjusted" = "Unadjusted (N=55/66, u=1)"
    )
  ) +
  scale_y_discrete(labels = y_axis_labels) +  # Customize y-axis labels and order
  theme_minimal() +  # Set a minimal theme
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Use linewidth instead of size
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

############# 4. Conduct post-matching analyses#############
##only for matching analyses where n(control)>30 and a maximum of 3 unbalanced covariates (which must then be included as covariates in the model)

###Full matching with mahalanobis distance performed best --> main analysis
boot_fun <- function(data, i) {
  boot_data <- data[i,]
  
  #Do nm mahalanobis matching with replacement
  m <- matchit(farm_type ~  suitability_index,
               data = yield13,
               method = "full", distance = "mahalanobis", estimand = "ATT", pop.size = 1000) ###change for the different matching runs accordingly
  
  #Extract matched dataset
  md <- match.data(m, data = boot_data)
  
  #Fit outcome model
  fit <- lmer(profit_point_acre ~ farm_type + (1|Square.ID),
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
  
  #Risk ratio (RR) - Ep1 / Ep0 - use to get % change in profit due to ZBNF programme  
  #change to Ep1 - Ep0 to get RD - risk difference, i.e. measure in outcome unit (yield)
  ##both are measures of marginal effect
  ##return Ep1 to get estimated yield for ZBNF, and Ep0 to get estimated yield for chemical farming
  return(Ep1 / Ep0)
}

set.seed(54321)

boot_out <- boot(yield13, boot_fun, R = 999) ##change to 999 minimum
boot_out

boot_ci <- boot.ci(boot_out, type = "bca") 
boot_ci

