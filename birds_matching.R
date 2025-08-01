###Code to conduct matching and obtain matching weights for the analyses on birds

###Download the data from 10.5281/zenodo.16687021

##load packages
library(dplyr)
library(car)
library(MatchIt)
library(cobalt)

##get covariates
covar <- read.csv("covariates_combined.csv")
covar2 <- covar%>% dplyr::rename(Point.ID = point_ID) 
###get data on the % vegetation patches in each square
wfh.square <- read.csv("wfh_squares_final.csv")
covar3 <- full_join(covar2, wfh.square, by = "site_ID")
load("density_data_fitted_only.rda")
birds17 <- birds.counts6%>% dplyr::rename(site_ID = site_ID.x) 
covar4 <- covar3%>% dplyr::rename(point_ID = Point.ID) 
covar5 <- covar4 %>% select(c(point_ID, rain, temp, elevation, perc_quality_wfh_100m,
                              perc_quality))
birds18 <- plyr::join_all(list(birds17, covar5),
                          by = c("point_ID"), type = 'full')
birds19 <- birds18 %>% select(c(point_ID, name_latin, site_ID, site_type.x, farm_type,
                                rain...83, temp...84, 
                                elevation...85, perc_quality_wfh_100m...86,
                                perc_quality...87, TotalClustersVisit, TotalClusterPoint,
                                effArea2))
birds20<- birds19%>% dplyr::rename(rain = rain...83, temp = temp...84, 
                                   elevation =elevation...85, perc_quality_100 =perc_quality_wfh_100m...86,
                                   perc_quality = perc_quality...87, ) 
#####collinearity checks
modt <- car::vif(lm(TotalClusterPoint ~ site_type.x + perc_quality + elevation+ temp + rain,
                    data = birds20))

modt ##all GVIF^(1/(2*Df)) should be below 2; rainfall is greater so exclude

modt <- car::vif(lm(TotalClusterPoint ~ site_type.x + perc_quality + elevation+ temp,
                    data = birds20))
modt

########################MATCHING ############################
birds21 <- birds20 %>%
  dplyr::mutate(farm_type = case_when(site_type.x == "ZBNF" ~ 1,site_type.x == "chemical" ~0)) %>%
  dplyr::group_by(farm_type) %>% ungroup()

birds22 <- birds21 %>%
  group_by(site_ID, point_ID) %>%
  mutate(elevation2 = sample(60:260, 1)) %>%
  distinct(point_ID, .keep_all = TRUE) %>%
  ungroup()
#####matching done at the point level (as no matching covariates are at the visit level)
###check covariate balance prior to matching
mod1pf <- matchit(farm_type ~ perc_quality + elevation2+ temp, ##y needs to be in 0 and 1
                  data = birds22,
                  method = NULL, distance = "mahalanobis", estimand = "ATT")
mod1pf
summary(mod1pf)
###evaluate balanace using cobalt
wfarm <- bal.tab(mod1pf, disp = c("means", "sds"), un = TRUE, 
                 stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                 thresholds = c(m = .25),
                 disp.means = TRUE, 
                 s.d.denom = "treated") 

wfarm
love.plot(wfarm, binary = "std", thresholds = c(m = .25)) ###suitability and % WFH are unmatched, rest fine
######extract SDMs for comparison between matching methods
colSDM <- summary(mod1pf)[["sum.all"]][, "Std. Mean Diff."]
SDMun <- as.data.frame(colSDM)
SDMun$covariate <- rownames(SDMun)
SDMun$method <- "unadjusted"

###Nearest neighbour, Mahalanobis distance, with replacement
mod_near_mahal <- matchit(farm_type ~ perc_quality + elevation2+ temp, ##y needs to be in 0 and 1
                          data = birds22,
                          method = "nearest", distance = "mahalanobis", estimand = "ATT", replace = TRUE, pop.size = 1000)
mod_near_mahal
summary(mod_near_mahal)
w_near_mahal <- bal.tab(mod_near_mahal, disp = c("means", "sds"), un = TRUE, 
                        stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                        thresholds = c(m = .25),
                        disp.means = TRUE, 
                        s.d.denom = "treated")
w_near_mahal

love.plot(mod_near_mahal, binary = "std", thresholds = c(m = .25)) ##%WFH and suitability still unmatched
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_near_mahal)[["sum.matched"]][, "Std. Mean Diff."]
SDMnm <- as.data.frame(colSDM)
SDMnm$covariate <- rownames(SDMnm)
SDMnm$method <- "near_mahal_repl"

###Nearest neighbour, Mahalanobis distance, without replacement
mod_near_mahal2 <- matchit(farm_type ~ perc_quality + elevation2+ temp, ##y needs to be in 0 and 1
                           data = birds22,
                           method = "nearest", distance = "mahalanobis", estimand = "ATT", replace = FALSE, pop.size = 1000)
mod_near_mahal2
summary(mod_near_mahal2)
w_near_mahal2 <- bal.tab(mod_near_mahal2, disp = c("means", "sds"), un = TRUE, 
                         stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                         thresholds = c(m = .25),
                         disp.means = TRUE, 
                         s.d.denom = "treated")
w_near_mahal2

love.plot(mod_near_mahal2, binary = "std", thresholds = c(m = .25)) ##%WFH and suitability still unmatched
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_near_mahal2)[["sum.matched"]][, "Std. Mean Diff."]
SDMnm2 <- as.data.frame(colSDM)
SDMnm2$covariate <- rownames(SDMnm2)
SDMnm2$method <- "near_mahal_no_repl"

###Nearest neighbour, Propensity score distance, with replacement
mod_near_prop <- matchit(farm_type ~ perc_quality + elevation2+ temp, ##y needs to be in 0 and 1
                         data = birds22,
                         method = "nearest", distance = "glm", estimand = "ATT", replace = TRUE, pop.size = 1000)
mod_near_prop
summary(mod_near_prop)
w_near_prop <- bal.tab(mod_near_prop, disp = c("means", "sds"), un = TRUE, 
                       stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                       thresholds = c(m = .25),
                       disp.means = TRUE, 
                       s.d.denom = "treated") 
w_near_prop
love.plot(mod_near_prop, binary = "std", thresholds = c(m = .25)) ##still as above, even worse than original
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_near_prop)[["sum.matched"]][, "Std. Mean Diff."]
SDMnp <- as.data.frame(colSDM)
SDMnp$covariate <- rownames(SDMnp)
SDMnp$method <- "near_prop_repl"
###Nearest neighbour, Propensity score distance, without replacement
mod_near_prop2 <- matchit(farm_type ~ perc_quality + elevation2+ temp, ##y needs to be in 0 and 1
                          data = birds22,
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
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_near_prop2)[["sum.matched"]][, "Std. Mean Diff."]
SDMnp2 <- as.data.frame(colSDM)
SDMnp2$covariate <- rownames(SDMnp2)
SDMnp2$method <- "near_prop_no_repl"

####genetic matching, with replacement
mod_genet <- matchit(farm_type ~ perc_quality + elevation2+ temp, ##y needs to be in 0 and 1
                     data = birds22, pop.size=1000,
                     method = "genetic", estimand = "ATT", replace = TRUE)
mod_genet
summary(mod_genet)

w_genet <- bal.tab(mod_genet, disp = c("means", "sds"), un = TRUE, 
                   stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                   thresholds = c(m = .25),
                   disp.means = TRUE, 
                   s.d.denom = "treated") 
w_genet
love.plot(mod_genet, binary = "std", thresholds = c(m = .25)) ##everything apart from distance is balanced, ESS is only 20
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_genet)[["sum.matched"]][, "Std. Mean Diff."]
SDMg <- as.data.frame(colSDM)
SDMg$covariate <- rownames(SDMg)
SDMg$method <- "genet_repl"
####genetic matching, without replacement
mod_genet2 <- matchit(farm_type ~ perc_quality + elevation2+ temp, ##y needs to be in 0 and 1
                      data = birds22, pop.size=1000,
                      method = "genetic", estimand = "ATT", replace = FALSE)
mod_genet2
summary(mod_genet2)
w_genet2 <- bal.tab(mod_genet2, disp = c("means", "sds"), un = TRUE, 
                    stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                    thresholds = c(m = .25),
                    disp.means = TRUE, 
                    s.d.denom = "treated") 
w_genet2
love.plot(mod_genet2, binary = "std", thresholds = c(m = .25)) ##everything apart from distance is balanced, ESS is only 20
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_genet2)[["sum.matched"]][, "Std. Mean Diff."]
SDMg2 <- as.data.frame(colSDM)
SDMg2$covariate <- rownames(SDMg2)
SDMg2$method <- "genet_no_repl"

#######genetic matching, with replacement, also matched on propensity score
ps_model <- glm(farm_type ~ perc_quality + elevation2+ temp, ##y needs to be in 0 and 1
                data = birds22, family = "binomial")
propensity_scores <- predict(ps_model, type = "response")

mod_genet_prop <- matchit(farm_type ~ perc_quality + elevation2+ temp, ##y needs to be in 0 and 1
                          data = birds22, pop.size=1000,
                          method = "genetic", estimand = "ATT", replace = TRUE)
mod_genet_prop
summary(mod_genet_prop)
w_genet_prop <- bal.tab(mod_genet_prop, disp = c("means", "sds"), un = TRUE, 
                        stats = c("mean.diffs", "variance.ratios"), binary = "std", continuous = "std", ##std for standardised (not raw) mean difference
                        thresholds = c(m = .25),
                        disp.means = TRUE, 
                        s.d.denom = "treated") 
w_genet_prop
love.plot(mod_genet_prop, binary = "std", thresholds = c(m = .25)) ##everything apart from distance is balanced, ESS is only 20
######extract SDMs for comparison between matching methods
colSDM <- summary(mod_genet_prop)[["sum.matched"]][, "Std. Mean Diff."]
SDMg3 <- as.data.frame(colSDM)
SDMg3$covariate <- rownames(SDMg3)
SDMg3$method <- "genet_repl_prop"

##Full matching, mahalanobis 
mod_full_mahal <- matchit(farm_type ~ perc_quality + elevation2+ temp, ##y needs to be in 0 and 1
                          data = birds22,
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

##Full matching,propensity score
mod_full_prop <- matchit(farm_type ~ perc_quality + elevation2+ temp, ##y needs to be in 0 and 1
                         data = birds22,
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


##crate dataset containing SDMs from all matching methods######

sdms.data <- purrr::reduce(
  list(SDMun, SDMnm, SDMnm2, SDMnp, SDMnp2, SDMg, SDMg2, SDMg3, SDMfm, SDMfp),
  function(left, right) {
    dplyr::full_join(left, right)
  }
)

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
  perc_quality = "Proportion of habitat patches",
  elevation2 ="Elevation",
  temp = "Annual average temperature"
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
      "full_mahal" = "Full matching, Mahalanobis (N=52/52, u=0)", ##ess=28 ######choose
      "full_prop" = "Full matching, Propensity Score (N=52/52, u=0)", ##ess=26
      "genet_no_repl" = "Genetic matching, without replacement (N=52/52, u=2)", ##
      "genet_repl" = "Genetic matching, with replacement (N=29/52, u=1)", ##ess 17
      "genet_repl_prop" = "Genetic matching, with replacement and propensity score (N=31/52, u=0)", ##ess 17
      "near_mahal_no_repl" = "Nearest Neighbour, Mahalanobis, without replacement (N=52/52, u=2)", ##
      "near_mahal_repl" = "Nearest Neighbour, Mahalanobis, with replacement (N=27/52, u=0)",##17 ess
      "near_prop_no_repl" = "Nearest Neighbour, Propensity Score, without replacement (N=52/52, u=2)",##ess 20
      "near_prop_repl" = "Nearest Neighbour, Propensity Score, with replacement (N=24/52, u=5)", ##ess 20
      "unadjusted" = "Unadjusted (N=52/52, u=2)"
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


###full matching mahalanobis is best - extract
#Extract matched dataset
md <- match.data(mod_full_mahal)
write.csv(md, file ="match_weights_point_farm_birds.csv")
