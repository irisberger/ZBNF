###Code to estimate the effect of increasing yield on birds (at the species and niche level) and to construct density-yield curves
##(this multi-species count model is described in Formula 2 in the main manuscript)
##The code below is for the yield model, replace with profit to estimate the effect
##of increasing profit and to construct density-profit curves

###Download the data from 10.5281/zenodo.16687021

##load packages
library(dplyr)
library(brms)
library(tidybayes)
library(tidyverse)

###import data
load("data_for_big_brm_.rda")

##for profit analyses also import profit data
#profit <- read.csv("profit_squares.csv")
#prof <- merge(profit1, birds.p, by = "square_ID", all.y  = T) ##
#prof$prof_square_scaled <- scale(prof$profit_square)

#### 1. Fit the model#####
##scale yield data
##zero yield is -0.8201702 when scaled; zero profit is -0.5085231 when scaled
yield_0_scaled <- -0.8201702

birds1 <- birds.p %>% ##replace below with prof_square_scaled as appropriate
  mutate(zbnf = ifelse(site_type == "ZBNF", "1", "0"),
         chem = ifelse(site_type == "chemical", "1", "0"),
         zbnfyield = ifelse(site_type == "ZBNF", yield_square_scaled, yield_0_scaled),
         chemyield = ifelse(site_type == "chemical", yield_square_scaled, yield_0_scaled),)

###set priors
my_priors <- c(
  prior(normal(0, 10), class = "Intercept"),
  prior(normal(0, 10), class = "b")
)

##brms model
my_brm_forest <- brm(
  bf(CountVisit ~ 1+ Trophic.Niche * (zbnf + chem + zbnfyield + chemyield)+ elevation_scaled + temp_scaled +
       (1 + zbnf + chem + zbnfyield + chemyield | name_latin) 
     + (1 | square_ID/point_ID) + (offset(effArea2)),
     zi ~ zbnf + chem + (1 | name_latin)
  ),
  data = birds1,
  family = zero_inflated_poisson(),
  chains = 4,
  cores = 4,
  iter = 3000,
  warmup = 1000,
  prior = my_priors
)


######## 2. Extract density-yield/profit curves (at the guild-level) #################
## 2.1 Agrichemical farming 
yield_square_sd <- sd(birds1$yield_square)
yield_square_mean <- mean(birds1$yield_square)
##new 
new_chem_dat <- expand.grid(
  chem = 1,
  zbnf = 0,
  Trophic.Niche = unique(birds1$Trophic.Niche),
  chemyield = seq(min(birds1$chemyield), max(birds1$chemyield), length.out = 100),
  zbnfyield = -0.8201702, ##which is 0 yield when unscaled
  effArea2 = 1,
  elevation_scaled = 0
) %>% mutate(
  yield_chem_backtransformed = chemyield * yield_square_sd + yield_square_mean)

preds.chem.niche.DY <- tidybayes::add_epred_draws(newdata = new_chem_dat, 
                                                  my_brm_forest4, dpar = "mu", ndraws = 300, 
                                                  re_formula = NA) %>%
  group_by(Trophic.Niche, chemyield) %>%
  median_hdci(.epred, .width = .9)
###predicted value is the number of individuals per point count (10min point count, r=100m)
###to express counts per 500mx500m square instead:
pred.chem.DY <- preds.chem.niche.DY %>%
  mutate(
    # Calculate the area of the circle
    circle_area = pi * 100^2,
    # Calculate the area of the square
    square_area = 500^2,
    # Convert median estimate
    chem_count_square_med = (.epred / circle_area) * square_area,
    # Convert lower bound estimate
    chem_count_square_lower = (.lower / circle_area) * square_area,
    # Convert upper bound estimate
    chem_count_square_upper = (.upper / circle_area) * square_area
  )
###back-transform yield data
yield_square_sd <- sd(birds1$yield_square)
yield_square_mean <- mean(birds1$yield_square)
pred.chem.DY2 <- pred.chem.DY %>%
  mutate(yield_square_backtransformed = chemyield * yield_square_sd + yield_square_mean)

###Plot DY curves for each trophic niche
trophic_niches <- unique(pred.chem.DY2$Trophic.Niche)
for (niche in trophic_niches) {
  niche_data <- pred.chem.DY2 %>% filter(Trophic.Niche == niche)
  p <- ggplot(niche_data, aes(yield_square_backtransformed, chem_count_square_med)) +
    geom_line(size = 4) +
    geom_ribbon(aes(ymin = chem_count_square_lower, ymax = chem_count_square_upper), fill = "slateblue4", alpha = 0.5) +
    ggtitle(paste("Trophic Niche:", niche)) +
    theme_minimal() +
    labs(x = "Yield Square", y = "Individuals per square") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )
  print(p)
}

## 2.2 ZBNF
yield_square_sd <- sd(birds1$yield_square)
yield_square_mean <- mean(birds1$yield_square)
new_zbnf_dat <- expand.grid(
  chem = 0,
  zbnf = 1,
  Trophic.Niche = unique(birds1$Trophic.Niche),
  zbnfyield = seq(min(birds1$zbnfyield), max(birds1$zbnfyield), length.out = 100),
  chemyield = -0.8201702, ##which is 0 yield when unscaled
  effArea2 = 1,
  elevation_scaled = 0
) %>% mutate(
  yield_chem_backtransformed = zbnfyield * yield_square_sd + yield_square_mean)

preds.zbnf.niche.DY <- tidybayes::add_epred_draws(newdata = new_zbnf_dat, 
                                                  my_brm_forest4, dpar = "mu", ndraws = 300, 
                                                  re_formula = NA) %>%
  group_by(Trophic.Niche, zbnfyield) %>%
  median_hdci(.epred, .width = .9)
###predicted value is the number of individulas per point count (10min point count, r=100m)
###to express counts per 500mx500m square instead:
pred.zbnf.DY <- preds.zbnf.niche.DY %>%
  mutate(
    # Calculate the area of the circle
    circle_area = pi * 100^2,
    # Calculate the area of the square
    square_area = 500^2,
    # Convert median estimate
    zbnf_count_square_med = (.epred / circle_area) * square_area,
    # Convert lower bound estimate
    zbnf_count_square_lower = (.lower / circle_area) * square_area,
    # Convert upper bound estimate
    zbnf_count_square_upper = (.upper / circle_area) * square_area
  )

###back-transform yield data
yield_square_sd <- sd(birds1$yield_square)
yield_square_mean <- mean(birds1$yield_square)
pred.zbnf.DY2 <- pred.zbnf.DY %>%
  mutate(yield_square_backtransformed = zbnfyield * yield_square_sd + yield_square_mean)

##get DY curve for each trophic niche
trophic_niches <- unique(pred.zbnf.DY2$Trophic.Niche)
for (niche in trophic_niches) {
  niche_data2 <- pred.zbnf.DY2 %>% filter(Trophic.Niche == niche)
  p <- ggplot(niche_data2, aes(yield_square_backtransformed, zbnf_count_square_med)) +
    geom_line(size = 4) +
    geom_ribbon(aes(ymin = zbnf_count_square_lower, ymax = zbnf_count_square_upper), fill = "khaki", alpha = 0.5) +
    ggtitle(paste("Trophic Niche:", niche)) +
    theme_minimal() +
    labs(x = "Yield Square", y = "Individuals per square") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )
  print(p)
}

##### 3. Estimate effect of yield/profit for each farming system
#### Species-wise effects of yield in ZBNF systems
sp1.zbnf <- birds1 %>% group_by(Trophic.Niche, name_latin, zbnf, chem) %>% 
  summarise(zbnfyield = 1, chemyield = 0, effArea2 = 1, elevation_scaled = 0) ###should elevation be 1 or 0?

sp0.zbnf <- birds1 %>% group_by(Trophic.Niche, name_latin, zbnf, chem) %>% 
  summarise(zbnfyield = 0, chemyield = 0, effArea2 = 1, elevation_scaled = 0)

preds.zbnf1 <- tidybayes::add_linpred_draws(newdata = filter(sp1.zbnf, zbnf == "1"), 
                                            my_brm_forest, dpar = "mu",
                                            re_formula = ~(1 + zbnf + chem + zbnfyield + chemyield | name_latin)) 

preds.zbnf0 <- tidybayes::add_linpred_draws(newdata = filter(sp0.zbnf, zbnf == "1"), 
                                            my_brm_forest, dpar = "mu",
                                            re_formula = ~(1 + zbnf + chem + zbnfyield + chemyield | name_latin)) 
yield_percent_zbnf <- preds.zbnf1 %>% ungroup() %>%
  mutate(baseline0 = preds.zbnf0$mu, diff = mu - baseline0,
         perc_diff = (exp(diff) - 1)*100) %>%
  group_by(Trophic.Niche, name_latin) %>%
  mutate(PD = (sum(sign(perc_diff) == sign(median(perc_diff)))/n())*100) %>%
  group_by(name_latin, Trophic.Niche, PD) %>% median_hdci(perc_diff, .width = .9) %>%
  mutate(Interpretation = case_when(PD > 95 & PD < 97.5 & perc_diff < 0 ~ "PD > 95%, -ve trend",
                                    PD >= 97.5 & perc_diff < 0 ~ "PD > 97.5%, -ve trend",
                                    PD > 95 & PD < 97.5 & perc_diff > 0 ~ "PD > 95%, +ve trend",
                                    PD >= 97.5 & perc_diff > 0 ~ "PD > 97.5%, +ve trend",
                                    PD < 95 ~ "Uncertain"))
ggplot(yield_percent_zbnf, aes(perc_diff, reorder(name_latin, -perc_diff), colour = Interpretation)) +
  geom_point() +
  geom_errorbar(aes(xmin = .lower, xmax = .upper)) +
  scale_color_manual(values = c("orange", "skyblue", "tomato", "dodgerblue", "grey25")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #facet_wrap(~ site_type, scales = "free") +
  coord_cartesian(xlim = c(-100, 500))+
  theme_bw() +
  labs(x = "Percentage change in abundance with increasing yield", y = "Trophic Niche", title = "Effect of yield by farming practice")

######### Species-wise effects of yield in agrichemical systems
sp1.chem <- birds1 %>% group_by(Trophic.Niche, name_latin, zbnf, chem) %>% 
  summarise(zbnfyield = 0, chemyield = 1, effArea2 = 1, elevation_scaled = 0) ###should elevation be 1 or 0?

sp0.chem <- birds1 %>% group_by(Trophic.Niche, name_latin, zbnf, chem) %>% 
  summarise(zbnfyield = 0, chemyield = 0, effArea2 = 1, elevation_scaled = 0)

preds.chem1 <- tidybayes::add_linpred_draws(newdata = filter(sp1.chem, chem == "1"), 
                                            my_brm_forest, dpar = "mu",
                                            re_formula = ~(1 + zbnf + chem + zbnfyield + chemyield | name_latin)) 

preds.chem0 <- tidybayes::add_linpred_draws(newdata = filter(sp0.chem, chem == "1"), 
                                            my_brm_forest, dpar = "mu",
                                            re_formula = ~(1 + zbnf + chem + zbnfyield + chemyield | name_latin)) 
yield_percent_chem <- preds.chem1 %>% ungroup() %>%
  mutate(baseline0 = preds.chem0$mu, diff = mu - baseline0,
         perc_diff = (exp(diff) - 1)*100) %>%
  group_by(Trophic.Niche, name_latin) %>%
  mutate(PD = (sum(sign(perc_diff) == sign(median(perc_diff)))/n())*100) %>%
  group_by(name_latin, Trophic.Niche, PD) %>% median_hdci(perc_diff, .width = .9) %>%
  mutate(Interpretation = case_when(PD > 95 & PD < 97.5 & perc_diff < 0 ~ "PD > 95%, -ve trend",
                                    PD >= 97.5 & perc_diff < 0 ~ "PD > 97.5%, -ve trend",
                                    PD > 95 & PD < 97.5 & perc_diff > 0 ~ "PD > 95%, +ve trend",
                                    PD >= 97.5 & perc_diff > 0 ~ "PD > 97.5%, +ve trend",
                                    PD < 95 ~ "Uncertain"))
ggplot(yield_percent_chem, aes(perc_diff, reorder(name_latin, -perc_diff), colour = Interpretation)) +
  geom_point() +
  geom_errorbar(aes(xmin = .lower, xmax = .upper)) +
  scale_color_manual(values = c("orange", "tomato", "dodgerblue", "grey25")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #facet_wrap(~ site_type, scales = "free") +
  coord_cartesian(xlim = c(-100, 500))+
  theme_bw() +
  labs(x = "Percentage change in abundance with increasing yield", y = "Trophic Niche", title = "Effect of yield by farming practice")

#### Niche-wise effects of yield in ZBNF systems##########
sp1.zbnf.niche <- birds1 %>% group_by(Trophic.Niche, zbnf, chem) %>% 
  summarise(zbnfyield = 1, chemyield = 0, effArea2 = 1, elevation_scaled = 0) ###should elevation be 1 or 0?

sp0.zbnf.niche <- birds1 %>% group_by(Trophic.Niche, zbnf, chem) %>% 
  summarise(zbnfyield = 0, chemyield = 0, effArea2 = 1, elevation_scaled = 0)

preds.zbnf1.niche <- tidybayes::add_linpred_draws(newdata = filter(sp1.zbnf.niche, zbnf == "1"), 
                                                  my_brm_forest, dpar = "mu", re_formula =  NA) 

preds.zbnf0.niche <- tidybayes::add_linpred_draws(newdata = filter(sp0.zbnf.niche, zbnf == "1"), 
                                                  my_brm_forest, dpar = "mu",
                                                  re_formula = NA) 
yield_percent_zbnf_niche <- preds.zbnf1.niche %>% ungroup() %>%
  mutate(baseline0 = preds.zbnf0.niche$mu, diff = mu - baseline0,
         perc_diff = (exp(diff) - 1)*100) %>%
  group_by(Trophic.Niche) %>%
  mutate(PD = (sum(sign(perc_diff) == sign(median(perc_diff)))/n())*100) %>%
  group_by(Trophic.Niche, PD) %>% median_hdci(perc_diff, .width = .9) %>%
  mutate(Interpretation = case_when(PD > 95 & PD < 97.5 & perc_diff < 0 ~ "PD > 95%, -ve trend",
                                    PD >= 97.5 & perc_diff < 0 ~ "PD > 97.5%, -ve trend",
                                    PD > 95 & PD < 97.5 & perc_diff > 0 ~ "PD > 95%, +ve trend",
                                    PD >= 97.5 & perc_diff > 0 ~ "PD > 97.5%, +ve trend",
                                    PD < 95 ~ "Uncertain"))
ggplot(yield_percent_zbnf_niche, aes(perc_diff, reorder(Trophic.Niche, -perc_diff), colour = Interpretation)) +
  geom_point() +
  geom_errorbar(aes(xmin = .lower, xmax = .upper)) +
  scale_color_manual(values = c("orange", "skyblue", "tomato", "dodgerblue", "grey25")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #facet_wrap(~ site_type, scales = "free") +
  coord_cartesian(xlim = c(-100, 500))+
  theme_bw() +
  labs(x = "Percentage change in abundance with increasing yield", y = "Trophic Niche", title = "Effect of yield by farming practice")

##### #### Niche-wise effects of yield in agrichemical systems##########
sp1.chem.niche <- birds1 %>% group_by(Trophic.Niche, zbnf, chem) %>% 
  summarise(zbnfyield = 0, chemyield = 1, effArea2 = 1, elevation_scaled = 0) ###should elevation be 1 or 0?

sp0.chem.niche <- birds1 %>% group_by(Trophic.Niche, zbnf, chem) %>% 
  summarise(zbnfyield = 0, chemyield = 0, effArea2 = 1, elevation_scaled = 0)

preds.chem1.niche <- tidybayes::add_linpred_draws(newdata = filter(sp1.chem.niche, chem == "1"), 
                                                  my_brm_forest, dpar = "mu",
                                                  re_formula = NA) 

preds.chem0.niche <- tidybayes::add_linpred_draws(newdata = filter(sp0.chem.niche, chem == "1"), 
                                                  my_brm_forest, dpar = "mu",
                                                  re_formula = NA) 
## percentage change in abnundance for each 1 unit increase in yield
yield_percent_chem_niche <- preds.chem1.niche %>% ungroup() %>%
  mutate(baseline0 = preds.chem0.niche$mu, diff = mu - baseline0,
         perc_diff = (exp(diff) - 1)*100) %>%
  group_by(Trophic.Niche) %>%
  mutate(PD = (sum(sign(perc_diff) == sign(median(perc_diff)))/n())*100) %>%
  group_by(Trophic.Niche, PD) %>% median_hdci(perc_diff, .width = .9) %>%
  mutate(Interpretation = case_when(PD > 95 & PD < 97.5 & perc_diff < 0 ~ "PD > 95%, -ve trend",
                                    PD >= 97.5 & perc_diff < 0 ~ "PD > 97.5%, -ve trend",
                                    PD > 95 & PD < 97.5 & perc_diff > 0 ~ "PD > 95%, +ve trend",
                                    PD >= 97.5 & perc_diff > 0 ~ "PD > 97.5%, +ve trend",
                                    PD < 95 ~ "Uncertain"))
ggplot(yield_percent_chem_niche, aes(perc_diff, reorder(Trophic.Niche, -perc_diff), colour = Interpretation)) +
  geom_point() +
  geom_errorbar(aes(xmin = .lower, xmax = .upper)) +
  scale_color_manual(values = c("tomato", "grey25")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #facet_wrap(~ site_type, scales = "free") +
  coord_cartesian(xlim = c(-100, 500))+
  theme_bw() +
  labs(x = "Percentage change in abundance with increasing yield", y = "Trophic Niche", title = "Effect of yield by farming practice")

