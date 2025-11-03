###Code to estimate the impact of the ZBNF transition on birds 
##(this multi-species count model is described in Formula 1 in the main manuscript)

###Download the data from 10.5281/zenodo.16687021

##load packages
library(dplyr)
library(brms)
library(tidybayes)
library(tidyverse)

#### 1. Import and prepare data
##get bird data
load("density_data_fitted_only.rda")
##only use zbnf & chemical farming data for this model
birds1 <- data.sum4 %>% subset(site_type != "forest")
##calculate the total number of individulas per species per visit 
bird2 <- birds1 %>%
  group_by(point_ID, point_count_repeat_number, name_latin) %>%
  mutate(CountVisit = sum(size)) %>%
  distinct(point_ID, point_count_repeat_number, name_latin, .keep_all = TRUE)
###remove rows for species never recorded in agriculture
bird3 <- bird2 %>%
  group_by(name_latin) %>%
  mutate(CountTotal = sum(CountVisit)) 
bird4 <- bird3 %>% filter(CountTotal != 0)
birds5 <- bird4 %>% select(point_ID, name_latin, site_type, square_ID, CountVisit, effArea, point_count_repeat_number)

##get the effective area surveys for non-observations (i.e. when a species was not recorded during a given visit)
####ideally use information from the same species during a different visit to the same point,
##if unavailable use it from the same species in the same square, and if this is unavailable then from the same species
birds12 <- birds5 %>%
  dplyr::group_by(point_ID, name_latin) %>%
  dplyr::mutate(effArea2 = ifelse(is.na(effArea), mean(effArea, na.rm = TRUE), effArea)) %>%
  ungroup()
birds13 <- birds12 %>%
  dplyr::group_by(square_ID, name_latin) %>%
  dplyr::mutate(effArea2 = ifelse(is.na(effArea), mean(effArea, na.rm = TRUE), effArea)) %>%
  ungroup()
birds14 <- birds13 %>%
  dplyr::group_by(name_latin) %>%
  dplyr::mutate(effArea2 = ifelse(is.na(effArea), mean(effArea, na.rm = TRUE), effArea)) %>%
  ungroup()

###import matching weights
mweights <- read.csv("match_weights_point_farm_birds.csv")
mweights2 <- mweights %>% select(point_ID, weights)
###import data on species' forest dependency
fdepend <- read.csv("forest_dependency_BL.csv")
fdepend2 <- fdepend %>% select(name_latin, Forest.dependency)
fdepend2$name_latin <- gsub(" ", "_", fdepend2$name_latin)
#combine datasets
birds.a <- plyr::join_all(list(birds14, fdepend2),
                          by = c("name_latin"), type = 'full')
birds.b <- birds.a %>%
  filter(!is.na(site_type))
##combine with matching weights
birds.c <- plyr::join_all(list(birds.b, mweights2),
                          by = c("point_ID"), type = 'full')
###some forest dependency data missing - use medium for now
birds.c$Forest.dependency <- ifelse(is.na(birds.c$Forest.dependency) | birds.c$Forest.dependency == "", "Medium", birds.c$Forest.dependency)

###import trophic niche data
niche <- read.csv("BirdLife_traits.csv")
niche2 <- niche %>% select(Species1, Trophic.Niche)
niche2$name_latin <- gsub(" ", "_", niche2$Species1)
birds.d <- plyr::join_all(list(birds.c, niche2),
                          by = c("name_latin"), type = 'full')
birds.e <- birds.d %>%
  filter(!is.na(site_type))
##add 'aquatic predator' to the vertivore group
birds.f <- birds.e %>%
  mutate(Trophic.Niche = if_else(Trophic.Niche == "Aquatic predator", "Vertivore", Trophic.Niche))
unique(birds.f$name_latin[is.na(birds.f$Trophic.Niche)])
birds.g <- birds.f %>%
  mutate(Trophic.Niche = if_else(name_latin == "Argya_affinis" | name_latin == "Argya_striata", "Invertivore",
                                 if_else(name_latin == "Ortygornis_pondicerianus", "Omnivore", Trophic.Niche)))

birds.h <- birds.g %>% select(point_ID, name_latin, site_type, square_ID, CountVisit, effArea2, weights, Trophic.Niche, Forest.dependency)

############2. Fit hierarchical count model #####################
###set priors
my_priors <- c(
  prior(normal(0, 10), class = "Intercept"),
  prior(normal(0, 10), class = "b")
)
##fit brms model
my_brm <- brm(
  bf(CountVisit | weights(weights) ~ site_type * Trophic.Niche + ##site type is farming practice (ZBNF or agrichemical)
       (1 + site_type | name_latin) + ##name_latin is the species name
       (1 | square_ID / point_ID) +
       (offset(effArea2)), ##effArea2 is on the log scale
     zi ~ site_type + (1|name_latin)
  ),
  data = birds.h,
  family = zero_inflated_poisson(),
  chains = 4,
  iter = 3000,
  warmup = 1000,
  prior = my_priors
)


########### Extract species-level effects############
newdata.sp.niche1 <- birds.h %>% group_by(name_latin, site_type, Trophic.Niche) %>% tally() %>%
  mutate(effArea2 = 1)  

eff.sp.niche.ZBNF1 <- add_epred_draws(newdata = filter(newdata.sp.niche1, site_type == "ZBNF"), my_brm, 
                                      re_formula = ~(1 + site_type| name_latin))
eff.sp.niche.CH1 <- add_epred_draws(newdata = filter(newdata.sp.niche1, site_type == "chemical"), my_brm, 
                                    re_formula = ~(1 + site_type | name_latin))

sp.niche.perc.diff1 <- data.frame(name_latin = eff.sp.niche.ZBNF1$name_latin, 
                                  Trophic.Niche = eff.sp.niche.ZBNF1$Trophic.Niche,
                                  ZBNF = eff.sp.niche.ZBNF1$.epred,
                                  chemical = eff.sp.niche.CH1$.epred) %>%
  group_by(name_latin, Trophic.Niche) %>%
  mutate(perc = (ZBNF/chemical -1 )*100,
         PD = (sum(sign(perc) == sign(median(perc)))/n())*100)  %>%
  group_by(name_latin, Trophic.Niche, PD) %>% median_hdci(perc, .width = .9) %>%
  mutate(Interpretation = case_when(PD > 95 & PD < 97.5 & perc < 0 ~ "PD > 95%, -ve trend",
                                    PD >= 97.5 & perc < 0 ~ "PD > 97.5%, -ve trend",
                                    PD > 95 & PD < 97.5 & perc > 0 ~ "PD > 95%, +ve trend",
                                    PD >= 97.5 & perc > 0 ~ "PD > 97.5%, +ve trend",
                                    PD < 95 ~ "Uncertain"))

ggplot(sp.niche.perc.diff1, aes(perc, reorder(name_latin, -perc), colour = Interpretation)) +
  geom_point(size = 4) +  # Adjust the size here
  geom_errorbarh(aes(xmin = .lower, xmax= .upper), height = 0, size = 1) +  # Adjust the size here
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("orange", "skyblue", "tomato", "dodgerblue", "grey25")) +
  xlab("Percentage change in abundance\n(chemical farming to ZBNF transition)") +
  ylab("Species") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")


#### Extract niche-level effects#####
newdata.niche1 <- birds.h %>% group_by(site_type, Trophic.Niche) %>% tally() %>%
  mutate(effArea2 = 1)  
eff.niche.ZBNF1 <- add_epred_draws(newdata = filter(newdata.niche1, site_type == "ZBNF"), my_brm, 
                                   re_formula = NA)
eff.niche.CH1 <- add_epred_draws(newdata = filter(newdata.niche1, site_type == "chemical"), my_brm, 
                                 re_formula = NA)

niche.perc.diff1 <- data.frame(Trophic.Niche = eff.niche.ZBNF1$Trophic.Niche,
                               ZBNF = eff.niche.ZBNF1$.epred,
                               chemical = eff.niche.CH1$.epred) %>%
  group_by(Trophic.Niche) %>%
  mutate(perc = (ZBNF/chemical -1 )*100,
         PD = (sum(sign(perc) == sign(median(perc)))/n())*100)  %>%
  group_by(Trophic.Niche, PD) %>% median_hdci(perc, .width = .9) %>%
  mutate(Interpretation = case_when(PD > 95 & PD < 97.5 & perc < 0 ~ "PD > 95%, -ve trend",
                                    PD >= 97.5 & perc < 0 ~ "PD > 97.5%, -ve trend",
                                    PD > 95 & PD < 97.5 & perc > 0 ~ "PD > 95%, +ve trend",
                                    PD >= 97.5 & perc > 0 ~ "PD > 97.5%, +ve trend",
                                    PD < 95 ~ "Uncertain"))

ggplot(niche.perc.diff1, aes(perc, reorder(Trophic.Niche, -perc), colour = Interpretation)) +
  geom_point(size = 4) +  # Adjust the size here
  geom_errorbarh(aes(xmin = .lower, xmax= .upper), height = 0, size = 1) +  # Adjust the size here
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("dodgerblue", "grey25")) +
  xlab("Percentage change in abundance\n(chemical farming to ZBNF transition)") +
  ylab("Trophic niche") +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        legend.position = "bottom")
