###Code to claculate the Bray Cutis index between each forest point and each
##agricultural points & to run the model described in Formula 4
##of the main manuscript & to then extract effect sizes

###Download the data from 10.5281/zenodo.16687021

######### 1. Calculate Bray Curtis index ############
##load packages
library(tidyr)
library(vegan)
library(dplyr)
library(brms)
library(tidybayes)
library(ggplot2)


load("data_for_big_brm_.rda")

# Fill missing combinations (ie no records of a given spp at a given point) before reshaping to avoid issues with missing values
birds_complete <- tidyr::complete(birds.p, point_ID, name_latin, fill = list(CountVisit = 0))
###add observations from each visit get total sum for each point
birds_complete2 <- birds_complete %>%
  dplyr::group_by(point_ID, name_latin, site_type, square_ID) %>%  # Group by point_ID, name_latin, site_type, and square_ID
  dplyr::summarise(CountPoint = sum(CountVisit), .groups = 'drop') # Sum CountVisit and drop grouping
##change fromat of data needed to calculate the index
birds_wide <- tidyr::pivot_wider(
  data = birds_complete2, 
  id_cols = c(point_ID, site_type, square_ID),  # Include the columns I want to keep
  names_from = name_latin, 
  values_from = CountPoint
)
# Check the structure of the new dataframe to ensure it’s been reshaped correctly
str(birds_wide)
head(birds_wide)

#####1.1 Each ZBNF point compared with each forest point######
ZBNF_data <- birds_wide %>% filter(site_type == "ZBNF")
forest_data <- birds_wide %>% filter(site_type == "forest")
# Convert character counts to numeric
ZBNF_counts <- apply(ZBNF_data[, c(-1, -2, -3)], 2, as.numeric)
forest_counts <- apply(forest_data[, c(-1, -2, -3)], 2, as.numeric)
# Create an empty matrix to store results
bray_curtis_results <- matrix(NA, nrow = nrow(ZBNF_counts), ncol = nrow(forest_counts))
# Create vectors to store row and column names
ZBNF_point_IDs <- ZBNF_data$point_ID
forest_point_IDs <- forest_data$point_ID
# Iterate over each ZBNF point and each forest point to calculate Bray-Curtis dissimilarity
for (i in 1:nrow(ZBNF_counts)) {
  for (j in 1:nrow(forest_counts)) {
    bray_curtis_results[i, j] <- vegdist(rbind(ZBNF_counts[i,], forest_counts[j,]), method = "bray")
  }
}
# Convert the results to a dataframe
bray_curtis_results_df <- as.data.frame(bray_curtis_results)
rownames(bray_curtis_results_df) <- ZBNF_point_IDs
colnames(bray_curtis_results_df) <- forest_point_IDs

bray_curtis_results_df$ZBNF_point_ID <- ZBNF_point_IDs

##get long format
bc_zbnf_for <- bray_curtis_results_df %>%
  pivot_longer(!ZBNF_point_ID, names_to = "forest_point", values_to = "bray_index")


######1.2 Compare each agrchemical point to each forest point
chemical_data <- birds_wide %>% subset(site_type == "chemical")
# Convert character counts to numeric
chemical_counts <- apply(chemical_data[, c(-1, -2, -3)], 2, as.numeric) 
# Create an empty matrix to store results
bray_curtis_results <- matrix(NA, nrow = nrow(chemical_counts), ncol = nrow(forest_counts))
# Create vectors to store row and column names
chemical_point_IDs <- chemical_data$point_ID
forest_point_IDs <- forest_data$point_ID
# Iterate over each chemical point and each forest point to calculate Bray-Curtis dissimilarity
for (i in 1:nrow(chemical_counts)) {
  for (j in 1:nrow(forest_counts)) {
    bray_curtis_results[i, j] <- vegdist(rbind(chemical_counts[i,], forest_counts[j,]), method = "bray")
  }
}
# Convert the results to a dataframe
bray_curtis_results_df <- as.data.frame(bray_curtis_results)
rownames(bray_curtis_results_df) <- chemical_point_IDs
colnames(bray_curtis_results_df) <- forest_point_IDs

bray_curtis_results_df$chemical_point_ID <- chemical_point_IDs

##get long format
bc_chemical_for <- bray_curtis_results_df %>%
  pivot_longer(!chemical_point_ID, names_to = "forest_point", values_to = "bray_index")


############ 2. Run brms model ################
load("final_farm_type_birds_model2.rda")
my_covars <- birds.h %>% select(point_ID, site_type, square_ID, weights)
my_covars2 <- unique(my_covars)
my_covars3 <- my_covars2 %>% subset(square_ID != "NACH")

zbnf_index <- bc_zbnf_for %>% select(ZBNF_point_ID, bray_index, forest_point) %>% dplyr::rename(point_ID = ZBNF_point_ID)
chem_index <- bc_chemical_for %>% select(chemical_point_ID, bray_index, forest_point) %>% dplyr::rename(point_ID = chemical_point_ID)

bray_index1 <- bind_rows(
  zbnf_index %>% select(point_ID, bray_index, forest_point),
  chem_index %>% select(point_ID, bray_index, forest_point)
)

bray_index3 <- plyr::join_all(list(bray_index1, my_covars3),
                              by = c("point_ID"), type = 'full')
my_priors <- c(
  prior(normal(0, 10), class = "Intercept"),
  prior(normal(0, 10), class = "b"))

my_brm_bray_hw <- brm(
  bf(
    bray_index | weights(weights) ~ site_type +
      (1 | square_ID / point_ID) + (1 | forest_point),
    phi ~ 1,
    zoi ~ site_type + (1 | square_ID / point_ID) + (1 | forest_point),
    coi ~ site_type + (1 | square_ID / point_ID) + (1 | forest_point)
  ),
  data = bray_index3,
  family = zero_one_inflated_beta(),
  cores = 4,
  chains = 4,
  iter = 3000,
  warmup = 1000,      
  prior = my_priors
)



##########3. Compare bary curtis index between farming practices ########
newdata<- bray_index3 %>% group_by(site_type) %>% tally() 
eff.ZBNF <- add_epred_draws(newdata = filter(newdata, site_type == "ZBNF"), my_brm_bray_hw,
                            re_formula = NA)
eff.CH <- add_epred_draws(newdata = filter(newdata, site_type == "chemical"), my_brm_bray_hw,
                          re_formula = NA)
perc.diff <- data.frame(ZBNF = eff.ZBNF$.epred,
                        chemical = eff.CH$.epred) %>%
  mutate(perc = (ZBNF/chemical -1 )*100,
         PD = (sum(sign(perc) == sign(median(perc)))/n())*100)  %>%
  group_by(PD) %>% median_hdci(perc, .width = .9) %>%
  mutate(Interpretation = case_when(PD > 95 & PD < 97.5 & perc < 0 ~ "PD > 95%, -ve trend",
                                    PD >= 97.5 & perc < 0 ~ "PD > 97.5%, -ve trend",
                                    PD > 95 & PD < 97.5 & perc > 0 ~ "PD > 95%, +ve trend",
                                    PD >= 97.5 & perc > 0 ~ "PD > 97.5%, +ve trend",
                                    PD < 95 ~ "Uncertain"))
ggplot(perc.diff, aes(perc, reorder(-perc), colour = Interpretation)) +
  geom_point()+
  geom_errorbarh(aes(xmin = .lower, xmax= .upper), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("orange", "skyblue", "tomato", "dodgerblue", "grey25")) +
  xlab("Percentage change in abundance\n(chemical to ZBNF transition)") +
  ylab("Taxon") +
  theme_bw() +
  theme(axis.text.y = element_blank(), legend.position = "bottom")

bray.zbnf <- data.frame(ZBNF = eff.ZBNF$.epred,
                        chemical = eff.CH$.epred)%>% median_hdci(ZBNF, .width = .9) 

bray.chem <- data.frame(ZBNF = eff.ZBNF$.epred,
                        chemical = eff.CH$.epred)%>% median_hdci(chemical, .width = .9) 
