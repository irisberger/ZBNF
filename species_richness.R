##Code to estimate the species richness of species of conservation importance in 
##ZBNF, agrichemical, and forest systems.
##The model is described in Formula 3 in the main manuscript

###Download the data from 10.5281/zenodo.16687021

##load packages
library(dplyr)
library(brms)
library(ggplot2)
library(tidybayes)

### 1. Prepare data####
##import bird data
load("density_data_fitted_only.rda")
bird2 <- data.sum4 %>%
  dplyr::group_by(point_ID, point_count_repeat_number, name_latin) %>%
  dplyr::mutate(CountVisit = sum(size)) %>%
  dplyr::distinct(point_ID, point_count_repeat_number, name_latin, .keep_all = TRUE)
bird3 <- bird2 %>%
  dplyr::group_by(name_latin) %>%
  dplyr::mutate(CountTotal = sum(CountVisit)) 
bird4 <- bird3 %>% filter(CountTotal != 0)
birds5 <- bird4 %>% select(point_ID, name_latin, site_type, square_ID, CountVisit, effArea, point_count_repeat_number)


####GET EFFECTIVE AREAS FOR NON-OBSERVATIONS
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

####import species' conservation status
cons.status <- read.csv("all_spp_info.csv")
cons.status$name_latin <- cons.status$Species
####combine 
birds.d <- plyr::join_all(list(birds14, cons.status), by = c("name_latin"),
                          type  = 'left')
###remove rows where a given species was not recorded
birds.f <- birds.d %>% subset(CountVisit != 0)

birds.f2 <- birds.f %>% filter(Cons.importance.new == "yes") ##keep only species of conservation importance

###get number of species (of conservation importance) per visit to each point
birds.g <- birds.f2 %>%
  dplyr::group_by(point_ID, point_count_repeat_number) %>%
  dplyr::mutate(spp_richness_cons_concern = n_distinct(name_latin)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(point_ID, point_count_repeat_number, .keep_all = TRUE)

##add rows with 0 species of cons concern
all_combinations <- data.sum4 %>%
  dplyr::distinct(point_ID, site_type, square_ID, point_count_repeat_number)

complete_data <- all_combinations %>%
  left_join(birds.g, by = c("point_ID", "point_count_repeat_number"))

# Fill NAs in spp_richness_cons_concern with 0
birds.h <- complete_data %>%
  mutate(spp_richness_cons_concern = ifelse(is.na(spp_richness_cons_concern), 0, spp_richness_cons_concern))

#####get the average effective area surveyed across all species per point
point_avg <- data.sum4 %>%
  group_by(point_ID) %>%
  summarize(avg_effArea = mean(effArea, na.rm = TRUE))

# Calculate average effArea for each square
square_avg <- data.sum4 %>%
  group_by(square_ID) %>%
  summarize(avg_effArea = mean(effArea, na.rm = TRUE))
final_avg <- point_avg %>%
  left_join(data.sum4 %>% select(point_ID, square_ID) %>% distinct(), by = "point_ID") %>%
  left_join(square_avg, by = "square_ID", suffix = c("_point", "_square")) %>%
  mutate(final_effArea = ifelse(is.na(avg_effArea_point), avg_effArea_square, avg_effArea_point)) %>%
  select(point_ID, final_effArea)

birds.i <- birds.h %>%
  left_join(final_avg, by = "point_ID")

#### 2.Fit model####
###set priors
my_priors <- c(
  prior(normal(0, 10), class = "Intercept"),
  prior(normal(0, 10), class = "b")
)
##brms model
my_brm_richness_new2 <- brm(
  bf(spp_richness_cons_concern ~ site_type.x  +
       (1 | square_ID.x / point_ID)+
       (offset(final_effArea)), ## final_effArea is on the log scale
     zi ~ site_type.x
  ),
  data = birds.i,
  family = zero_inflated_poisson(),
  chains = 4,
  cores = 4,
  iter = 3000,
  warmup = 1000,
  prior = my_priors
)

####3. Extract effects
# Summarize the model predictions
newdata <- birds.i %>%
  distinct(site_type.x) %>%
  mutate(final_effArea = 1)

# Add predicted draws
predictions <- add_epred_draws(newdata, my_brm_richness_new2, re_formula = NA)

# Calculate median and credible intervals
summary_predictions <- predictions %>%
  mutate(site_type = factor(site_type.x, levels = c("forest", "ZBNF", "chemical"))) %>%
  group_by(site_type) %>%
  summarise(
    median = median(.epred),
    lower = quantile(.epred, 0.05),
    upper = quantile(.epred, 0.95)
  )

# Define colors
site_colors <- c("forest" = "palegreen3", "ZBNF" =  "khaki", "chemical" = "slateblue4")

# Plot the results
ggplot(summary_predictions, aes(x = site_type, y = median)) +
  geom_point(size = 8, aes(color = site_type), fill = "white") +  # Points with outline and fill color
  geom_errorbar(aes(ymin = lower, ymax = upper, color = site_type), width = 0.4, size = 3) +
  labs(
    title = "Estimated Median Species Richness by Site Type",
    x = "Site Type",
    y = "Estimated Median Species Richness"
  ) +
  scale_color_manual(values = site_colors) +  # Set manual color scale
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Rotate x-axis labels
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    axis.title = element_text(size = 12),  # Adjust axis title size
    axis.line = element_line(size = 0.5, colour = "black"),  # Axis line color and size
    axis.ticks = element_line(size = 0.5, colour = "black"),  # Axis tick color and size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line.y = element_line(colour = "black"),  # Y-axis line color
    axis.line.x = element_line(colour = "black")  # X-axis line color

  )
