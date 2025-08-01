# ZBNF
This R code is associated with the paper:
Berger, I., Kamble, A., Morton, O., Raj, V., Nair, S.R., Edwards, D.P., Wauchope, H.S., Joshi, V., Basu, P., Smith, B., and Dicks, L.V. (2025) "India's agroecology programme,'Zero Budget Natural Farming', delivers biodiversity and economic benefits without lowering yields" Nature Ecology & Evolution
It was written by  Iris Berger. Please do not hesitate to direct any queries to irisberger1996@gmail.com
The data used in the scripts is available via Zenodo: 10.5281/zenodo.16687021

yield_analyses.R - This script cleans and prepares the harvest-level yield data, conducts all matching runs, and conducts g-computation to quantify the effect of the ZBNF programme on harvest-level yield

profit_analyses.R - This script cleans and prepares the field-level profit data, conducts all matching runs, and conducts g-computation to quantify the effect of the ZBNF programme on field-level economic profit

square_yield_profit.R - This script calculates square-level yield and economic profit

detection_functions.R - This script fits detection functions, identifies the best-fitting detection function for each bird species at each point, and calculates the effective area surveyed for each species at each point (which will be included as an offset in subsequent models)

birds_matching.R - This script conducts all matching runs and obtains the matching weights from the best-performing run. The matching weights will be included in the models described in count_mod1.R and bray_curtis_index.R

count_mod1.R - This script prepares the bird data and fits a hierarchical count model (see Model Formula 1 in the main manuscript). It then extracts the estimated effects of the impact of the ZBNF programme on bird densities (at the species and trophic guild levels).

count_mod2.R - This script fits a hierarchical count model (see Model Formula 2 in the main manuscript), extracts density-yield/profit curves and estimates effects of increasing yield/profit on bird abundances (at the species and guild level) for each farming system

species_richness.R - This script estimates the species richness of species of conservation importance in forest, ZBNF, and agrichemical systems by fitting the model described in Formula 3 of the main manuscript

bray_curtis_index.R - This script calculates the Bray-Curtis dissimilarity index between each forest point and each agricultural point and fits the model described in Formula 4 of the main manuscript to estimate the effect of the ZBNF programme on bird community integrity
