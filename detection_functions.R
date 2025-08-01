###### Code to fit detection functions, identify the best detection function for each species 
##at each point, and estimate the effective area surveyed

###Download the data from 10.5281/zenodo.16687021

##load packages
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(mrds)
library(reshape)
library(doMC)
registerDoMC()
library(purrr)
library(dplyr)
library(ggplot2)
library(MuMIn)
library(tidyr)
library(Distance)


# Load the required data and functions
rm(list = ls())
load("DataForDetectionFunctionFittingFinal.rda")
load("DetectionFunction_functions.rda")

# Check that cluster size has no effect on detection
cluster.size.models <- dlply(single.species.data, "Own", function(x) lm(x$size~x$distance))
cluster.size.summaries <- llply(cluster.size.models, summary)
p.vals <- ldply(cluster.size.summaries, function(x) x$coefficients[2,4])
sig.p.vals <- subset(p.vals, V1 <= 0.05)
sig.cs.spp <- subset(single.species.data, name_latin %in% sig.p.vals$.id)
plot.list <- list()
for(i in 1:length(sig.p.vals$.id)) {
  species <- sig.p.vals$.id[i]
  data = subset(single.species.data, name_latin == species)
  plot.list[[paste(species)]] <- plot(data$distance, data$size, main=paste(species))
}

cs.models.multi <- dlply(multi.species.data, "Group", function(x) lm(x$size~x$distance))
cs.summaries.multi <- llply(cs.models.multi, summary)
p.vals.multi <- ldply(cs.summaries.multi, function(x) x$coefficients[2,4])
sig.p.vals <- subset(p.vals.multi, V1 <= 0.05)

sig.cs.spp <- subset(multi.species.data, name_latin %in% sig.p.vals$.id)
plot.list <- list()
for(i in 1:length(sig.p.vals$.id)) {
  group <- sig.p.vals$.id[i]
  data = subset(multi.species.data, Group == group)
  plot.list[[paste(group)]] <- plot(data$distance, data$size, main=paste(group))
}
##### all ok
rm(data, group, sig.cs.spp, sig.p.vals, p.vals, p.vals.multi, plot.list, cs.models.multi, cs.summaries.multi, i, species)


###### SINGLE SPECIES DETECTION FUNCTIONS#########
#####1) One 'global' detection function
# half-normal, no adjustment terms
models.hn.single.no.adj <- dlply(single.species.data, .variables = "Own", 
                                 .fun = possibly(hn.single.no.adj, otherwise = NULL), 
                                 .parallel = TRUE)
failed <- names(which(vapply(models.hn.single.no.adj, is.null, logical(1)))) 
failed.single.sp <- c(paste(failed, "hn.single.no.adj", sep="."))
names(models.hn.single.no.adj) <- c(paste(names(models.hn.single.no.adj),"hn.single.no.adj", sep="."))
models.hn.single.no.adj <- models.hn.single.no.adj[!(names(models.hn.single.no.adj) %in% failed.single.sp)]
##### half-normal, cosine 2
models.hn.single.cos.2 <- dlply(single.species.data, .variables = "Own", 
                                .fun = possibly(hn.single.cos.2, otherwise = NULL), 
                                .parallel = TRUE)
names(models.hn.single.cos.2) <- c(paste(names(models.hn.single.cos.2), "hn.single.cos.2", sep="."))
failed <- names(which(vapply(models.hn.single.cos.2, is.null, logical(1)))) 
failed.single.sp <- c(failed.single.sp, failed)
models.hn.single.cos.2 <- models.hn.single.cos.2[!(names(models.hn.single.cos.2) %in% failed)]
##### half-normal, cosine 3
models.hn.single.cos.3 <- dlply(single.species.data, .variables = "Own", 
                                .fun = possibly(hn.single.cos.3, otherwise = NULL), 
                                .parallel = TRUE)
names(models.hn.single.cos.3) <- c(paste(names(models.hn.single.cos.3), "hn.single.cos.3", sep="."))
failed <- names(which(vapply(models.hn.single.cos.3, is.null, logical(1)))) 
failed.single.sp <- c(failed.single.sp, failed)
models.hn.single.cos.3 <- models.hn.single.cos.3[!(names(models.hn.single.cos.3) %in% failed)]
#####half-normal, hermitic polynomial 4
models.hn.single.herm.4 <- dlply(single.species.data, .variables = "Own", 
                                 .fun = possibly(hn.single.herm.4, otherwise = NULL), 
                                 .parallel = TRUE)
names(models.hn.single.herm.4) <- c(paste(names(models.hn.single.herm.4), "hn.single.herm.4", sep="."))
failed <- names(which(vapply(models.hn.single.herm.4, is.null, logical(1))))
models.hn.single.herm.4 <- models.hn.single.herm.4[!(names(models.hn.single.herm.4) %in% failed)]
failed.single.sp <- c(failed.single.sp, failed)
#### half-normal, hermitic polynomial 6
models.hn.single.herm.6 <- dlply(single.species.data, .variables = "Own", 
                                 .fun = possibly(hn.single.herm.6, otherwise = NULL), 
                                 .parallel = TRUE)
names(models.hn.single.herm.6) <- c(paste(names(models.hn.single.herm.6), "hn.single.herm.6", sep="."))
failed <- names(which(vapply(models.hn.single.herm.6, is.null, logical(1))))
models.hn.single.herm.6 <- models.hn.single.herm.6[!(names(models.hn.single.herm.6) %in% failed)]
failed.single.sp <- c(failed.single.sp, failed)
#####hazard rate, no adjustment terms
models.hr.single.no.adj <- dlply(single.species.data, .variables = "Own", 
                                 .fun = possibly(hr.single.no.adj, otherwise = NULL), 
                                 .parallel = TRUE)
names(models.hr.single.no.adj) <- c(paste(names(models.hr.single.no.adj), "hr.single.no.adj", sep="."))
failed <- names(which(vapply(models.hr.single.no.adj, is.null, logical(1))))
models.hr.single.no.adj <- models.hr.single.no.adj[!(names(models.hr.single.no.adj) %in% failed)]
failed.single.sp <- c(failed.single.sp, failed)
#####hazard-rate, cosine 2
models.hr.single.cos.2 <- dlply(single.species.data, .variables = "Own", 
                                .fun = possibly(hr.single.cos.2, otherwise = NULL), 
                                .parallel = TRUE)
names(models.hr.single.cos.2) <- c(paste(names(models.hr.single.cos.2), "single.hr.cos.2", sep="."))
failed <- names(which(vapply(models.hr.single.cos.2, is.null, logical(1))))
models.hr.single.cos.2 <- models.hr.single.cos.2[!(names(models.hr.single.cos.2) %in% failed)]
failed.single.sp <- c(failed.single.sp, failed)
##### hazard-rate, cosine 3
models.hr.single.cos.3 <- dlply(single.species.data, .variables = "Own", 
                                .fun = possibly(hr.single.cos.3, otherwise = NULL), 
                                .parallel = TRUE)
names(models.hr.single.cos.3) <- c(paste(names(models.hr.single.cos.3),"hr.single.cos.3", sep="."))
failed <- names(which(vapply(models.hr.single.cos.3, is.null, logical(1))))
models.hr.single.cos.3 <- models.hr.single.cos.3[!(names(models.hr.single.cos.3) %in% failed)]
failed.single.sp <- c(failed.single.sp, failed)
#####hazard-rate, polynomial 2
models.hr.single.poly.2 <- dlply(single.species.data, .variables = "Own", 
                                 .fun = possibly(hr.single.poly.2, otherwise = NULL), 
                                 .parallel = TRUE)
names(models.hr.single.poly.2) <- c(paste(names(models.hr.single.poly.2),"hr.single.poly.2", sep="."))
failed <- names(which(vapply(models.hr.single.poly.2, is.null, logical(1))))
models.hr.single.poly.2 <- models.hr.single.poly.2[!(names(models.hr.single.poly.2) %in% failed)]
failed.single.sp <- c(failed.single.sp, failed)
#####hazard-rate, polynomial 4
models.hr.single.poly.4 <- dlply(single.species.data, .variables = "Own", 
                                 .fun = possibly(hr.single.poly.4, otherwise = NULL), 
                                 .parallel = TRUE)
names(models.hr.single.poly.4) <- c(paste(names(models.hr.single.poly.4),"hr.single.poly.4", sep="."))
failed <- names(which(vapply(models.hr.single.poly.4, is.null, logical(1))))
models.hr.single.poly.4 <- models.hr.single.poly.4[!(names(models.hr.single.poly.4) %in% failed)]
failed.single.sp <- c(failed.single.sp, failed)

##########################
## (2.) Single detection function with proportion closed habitat as covariate
# half-normal, no adjustment terms
models.hn.single.closed.no.adj <- dlply(single.species.data, .variables = "Own", 
                                        .fun = possibly(hn.single.closed.no.adj, otherwise = NULL), 
                                        .parallel = TRUE)
failed <- names(which(vapply(models.hn.single.closed.no.adj, is.null, logical(1)))) 
failed.single.closed.sp <- c(paste(failed, "hn.single.closed.no.adj", sep="."))
names(models.hn.single.closed.no.adj) <- c(paste(names(models.hn.single.closed.no.adj),"hn.single.closed.no.adj", sep="."))
models.hn.single.closed.no.adj <- models.hn.single.closed.no.adj[!(names(models.hn.single.closed.no.adj) %in% failed.single.closed.sp)]
##### half-normal, cosine 2
models.hn.single.closed.cos.2 <- dlply(single.species.data, .variables = "Own", 
                                       .fun = possibly(hn.single.closed.cos.2, otherwise = NULL), 
                                       .parallel = TRUE)
names(models.hn.single.closed.cos.2) <- c(paste(names(models.hn.single.closed.cos.2), "hn.single.closed.cos.2", sep="."))
failed <- names(which(vapply(models.hn.single.closed.cos.2, is.null, logical(1)))) 
failed.single.closed.sp <- c(failed.single.closed.sp, failed)
models.hn.single.closed.cos.2 <- models.hn.single.closed.cos.2[!(names(models.hn.single.closed.cos.2) %in% failed)]
##### half-normal, cosine 3
models.hn.single.closed.cos.3 <- dlply(single.species.data, .variables = "Own", 
                                       .fun = possibly(hn.single.closed.cos.3, otherwise = NULL), 
                                       .parallel = TRUE)
names(models.hn.single.closed.cos.3) <- c(paste(names(models.hn.single.closed.cos.3), "hn.single.closed.cos.3", sep="."))
failed <- names(which(vapply(models.hn.single.closed.cos.3, is.null, logical(1)))) 
failed.single.closed.sp <- c(failed.single.closed.sp, failed)
models.hn.single.closed.cos.3 <- models.hn.single.closed.cos.3[!(names(models.hn.single.closed.cos.3) %in% failed)]
#####half-normal, hermitic polynomial 4 
models.hn.single.closed.herm.4 <- dlply(single.species.data, .variables = "Own", 
                                        .fun = possibly(hn.single.closed.herm.4, otherwise = NULL), 
                                        .parallel = TRUE)
names(models.hn.single.closed.herm.4) <- c(paste(names(models.hn.single.closed.herm.4), "hn.single.closed.herm.4", sep="."))
failed <- names(which(vapply(models.hn.single.closed.herm.4, is.null, logical(1)))) 
failed.single.closed.sp <- c(failed.single.closed.sp, failed)
models.hn.single.closed.herm.4 <- models.hn.single.closed.herm.4[!(names(models.hn.single.closed.herm.4) %in% failed)]
#####half-normal, hermitic polynomial 6
models.hn.single.closed.herm.6 <- dlply(single.species.data, .variables = "Own", 
                                        .fun = possibly(hn.single.closed.herm.6, otherwise = NULL), 
                                        .parallel = TRUE)
names(models.hn.single.closed.herm.6) <- c(paste(names(models.hn.single.closed.herm.6), "hn.single.closed.herm.6", sep="."))
failed <- names(which(vapply(models.hn.single.closed.herm.6, is.null, logical(1)))) 
failed.single.closed.sp <- c(failed.single.closed.sp, failed)
models.hn.single.closed.herm.6 <- models.hn.single.closed.herm.6[!(names(models.hn.single.closed.herm.6) %in% failed)]
#####hazard rate, no adjustment terms
models.hr.single.closed.no.adj <- dlply(single.species.data, .variables = "Own", 
                                        .fun = possibly(hr.single.closed.no.adj, otherwise = NULL), 
                                        .parallel = TRUE)
failed <- names(which(vapply(models.hr.single.closed.no.adj, is.null, logical(1)))) 
failed.single.closed.sp <- c(paste(failed, "hr.single.closed.no.adj", sep="."))
names(models.hr.single.closed.no.adj) <- c(paste(names(models.hr.single.closed.no.adj),"hr.single.closed.no.adj", sep="."))
models.hr.single.closed.no.adj <- models.hr.single.closed.no.adj[!(names(models.hr.single.closed.no.adj) %in% failed.single.closed.sp)]
#####hazard-rate, cosine 2
models.hr.single.closed.cos.2 <- dlply(single.species.data, .variables = "Own", 
                                       .fun = possibly(hr.single.closed.cos.2, otherwise = NULL), 
                                       .parallel = TRUE)
failed <- names(which(vapply(models.hr.single.closed.cos.2, is.null, logical(1)))) 
failed.single.closed.sp <- c(paste(failed, "hr.single.closed.cos.2", sep="."))
names(models.hr.single.closed.cos.2) <- c(paste(names(models.hr.single.closed.cos.2),"hr.single.closed.cos.2", sep="."))
models.hr.single.closed.cos.2 <- models.hr.single.closed.cos.2[!(names(models.hr.single.closed.cos.2) %in% failed.single.closed.sp)]
#####hazard-rate, cosine 3
models.hr.single.closed.cos.3 <- dlply(single.species.data, .variables = "Own", 
                                       .fun = possibly(hr.single.closed.cos.3, otherwise = NULL), 
                                       .parallel = TRUE)
failed <- names(which(vapply(models.hr.single.closed.cos.3, is.null, logical(1)))) 
failed.single.closed.sp <- c(paste(failed, "hr.single.closed.cos.3", sep="."))
names(models.hr.single.closed.cos.3) <- c(paste(names(models.hr.single.closed.cos.3),"hr.single.closed.cos.3", sep="."))
models.hr.single.closed.cos.3 <- models.hr.single.closed.cos.3[!(names(models.hr.single.closed.cos.3) %in% failed.single.closed.sp)]
####hazard-rate, polynomial 2
models.hr.single.closed.poly.2 <- dlply(single.species.data, .variables = "Own", 
                                        .fun = possibly(hr.single.closed.poly.2, otherwise = NULL), 
                                        .parallel = TRUE)
failed <- names(which(vapply(models.hr.single.closed.poly.2, is.null, logical(1)))) 
failed.single.closed.sp <- c(paste(failed, "hr.single.closed.poly.2", sep="."))
names(models.hr.single.closed.poly.2) <- c(paste(names(models.hr.single.closed.poly.2),"hr.single.closed.poly.2", sep="."))
models.hr.single.closed.poly.2 <- models.hr.single.closed.poly.2[!(names(models.hr.single.closed.poly.2) %in% failed.single.closed.sp)]
#####hazard-rate, polynomial 4
models.hr.single.closed.poly.4 <- dlply(single.species.data, .variables = "Own", 
                                        .fun = possibly(hr.single.closed.poly.4, otherwise = NULL), 
                                        .parallel = TRUE)
failed <- names(which(vapply(models.hr.single.closed.poly.4, is.null, logical(1)))) 
failed.single.closed.sp <- c(paste(failed, "hr.single.closed.poly.4", sep="."))
names(models.hr.single.closed.poly.4) <- c(paste(names(models.hr.single.closed.poly.4),"hr.single.closed.poly.4", sep="."))
models.hr.single.closed.poly.4 <- models.hr.single.closed.poly.4[!(names(models.hr.single.closed.poly.4) %in% failed.single.closed.sp)]

#################################################
## 3.) FOREST only detection functions #### 
# half-normal, no adjustment terms
models.hn.single.forest.no.adj <- dlply(single.species.data.forest, .variables = "Own", 
                                        .fun = possibly(hn.single.no.adj, otherwise = NULL), 
                                        .parallel = TRUE)
failed <- names(which(vapply(models.hn.single.forest.no.adj, is.null, logical(1)))) 
failed.single.forest.sp <- c(paste(failed, "hn.single.forest.no.adj", sep="."))
names(models.hn.single.forest.no.adj) <- c(paste(names(models.hn.single.forest.no.adj),"hn.single.forest.no.adj", sep="."))
models.hn.single.forest.no.adj <- models.hn.single.forest.no.adj[!(names(models.hn.single.forest.no.adj) %in% failed.single.forest.sp)]
#####half-normal, cosine 2
models.hn.single.forest.cos.2 <- dlply(single.species.data.forest, .variables = "Own", 
                                       .fun = possibly(hn.single.cos.2, otherwise = NULL), 
                                       .parallel = TRUE)
names(models.hn.single.forest.cos.2) <- c(paste(names(models.hn.single.forest.cos.2), "hn.single.forest.cos.2", sep="."))
failed <- names(which(vapply(models.hn.single.forest.cos.2, is.null, logical(1)))) 
failed.single.forest.sp <- c(failed.single.forest.sp, failed)
models.hn.single.forest.cos.2 <- models.hn.single.forest.cos.2[!(names(models.hn.single.forest.cos.2) %in% failed)]
#####half-normal, cosine 3
models.hn.single.forest.cos.3 <- dlply(single.species.data.forest, .variables = "Own", 
                                       .fun = possibly(hn.single.cos.3, otherwise = NULL), 
                                       .parallel = TRUE)
names(models.hn.single.forest.cos.3) <- c(paste(names(models.hn.single.forest.cos.3), "hn.single.forest.cos.3", sep="."))
failed <- names(which(vapply(models.hn.single.forest.cos.3, is.null, logical(1)))) 
failed.single.forest.sp <- c(failed.single.forest.sp, failed)
models.hn.single.forest.cos.3 <- models.hn.single.forest.cos.3[!(names(models.hn.single.forest.cos.3) %in% failed)]
#####half-normal, hermitic polynomial 4 
models.hn.single.forest.herm.4 <- dlply(single.species.data.forest, .variables = "Own", 
                                        .fun = possibly(hn.single.herm.4, otherwise = NULL), 
                                        .parallel = TRUE)
names(models.hn.single.forest.herm.4) <- c(paste(names(models.hn.single.forest.herm.4), "hn.single.forest.herm.4", sep="."))
failed <- names(which(vapply(models.hn.single.forest.herm.4, is.null, logical(1)))) 
failed.single.forest.sp <- c(failed.single.forest.sp, failed)
models.hn.single.forest.herm.4 <- models.hn.single.forest.herm.4[!(names(models.hn.single.forest.herm.4) %in% failed)]
#####half-normal, hermitic polynomial 6
models.hn.single.forest.herm.6 <- dlply(single.species.data.forest, .variables = "Own", 
                                        .fun = possibly(hn.single.herm.6, otherwise = NULL), 
                                        .parallel = TRUE)
names(models.hn.single.forest.herm.6) <- c(paste(names(models.hn.single.forest.herm.6), "hn.single.forest.herm.6", sep="."))
failed <- names(which(vapply(models.hn.single.forest.herm.6, is.null, logical(1)))) 
failed.single.forest.sp <- c(failed.single.forest.sp, failed)
models.hn.single.forest.herm.6 <- models.hn.single.forest.herm.6[!(names(models.hn.single.forest.herm.6) %in% failed)]
##### hazard rate, no adjustment terms
models.hr.single.forest.no.adj <- dlply(single.species.data.forest, .variables = "Own", 
                                        .fun = possibly(hr.single.no.adj, otherwise = NULL), 
                                        .parallel = TRUE)
failed <- names(which(vapply(models.hr.single.forest.no.adj, is.null, logical(1)))) 
failed.single.forest.sp <- c(paste(failed, "hr.single.forest.no.adj", sep="."))
names(models.hr.single.forest.no.adj) <- c(paste(names(models.hr.single.forest.no.adj),"hr.single.forest.no.adj", sep="."))
models.hr.single.forest.no.adj <- models.hr.single.forest.no.adj[!(names(models.hr.single.forest.no.adj) %in% failed.single.forest.sp)]
#####hazard-rate, cosine 2
models.hr.single.forest.cos.2 <- dlply(single.species.data.forest, .variables = "Own", 
                                       .fun = possibly(hr.single.cos.2, otherwise = NULL), 
                                       .parallel = TRUE)
failed <- names(which(vapply(models.hr.single.forest.cos.2, is.null, logical(1)))) 
failed.single.forest.sp <- c(paste(failed, "hr.single.forest.cos.2", sep="."))
names(models.hr.single.forest.cos.2) <- c(paste(names(models.hr.single.forest.cos.2),"hr.single.forest.cos.2", sep="."))
models.hr.single.forest.cos.2 <- models.hr.single.forest.cos.2[!(names(models.hr.single.forest.cos.2) %in% failed.single.forest.sp)]
##### hazard-rate, cosine 3
models.hr.single.forest.cos.3 <- dlply(single.species.data.forest, .variables = "Own", 
                                       .fun = possibly(hr.single.cos.3, otherwise = NULL), 
                                       .parallel = TRUE)
failed <- names(which(vapply(models.hr.single.forest.cos.3, is.null, logical(1)))) 
failed.single.forest.sp <- c(paste(failed, "hr.single.forest.cos.3", sep="."))
names(models.hr.single.forest.cos.3) <- c(paste(names(models.hr.single.forest.cos.3),"hr.single.forest.cos.3", sep="."))
models.hr.single.forest.cos.3 <- models.hr.single.forest.cos.3[!(names(models.hr.single.forest.cos.3) %in% failed.single.forest.sp)]
##hazard-rate, polynomial 2
models.hr.single.forest.poly.2 <- dlply(single.species.data.forest, .variables = "Own", 
                                        .fun = possibly(hr.single.poly.2, otherwise = NULL), 
                                        .parallel = TRUE)
failed <- names(which(vapply(models.hr.single.forest.poly.2, is.null, logical(1)))) 
failed.single.forest.sp <- c(paste(failed, "hr.single.forest.poly.2", sep="."))
names(models.hr.single.forest.poly.2) <- c(paste(names(models.hr.single.forest.poly.2),"hr.single.forest.poly.2", sep="."))
models.hr.single.forest.poly.2 <- models.hr.single.forest.poly.2[!(names(models.hr.single.forest.poly.2) %in% failed.single.forest.sp)]
## hazard-rate, polynomial 4
models.hr.single.forest.poly.4 <- dlply(single.species.data.forest, .variables = "Own", 
                                        .fun = possibly(hr.single.poly.4, otherwise = NULL), 
                                        .parallel = TRUE)
failed <- names(which(vapply(models.hr.single.forest.poly.4, is.null, logical(1)))) 
failed.single.forest.sp <- c(paste(failed, "hr.single.forest.poly.4", sep="."))
names(models.hr.single.forest.poly.4) <- c(paste(names(models.hr.single.forest.poly.4),"hr.single.forest.poly.4", sep="."))
models.hr.single.forest.poly.4 <- models.hr.single.forest.poly.4[!(names(models.hr.single.forest.poly.4) %in% failed.single.forest.sp)]

#################################################
## 4.) FARMLAnd only detection functions #### 
# half-normal, no adjustment terms
models.hn.single.farm.no.adj <- dlply(single.species.data.farm, .variables = "Own", 
                                      .fun = possibly(hn.single.no.adj, otherwise = NULL), 
                                      .parallel = TRUE)
failed <- names(which(vapply(models.hn.single.farm.no.adj, is.null, logical(1)))) 
failed.single.farm.sp <- c(paste(failed, "hn.single.farm.no.adj", sep="."))
names(models.hn.single.farm.no.adj) <- c(paste(names(models.hn.single.farm.no.adj),"hn.single.farm.no.adj", sep="."))
models.hn.single.farm.no.adj <- models.hn.single.farm.no.adj[!(names(models.hn.single.farm.no.adj) %in% failed.single.farm.sp)]
#####half-normal, cosine 2
models.hn.single.farm.cos.2 <- dlply(single.species.data.farm, .variables = "Own", 
                                     .fun = possibly(hn.single.cos.2, otherwise = NULL), 
                                     .parallel = TRUE)
names(models.hn.single.farm.cos.2) <- c(paste(names(models.hn.single.farm.cos.2), "hn.single.farm.cos.2", sep="."))
failed <- names(which(vapply(models.hn.single.farm.cos.2, is.null, logical(1)))) 
failed.single.farm.sp <- c(failed.single.farm.sp, failed)
models.hn.single.farm.cos.2 <- models.hn.single.farm.cos.2[!(names(models.hn.single.farm.cos.2) %in% failed)]
#####half-normal, cosine 3
models.hn.single.farm.cos.3 <- dlply(single.species.data.farm, .variables = "Own", 
                                     .fun = possibly(hn.single.cos.3, otherwise = NULL), 
                                     .parallel = TRUE)
names(models.hn.single.farm.cos.3) <- c(paste(names(models.hn.single.farm.cos.3), "hn.single.farm.cos.3", sep="."))
failed <- names(which(vapply(models.hn.single.farm.cos.3, is.null, logical(1)))) 
failed.single.farm.sp <- c(failed.single.farm.sp, failed)
models.hn.single.farm.cos.3 <- models.hn.single.farm.cos.3[!(names(models.hn.single.farm.cos.3) %in% failed)]
# half-normal, hermitic polynomial 4 
models.hn.single.farm.herm.4 <- dlply(single.species.data.farm, .variables = "Own", 
                                      .fun = possibly(hn.single.herm.4, otherwise = NULL), 
                                      .parallel = TRUE)
names(models.hn.single.farm.herm.4) <- c(paste(names(models.hn.single.farm.herm.4), "hn.single.farm.herm.4", sep="."))
failed <- names(which(vapply(models.hn.single.farm.herm.4, is.null, logical(1)))) 
failed.single.farm.sp <- c(failed.single.farm.sp, failed)
models.hn.single.farm.herm.4 <- models.hn.single.farm.herm.4[!(names(models.hn.single.farm.herm.4) %in% failed)]
####half-normal, hermitic polynomial 6
models.hn.single.farm.herm.6 <- dlply(single.species.data.farm, .variables = "Own", 
                                      .fun = possibly(hn.single.herm.6, otherwise = NULL), 
                                      .parallel = TRUE)
names(models.hn.single.farm.herm.6) <- c(paste(names(models.hn.single.farm.herm.6), "hn.single.farm.herm.6", sep="."))
failed <- names(which(vapply(models.hn.single.farm.herm.6, is.null, logical(1)))) 
failed.single.farm.sp <- c(failed.single.farm.sp, failed)
models.hn.single.farm.herm.6 <- models.hn.single.farm.herm.6[!(names(models.hn.single.farm.herm.6) %in% failed)]
#####hazard rate, no adjustment terms
models.hr.single.farm.no.adj <- dlply(single.species.data.farm, .variables = "Own", 
                                      .fun = possibly(hr.single.no.adj, otherwise = NULL), 
                                      .parallel = TRUE)
failed <- names(which(vapply(models.hr.single.farm.no.adj, is.null, logical(1)))) 
failed.single.farm.sp <- c(paste(failed, "hr.single.farm.no.adj", sep="."))
names(models.hr.single.farm.no.adj) <- c(paste(names(models.hr.single.farm.no.adj),"hr.single.farm.no.adj", sep="."))
models.hr.single.farm.no.adj <- models.hr.single.farm.no.adj[!(names(models.hr.single.farm.no.adj) %in% failed.single.farm.sp)]
###hazard-rate, cosine 2
models.hr.single.farm.cos.2 <- dlply(single.species.data.farm, .variables = "Own", 
                                     .fun = possibly(hr.single.cos.2, otherwise = NULL), 
                                     .parallel = TRUE)
failed <- names(which(vapply(models.hr.single.farm.cos.2, is.null, logical(1)))) 
failed.single.farm.sp <- c(paste(failed, "hr.single.farm.cos.2", sep="."))
names(models.hr.single.farm.cos.2) <- c(paste(names(models.hr.single.farm.cos.2),"hr.single.farm.cos.2", sep="."))
models.hr.single.farm.cos.2 <- models.hr.single.farm.cos.2[!(names(models.hr.single.farm.cos.2) %in% failed.single.farm.sp)]
#####hazard-rate, cosine 3
models.hr.single.farm.cos.3 <- dlply(single.species.data.farm, .variables = "Own", 
                                     .fun = possibly(hr.single.cos.3, otherwise = NULL), 
                                     .parallel = TRUE)
failed <- names(which(vapply(models.hr.single.farm.cos.3, is.null, logical(1)))) 
failed.single.farm.sp <- c(paste(failed, "hr.single.farm.cos.3", sep="."))
names(models.hr.single.farm.cos.3) <- c(paste(names(models.hr.single.farm.cos.3),"hr.single.farm.cos.3", sep="."))
models.hr.single.farm.cos.3 <- models.hr.single.farm.cos.3[!(names(models.hr.single.farm.cos.3) %in% failed.single.farm.sp)]
#####hazard-rate, polynomial 2
models.hr.single.farm.poly.2 <- dlply(single.species.data.farm, .variables = "Own", 
                                      .fun = possibly(hr.single.poly.2, otherwise = NULL), 
                                      .parallel = TRUE)
failed <- names(which(vapply(models.hr.single.farm.poly.2, is.null, logical(1)))) 
failed.single.farm.sp <- c(paste(failed, "hr.single.farm.poly.2", sep="."))
names(models.hr.single.farm.poly.2) <- c(paste(names(models.hr.single.farm.poly.2),"hr.single.farm.poly.2", sep="."))
models.hr.single.farm.poly.2 <- models.hr.single.farm.poly.2[!(names(models.hr.single.farm.poly.2) %in% failed.single.farm.sp)]
#####hazard-rate, polynomial 4
models.hr.single.farm.poly.4 <- dlply(single.species.data.farm, .variables = "Own", 
                                      .fun = possibly(hr.single.poly.4, otherwise = NULL), 
                                      .parallel = TRUE)
failed <- names(which(vapply(models.hr.single.farm.poly.4, is.null, logical(1)))) 
failed.single.farm.sp <- c(paste(failed, "hr.single.farm.poly.4", sep="."))
names(models.hr.single.farm.poly.4) <- c(paste(names(models.hr.single.farm.poly.4),"hr.single.farm.poly.4", sep="."))
models.hr.single.farm.poly.4 <- models.hr.single.farm.poly.4[!(names(models.hr.single.farm.poly.4) %in% failed.single.farm.sp)]

################################################## save 
# Put all the single spp models into a single list, save, clean up
#####
single.sp.models <- c(models.hn.single.cos.2, models.hn.single.cos.3, models.hn.single.herm.4, models.hn.single.herm.6, 
                      models.hn.single.no.adj, models.hr.single.cos.2, models.hr.single.cos.3, models.hr.single.no.adj,
                      models.hr.single.poly.2, models.hr.single.poly.4,
                      models.hn.single.closed.cos.2, models.hn.single.closed.cos.3, models.hn.single.closed.herm.4, models.hn.single.closed.herm.6, 
                      models.hn.single.closed.no.adj, models.hr.single.closed.cos.2, models.hr.single.closed.cos.3, models.hr.single.closed.no.adj,
                      models.hr.single.closed.poly.2, models.hr.single.closed.poly.4,
                      models.hn.single.forest.cos.2, models.hn.single.forest.cos.3, models.hn.single.forest.herm.4, models.hn.single.forest.herm.6, 
                      models.hn.single.forest.no.adj, models.hr.single.forest.cos.2, models.hr.single.forest.cos.3, models.hr.single.forest.no.adj,
                      models.hr.single.forest.poly.2, models.hr.single.forest.poly.4,
                      models.hn.single.farm.cos.2, models.hn.single.farm.cos.3, models.hn.single.farm.herm.4, models.hn.single.farm.herm.6, 
                      models.hn.single.farm.no.adj, models.hr.single.farm.cos.2, models.hr.single.farm.cos.3, models.hr.single.farm.no.adj,
                      models.hr.single.farm.poly.2, models.hr.single.farm.poly.4)

single.sp.models.habitats <- c(
  models.hn.single.forest.cos.2, models.hn.single.forest.cos.3, models.hn.single.forest.herm.4, models.hn.single.forest.herm.6, 
  models.hn.single.forest.no.adj, models.hr.single.forest.cos.2, models.hr.single.forest.cos.3, models.hr.single.forest.no.adj,
  models.hr.single.forest.poly.2, models.hr.single.forest.poly.4,
  models.hn.single.farm.cos.2, models.hn.single.farm.cos.3, models.hn.single.farm.herm.4, models.hn.single.farm.herm.6, 
  models.hn.single.farm.no.adj, models.hr.single.farm.cos.2, models.hr.single.farm.cos.3, models.hr.single.farm.no.adj,
  models.hr.single.farm.poly.2, models.hr.single.farm.poly.4)

save(single.sp.models.habitats, file = "SingleSpModsHabitats.rda")

########################MULTI SPECIES MODELS (DETECTION GROUPS)##########################
#########1) 1 'global' detection function
# half-normal, no adjustment terms
models.hn.multi.no.adj <- dlply(multi.species.data, .variables = "Group", 
                                .fun = possibly(hn.multi.no.adj, otherwise = NULL), 
                                .parallel = TRUE)
failed <- names(which(vapply(models.hn.multi.no.adj, is.null, logical(1)))) ##all ok
failed.multi.sp <- c(paste(failed, "hn.multi.no.adj", sep="."))
names(models.hn.multi.no.adj) <- c(paste(names(models.hn.multi.no.adj),"hn.multi.no.adj", sep="."))
models.hn.multi.no.adj <- models.hn.multi.no.adj[!(names(models.hn.multi.no.adj) %in% failed.multi.sp)]
#####half-normal, cosine 2
models.hn.multi.cos.2 <- dlply(multi.species.data, .variables = "Group", 
                               .fun = possibly(hn.multi.cos.2, otherwise = NULL), 
                               .parallel = TRUE)
names(models.hn.multi.cos.2) <- c(paste(names(models.hn.multi.cos.2), "hn.multi.cos.2", sep="."))
failed <- names(which(vapply(models.hn.multi.cos.2, is.null, logical(1)))) 
failed.multi.sp <- c(failed.multi.sp, failed)
models.hn.multi.cos.2 <- models.hn.multi.cos.2[!(names(models.hn.multi.cos.2) %in% failed)]
##### half-normal, cosine 3
models.hn.multi.cos.3 <- dlply(multi.species.data, .variables = "Group", 
                               .fun = possibly(hn.multi.cos.3, otherwise = NULL), 
                               .parallel = TRUE)
names(models.hn.multi.cos.3) <- c(paste(names(models.hn.multi.cos.3), "hn.multi.cos.3", sep="."))
failed <- names(which(vapply(models.hn.multi.cos.3, is.null, logical(1)))) 
failed.multi.sp <- c(failed.multi.sp, failed)
models.hn.multi.cos.3 <- models.hn.multi.cos.3[!(names(models.hn.multi.cos.3) %in% failed)]
#####half-normal, hermitic polynomial 4 - THESE DON'T WORK
models.hn.multi.herm.4 <- dlply(multi.species.data, .variables = "Group", 
                                .fun = possibly(hn.multi.herm.4, otherwise = NULL), 
                                .parallel = TRUE)
names(models.hn.multi.herm.4) <- c(paste(names(models.hn.multi.herm.4), "hn.multi.herm.4", sep="."))
failed <- names(which(vapply(models.hn.multi.herm.4, is.null, logical(1))))
models.hn.multi.herm.4 <- models.hn.multi.herm.4[!(names(models.hn.multi.herm.4) %in% failed)]
failed.multi.sp <- c(failed.multi.sp, failed)
####half-normal, hermitic polynomial 6
models.hn.multi.herm.6 <- dlply(multi.species.data, .variables = "Group", 
                                .fun = possibly(hn.multi.herm.6, otherwise = NULL), 
                                .parallel = TRUE)
names(models.hn.multi.herm.6) <- c(paste(names(models.hn.multi.herm.6), "hn.multi.herm.6", sep="."))
failed <- names(which(vapply(models.hn.multi.herm.6, is.null, logical(1))))
models.hn.multi.herm.6 <- models.hn.multi.herm.6[!(names(models.hn.multi.herm.6) %in% failed)]
failed.multi.sp <- c(failed.multi.sp, failed)
###hazard rate, no adjustment terms
models.hr.multi.no.adj <- dlply(multi.species.data, .variables = "Group", 
                                .fun = possibly(hr.multi.no.adj, otherwise = NULL), 
                                .parallel = TRUE)
names(models.hr.multi.no.adj) <- c(paste(names(models.hr.multi.no.adj), "hr.multi.no.adj", sep="."))
failed <- names(which(vapply(models.hr.multi.no.adj, is.null, logical(1))))
models.hr.multi.no.adj <- models.hr.multi.no.adj[!(names(models.hr.multi.no.adj) %in% failed)]
failed.multi.sp <- c(failed.multi.sp, failed)
####hazard-rate, cosine 2
models.hr.multi.cos.2 <- dlply(multi.species.data, .variables = "Group", 
                               .fun = possibly(hr.multi.cos.2, otherwise = NULL), 
                               .parallel = TRUE)
names(models.hr.multi.cos.2) <- c(paste(names(models.hr.multi.cos.2), "multi.hr.cos.2", sep="."))
failed <- names(which(vapply(models.hr.multi.cos.2, is.null, logical(1))))
models.hr.multi.cos.2 <- models.hr.multi.cos.2[!(names(models.hr.multi.cos.2) %in% failed)]
failed.multi.sp <- c(failed.multi.sp, failed)
####hazard-rate, cosine 3
models.hr.multi.cos.3 <- dlply(multi.species.data, .variables = "Group", 
                               .fun = possibly(hr.multi.cos.3, otherwise = NULL), 
                               .parallel = TRUE)
names(models.hr.multi.cos.3) <- c(paste(names(models.hr.multi.cos.3),"hr.multi.cos.3", sep="."))
failed <- names(which(vapply(models.hr.multi.cos.3, is.null, logical(1))))
models.hr.multi.cos.3 <- models.hr.multi.cos.3[!(names(models.hr.multi.cos.3) %in% failed)]
failed.multi.sp <- c(failed.multi.sp, failed)
#####hazard-rate, polynomial 2
models.hr.multi.poly.2 <- dlply(multi.species.data, .variables = "Group", 
                                .fun = possibly(hr.multi.poly.2, otherwise = NULL), 
                                .parallel = TRUE)
names(models.hr.multi.poly.2) <- c(paste(names(models.hr.multi.poly.2),"hr.multi.poly.2", sep="."))
failed <- names(which(vapply(models.hr.multi.poly.2, is.null, logical(1))))
models.hr.multi.poly.2 <- models.hr.multi.poly.2[!(names(models.hr.multi.poly.2) %in% failed)]
failed.multi.sp <- c(failed.multi.sp, failed)
####hazard-rate, polynomial 4
models.hr.multi.poly.4 <- dlply(multi.species.data, .variables = "Group", 
                                .fun = possibly(hr.multi.poly.4, otherwise = NULL), 
                                .parallel = TRUE)
names(models.hr.multi.poly.4) <- c(paste(names(models.hr.multi.poly.4),"hr.multi.poly.4", sep="."))
failed <- names(which(vapply(models.hr.multi.poly.4, is.null, logical(1))))
models.hr.multi.poly.4 <- models.hr.multi.poly.4[!(names(models.hr.multi.poly.4) %in% failed)]
failed.multi.sp <- c(failed.multi.sp, failed)
########################################
## (2.) multi detection function with proportion closed habitat as covariate
# half-normal, no adjustment terms
models.hn.multi.closed.no.adj <- dlply(multi.species.data, .variables = "Group", 
                                       .fun = possibly(hn.multi.closed.no.adj, otherwise = NULL), 
                                       .parallel = TRUE)
failed <- names(which(vapply(models.hn.multi.closed.no.adj, is.null, logical(1)))) 
failed.multi.closed.sp <- c(paste(failed, "hn.multi.closed.no.adj", sep="."))
names(models.hn.multi.closed.no.adj) <- c(paste(names(models.hn.multi.closed.no.adj),"hn.multi.closed.no.adj", sep="."))
models.hn.multi.closed.no.adj <- models.hn.multi.closed.no.adj[!(names(models.hn.multi.closed.no.adj) %in% failed.multi.closed.sp)]
#####half-normal, cosine 2
models.hn.multi.closed.cos.2 <- dlply(multi.species.data, .variables = "Group", 
                                      .fun = possibly(hn.multi.closed.cos.2, otherwise = NULL), 
                                      .parallel = TRUE)
names(models.hn.multi.closed.cos.2) <- c(paste(names(models.hn.multi.closed.cos.2), "hn.multi.closed.cos.2", sep="."))
failed <- names(which(vapply(models.hn.multi.closed.cos.2, is.null, logical(1)))) 
failed.multi.closed.sp <- c(failed.multi.closed.sp, failed)
models.hn.multi.closed.cos.2 <- models.hn.multi.closed.cos.2[!(names(models.hn.multi.closed.cos.2) %in% failed)]
#####half-normal, cosine 3
models.hn.multi.closed.cos.3 <- dlply(multi.species.data, .variables = "Group", 
                                      .fun = possibly(hn.multi.closed.cos.3, otherwise = NULL), 
                                      .parallel = TRUE)
names(models.hn.multi.closed.cos.3) <- c(paste(names(models.hn.multi.closed.cos.3), "hn.multi.closed.cos.3", sep="."))
failed <- names(which(vapply(models.hn.multi.closed.cos.3, is.null, logical(1)))) 
failed.multi.closed.sp <- c(failed.multi.closed.sp, failed)
models.hn.multi.closed.cos.3 <- models.hn.multi.closed.cos.3[!(names(models.hn.multi.closed.cos.3) %in% failed)]
# half-normal, hermitic polynomial 4 
models.hn.multi.closed.herm.4 <- dlply(multi.species.data, .variables = "Group", 
                                       .fun = possibly(hn.multi.closed.herm.4, otherwise = NULL), 
                                       .parallel = TRUE)
names(models.hn.multi.closed.herm.4) <- c(paste(names(models.hn.multi.closed.herm.4), "hn.multi.closed.herm.4", sep="."))
failed <- names(which(vapply(models.hn.multi.closed.herm.4, is.null, logical(1)))) 
failed.multi.closed.sp <- c(failed.multi.closed.sp, failed)
models.hn.multi.closed.herm.4 <- models.hn.multi.closed.herm.4[!(names(models.hn.multi.closed.herm.4) %in% failed)]
#####half-normal, hermitic polynomial 6
models.hn.multi.closed.herm.6 <- dlply(multi.species.data, .variables = "Group", 
                                       .fun = possibly(hn.multi.closed.herm.6, otherwise = NULL), 
                                       .parallel = TRUE)
names(models.hn.multi.closed.herm.6) <- c(paste(names(models.hn.multi.closed.herm.6), "hn.multi.closed.herm.6", sep="."))
failed <- names(which(vapply(models.hn.multi.closed.herm.6, is.null, logical(1)))) 
failed.multi.closed.sp <- c(failed.multi.closed.sp, failed)
models.hn.multi.closed.herm.6 <- models.hn.multi.closed.herm.6[!(names(models.hn.multi.closed.herm.6) %in% failed)]
#####hazard rate, no adjustment terms
models.hr.multi.closed.no.adj <- dlply(multi.species.data, .variables = "Group", 
                                       .fun = possibly(hr.multi.closed.no.adj, otherwise = NULL), 
                                       .parallel = TRUE)
failed <- names(which(vapply(models.hr.multi.closed.no.adj, is.null, logical(1)))) 
failed.multi.closed.sp <- c(paste(failed, "hr.multi.closed.no.adj", sep="."))
names(models.hr.multi.closed.no.adj) <- c(paste(names(models.hr.multi.closed.no.adj),"hr.multi.closed.no.adj", sep="."))
models.hr.multi.closed.no.adj <- models.hr.multi.closed.no.adj[!(names(models.hr.multi.closed.no.adj) %in% failed.multi.closed.sp)]
#####hazard-rate, cosine 2
models.hr.multi.closed.cos.2 <- dlply(multi.species.data, .variables = "Group", 
                                      .fun = possibly(hr.multi.closed.cos.2, otherwise = NULL), 
                                      .parallel = TRUE)
failed <- names(which(vapply(models.hr.multi.closed.cos.2, is.null, logical(1)))) 
failed.multi.closed.sp <- c(paste(failed, "hr.multi.closed.cos.2", sep="."))
names(models.hr.multi.closed.cos.2) <- c(paste(names(models.hr.multi.closed.cos.2),"hr.multi.closed.cos.2", sep="."))
models.hr.multi.closed.cos.2 <- models.hr.multi.closed.cos.2[!(names(models.hr.multi.closed.cos.2) %in% failed.multi.closed.sp)]
#####hazard-rate, cosine 3
models.hr.multi.closed.cos.3 <- dlply(multi.species.data, .variables = "Group", 
                                      .fun = possibly(hr.multi.closed.cos.3, otherwise = NULL), 
                                      .parallel = TRUE)
failed <- names(which(vapply(models.hr.multi.closed.cos.3, is.null, logical(1)))) 
failed.multi.closed.sp <- c(paste(failed, "hr.multi.closed.cos.3", sep="."))
names(models.hr.multi.closed.cos.3) <- c(paste(names(models.hr.multi.closed.cos.3),"hr.multi.closed.cos.3", sep="."))
models.hr.multi.closed.cos.3 <- models.hr.multi.closed.cos.3[!(names(models.hr.multi.closed.cos.3) %in% failed.multi.closed.sp)]
####hazard-rate, polynomial 2
models.hr.multi.closed.poly.2 <- dlply(multi.species.data, .variables = "Group", 
                                       .fun = possibly(hr.multi.closed.poly.2, otherwise = NULL), 
                                       .parallel = TRUE)
failed <- names(which(vapply(models.hr.multi.closed.poly.2, is.null, logical(1)))) 
failed.multi.closed.sp <- c(paste(failed, "hr.multi.closed.poly.2", sep="."))
names(models.hr.multi.closed.poly.2) <- c(paste(names(models.hr.multi.closed.poly.2),"hr.multi.closed.poly.2", sep="."))
models.hr.multi.closed.poly.2 <- models.hr.multi.closed.poly.2[!(names(models.hr.multi.closed.poly.2) %in% failed.multi.closed.sp)]
#####hazard-rate, polynomial 4
models.hr.multi.closed.poly.4 <- dlply(multi.species.data, .variables = "Group", 
                                       .fun = possibly(hr.multi.closed.poly.4, otherwise = NULL), 
                                       .parallel = TRUE)
failed <- names(which(vapply(models.hr.multi.closed.poly.4, is.null, logical(1)))) 
failed.multi.closed.sp <- c(paste(failed, "hr.multi.closed.poly.4", sep="."))
names(models.hr.multi.closed.poly.4) <- c(paste(names(models.hr.multi.closed.poly.4),"hr.multi.closed.poly.4", sep="."))
models.hr.multi.closed.poly.4 <- models.hr.multi.closed.poly.4[!(names(models.hr.multi.closed.poly.4) %in% failed.multi.closed.sp)]

#################################################
## 3.) FOREST only detection functions #### 
# half-normal, no adjustment terms
models.hn.multi.forest.no.adj <- dlply(multi.species.data.forest, .variables = "Group", 
                                       .fun = possibly(hn.multi.no.adj, otherwise = NULL), 
                                       .parallel = TRUE)
failed <- names(which(vapply(models.hn.multi.forest.no.adj, is.null, logical(1)))) 
failed.multi.forest.sp <- c(paste(failed, "hn.multi.forest.no.adj", sep="."))
names(models.hn.multi.forest.no.adj) <- c(paste(names(models.hn.multi.forest.no.adj),"hn.multi.forest.no.adj", sep="."))
models.hn.multi.forest.no.adj <- models.hn.multi.forest.no.adj[!(names(models.hn.multi.forest.no.adj) %in% failed.multi.forest.sp)]
#####half-normal, cosine 2
models.hn.multi.forest.cos.2 <- dlply(multi.species.data.forest, .variables = "Group", 
                                      .fun = possibly(hn.multi.cos.2, otherwise = NULL), 
                                      .parallel = TRUE)
names(models.hn.multi.forest.cos.2) <- c(paste(names(models.hn.multi.forest.cos.2), "hn.multi.forest.cos.2", sep="."))
failed <- names(which(vapply(models.hn.multi.forest.cos.2, is.null, logical(1)))) 
failed.multi.forest.sp <- c(failed.multi.forest.sp, failed)
models.hn.multi.forest.cos.2 <- models.hn.multi.forest.cos.2[!(names(models.hn.multi.forest.cos.2) %in% failed)]
#####half-normal, cosine 3
models.hn.multi.forest.cos.3 <- dlply(multi.species.data.forest, .variables = "Group", 
                                      .fun = possibly(hn.multi.cos.3, otherwise = NULL), 
                                      .parallel = TRUE)
names(models.hn.multi.forest.cos.3) <- c(paste(names(models.hn.multi.forest.cos.3), "hn.multi.forest.cos.3", sep="."))
failed <- names(which(vapply(models.hn.multi.forest.cos.3, is.null, logical(1)))) 
failed.multi.forest.sp <- c(failed.multi.forest.sp, failed)
models.hn.multi.forest.cos.3 <- models.hn.multi.forest.cos.3[!(names(models.hn.multi.forest.cos.3) %in% failed)]
#####half-normal, hermitic polynomial 4 
models.hn.multi.forest.herm.4 <- dlply(multi.species.data.forest, .variables = "Group", 
                                       .fun = possibly(hn.multi.herm.4, otherwise = NULL), 
                                       .parallel = TRUE)
names(models.hn.multi.forest.herm.4) <- c(paste(names(models.hn.multi.forest.herm.4), "hn.multi.forest.herm.4", sep="."))
failed <- names(which(vapply(models.hn.multi.forest.herm.4, is.null, logical(1)))) 
failed.multi.forest.sp <- c(failed.multi.forest.sp, failed)
models.hn.multi.forest.herm.4 <- models.hn.multi.forest.herm.4[!(names(models.hn.multi.forest.herm.4) %in% failed)]
#####half-normal, hermitic polynomial 6
models.hn.multi.forest.herm.6 <- dlply(multi.species.data.forest, .variables = "Group", 
                                       .fun = possibly(hn.multi.herm.6, otherwise = NULL), 
                                       .parallel = TRUE)
names(models.hn.multi.forest.herm.6) <- c(paste(names(models.hn.multi.forest.herm.6), "hn.multi.forest.herm.6", sep="."))
failed <- names(which(vapply(models.hn.multi.forest.herm.6, is.null, logical(1)))) 
failed.multi.forest.sp <- c(failed.multi.forest.sp, failed)
models.hn.multi.forest.herm.6 <- models.hn.multi.forest.herm.6[!(names(models.hn.multi.forest.herm.6) %in% failed)]
#####hazard rate, no adjustment terms
models.hr.multi.forest.no.adj <- dlply(multi.species.data.forest, .variables = "Group", 
                                       .fun = possibly(hr.multi.no.adj, otherwise = NULL), 
                                       .parallel = TRUE)
failed <- names(which(vapply(models.hr.multi.forest.no.adj, is.null, logical(1)))) 
failed.multi.forest.sp <- c(paste(failed, "hr.multi.forest.no.adj", sep="."))
names(models.hr.multi.forest.no.adj) <- c(paste(names(models.hr.multi.forest.no.adj),"hr.multi.forest.no.adj", sep="."))
models.hr.multi.forest.no.adj <- models.hr.multi.forest.no.adj[!(names(models.hr.multi.forest.no.adj) %in% failed.multi.forest.sp)]
#####hazard-rate, cosine 2
models.hr.multi.forest.cos.2 <- dlply(multi.species.data.forest, .variables = "Group", 
                                      .fun = possibly(hr.multi.cos.2, otherwise = NULL), 
                                      .parallel = TRUE)
failed <- names(which(vapply(models.hr.multi.forest.cos.2, is.null, logical(1)))) 
failed.multi.forest.sp <- c(paste(failed, "hr.multi.forest.cos.2", sep="."))
names(models.hr.multi.forest.cos.2) <- c(paste(names(models.hr.multi.forest.cos.2),"hr.multi.forest.cos.2", sep="."))
models.hr.multi.forest.cos.2 <- models.hr.multi.forest.cos.2[!(names(models.hr.multi.forest.cos.2) %in% failed.multi.forest.sp)]
#####hazard-rate, cosine 3
models.hr.multi.forest.cos.3 <- dlply(multi.species.data.forest, .variables = "Group", 
                                      .fun = possibly(hr.multi.cos.3, otherwise = NULL), 
                                      .parallel = TRUE)
failed <- names(which(vapply(models.hr.multi.forest.cos.3, is.null, logical(1)))) 
failed.multi.forest.sp <- c(paste(failed, "hr.multi.forest.cos.3", sep="."))
names(models.hr.multi.forest.cos.3) <- c(paste(names(models.hr.multi.forest.cos.3),"hr.multi.forest.cos.3", sep="."))
models.hr.multi.forest.cos.3 <- models.hr.multi.forest.cos.3[!(names(models.hr.multi.forest.cos.3) %in% failed.multi.forest.sp)]
#####hazard-rate, polynomial 2
models.hr.multi.forest.poly.2 <- dlply(multi.species.data.forest, .variables = "Group", 
                                       .fun = possibly(hr.multi.poly.2, otherwise = NULL), 
                                       .parallel = TRUE)
failed <- names(which(vapply(models.hr.multi.forest.poly.2, is.null, logical(1)))) 
failed.multi.forest.sp <- c(paste(failed, "hr.multi.forest.poly.2", sep="."))
names(models.hr.multi.forest.poly.2) <- c(paste(names(models.hr.multi.forest.poly.2),"hr.multi.forest.poly.2", sep="."))
models.hr.multi.forest.poly.2 <- models.hr.multi.forest.poly.2[!(names(models.hr.multi.forest.poly.2) %in% failed.multi.forest.sp)]
#####hazard-rate, polynomial 4
models.hr.multi.forest.poly.4 <- dlply(multi.species.data.forest, .variables = "Group", 
                                       .fun = possibly(hr.multi.poly.4, otherwise = NULL), 
                                       .parallel = TRUE)
failed <- names(which(vapply(models.hr.multi.forest.poly.4, is.null, logical(1)))) 
failed.multi.forest.sp <- c(paste(failed, "hr.multi.forest.poly.4", sep="."))
names(models.hr.multi.forest.poly.4) <- c(paste(names(models.hr.multi.forest.poly.4),"hr.multi.forest.poly.4", sep="."))
models.hr.multi.forest.poly.4 <- models.hr.multi.forest.poly.4[!(names(models.hr.multi.forest.poly.4) %in% failed.multi.forest.sp)]

#################################################
## 4.) FARMLAnd only detection functions #### 
# half-normal, no adjustment terms
models.hn.multi.farm.no.adj <- dlply(multi.species.data.farm, .variables = "Group", 
                                     .fun = possibly(hn.multi.no.adj, otherwise = NULL), 
                                     .parallel = TRUE)
failed <- names(which(vapply(models.hn.multi.farm.no.adj, is.null, logical(1)))) 
failed.multi.farm.sp <- c(paste(failed, "hn.multi.farm.no.adj", sep="."))
names(models.hn.multi.farm.no.adj) <- c(paste(names(models.hn.multi.farm.no.adj),"hn.multi.farm.no.adj", sep="."))
models.hn.multi.farm.no.adj <- models.hn.multi.farm.no.adj[!(names(models.hn.multi.farm.no.adj) %in% failed.multi.farm.sp)]
##### half-normal, cosine 2
models.hn.multi.farm.cos.2 <- dlply(multi.species.data.farm, .variables = "Group", 
                                    .fun = possibly(hn.multi.cos.2, otherwise = NULL), 
                                    .parallel = TRUE)
names(models.hn.multi.farm.cos.2) <- c(paste(names(models.hn.multi.farm.cos.2), "hn.multi.farm.cos.2", sep="."))
failed <- names(which(vapply(models.hn.multi.farm.cos.2, is.null, logical(1)))) 
failed.multi.farm.sp <- c(failed.multi.farm.sp, failed)
models.hn.multi.farm.cos.2 <- models.hn.multi.farm.cos.2[!(names(models.hn.multi.farm.cos.2) %in% failed)]
#####half-normal, cosine 3
models.hn.multi.farm.cos.3 <- dlply(multi.species.data.farm, .variables = "Group", 
                                    .fun = possibly(hn.multi.cos.3, otherwise = NULL), 
                                    .parallel = TRUE)
names(models.hn.multi.farm.cos.3) <- c(paste(names(models.hn.multi.farm.cos.3), "hn.multi.farm.cos.3", sep="."))
failed <- names(which(vapply(models.hn.multi.farm.cos.3, is.null, logical(1)))) 
failed.multi.farm.sp <- c(failed.multi.farm.sp, failed)
models.hn.multi.farm.cos.3 <- models.hn.multi.farm.cos.3[!(names(models.hn.multi.farm.cos.3) %in% failed)]
#####half-normal, hermitic polynomial 4 
models.hn.multi.farm.herm.4 <- dlply(multi.species.data.farm, .variables = "Group", 
                                     .fun = possibly(hn.multi.herm.4, otherwise = NULL), 
                                     .parallel = TRUE)
names(models.hn.multi.farm.herm.4) <- c(paste(names(models.hn.multi.farm.herm.4), "hn.multi.farm.herm.4", sep="."))
failed <- names(which(vapply(models.hn.multi.farm.herm.4, is.null, logical(1)))) 
failed.multi.farm.sp <- c(failed.multi.farm.sp, failed)
models.hn.multi.farm.herm.4 <- models.hn.multi.farm.herm.4[!(names(models.hn.multi.farm.herm.4) %in% failed)]
#####half-normal, hermitic polynomial 6
models.hn.multi.farm.herm.6 <- dlply(multi.species.data.farm, .variables = "Group", 
                                     .fun = possibly(hn.multi.herm.6, otherwise = NULL), 
                                     .parallel = TRUE)
names(models.hn.multi.farm.herm.6) <- c(paste(names(models.hn.multi.farm.herm.6), "hn.multi.farm.herm.6", sep="."))
failed <- names(which(vapply(models.hn.multi.farm.herm.6, is.null, logical(1)))) 
failed.multi.farm.sp <- c(failed.multi.farm.sp, failed)
models.hn.multi.farm.herm.6 <- models.hn.multi.farm.herm.6[!(names(models.hn.multi.farm.herm.6) %in% failed)]
#####hazard rate, no adjustment terms
models.hr.multi.farm.no.adj <- dlply(multi.species.data.farm, .variables = "Group", 
                                     .fun = possibly(hr.multi.no.adj, otherwise = NULL), 
                                     .parallel = TRUE)
failed <- names(which(vapply(models.hr.multi.farm.no.adj, is.null, logical(1)))) 
failed.multi.farm.sp <- c(paste(failed, "hr.multi.farm.no.adj", sep="."))
names(models.hr.multi.farm.no.adj) <- c(paste(names(models.hr.multi.farm.no.adj),"hr.multi.farm.no.adj", sep="."))
models.hr.multi.farm.no.adj <- models.hr.multi.farm.no.adj[!(names(models.hr.multi.farm.no.adj) %in% failed.multi.farm.sp)]
#####hazard-rate, cosine 2
models.hr.multi.farm.cos.2 <- dlply(multi.species.data.farm, .variables = "Group", 
                                    .fun = possibly(hr.multi.cos.2, otherwise = NULL), 
                                    .parallel = TRUE)
failed <- names(which(vapply(models.hr.multi.farm.cos.2, is.null, logical(1)))) 
failed.multi.farm.sp <- c(paste(failed, "hr.multi.farm.cos.2", sep="."))
names(models.hr.multi.farm.cos.2) <- c(paste(names(models.hr.multi.farm.cos.2),"hr.multi.farm.cos.2", sep="."))
models.hr.multi.farm.cos.2 <- models.hr.multi.farm.cos.2[!(names(models.hr.multi.farm.cos.2) %in% failed.multi.farm.sp)]
##### hazard-rate, cosine 3
models.hr.multi.farm.cos.3 <- dlply(multi.species.data.farm, .variables = "Group", 
                                    .fun = possibly(hr.multi.cos.3, otherwise = NULL), 
                                    .parallel = TRUE)
failed <- names(which(vapply(models.hr.multi.farm.cos.3, is.null, logical(1)))) 
failed.multi.farm.sp <- c(paste(failed, "hr.multi.farm.cos.3", sep="."))
names(models.hr.multi.farm.cos.3) <- c(paste(names(models.hr.multi.farm.cos.3),"hr.multi.farm.cos.3", sep="."))
models.hr.multi.farm.cos.3 <- models.hr.multi.farm.cos.3[!(names(models.hr.multi.farm.cos.3) %in% failed.multi.farm.sp)]
#####hazard-rate, polynomial 2
models.hr.multi.farm.poly.2 <- dlply(multi.species.data.farm, .variables = "Group", 
                                     .fun = possibly(hr.multi.poly.2, otherwise = NULL), 
                                     .parallel = TRUE)
failed <- names(which(vapply(models.hr.multi.farm.poly.2, is.null, logical(1)))) 
failed.multi.farm.sp <- c(paste(failed, "hr.multi.farm.poly.2", sep="."))
names(models.hr.multi.farm.poly.2) <- c(paste(names(models.hr.multi.farm.poly.2),"hr.multi.farm.poly.2", sep="."))
models.hr.multi.farm.poly.2 <- models.hr.multi.farm.poly.2[!(names(models.hr.multi.farm.poly.2) %in% failed.multi.farm.sp)]
#####hazard-rate, polynomial 4
models.hr.multi.farm.poly.4 <- dlply(multi.species.data.farm, .variables = "Group", 
                                     .fun = possibly(hr.multi.poly.4, otherwise = NULL), 
                                     .parallel = TRUE)
failed <- names(which(vapply(models.hr.multi.farm.poly.4, is.null, logical(1)))) 
failed.multi.farm.sp <- c(paste(failed, "hr.multi.farm.poly.4", sep="."))
names(models.hr.multi.farm.poly.4) <- c(paste(names(models.hr.multi.farm.poly.4),"hr.multi.farm.poly.4", sep="."))
models.hr.multi.farm.poly.4 <- models.hr.multi.farm.poly.4[!(names(models.hr.multi.farm.poly.4) %in% failed.multi.farm.sp)]

################################################## save 
# Put all the multi spp models into a single list, save, clean up
#####
multi.sp.models <- c(models.hn.multi.cos.2, models.hn.multi.cos.3, models.hn.multi.herm.4, models.hn.multi.herm.6, 
                     models.hn.multi.no.adj, models.hr.multi.cos.2, models.hr.multi.cos.3, models.hr.multi.no.adj,
                     models.hr.multi.poly.2, models.hr.multi.poly.4,
                     models.hn.multi.closed.cos.2, models.hn.multi.closed.cos.3, models.hn.multi.closed.herm.4, models.hn.multi.closed.herm.6, 
                     models.hn.multi.closed.no.adj, models.hr.multi.closed.cos.2, models.hr.multi.closed.cos.3, models.hr.multi.closed.no.adj,
                     models.hr.multi.closed.poly.2, models.hr.multi.closed.poly.4,
                     models.hn.multi.forest.cos.2, models.hn.multi.forest.cos.3, models.hn.multi.forest.herm.4, models.hn.multi.forest.herm.6, 
                     models.hn.multi.forest.no.adj, models.hr.multi.forest.cos.2, models.hr.multi.forest.cos.3, models.hr.multi.forest.no.adj,
                     models.hr.multi.forest.poly.2, models.hr.multi.forest.poly.4,
                     models.hn.multi.farm.cos.2, models.hn.multi.farm.cos.3, models.hn.multi.farm.herm.4, models.hn.multi.farm.herm.6, 
                     models.hn.multi.farm.no.adj, models.hr.multi.farm.cos.2, models.hr.multi.farm.cos.3, models.hr.multi.farm.no.adj,
                     models.hr.multi.farm.poly.2, models.hr.multi.farm.poly.4)
save(multi.sp.models, file = "MultiSpMods.rda")

save(multi.sp.models, file = "MultiSpMods.rda")

multi.sp.models.habitats <- c(
  models.hn.multi.forest.cos.2, models.hn.multi.forest.cos.3, models.hn.multi.forest.herm.4, models.hn.multi.forest.herm.6, 
  models.hn.multi.forest.no.adj, models.hr.multi.forest.cos.2, models.hr.multi.forest.cos.3, models.hr.multi.forest.no.adj,
  models.hr.multi.forest.poly.2, models.hr.multi.forest.poly.4,
  models.hn.multi.farm.cos.2, models.hn.multi.farm.cos.3, models.hn.multi.farm.herm.4, models.hn.multi.farm.herm.6, 
  models.hn.multi.farm.no.adj, models.hr.multi.farm.cos.2, models.hr.multi.farm.cos.3, models.hr.multi.farm.no.adj,
  models.hr.multi.farm.poly.2, models.hr.multi.farm.poly.4)

save(multi.sp.models.habitats, distance.data, distance.data.orig, failed.single.sp, failed.multi.sp, bird.counts6,
     file = "MultiSpModsHabitats.rda")

save(single.sp.models.habitats, single.species.data, region.table, sample.table, 
     multi.sp.models.habitats, multi.species.data,
     distance.data, distance.data.orig, failed.single.sp, failed.multi.sp, bird.counts6, file = "AllFittedDetectionFunctionsHabitats.rda")


###################################### 2. EVALUATE DETECTION FUNCTIONS
##############Code runs goodness of fit (GOF) tests, rejecting models where all three are significant,
#looks at the detection functions themselves, checks that each species has at least one model (either a single species one, or a group one)
   

# GOF tests
##if significant then the observed records/data are significantly different from the model
gof4 <- function(x) {
  ddf.gof(x, main=x)
}
gof.list.single4 <- llply(single.sp.models, purrr::possibly(gof4, otherwise = NULL))
which(vapply(gof.list.single4, is.null, logical(1))) ##all are ok

gof.list.multi4 <- llply(multi.sp.models, purrr::possibly(gof4, otherwise = NULL))
which(vapply(gof.list.multi4, is.null, logical(1))) ##all are ok

#########forest & farmland specific ones
gof.list.single4h <- llply(single.sp.models.habitats, purrr::possibly(gof4, otherwise = NULL))
which(vapply(gof.list.single4h, is.null, logical(1))) ##all are ok

gof.list.multi4h <- llply(multi.sp.models.habitats, purrr::possibly(gof4, otherwise = NULL))
which(vapply(gof.list.multi4h, is.null, logical(1))) ##all are ok

####Extract GOF p-values
chi.p <- function(x) x$chisquare$chi1$p
KS <- function(x) x$dsgof$ks$p
CvM <- function(x) x$dsgof$CvM$p
# Single species models
gof.table.single <- data.frame(model.name=names(gof.list.single4), chi.p=ldply(gof.list.single4,chi.p)[,2], 
                               KS=ldply(gof.list.single4, KS)[,2], CvM=ldply(gof.list.single4, CvM)[,2])
gof.table.single$ID <- 1:length(gof.table.single$model.name)                      
###multi spp
gof.table.multi <- data.frame(model.name=names(gof.list.multi4), chi.p=ldply(gof.list.multi4,chi.p)[,2], 
                              KS=ldply(gof.list.multi4, KS)[,2], CvM=ldply(gof.list.multi4, CvM)[,2])
gof.table.multi$ID <- 1:length(gof.table.multi$model.name)                      
##farm & forest
# Single species models
gof.table.single.h <- data.frame(model.name=names(gof.list.single4h), chi.p=ldply(gof.list.single4h,chi.p)[,2], 
                                 KS=ldply(gof.list.single4h, KS)[,2], CvM=ldply(gof.list.single4h, CvM)[,2])
gof.table.single.h$ID <- 1:length(gof.table.single.h$model.name)                      
###multi spp
gof.table.multi.h <- data.frame(model.name=names(gof.list.multi4h), chi.p=ldply(gof.list.multi4h,chi.p)[,2], 
                                KS=ldply(gof.list.multi4h, KS)[,2], CvM=ldply(gof.list.multi4h, CvM)[,2])
gof.table.multi.h$ID <- 1:length(gof.table.multi.h$model.name)   
#####Reject models with all GOF tests giving significant values
reject.gof.table.single <- subset(gof.table.single, chi.p < 0.05 & KS < 0.05 & CvM < 0.05)
kept.single.sp.models <- single.sp.models[!(names(single.sp.models) %in% reject.gof.table.single$model.name)]

reject.gof.table.multi <- subset(gof.table.multi, chi.p < 0.05 & KS < 0.05 & CvM < 0.05)
kept.multi.sp.models <- multi.sp.models[!(names(multi.sp.models) %in% reject.gof.table.multi$model.name)]

gof.rejects.single <- reject.gof.table.single$model.name ##none rejected based on GOF values
gof.rejects.multi <- reject.gof.table.multi$model.name
#####
reject.gof.table.single.h <- subset(gof.table.single.h, chi.p < 0.05 & KS < 0.05 & CvM < 0.05)
kept.single.sp.models.h <- single.sp.models.habitats[!(names(single.sp.models.habitats) %in% reject.gof.table.single.h$model.name)]

reject.gof.table.multi.h <- subset(gof.table.multi.h, chi.p < 0.05 & KS < 0.05 & CvM < 0.05)
kept.multi.sp.models.h <- multi.sp.models.habitats[!(names(multi.sp.models.habitats) %in% reject.gof.table.multi.h$model.name)]

gof.rejects.single.h <- reject.gof.table.single.h$model.name ##none rejected based on GOF values
gof.rejects.multi.h <- reject.gof.table.multi.h$model.name

##################################################################
# Cutting bad detection functions
###remove functions where g(x) > 0 at any point 
##### 1. SINGLE
bad.det.funs.single <- c("Acridotheres_tristis.hn.single.closed.cos.2",
                         "Acridotheres_tristis.hn.single.closed.cos.3",
                         "Acridotheres_tristis.hn.single.cos.2",
                         "Acridotheres_tristis.hn.single.cos.3",
                         "Acridotheres_tristis.hr.single.closed.cos.2",
                         "Acridotheres_tristis.hr.single.closed.cos.3",
                         "Acridotheres_tristis.hr.single.closed.poly.4",
                         "Acridotheres_tristis.hr.single.cos.3",
                         "Acridotheres_tristis.hr.single.poly.4",
                         "Acridotheres_tristis.single.hr.cos.2",
                         "Acrocephalus_dumetorum.single.hr.cos.2",
                         "Alauda_gulgula.hn.single.cos.3",
                         "Alauda_gulgula.hn.single.closed.cos.3",
                         "Alauda_gulgula.hr.single.closed.cos.3",
                         "Alauda_gulgula.hr.single.closed.poly.4",
                         "Alauda_gulgula.hr.single.cos.3",
                         "Alauda_gulgula.hr.single.poly.4",
                         "Alexandrinus_krameri.hn.single.cos.3",
                         "Alexandrinus_krameri.hr.single.closed.cos.3",
                         "Anthus_rufulus.hn.single.cos.3",
                         "Anthus_rufulus.hn.single.closed.cos.3",
                         "Alexandrinus_krameri.single.hr.cos.2",
                         "Anthus_rufulus.hr.single.closed.cos.2",
                         "Anthus_rufulus.hr.single.closed.cos.3",
                         "Anthus_rufulus.hr.single.cos.3",
                         "Anthus_rufulus.single.hr.cos.2",
                         "Ardeola_grayii.hn.single.closed.cos.2",
                         "Ardeola_grayii.hn.single.closed.cos.3",
                         "Ardeola_grayii.hn.single.closed.herm.4",
                         "Ardeola_grayii.hn.single.closed.herm.6",
                         "Ardeola_grayii.hn.single.closed.no.adj",
                         "Ardeola_grayii.hn.single.cos.2",
                         "Ardeola_grayii.hn.single.cos.3",
                         "Ardeola_grayii.hn.single.herm.4",
                         "Ardeola_grayii.hn.single.herm.6",
                         "Ardeola_grayii.hr.single.closed.cos.2",
                         "Ardeola_grayii.hr.single.closed.poly.4",
                         "Ardeola_grayii.hr.single.poly.4",
                         "Ardeola_grayii.single.hr.cos.2",
                         "Argya_striata.hn.single.closed.cos.2",
                         "Argya_striata.hn.single.closed.cos.3",
                         "Argya_striata.hn.single.closed.herm.4",
                         "Argya_striata.hn.single.closed.herm.6",
                         "Argya_striata.hn.single.closed.no.adj",
                         "Argya_striata.hn.single.cos.2",
                         "Argya_striata.hn.single.cos.3",
                         "Argya_striata.hr.single.closed.cos.3",
                         "Argya_striata.hr.single.closed.poly.4",
                         "Argya_striata.hr.single.cos.3",
                         "Argya_striata.single.hr.cos.2",
                         "Bubulcus_ibis.hn.single.closed.cos.2",
                         "Bubulcus_ibis.hn.single.closed.cos.3",
                         "Bubulcus_ibis.hn.single.closed.herm.4",
                         "Bubulcus_ibis.hn.single.closed.herm.6",
                         "Bubulcus_ibis.hn.single.cos.2",
                         "Bubulcus_ibis.hn.single.cos.3",
                         "Bubulcus_ibis.hn.single.herm.4",
                         "Bubulcus_ibis.hn.single.herm.6",
                         "Bubulcus_ibis.hr.single.closed.cos.2",
                         "Bubulcus_ibis.hr.single.closed.cos.3",
                         "Bubulcus_ibis.hr.single.closed.poly.4",
                         "Bubulcus_ibis.hr.single.cos.3",
                         "Bubulcus_ibis.hr.single.poly.4",
                         "Bubulcus_ibis.single.hr.cos.2",
                         "Cacomantis_sonneratii.hn.single.cos.2",
                         "Cacomantis_sonneratii.hn.single.cos.3",
                         "Cacomantis_sonneratii.hn.single.herm.4",
                         "Cacomantis_sonneratii.hn.single.herm.6",
                         "Cacomantis_sonneratii.hn.single.no.adj",
                         "Cacomantis_sonneratii.single.hr.cos.2",
                         "Centropus_sinensis.hn.single.closed.cos.3",
                         "Centropus_sinensis.hn.single.closed.herm.4",
                         "Centropus_sinensis.hn.single.closed.herm.6",
                         "Centropus_sinensis.hn.single.closed.no.adj",
                         "Centropus_sinensis.hn.single.cos.3",
                         "Centropus_sinensis.hr.single.closed.cos.2",
                         "Centropus_sinensis.hr.single.closed.cos.3",
                         "Centropus_sinensis.hr.single.closed.no.adj",
                         "Chloropsis_aurifrons.hr.single.cos.3",
                         "Chloropsis_aurifrons.hr.single.no.adj",
                         "Chloropsis_aurifrons.hr.single.poly.4",
                         "Chloropsis_aurifrons.single.hr.cos.2",
                         "Chloropsis_jerdoni.hn.single.closed.cos.3",
                         "Chloropsis_jerdoni.hn.single.cos.3",
                         "Chloropsis_jerdoni.hr.single.closed.cos.3",
                         "Chloropsis_jerdoni.hr.single.cos.3",
                         "Cinnyris_asiaticus.hr.single.closed.cos.2",
                         "Cinnyris_asiaticus.hr.single.closed.cos.3",
                         "Cinnyris_asiaticus.hr.single.cos.3",
                         "Cinnyris_asiaticus.single.hr.cos.2",
                         "Copsychus_saularis.hn.single.closed.cos.3",
                         "Copsychus_saularis.hn.single.cos.3",
                         "Copsychus_saularis.hr.single.closed.cos.3",
                         "Copsychus_saularis.hr.single.closed.poly.4",
                         "Copsychus_saularis.hr.single.cos.3",
                         "Copsychus_saularis.hr.single.poly.4",
                         "Coracias_benghalensis.hn.single.closed.cos.2",
                         "Coracias_benghalensis.hn.single.closed.cos.3",
                         "Coracias_benghalensis.hn.single.cos.2",
                         "Coracias_benghalensis.hn.single.cos.3",
                         "Coracias_benghalensis.hr.single.closed.cos.2",
                         "Coracias_benghalensis.hr.single.closed.cos.3",
                         "Coracias_benghalensis.hr.single.closed.poly.4",
                         "Coracias_benghalensis.hr.single.cos.3",
                         "Coracias_benghalensis.hr.single.poly.4",
                         "Coracias_benghalensis.single.hr.cos.2",
                         "Corvus_macrorhynchos.hn.single.closed.cos.3",
                         "Corvus_macrorhynchos.hn.single.cos.3",
                         "Corvus_macrorhynchos.hr.single.closed.cos.2",
                         "Corvus_macrorhynchos.hr.single.closed.cos.3",
                         "Corvus_macrorhynchos.hr.single.cos.3",
                         "Corvus_macrorhynchos.hr.single.no.adj",
                         "Corvus_macrorhynchos.hr.single.poly.4",
                         "Cyornis_tickelliae.single.hr.cos.2",
                         "Dicaeum_agile.single.hr.cos.2",
                         "Dicaeum_erythrorhynchos.hr.single.cos.3",
                         "Dicaeum_erythrorhynchos.single.hr.cos.2",
                         "Dicrurus_aeneus.hn.single.cos.3",
                         "Dicrurus_aeneus.hr.single.cos.3",
                         "Dicrurus_aeneus.hr.single.poly.4",
                         "Dicrurus_macrocercus.hr.single.closed.cos.3",
                         "Dinopium_benghalense.hn.single.closed.cos.3",
                         "Dinopium_benghalense.hn.single.cos.3",
                         "Dinopium_benghalense.hr.single.closed.cos.2",
                         "Dinopium_benghalense.hr.single.closed.cos.3",
                         "Dinopium_benghalense.hr.single.no.adj",
                         "Dinopium_benghalense.hr.single.poly.4",
                         "Dinopium_benghalense.single.hr.cos.2",
                         "Dinopium_benghalense.hr.single.cos.3",
                         "Gallus_gallus.hr.single.poly.4",
                         "Gallus_gallus.single.hr.cos.2",
                         "Halcyon_smyrnensis.hn.single.closed.cos.3",
                         "Halcyon_smyrnensis.hn.single.cos.3",
                         "Halcyon_smyrnensis.hr.single.closed.cos.2",
                         "Halcyon_smyrnensis.hr.single.closed.poly.4",
                         "Halcyon_smyrnensis.hr.single.cos.3",
                         "Halcyon_smyrnensis.hr.single.poly.4",
                         "Hypothymis_azurea.hn.single.cos.3",
                         "Hypothymis_azurea.hn.single.herm.4",
                         "Hypothymis_azurea.hn.single.herm.6",
                         "Hypothymis_azurea.hn.single.no.adj",
                         "Hypothymis_azurea.hr.single.cos.3",
                         "Hypothymis_azurea.hr.single.no.adj",
                         "Hypothymis_azurea.hr.single.poly.4",
                         "Leptocoma_zeylonica.hn.single.closed.cos.3",
                         "Leptocoma_zeylonica.hn.single.cos.3",
                         "Leptocoma_zeylonica.hr.single.closed.cos.3",
                         "Leptocoma_zeylonica.hr.single.cos.3",
                         "Loriculus_vernalis.hn.single.closed.herm.4",
                         "Loriculus_vernalis.hn.single.closed.herm.6",
                         "Loriculus_vernalis.hn.single.closed.no.adj",
                         "Loriculus_vernalis.hn.single.herm.4",
                         "Loriculus_vernalis.hn.single.herm.6",
                         "Loriculus_vernalis.hn.single.no.adj",
                         "Machlolophus_xanthogenys.hn.single.closed.cos.3",
                         "Machlolophus_xanthogenys.hn.single.closed.no.adj",
                         "Machlolophus_xanthogenys.hn.single.cos.3",
                         "Machlolophus_xanthogenys.hr.single.closed.cos.3",
                         "Machlolophus_xanthogenys.hr.single.closed.no.adj",
                         "Machlolophus_xanthogenys.hr.single.cos.3",
                         "Mixornis_gularis.hr.single.cos.3",
                         "Mixornis_gularis.hr.single.closed.poly.4",
                         "Mixornis_gularis.hr.single.closed.cos.3",
                         "Mixornis_gularis.hr.single.poly.4",
                         "Orthotomus_sutorius.hr.single.closed.cos.2",
                         "Orthotomus_sutorius.hr.single.closed.cos.3",
                         "Orthotomus_sutorius.hr.single.cos.3",
                         "Orthotomus_sutorius.single.hr.cos.2",
                         "Parus_major.hn.single.cos.3",
                         "Parus_major.hr.single.cos.3",
                         "Pellorneum_ruficeps.hn.single.cos.3",
                         "Pellorneum_ruficeps.hr.single.cos.3",
                         "Pellorneum_ruficeps.hr.single.no.adj",
                         "Pellorneum_ruficeps.hr.single.poly.4",
                         "Pellorneum_ruficeps.single.hr.cos.2",
                         "Pericrocotus_flammeus.hn.single.cos.3",
                         "Pericrocotus_flammeus.hr.single.cos.3",
                         "Phylloscopus_trochiloides.hn.single.closed.herm.4",
                         "Phylloscopus_trochiloides.hn.single.closed.herm.6",
                         "Phylloscopus_trochiloides.hn.single.closed.no.adj",
                         "Phylloscopus_trochiloides.hn.single.herm.4",
                         "Phylloscopus_trochiloides.hn.single.herm.6",
                         "Phylloscopus_trochiloides.hn.single.no.adj",
                         "Pomatorhinus_horsfieldii.hn.single.cos.2",
                         "Pomatorhinus_horsfieldii.hn.single.cos.3",
                         "Pomatorhinus_horsfieldii.hn.single.herm.4",
                         "Pomatorhinus_horsfieldii.hn.single.herm.6",
                         "Pomatorhinus_horsfieldii.hn.single.no.adj",
                         "Pomatorhinus_horsfieldii.hr.single.cos.3",
                         "Pomatorhinus_horsfieldii.hr.single.no.adj",
                         "Pomatorhinus_horsfieldii.hr.single.poly.4",
                         "Prinia_socialis.hn.single.closed.cos.2",
                         "Prinia_socialis.hn.single.closed.cos.3",
                         "Prinia_socialis.hr.single.closed.cos.2",
                         "Prinia_socialis.hr.single.no.adj",
                         "Psilopogon_haemacephalus.hn.single.closed.cos.2",
                         "Psilopogon_haemacephalus.hn.single.closed.cos.3",
                         "Psilopogon_haemacephalus.hn.single.closed.herm.4",
                         "Psilopogon_haemacephalus.hn.single.closed.herm.6",
                         "Psilopogon_haemacephalus.hn.single.closed.no.adj",
                         "Psilopogon_haemacephalus.hn.single.cos.3",
                         "Psilopogon_haemacephalus.hn.single.herm.4",
                         "Psilopogon_haemacephalus.hn.single.herm.6",
                         "Psilopogon_haemacephalus.hn.single.no.adj",
                         "Psilopogon_haemacephalus.hr.single.closed.cos.2",
                         "Psilopogon_haemacephalus.hr.single.closed.cos.3",
                         "Psilopogon_haemacephalus.hr.single.closed.no.adj",
                         "Psilopogon_haemacephalus.hr.single.closed.poly.4",
                         "Psilopogon_haemacephalus.hr.single.cos.3",
                         "Psilopogon_haemacephalus.hr.single.no.adj",
                         "Psilopogon_haemacephalus.hr.single.poly.4",
                         "Psilopogon_zeylanicus.hn.single.closed.cos.3",
                         "Psilopogon_zeylanicus.hn.single.closed.herm.4",
                         "Psilopogon_zeylanicus.hn.single.closed.herm.6",
                         "Psilopogon_zeylanicus.hn.single.cos.3",
                         "Psilopogon_zeylanicus.hn.single.herm.4",
                         "Psilopogon_zeylanicus.hn.single.herm.6",
                         "Psilopogon_zeylanicus.hn.single.no.adj",
                         "Psilopogon_zeylanicus.hr.single.closed.cos.3",
                         "Psilopogon_zeylanicus.hr.single.closed.no.adj",
                         "Psilopogon_zeylanicus.hr.single.closed.poly.4",
                         "Psilopogon_zeylanicus.hr.single.cos.3",
                         "Psilopogon_zeylanicus.hr.single.no.adj",
                         "Psilopogon_zeylanicus.hr.single.poly.4",
                         "Pycnonotus_cafer.hr.single.closed.cos.2",
                         "Pycnonotus_cafer.single.hr.cos.2",
                         "Pycnonotus_jocosus.hn.single.cos.2",
                         "Pycnonotus_jocosus.hr.single.closed.cos.2",
                         "Pycnonotus_jocosus.single.hr.cos.2",
                         "Pycnonotus_luteolus.hr.single.closed.cos.3",
                         "Pycnonotus_luteolus.hr.single.cos.3",
                         "Rubigula_flaviventris.hn.single.closed.herm.4",
                         "Rubigula_flaviventris.hn.single.closed.herm.6",
                         "Rubigula_flaviventris.hn.single.closed.no.adj",
                         "Rubigula_flaviventris.hr.single.cos.3",
                         "Rubigula_flaviventris.hr.single.closed.poly.4",
                         "Rubigula_flaviventris.hr.single.closed.no.adj",
                         "Rubigula_flaviventris.hr.single.closed.cos.3",
                         "Rubigula_flaviventris.hr.single.closed.cos.2",
                         "Rubigula_flaviventris.hn.single.no.adj",
                         "Rubigula_flaviventris.hn.single.herm.6",
                         "Rubigula_flaviventris.hn.single.herm.4",
                         "Rubigula_flaviventris.hr.single.no.adj",
                         "Rubigula_flaviventris.hr.single.poly.4",
                         "Rubigula_flaviventris.single.hr.cos.2",
                         "Sitta_frontalis.hr.single.cos.3",
                         "Sitta_frontalis.hr.single.no.adj",
                         "Sitta_frontalis.hr.single.poly.4",
                         "Sitta_frontalis.single.hr.cos.2",
                         "Spilopelia_chinensis.hn.single.closed.cos.3",
                         "Spilopelia_chinensis.hn.single.cos.3",
                         "Spilopelia_chinensis.hr.single.closed.cos.3",
                         "Spilopelia_chinensis.single.hr.cos.2",
                         "Spilopelia_suratensis.hn.single.closed.cos.3",
                         "Spilopelia_suratensis.hn.single.cos.3",
                         "Spilopelia_suratensis.hr.single.closed.cos.2",
                         "Spilopelia_suratensis.hr.single.closed.cos.3",
                         "Spilopelia_suratensis.hr.single.cos.3",
                         "Spilopelia_suratensis.single.hr.cos.2",
                         "Vanellus_indicus.hn.single.closed.cos.3",
                         "Vanellus_indicus.hn.single.cos.3",
                         "Vanellus_indicus.hr.single.closed.cos.2",
                         "Vanellus_indicus.hr.single.closed.cos.3",
                         "Vanellus_indicus.hr.single.closed.no.adj",
                         "Vanellus_indicus.hr.single.closed.poly.4",
                         "Vanellus_indicus.hr.single.cos.3",
                         "Zosterops_palpebrosus.hr.single.closed.cos.3",
                         "Zosterops_palpebrosus.hr.single.closed.poly.4",
                         "Zosterops_palpebrosus.hr.single.cos.3",
                         "Zosterops_palpebrosus.hr.single.poly.4")
######multi
bad.det.funs.multi <- c("arboreal_invertivore.hn.multi.closed.no.adj",
                        "arboreal_invertivore.hr.multi.closed.cos.3",
                        "arboreal_invertivore.hr.multi.cos.3",
                        "arboreal_invertivore.hr.multi.closed.poly.4",
                        "babbler.hr.multi.closed.cos.3",
                        "babbler.hr.multi.cos.3",
                        "crow.hn.multi.closed.cos.3",
                        "crow.hn.multi.cos.3",
                        "crow.hr.multi.closed.cos.3",
                        "crow.hr.multi.cos.3",
                        "doves.hn.multi.closed.cos.3",
                        "doves.hn.multi.cos.3",
                        "doves.hr.multi.closed.cos.3",
                        "doves.hr.multi.cos.3",
                        "flycatcher.hr.multi.closed.cos.3",
                        "flycatcher.hr.multi.cos.3",
                        "ground_dweller_large.hn.multi.closed.cos.3",
                        "ground_dweller_large.hn.multi.cos.3",
                        "ground_feeding_passerines.hn.multi.closed.cos.3",
                        "ground_feeding_passerines.hn.multi.cos.3",
                        "kingfisher.hr.multi.closed.cos.2",
                        "large_wader.hn.multi.closed.cos.2",
                        "large_wader.hn.multi.closed.herm.4",
                        "large_wader.hn.multi.closed.herm.6",
                        "large_wader.hn.multi.cos.2",
                        "large_wader.hn.multi.herm.4",
                        "large_wader.hn.multi.herm.6",
                        "large_wader.hr.multi.closed.cos.2",
                        "large_wader.hr.multi.closed.cos.3",
                        "large_wader.hr.multi.closed.poly.4",
                        "large_wader.hr.multi.cos.3",
                        "large_wader.hr.multi.poly.4",
                        "large_wader.multi.hr.cos.2",
                        "myna.hn.multi.closed.cos.2",
                        "myna.hn.multi.cos.2",
                        "myna.hn.multi.cos.3",
                        "myna.hr.multi.closed.cos.2",
                        "myna.hr.multi.closed.cos.3",
                        "myna.hr.multi.closed.poly.4",
                        "myna.hr.multi.cos.3",
                        "myna.hr.multi.poly.4",
                        "pipit.hn.multi.closed.cos.3",
                        "pipit.hn.multi.cos.3",
                        "pipit.hr.multi.closed.cos.2",
                        "pipit.multi.hr.cos.2",
                        "roller_weaver.hn.multi.closed.cos.2",
                        "roller_weaver.hn.multi.closed.cos.3",
                        "roller_weaver.hn.multi.cos.2",
                        "roller_weaver.hn.multi.cos.3",
                        "roller_weaver.hn.multi.farm.cos.2",
                        "roller_weaver.hr.multi.closed.cos.2",
                        "roller_weaver.hr.multi.closed.cos.3",
                        "roller_weaver.hr.multi.cos.3",
                        "roller_weaver.hr.multi.poly.4",
                        "roller_weaver.multi.hr.cos.2",
                        "sunbirds.hr.multi.closed.cos.3",
                        "sunbirds.hr.multi.cos.3",
                        "warblers.hr.multi.closed.cos.3",
                        "warblers.hr.multi.cos.3",
                        "woodpecker.hr.multi.closed.cos.2",
                        "woodpecker.hr.multi.closed.no.adj",
                        "woodpecker.hr.multi.no.adj",
                        "woodpecker.hr.multi.poly.4")



bad.det.funs.single.habitat <- c("Acridotheres_tristis.hn.single.farm.cos.3",
                                 "Acridotheres_tristis.hn.single.farm.cos.2",
                                 "Acridotheres_tristis.hr.single.farm.cos.2",
                                 "Acridotheres_tristis.hr.single.farm.cos.3",
                                 "Acridotheres_tristis.hr.single.farm.poly.4",
                                 "Alauda_gulgula.hn.single.farm.cos.3",
                                 "Alauda_gulgula.hr.single.farm.cos.3",
                                 "Alauda_gulgula.hr.single.farm.poly.4",
                                 "Alcippe_poioicephala.hn.single.forest.cos.3",
                                 "Alcippe_poioicephala.hn.single.forest.herm.4",
                                 "Alcippe_poioicephala.hn.single.forest.herm.6",
                                 "Alcippe_poioicephala.hn.single.forest.no.adj",
                                 "Alexandrinus_krameri.hn.single.farm.cos.3",
                                 "Alexandrinus_krameri.hr.single.farm.cos.2",
                                 "Alexandrinus_krameri.hr.single.farm.cos.3",
                                 "Anthus_rufulus.hn.single.farm.cos.3",
                                 "Anthus_rufulus.hr.single.farm.cos.2",
                                 "Anthus_rufulus.hr.single.farm.cos.3",
                                 "Ardeola_grayii.hn.single.farm.cos.2",
                                 "Ardeola_grayii.hn.single.farm.cos.3",
                                 "Ardeola_grayii.hn.single.farm.herm.4",
                                 "Ardeola_grayii.hn.single.farm.herm.6",
                                 "Ardeola_grayii.hr.single.farm.cos.2",
                                 "Ardeola_grayii.hr.single.farm.poly.4",
                                 "Bubulcus_ibis.hn.single.farm.cos.2",
                                 "Bubulcus_ibis.hn.single.farm.cos.3",
                                 "Bubulcus_ibis.hn.single.farm.herm.4",
                                 "Bubulcus_ibis.hn.single.farm.herm.6",
                                 "Bubulcus_ibis.hr.single.farm.cos.2",
                                 "Bubulcus_ibis.hr.single.farm.cos.3",
                                 "Bubulcus_ibis.hr.single.farm.poly.4",
                                 "Cacomantis_sonneratii.hr.single.forest.cos.2",
                                 "Centropus_sinensis.hn.single.farm.cos.3",
                                 "Centropus_sinensis.hr.single.farm.cos.3",
                                 "Centropus_sinensis.hr.single.farm.poly.4",
                                 "Chloropsis_aurifrons.hr.single.forest.cos.2",
                                 "Chloropsis_aurifrons.hr.single.forest.cos.3",
                                 "Chloropsis_aurifrons.hr.single.forest.no.adj",
                                 "Chloropsis_aurifrons.hr.single.forest.poly.4",
                                 "Chloropsis_jerdoni.hn.single.forest.cos.3",
                                 "Chloropsis_jerdoni.hr.single.forest.cos.3",
                                 "Cinnyris_asiaticus.hr.single.forest.cos.2",
                                 "Cinnyris_asiaticus.hr.single.forest.cos.3",
                                 "Copsychus_saularis.hn.single.forest.cos.2",
                                 "Copsychus_saularis.hn.single.forest.cos.3",
                                 "Copsychus_saularis.hr.single.forest.cos.3",
                                 "Copsychus_saularis.hr.single.forest.poly.4",
                                 "Coracias_benghalensis.hn.single.farm.cos.2",
                                 "Coracias_benghalensis.hn.single.farm.cos.3",
                                 "Coracias_benghalensis.hr.single.farm.cos.2",
                                 "Coracias_benghalensis.hr.single.farm.cos.3",
                                 "Coracias_benghalensis.hr.single.farm.poly.4",
                                 "Corvus_macrorhynchos.hn.single.farm.cos.3",
                                 "Corvus_macrorhynchos.hr.single.farm.cos.3",
                                 "Cyornis_tickelliae.hr.single.forest.cos.2",
                                 "Dicaeum_agile.hr.single.forest.cos.2",
                                 "Dicaeum_erythrorhynchos.hr.single.forest.cos.2",
                                 "Dicaeum_erythrorhynchos.hr.single.forest.cos.3",
                                 "Dicrurus_aeneus.hn.single.forest.cos.3",
                                 "Dicrurus_aeneus.hr.single.forest.cos.3",
                                 "Dinopium_benghalense.hn.single.forest.cos.2",
                                 "Dinopium_benghalense.hn.single.forest.cos.3",
                                 "Dinopium_benghalense.hr.single.forest.cos.2",
                                 "Dinopium_benghalense.hr.single.forest.cos.3",
                                 "Gallus_gallus.hr.single.forest.cos.2",
                                 "Gallus_gallus.hr.single.forest.poly.4",
                                 "Halcyon_smyrnensis.hn.single.farm.cos.2",
                                 "Halcyon_smyrnensis.hn.single.farm.herm.4",
                                 "Halcyon_smyrnensis.hn.single.farm.herm.6",
                                 "Halcyon_smyrnensis.hr.single.farm.cos.2",
                                 "Halcyon_smyrnensis.hr.single.farm.cos.3",
                                 "Halcyon_smyrnensis.hr.single.farm.poly.4",
                                 "Hypothymis_azurea.hr.single.forest.cos.3",
                                 "Hypothymis_azurea.hr.single.forest.no.adj",
                                 "Hypothymis_azurea.hr.single.forest.poly.4",
                                 "Leptocoma_zeylonica.hr.single.forest.cos.3",
                                 "Machlolophus_xanthogenys.hn.single.forest.cos.3",
                                 "Machlolophus_xanthogenys.hr.single.forest.cos.3",
                                 "Mixornis_gularis.hr.single.forest.cos.3",
                                 "Mixornis_gularis.hr.single.forest.poly.4",
                                 "Orthotomus_sutorius.hr.single.farm.cos.2",
                                 "Orthotomus_sutorius.hr.single.farm.cos.3",
                                 "Orthotomus_sutorius.hr.single.forest.cos.2",
                                 "Orthotomus_sutorius.hr.single.forest.cos.3",
                                 "Parus_major.hn.single.forest.cos.3",
                                 "Parus_major.hr.single.forest.cos.3",
                                 "Pellorneum_ruficeps.hn.single.forest.cos.3",
                                 "Pellorneum_ruficeps.hr.single.forest.cos.2",
                                 "Pellorneum_ruficeps.hr.single.forest.cos.3",
                                 "Pellorneum_ruficeps.hr.single.forest.no.adj",
                                 "Pellorneum_ruficeps.hr.single.forest.poly.4",
                                 "Pericrocotus_flammeus.hn.single.forest.cos.3",
                                 "Pericrocotus_flammeus.hr.single.forest.cos.3",
                                 "Pomatorhinus_horsfieldii.hn.single.forest.cos.3",
                                 "Pomatorhinus_horsfieldii.hn.single.forest.herm.4",
                                 "Pomatorhinus_horsfieldii.hn.single.forest.herm.6",
                                 "Pomatorhinus_horsfieldii.hn.single.forest.no.adj",
                                 "Pomatorhinus_horsfieldii.hr.single.forest.cos.3",
                                 "Pomatorhinus_horsfieldii.hr.single.forest.no.adj",
                                 "Pomatorhinus_horsfieldii.hr.single.forest.poly.4",
                                 "Prinia_socialis.hr.single.farm.cos.3",
                                 "Psilopogon_haemacephalus.hn.single.forest.cos.3",
                                 "Psilopogon_haemacephalus.hn.single.forest.herm.4",
                                 "Psilopogon_haemacephalus.hn.single.forest.herm.6",
                                 "Psilopogon_haemacephalus.hn.single.forest.no.adj",
                                 "Psilopogon_haemacephalus.hr.single.forest.cos.3",
                                 "Psilopogon_haemacephalus.hr.single.forest.no.adj",
                                 "Psilopogon_haemacephalus.hr.single.forest.poly.4",
                                 "Psilopogon_zeylanicus.hr.single.forest.cos.3",
                                 "Psilopogon_zeylanicus.hr.single.forest.no.adj",
                                 "Psilopogon_zeylanicus.hr.single.forest.poly.4",
                                 "Pycnonotus_cafer.hn.single.farm.cos.3",
                                 "Pycnonotus_cafer.hr.single.farm.cos.2",
                                 "Pycnonotus_cafer.hr.single.farm.cos.3",
                                 "Pycnonotus_cafer.hr.single.forest.cos.2",
                                 "Pycnonotus_jocosus.hr.single.forest.cos.2",
                                 "Rubigula_flaviventris.hr.single.forest.cos.2",
                                 "Rubigula_flaviventris.hr.single.forest.cos.3",
                                 "Rubigula_flaviventris.hr.single.forest.no.adj",
                                 "Rubigula_flaviventris.hr.single.forest.poly.4",
                                 "Sitta_frontalis.hr.single.forest.cos.2",
                                 "Sitta_frontalis.hr.single.forest.cos.3",
                                 "Sitta_frontalis.hr.single.forest.no.adj",
                                 "Sitta_frontalis.hr.single.forest.poly.4",
                                 "Spilopelia_chinensis.hn.single.farm.cos.2",
                                 "Spilopelia_chinensis.hn.single.farm.cos.3",
                                 "Spilopelia_chinensis.hr.single.farm.cos.3",
                                 "Spilopelia_chinensis.hr.single.farm.poly.4",
                                 "Spilopelia_chinensis.hr.single.forest.cos.3",
                                 "Vanellus_indicus.hn.single.farm.cos.3",
                                 "Vanellus_indicus.hr.single.farm.cos.3",
                                 "Zosterops_palpebrosus.hr.single.forest.cos.3")

bad.det.funs.multi.habitat <- c("arboreal_frugivore.hr.multi.farm.cos.3",
                                "arboreal_invertivore.hr.multi.farm.cos.3",
                                "babbler.hr.multi.farm.cos.2",
                                "babbler.hr.multi.farm.cos.3",
                                "crow.hn.multi.farm.cos.3",
                                "crow.hr.multi.farm.cos.3",
                                "doves.hn.multi.forest.cos.2",
                                "doves.hn.multi.forest.cos.3",
                                "doves.hr.multi.farm.cos.3",
                                "doves.hr.multi.forest.cos.3",
                                "flycatcher.hr.multi.farm.cos.2",
                                "flycatcher.hr.multi.farm.cos.3",
                                "ground_dweller_large.hn.multi.farm.cos.3",
                                "ground_dweller_large.hn.multi.forest.cos.2",
                                "ground_dweller_large.hn.multi.forest.cos.3",
                                "ground_dweller_large.hr.multi.farm.cos.3",
                                "ground_dweller_large.hr.multi.farm.poly.4",
                                "ground_dweller_large.hr.multi.forest.cos.2",
                                "ground_dweller_large.hr.multi.forest.poly.4",
                                "ground_feeding_passerines.hr.multi.forest.cos.3",
                                "kingfisher.hn.multi.farm.cos.2",
                                "kingfisher.hn.multi.farm.cos.3",
                                "kingfisher.hn.multi.farm.herm.4",
                                "kingfisher.hn.multi.farm.herm.6",
                                "kingfisher.hr.multi.farm.cos.2",
                                "kingfisher.hr.multi.farm.cos.3",
                                "kingfisher.hr.multi.farm.poly.4",
                                "large_wader.hn.multi.farm.cos.2",
                                "large_wader.hn.multi.farm.herm.4",
                                "large_wader.hn.multi.farm.herm.6",
                                "large_wader.hr.multi.farm.cos.2",
                                "large_wader.hr.multi.farm.cos.3",
                                "large_wader.hr.multi.farm.poly.4",
                                "myna.hn.multi.farm.cos.2",
                                "myna.hn.multi.farm.cos.3",
                                "myna.hr.multi.farm.cos.2",
                                "myna.hr.multi.farm.cos.3",
                                "myna.hr.multi.farm.poly.4",
                                "parakeet.hn.multi.farm.cos.3",
                                "parakeet.hr.multi.farm.cos.2",
                                "parakeet.hr.multi.farm.cos.3",
                                "pipit.hn.multi.farm.cos.3",
                                "pipit.hr.multi.farm.cos.2",
                                "roller_weaver.hn.multi.farm.cos.2",
                                "roller_weaver.hn.multi.farm.cos.3",
                                "roller_weaver.hr.multi.farm.cos.2",
                                "roller_weaver.hr.multi.farm.cos.3",
                                "roller_weaver.hr.multi.farm.poly.4",
                                "sunbirds.hn.multi.farm.cos.2",
                                "sunbirds.hr.multi.farm.cos.2",
                                "sunbirds.hr.multi.farm.cos.3",
                                "sunbirds.hr.multi.forest.cos.3",
                                "warblers.hr.multi.forest.cos.3",
                                "woodpecker.hr.multi.forest.no.adj",
                                "woodpecker.hr.multi.forest.poly.4")

kept.single.sp.models <- kept.single.sp.models[!(names(kept.single.sp.models) %in% bad.det.funs.single)]

kept.multi.sp.models <- kept.multi.sp.models[!(names(kept.multi.sp.models) %in% bad.det.funs.multi)]

kept.single.sp.models.h <- kept.single.sp.models.h[!(names(kept.single.sp.models.h) %in% bad.det.funs.single.habitat)]

kept.multi.sp.models.h <- kept.multi.sp.models.h[!(names(kept.multi.sp.models.h) %in% bad.det.funs.multi.habitat)]


##############pick the BEST forest/farmland specific detection function for each species ##################################
single.sp.no.model <- unique(single.species.data$name_latin[!(single.species.data$name_latin %in% sub("\\..*$", "", names(kept.single.sp.models)))])
single.sp.no.model

######Get FOREST specific detection functions
# Filter models with "forest" in the model name
kept.single.sp.models.forest <- kept.single.sp.models.h[grep("forest", names(kept.single.sp.models.h))]

separate_lists_forest <- list()
for (model_name in names(kept.single.sp.models.forest)) {
  # Extract the species name (model text before the first dot)
  group_name <- sub("\\..*$", "", model_name)
  # Check if the group exists, if not create a new list
  if (!(group_name %in% names(separate_lists_forest))) {
    separate_lists_forest[[group_name]] <- list()
  }
  # Append the model to the corresponding group list
  if (!is.null(kept.single.sp.models.forest[[model_name]])) {
    separate_lists_forest[[group_name]][[model_name]] <- kept.single.sp.models.forest[[model_name]]
  }
}
#######################get the best model for each species (in forests)
# Create a list to store the best models for each group
best_models_forest <- list()
# Iterate over each group in separate_lists
for (group_name in names(separate_lists_forest)) {
  # Get models in the current group
  models <- separate_lists_forest[[group_name]]
  # Calculate AICc for each model
  AICc_values <- sapply(models, AICc)
  # Find the model with the lowest AICc
  best_model_name <- names(models)[which.min(AICc_values)]
  # Extract the best model
  best_models_forest[[group_name]] <- models[[best_model_name]]
}

######Get farm specific detection functions
# Filter models with "farm" in the model name
kept.single.sp.models.farm <- kept.single.sp.models.h[grep("farm", names(kept.single.sp.models.h))]

separate_lists_farm <- list()
for (model_name in names(kept.single.sp.models.farm)) {
  # Extract the species name (model text before the first dot)
  group_name <- sub("\\..*$", "", model_name)
  # Check if the group exists, if not create a new list
  if (!(group_name %in% names(separate_lists_farm))) {
    separate_lists_farm[[group_name]] <- list()
  }
  # Append the model to the corresponding group list
  if (!is.null(kept.single.sp.models.farm[[model_name]])) {
    separate_lists_farm[[group_name]][[model_name]] <- kept.single.sp.models.farm[[model_name]]
  }
}
#######################get the best model for each species (in farms)
# Create a list to store the best models for each group
best_models_farm <- list()
# Iterate over each group in separate_lists
for (group_name in names(separate_lists_farm)) {
  # Get models in the current group
  models <- separate_lists_farm[[group_name]]
  # Calculate AICc for each model
  AICc_values <- sapply(models, AICc)
  # Find the model with the lowest AICc
  best_model_name <- names(models)[which.min(AICc_values)]
  # Extract the best model
  best_models_farm[[group_name]] <- models[[best_model_name]]
}

###look at the best detection functions to ensure they are ok
for (i in 1:length(best_models_forest)) {
  pdf(file=paste("SingleSpp/Habitats/BestModels/Forest/",names(best_models_forest[i]),".pdf",sep=""))
  plot(best_models_forest[[i]], showpoints=TRUE, main=paste(names(best_models_forest[i])))
  dev.off()
}

for (i in 1:length(best_models_farm)) {
  pdf(file=paste("SingleSpp/Habitats/BestModels/Farm/",names(best_models_farm[i]),".pdf",sep=""))
  plot(best_models_farm[[i]], showpoints=TRUE, main=paste(names(best_models_farm[i])))
  dev.off()
}


##############pick the BEST forest/farmland specific detection function for each species GROUP ##################################
######Get FOREST specific detection functions
# Filter models with "forest" in the model name
kept.multi.sp.models.forest <- kept.multi.sp.models.h[grep("forest", names(kept.multi.sp.models.h))]

separate_lists_multi_forest <- list()
for (model_name in names(kept.multi.sp.models.forest)) {
  # Extract the species name (model text before the first dot)
  group_name <- sub("\\..*$", "", model_name)
  # Check if the group exists, if not create a new list
  if (!(group_name %in% names(separate_lists_multi_forest))) {
    separate_lists_multi_forest[[group_name]] <- list()
  }
  # Append the model to the corresponding group list
  if (!is.null(kept.multi.sp.models.forest[[model_name]])) {
    separate_lists_multi_forest[[group_name]][[model_name]] <- kept.multi.sp.models.forest[[model_name]]
  }
}
#######################get the best model for each species (in forests)
# Create a list to store the best models for each group
best_models_multi_forest <- list()
# Iterate over each group in separate_lists
for (group_name in names(separate_lists_multi_forest)) {
  # Get models in the current group
  models <- separate_lists_multi_forest[[group_name]]
  # Calculate AICc for each model
  AICc_values <- sapply(models, AICc)
  # Find the model with the lowest AICc
  best_model_name <- names(models)[which.min(AICc_values)]
  # Extract the best model
  best_models_multi_forest[[group_name]] <- models[[best_model_name]]
}

############farmlands
kept.multi.sp.models.farm <- kept.multi.sp.models.h[grep("farm", names(kept.multi.sp.models.h))]

separate_lists_multi_farm <- list()
for (model_name in names(kept.multi.sp.models.farm)) {
  # Extract the species name (model text before the first dot)
  group_name <- sub("\\..*$", "", model_name)
  # Check if the group exists, if not create a new list
  if (!(group_name %in% names(separate_lists_multi_farm))) {
    separate_lists_multi_farm[[group_name]] <- list()
  }
  # Append the model to the corresponding group list
  if (!is.null(kept.multi.sp.models.farm[[model_name]])) {
    separate_lists_multi_farm[[group_name]][[model_name]] <- kept.multi.sp.models.farm[[model_name]]
  }
}
#######################get the best model for each species (in farms)
# Create a list to store the best models for each group
best_models_multi_farm <- list()
# Iterate over each group in separate_lists
for (group_name in names(separate_lists_multi_farm)) {
  # Get models in the current group
  models <- separate_lists_multi_farm[[group_name]]
  # Calculate AICc for each model
  AICc_values <- sapply(models, AICc)
  # Find the model with the lowest AICc
  best_model_name <- names(models)[which.min(AICc_values)]
  # Extract the best model
  best_models_multi_farm[[group_name]] <- models[[best_model_name]]
}

###look at the best detection functions to ensure they are ok
for (i in 1:length(best_models_multi_forest)) {
  pdf(file=paste("MultiSpp/Habitats/BestModels/Forest/",names(best_models_multi_forest[i]),".pdf",sep=""))
  plot(best_models_multi_forest[[i]], showpoints=TRUE, main=paste(names(best_models_multi_forest[i])))
  dev.off()
}

for (i in 1:length(best_models_multi_farm)) {
  pdf(file=paste("MultiSpp/Habitats/BestModels/Farm/",names(best_models_multi_farm[i]),".pdf",sep=""))
  plot(best_models_multi_farm[[i]], showpoints=TRUE, main=paste(names(best_models_multi_farm[i])))
  dev.off()
}
###all ok



##################global no covars & % closed habitat functions ################################################
######SINGLE SPECIES ###########################################
kept.single.sp.models.filtered <- kept.single.sp.models[!grepl("forest|farm", names(kept.single.sp.models))] ##remove the wrongly fitted farm & forest specific functions

separate_lists_global <- list()
for (model_name in names(kept.single.sp.models.filtered)) {
  # Extract the species name (model text before the first dot)
  group_name <- sub("\\..*$", "", model_name)
  # Check if the group exists, if not create a new list
  if (!(group_name %in% names(separate_lists_global))) {
    separate_lists_global[[group_name]] <- list()
  }
  # Append the model to the corresponding group list
  if (!is.null(kept.single.sp.models.filtered[[model_name]])) {
    separate_lists_global[[group_name]][[model_name]] <- kept.single.sp.models.filtered[[model_name]]
  }
}
#######################get the best model for each species (in farms)
# Create a list to store the best models for each group
best_models_global <- list()
# Iterate over each group in separate_lists
for (group_name in names(separate_lists_global)) {
  # Get models in the current group
  models <- separate_lists_global[[group_name]]
  # Calculate AICc for each model
  AICc_values <- sapply(models, AICc)
  # Find the model with the lowest AICc
  best_model_name <- names(models)[which.min(AICc_values)]
  # Extract the best model
  best_models_global[[group_name]] <- models[[best_model_name]]
}

###look at the best detection functions to ensure they are ok
for (i in 1:length(best_models_global)) {
  pdf(file=paste("DetectionFunctionCurves/SingleSpp/Global/",names(best_models_global[i]),".pdf",sep=""))
  plot(best_models_global[[i]], showpoints=TRUE, main=paste(names(best_models_global[i])))
  dev.off()
}

########MULTI SPECIES ###############################################################
kept.multi.sp.models.filtered <- kept.multi.sp.models[!grepl("forest|farm", names(kept.multi.sp.models))] ##remove the wrongly fitted farm & forest specific functions

separate_lists_global_multi <- list()
for (model_name in names(kept.multi.sp.models.filtered)) {
  # Extract the species name (model text before the first dot)
  group_name <- sub("\\..*$", "", model_name)
  # Check if the group exists, if not create a new list
  if (!(group_name %in% names(separate_lists_global_multi))) {
    separate_lists_global_multi[[group_name]] <- list()
  }
  # Append the model to the corresponding group list
  if (!is.null(kept.multi.sp.models.filtered[[model_name]])) {
    separate_lists_global_multi[[group_name]][[model_name]] <- kept.multi.sp.models.filtered[[model_name]]
  }
}
#######################get the best model for each species (in farms)
# Create a list to store the best models for each group
best_models_global_multi <- list()
# Iterate over each group in separate_lists
for (group_name in names(separate_lists_global_multi)) {
  # Get models in the current group
  models <- separate_lists_global_multi[[group_name]]
  # Calculate AICc for each model
  AICc_values <- sapply(models, AICc)
  # Find the model with the lowest AICc
  best_model_name <- names(models)[which.min(AICc_values)]
  # Extract the best model
  best_models_global_multi[[group_name]] <- models[[best_model_name]]
}

##############################################################################################
###look at the best detection functions to ensure they are ok
for (i in 1:length(best_models_global_multi)) {
  pdf(file=paste("DetectionFunctionCurves/MultiSpp/Global/",names(best_models_global_multi[i]),".pdf",sep=""))
  plot(best_models_global_multi[[i]], showpoints=TRUE, main=paste(names(best_models_global_multi[i])))
  dev.off()
}

#######save the best functions######
save(best_models_global_multi, best_models_global, best_models_multi_farm, best_models_multi_forest,
     best_models_forest, best_models_farm, multi.species.data, single.species.data, sample.table, 
     region.table, distance.data, distance.data.orig, bird.counts6, file="DetectionFunctionCurves/BestDetectionFunctions.rda")


############## 3. CALCULATE EFFECTIVE AREA SURVEYED ###########################

###import the data 
load("DataForDetectionFunctionFittingFinal.rda")
distance.data.use <- distance.data
single.species.data.farm.use <- single.species.data.farm
single.species.data.use <- single.species.data
single.species.data.forest.use <- single.species.data.forest
multi.species.data.use <- multi.species.data
multi.species.data.farm.use <- multi.species.data.farm
multi.species.data.forest.use <- multi.species.data.forest

####import all the best fitting detection functions
load("BestDetectionFunctions.rda")
unique(multi.species.data.use$name_latin) ##species name
unique(multi.species.data.use$detect.group) ##species group
multi.species.data.use$point_ID <- multi.species.data.use$Region.Label
multi.species.data.use$obs_code <- multi.species.data.use$object
# Create a template dataframe with all combinations of point_ID, name_latin, and point_count_repeat_number
template <- expand.grid(point_ID = unique(multi.species.data.use$point_ID),
                        name_latin = unique(multi.species.data.use$name_latin),
                        point_count_repeat_number = 1:4) ##4 visists to each point
# Merge template with bird.counts6 to fill in missing combinations
filled_data <- merge(template, multi.species.data.use, by = c("point_ID", "name_latin", "point_count_repeat_number"), all.x = TRUE)
# Fill missing number_individuals with 0
filled_data$size[is.na(filled_data$size)] <- 0
####information that needs to copied from the same square (points are nested in squares)
filled_data3 <- filled_data %>%
  dplyr::group_by(point_ID) %>%
  dplyr::mutate(
    site_type = na.omit(site_type)[1],
    square_ID = na.omit(square_ID)[1])
####information that needs to be copied for each species       
filled_data4 <- filled_data3 %>%
  dplyr::group_by(name_latin) %>%
  dplyr::mutate(
    detect.group = na.omit(detect.group)[1],
    no_records_total = na.omit(no_records_total)[1],
    no_records_forest = na.omit(no_records_forest)[1],
    no_records_farm = na.omit(no_records_farm)[1],
    obs_40_forest = na.omit(obs_40_forest)[1],
    obs_40_farm = na.omit(obs_40_farm)[1],
    Own = na.omit(Own)[1],
    Group = na.omit(Group)[1])  
####modelling will be conducted on the number of clusters rather than number of individuals
##so need number cluster sights and mean cluster size (per spp)
filled_data5 <- filled_data4 %>%
  dplyr::mutate(cluster_count = ifelse(size >= 1, 1, 0)) %>%
  dplyr::group_by(name_latin) %>%
  dplyr::mutate(mean_cluster_size = mean(size[size > 0], na.rm = TRUE))
####get total number of clusters per point (ie summed across visits)
data.sum <- filled_data5 %>%
  dplyr::group_by(point_ID, name_latin) %>%
  dplyr::mutate(TotalClustersPoint = sum(cluster_count))
# Object needs to be a unique ID
data.sum$object2 <- c(1:length(data.sum$distance))
all.data <- data.sum %>% dplyr::mutate(Effort = 4)
all.data <- all.data %>%
  dplyr::mutate(site_class = case_when(
    site_type == "forest" ~ "forest",
    site_type %in% c("ZBNF", "chemical") ~ "farm",
  ))

##get proportion of closed habitat data around each point 
prop.closed <- read.csv("points_prop_closed.csv")
prop.closed2 <- prop.closed %>% select(point_ID, prop_closed)

###get the proportion of vegetation patches in square
prop.habitat <- read.csv("wfh_squares_final.csv")
prop.habitat2 <- prop.habitat %>% select(site_ID, perc_quality)
prop.habitat3 <- dplyr::rename(prop.habitat2, square_ID = site_ID)
all.data2 <- plyr::join_all(list(all.data, prop.closed2),
                            by = c("point_ID"), type = 'full')
all.data2$prop_closed[is.na(all.data2$prop_closed)] <- 1 ##forst sites are 100% closed
all.data3 <- plyr::join_all(list(all.data2, prop.habitat3),
                            by = c("square_ID"), type = 'full')
all.data3 <- all.data3 %>%
  dplyr::filter(!is.na(name_latin))
unique(all.data3$name_latin)
unique(all.data3$detect.group)
all.data3$name_latin <- as.character(all.data3$name_latin)
all.data3$detect.group <- as.character(all.data3$detect.group)
all.data3.1 <- all.data3 %>% select(- p_closed)
all.data3.2 <- all.data3.1 %>% dplyr::rename(p_closed = prop_closed)


##3. CALCULATE THE EFFECTIVE AREA SURVEYD
##estimate species specific effective radii at each point
###test simple code first
effAreafn3 <- function(detfnobj, newdata, areaconv, truncation) {
  newdata$Pa <- NA
  k <- 0
  for (i in 1:dim(newdata)[1]) {
    if (newdata[i, "distance"] <= truncation & !is.na(newdata[i, "distance"])) {
      k <- k + 1
      newdata$Pa[i] <- detfnobj$fitted[k]
    }
  }
  newdata$effArea <- ifelse(is.na(newdata$Pa), NA, newdata$Pa * truncation^2 * pi / areaconv)
  result <- newdata[!duplicated(newdata$point_ID), ]
  result$effArea <- ifelse(result$cluster_count == 0,
                           predict(detfnobj, 
                                   newdata = data.frame(site_type = result$site_type),
                                   esw = TRUE)$fitted / areaconv,
                           result$effArea)
  return(result)
}

data_1spp <- all.data3.2 %>% dplyr::filter(name_latin == "Alexandrinus_krameri") ##test on one species
mody <- best_models_farm[["Alexandrinus_krameri"]]
testy2 <- effAreafn3(mody, data_1spp, 10000, 100)
testy3 <- effAreafn3(best_models_farm[["Alexandrinus_krameri"]], data_1spp, 10000, 100) ##works

##################check can loop through all species when using the same detection function for all
areaconv <- 10000
truncation <- 100
detfnobj <- best_models_farm[["Alexandrinus_krameri"]] ##use this to test

# Create an empty data frame to store the results
combined_results <- data.frame()
# Get unique values of all.data3$name_latin
unique_names <- unique(all.data3$name_latin)
# Iterate over each unique name
for (name in unique_names) {
  # Subset data based on the current unique name
  subset_data <- all.data3[all.data3$name_latin == name, ]
  # Apply effAreafn3 function to the subset
  subset_result <- effAreafn3(detfnobj, subset_data, areaconv, truncation)
  # Combine the results with the previous results
  combined_results <- rbind(combined_results, subset_result)
}

all.data4 <- all.data3.2 %>% dplyr::mutate(name_latin2 = as.character(name_latin)) 
all.data4$name_latin2 <- as.character(all.data4$name_latin2)
unique_names <- unique(all.data4$name_latin2)
unique_names

# Check if all names exist in best_models_farm
all_exist <- all(data_1spp$name_latin2 %in% names(best_models_farm))
print(all_exist)

all.data5 <- all.data4 %>%
  dplyr::filter(!is.na(name_latin2))

unique(all.data5$name_latin2)

name_latin2_exist <- sapply(all.data5$name_latin2, function(name) name %in% names(best_models_global)) ##single species global models
print(name_latin2_exist) ##works

#######see if can loop for single species models
funct_chosen <- lapply(all.data5$name_latin2, function(name) {
  if (name %in% names(best_models_farm)) {
    return(best_models_farm[[name]])
  } else {
    return(NULL)
  }
})
funct_chosen2 <- list()
for (name in all.data5$name_latin2) {
  if (name %in% names(best_models_farm)) {
    funct_chosen2[[name]] <- list(name = name, model = best_models_farm[[name]])
  }
}
funct_chosen3 <- list()
for (name in all.data5$name_latin2) {
  if (name %in% names(best_models_farm)) {
    funct_chosen3[[name]] <- list(name = name, model = best_models_farm[[name]], source = "best_models_farm")
  } else if (name %in% names(best_models_global)) {
    funct_chosen3[[name]] <- list(name = name, model = best_models_global[[name]], source = "best_models_global")
  }
} ###ok - all 52 single species have a model, model specified in 'source' of the large list

#######MULTI SPECIES GROUPS TESTS
###see if it works when selecting on the name of the detect.group
functm_chosen <- lapply(all.data5$detect.group, function(name) {
  if (name %in% names(best_models_multi_farm)) {
    return(best_models_multi_farm[[name]])
  } else {
    return(NULL)
  }
})
functm_chosen2 <- list()
for (name in all.data5$detect.group) {
  if (name %in% names(best_models_multi_farm)) {
    functm_chosen2[[name]] <- list(name = name, model = best_models_multi_farm[[name]])
  }
}
functm_chosen3 <- list()
for (name in all.data5$detect.group) {
  if (name %in% names(best_models_multi_farm)) {
    functm_chosen3[[name]] <- list(name = name, model = best_models_multi_farm[[name]], source = "best_models_multi_farm")
  } else if (name %in% names(best_models_global_multi)) {
    functm_chosen3[[name]] <- list(name = name, model = best_models_global_multi[[name]], source = "best_models_multi_global")
  }
} ##works fine

###code to pick the right detection function model for each species-site class combo
functm_chosen5 <- list()
unique_names <- unique(all.data5$name_latin2)
for (name in unique_names) {
  detect_group <- unique(all.data5$detect.group[all.data5$name_latin2 == name])
  own_status <- unique(all.data5$Own[all.data5$name_latin2 == name])
  site_class_farm <- all.data5$site_class[all.data5$name_latin2 == name & all.data5$site_class == "farm"]
  site_class_forest <- all.data5$site_class[all.data5$name_latin2 == name & all.data5$site_class == "forest"]
  # For farm site class
  if (length(site_class_farm) > 0) {
    if (own_status != "none") {
      # Check if the name of name_latin2 equals a model name in best_models_farm
      if (name %in% names(best_models_farm)) {
        functm_chosen5[[length(functm_chosen5) + 1]] <- list(name = name, 
                                                             model = best_models_farm[[name]], 
                                                             source = "best_models_farm", 
                                                             site_class = "farm")
      }
      # If not, use model from best_models_global
      else if (name %in% names(best_models_global)) {
        functm_chosen5[[length(functm_chosen5) + 1]] <- list(name = name, 
                                                             model = best_models_global[[name]], 
                                                             source = "best_models_global", 
                                                             site_class = "farm")
      }
    } else {
      # Check if the detect_group is available in best_models_multi_farm
      if (detect_group %in% names(best_models_multi_farm)) {
        functm_chosen5[[length(functm_chosen5) + 1]] <- list(name = name, 
                                                             model = best_models_multi_farm[[detect_group]], 
                                                             source = "best_models_multi_farm", 
                                                             site_class = "farm")
      } 
      # If not, use model from best_models_global_multi
      else if (detect_group %in% names(best_models_global_multi)) {
        functm_chosen5[[length(functm_chosen5) + 1]] <- list(name = name, 
                                                             model = best_models_global_multi[[detect_group]], 
                                                             source = "best_models_global_multi", 
                                                             site_class = "farm")
      }
    }
  }
  # For forest site class
  if (length(site_class_forest) > 0) {
    if (own_status != "none") {
      # Check if the name of name_latin2 equals a model name in best_models_forest
      if (name %in% names(best_models_forest)) {
        functm_chosen5[[length(functm_chosen5) + 1]] <- list(name = name, 
                                                             model = best_models_forest[[name]], 
                                                             source = "best_models_forest", 
                                                             site_class = "forest")
      }
      # If not, use model from best_models_global
      else if (name %in% names(best_models_global)) {
        functm_chosen5[[length(functm_chosen5) + 1]] <- list(name = name, 
                                                             model = best_models_global[[name]], 
                                                             source = "best_models_global", 
                                                             site_class = "forest")
      }
    } else {
      # Check if the detect_group is available in best_models_multi_forest
      if (detect_group %in% names(best_models_multi_forest)) {
        functm_chosen5[[length(functm_chosen5) + 1]] <- list(name = name, 
                                                             model = best_models_multi_forest[[detect_group]], 
                                                             source = "best_models_multi_forest", 
                                                             site_class = "forest")
      } 
      # If not, use model from best_models_global_multi
      else if (detect_group %in% names(best_models_global_multi)) {
        functm_chosen5[[length(functm_chosen5) + 1]] <- list(name = name, 
                                                             model = best_models_global_multi[[detect_group]], 
                                                             source = "best_models_global_multi", 
                                                             site_class = "forest")
      }
    }
  }
}

rename_list_elements <- function(lst) {
  for (i in seq_along(lst)) {
    name <- paste(lst[[i]][["name"]], lst[[i]][["site_class"]], sep = ".")
    names(lst)[i] <- name
  }
  return(lst)
}
functm_chosen6 <- rename_list_elements(functm_chosen5)
######works ###########


####CALCULATE EFFECTIVE AREA#####################
###Trial:
##get data from fitted ds model
af.farm.data <- (functm_chosen6[["Acridotheres_fuscus.farm"]][["model"]][["data"]])
af.farm.data
##get ds object
af.farm.mod <- functm_chosen6[["Acridotheres_fuscus.farm"]][["model"]]
af.farm.mod
class(af.farm.mod)
predict(af.farm.mod, newdata = NULL, esw = TRUE) ###works

# Create a template dataframe with all combinations of point_ID, name_latin, and point_count_repeat_number
template2 <- expand.grid(Region.Label = unique(af.farm.data$Region.Label),
                         name_latin = unique(af.farm.data$name_latin),
                         point_count_repeat_number = unique(af.farm.data$point_count_repeat_number))
# Merge template with bird.counts6 to fill in missing combinations
filled_data <- merge(template2, af.farm.data, by = c("Region.Label", "name_latin", "point_count_repeat_number"), all.x = TRUE)
# Fill missing number_individuals with 0
filled_data$size[is.na(filled_data$size)] <- 0
####information that needs to copied from the same square (points are nested in squares)
filled_data3 <- filled_data %>%
  dplyr::group_by(Region.Label) %>%
  dplyr::mutate(
    site_type = na.omit(site_type)[1],
    square_ID = na.omit(square_ID)[1])
####information that needs to be copied for each species       
filled_data4 <- filled_data3 %>%
  dplyr::group_by(name_latin) %>%
  dplyr::mutate(
    detect.group = na.omit(detect.group)[1],
    no_records_total = na.omit(no_records_total)[1],
    no_records_forest = na.omit(no_records_forest)[1],
    no_records_farm = na.omit(no_records_farm)[1],
    obs_40_forest = na.omit(obs_40_forest)[1],
    obs_40_farm = na.omit(obs_40_farm)[1],
    Own = na.omit(Own)[1],
    Group = na.omit(Group)[1])  

predict(af.farm.mod, newdata = NULL, esw = TRUE)

filled_data4$effArea <- predict(af.farm.mod, 
                                newdata=data.frame(name_latin=filled_data4$name_latin),
                                esw = TRUE)$fitted / 10000 


##tets
ardea.farm.data <- functm_chosen6[["Ardea_alba.farm"]][["model"]][["data"]] 
### need way of cycling through functm6_chosen -- i in length Name ?
ardea.farm.mod <- functm_chosen6[["Ardea_alba.farm"]][["model"]]

# Function to get Pa for a given fitted value
spf <- function(sphab, Code) {
  pa <- functm_chosen6[[sphab]][["model"]][["fitted"]][[as.character(Code)]] 
  print(pa)
}

spf("Ardea_alba.farm", 29) ###test

dffunct <- function(df){
  df$pa <- rep(NA,length(df$object))
  
  for (i in 1:length(df$object)) {
    
    sphab <- "Ardea_alba.farm" ##change to model name ## i in length functm6_chosen?
    
    #sphab <- paste0(df$name_latin[i],".",df$site_class[i])
    
    Code <- df$object[i]
    
    x <- spf(sphab=sphab,Code = Code)
    df$pa[i] <- x
    
    y <- x * 100^2 * pi / 10000 
    df$effArea[i] <- y
    
    
  }
  return(df)
}

new.ardea <- dffunct(ardea.farm.data) 
###yes - corresponds to fitted Pa values


#############################
# Ensure functm_chosen6 is populated correctly
print(functm_chosen6)

# Check the structure of the first element in functm_chosen6
print(str(functm_chosen6[[1]]))

# Define the rspdata function
rspdata <- function(sphab) {
  ds.data <- functm_chosen6[[sphab]][["model"]][["data"]]
  dataset_name <- paste0(sphab)
  return(ds.data)  # Return the dataset
}

# Test the rspdata function with a specific model
y <- rspdata("Amaurornis_phoenicurus.farm")
print(y)  # Print the dataset

# Create an empty list to store datasets
all_datasets <- list()

# Iterate over each model (Name) in functm_chosen6
for (sphab in names(functm_chosen6)) {
  ds_data <- rspdata(sphab)
  dataset_name <- paste0(sphab)
  all_datasets[[dataset_name]] <- ds_data
}

# Check the list of datasets
print(all_datasets)

###works - all_datasets contains the data that was used to fit the detection functions


####GET PAs & EFFECTIVE AREAS (fitted values)
spf <- function(sphab, Code) {
  pa <- functm_chosen6[[sphab]][["model"]][["fitted"]][[as.character(Code)]] 
  return(pa)
}

dffunct <- function(all_datasets) {
  for (i in names(all_datasets)) {
    df <- all_datasets[[i]]
    df$pa <- rep(NA, length(df$object))
    df$effArea <- rep(NA, length(df$object))
    
    for (j in 1:length(df$object)) {
      sphab <- i
      Code <- df$object[j]
      
      x <- spf(sphab = sphab, Code = Code)
      df$pa[j] <- x
      
      y <- x * 100^2 * pi / 10000 ### divided by 10k to get per ha density estimates
      df$effArea[j] <- y
    }
    all_datasets[[i]] <- df
  }
  return(all_datasets)
}

all_effareas <- dffunct(all_datasets) 
###all works

filter_dataframe <- function(df_name, df) {
  dataframe_name_clean <- gsub("^.*\\[\\[\"(.*)\"\\]\\]$", "\\1", df_name)
  
  species <- gsub("\\..*", "", dataframe_name_clean) # Extract species name from dataframe name
  habitat <- gsub(".*\\.", "", dataframe_name_clean) # Extract habitat type from dataframe name
  
  filtered_df <- df[df$name_latin == species & df$site_class == habitat, ]
  
  return(filtered_df)
}

filtered_all_effareas <- list()
for (name in names(all_effareas_with_siteclass)) {
  filtered_all_effareas[[name]] <- filter_dataframe(name, all_effareas_with_siteclass[[name]])
}

# Combine all data frames row-wise by matching column names
all_obs_fitted <- bind_rows(filtered_all_effareas, .id = "df_name")

###check that no observations were lost
missing_names <- setdiff(unique(bird.counts6$name_latin), unique(all_obs_fitted$name_latin))
print(missing_names) 
unique(all_obs_fitted$name_latin) ##199 species - all fine

###ensure that there is a row for each spp at each point (ie create rows for no detections)
template <- expand.grid(point_ID = unique(multi.species.data.use$point_ID),
                        name_latin = unique(multi.species.data.use$name_latin),
                        point_count_repeat_number = 1:4) ##4 visists to each point

all_obs_fitted2 <- all_obs_fitted %>%
  dplyr::rename(point_ID = Region.Label)

# Merge template with bird.counts6 to fill in missing combinations
filled_data <- merge(template, all_obs_fitted2, by = c("point_ID", "name_latin", "point_count_repeat_number"), all.x = TRUE)
# Fill missing number_individuals with 0
filled_data$size[is.na(filled_data$size)] <- 0
####information that needs to copied from the same square (points are nested in squares)
filled_data3 <- filled_data %>%
  dplyr::group_by(point_ID) %>%
  dplyr::mutate(
    site_type = na.omit(site_type)[1],
    square_ID = na.omit(square_ID)[1])
####information that needs to be copied for each species       
filled_data4 <- filled_data3 %>%
  dplyr::group_by(name_latin) %>%
  dplyr::mutate(
    detect.group = na.omit(detect.group)[1],
    no_records_total = na.omit(no_records_total)[1],
    no_records_forest = na.omit(no_records_forest)[1],
    no_records_farm = na.omit(no_records_farm)[1],
    obs_40_forest = na.omit(obs_40_forest)[1],
    obs_40_farm = na.omit(obs_40_farm)[1],
    Own = na.omit(Own)[1],
    Group = na.omit(Group)[1])  


filled_data5 <- filled_data4 %>%
  dplyr::mutate(cluster_count = ifelse(size >= 1, 1, 0)) %>%
  dplyr::group_by(name_latin) %>%
  dplyr::mutate(mean_cluster_size = mean(size[size > 0], na.rm = TRUE))

data.sum <- filled_data5 %>%
  dplyr::group_by(point_ID, point_count_repeat_number, name_latin) %>%
  dplyr::mutate(TotalClustersVisit = sum(cluster_count))

data.sum2 <- data.sum %>%
  dplyr::group_by(point_ID, point_count_repeat_number, name_latin) %>%
  dplyr::mutate(ClusterDensity = TotalClustersVisit/effArea)

data.sum3 <- data.sum2 %>%
  dplyr::group_by(point_ID, point_count_repeat_number, name_latin) %>%
  dplyr::mutate(PopDensity = ClusterDensity*mean_cluster_size)

unique(data.sum2$name_latin)
unique(data.sum2$point_ID) ##all fine

save(data.sum3, file = "density_data_fitted_only.csv")
