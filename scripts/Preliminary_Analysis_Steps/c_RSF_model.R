#---
#"Resource Selection Function for wolves in Washington State"
#"Lisanne Petracca"
#"November 2023"
#---

library(glmmTMB)
library(corrplot)
library(dotwhisker)
library(AICcmodavg)
library(data.table)
library(broom)
library(lme4)
library(GGally)
library(raster)
library(broom.mixed)
library(here)
library(data.table)
library(tidyverse)
library(nimble)

here()

######------ RUNNING AND PLOTTING RSF ------######

#testing for correlations among input
data <- as.data.frame(fread("C:/My_Directory/RSF_points_final.csv"))

#get rid of any cells that have NAs
data <- data[complete.cases(data), ]

#create truncated version of focal population density
hist(data$popdens_focal)
max(data$popdens_focal)
quantile(data$popdens_focal, probs = seq(0, 1, 0.01))

data$popdens_focal_truncated <- data$popdens_focal
data$popdens_focal_truncated[data$popdens_focal_truncated>0.09065171] <- 0.09065171

#logging distance to highway
data$logdisthwy <- log(data$disthwy+1)

#gather continuous variables and scale them
continuous_var <- data %>% dplyr::select(agriculture_focal, alldeer_final, developed_focal, logdisthwy, 
                                        elev_focal, forest_focal, grassland_focal, popdens_focal_truncated,
                                        roaddensity_final, rugged_focal_truncated, shrubland_focal)

scaled <- as.data.frame(scale(continuous_var))

#use corrplot to get correlations of continuous variables
M <- cor(scaled)
corrplot(M, method="number", type = "upper")

#developed_focal and popdens_focal cannot go together (r=0.66)

#get pack info
data$pack <- gsub("\\_.*","", data$packyear )
length(unique(data$pack))

#get continuous and dummy covariates together
covs <-cbind(scaled, data$allotment_01, data$elk_01)
head(covs)
covs <- covs %>% rename(allotment_01='data$allotment_01', elk_01='data$elk_01')

#get dataset in correct format for glmmTMB
data_glmmTMB <- cbind(data[,c("used", "pack")], covs)
head(data_glmmTMB)

#remove developed_focal from data_glmmTMB bc of correlation w popdens_focal
data_glmmTMB <- data_glmmTMB %>% dplyr::select(-c(developed_focal))
head(data_glmmTMB)

######---- LETS DO RSF IN glmmTMB ----#####

#using script based on Muff et al. (2020)

#assigns weights of 1000 to avail, 1 to used
data_glmmTMB$weight <- 1000^(1 - data_glmmTMB$used)

#'we recommend to manually fix the variance of the random intercept at a large value. 
#'This can be done in glmmTMB() by first setting up the model, but do not yet fit it:
wolf.tmp <- glmmTMB(used ~ agriculture_focal + allotment_01 + alldeer_final + logdisthwy + elev_focal + I(elev_focal^2) + 
                      forest_focal + grassland_focal + popdens_focal_truncated + roaddensity_final + 
                      rugged_focal_truncated + I(rugged_focal_truncated^2) + shrubland_focal + 
                      #allotment_01 
                      #elk_01 + 
                      (1|pack), family=binomial(), data = data_glmmTMB,
                      doFit=FALSE, weights = weight)

#' Then fix the standard deviation of the first random term, which is the `(1|pack)` component  
#' in the above model equation. We use $\sigma=10^3$, which corresponds to a variance of $10^6$:
#+ echo=TRUE, message=FALSE,cache=TRUE
wolf.tmp$parameters$theta[1] <- log(1e3)

#' We need to tell `glmmTMB` not to change the first entry of the vector of variances, and give all other variances another indicator to make sure they can be freely estimated:
#+ echo=TRUE, message=FALSE,cache=TRUE
#here we only have one random part, so unsure if this is doing anything
wolf.tmp$mapArg <- list(theta=factor(c(NA)))

#' Then fit the model and look at the results:
#+ echo=TRUE, message=FALSE, cache=TRUE 
wolf.rsf <- glmmTMB:::fitTMB(wolf.tmp)
summary(wolf.rsf)

saveRDS(wolf.rsf, "C:/My_Directory/wolf_rsf_wallotment_wtruncatedpopdens.rds")

wolf.rsf <- readRDS("C:/My_Directory/wolf_rsf_wallotment_wtruncatedpopdens.rds")
summary(wolf.rsf)

######---- LETS PLOT COEFFICIENTS ----#####

library(glmmTMB)
library(nloptr)

#plotting for glmmTMB
library(broom.mixed)
coef <- tidy(wolf.rsf, conf.int = TRUE)
coef <- coef[-1,]
coef$term <- c("agricultural cover", "allotment presence", "relative deer abundance",
               "distance to state highway", "elevation", "elevation^2", "forest cover", 
               "grassland cover", "human population density", "road density", "terrain ruggedness",
               "terrain ruggedness^2", "shrubland cover")
coef$term <- factor(coef$term, levels = c("allotment presence", "distance to state highway", "relative deer abundance", "elevation^2", 
                                          "elevation",
                                          "shrubland cover", "forest cover", "grassland cover",
                                          "road density", "terrain ruggedness^2", "terrain ruggedness", 
                                          "agricultural cover", "human population density"))
               
ggplot(coef, aes(estimate, term)) +
  geom_point(size=2) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),height=0.5) +
  geom_vline(xintercept=0, lty=2)
ggsave("C:/My_Directory/Chapter1_AppendixS4_FigureS2.eps")


#####----- LET'S SCALE ORIGINAL RASTERS -----#####

#read in rasters
rastlist <- list.files(path = "C:/My_Directory/a_Final_Rasters", pattern='.tif$', all.files=TRUE, full.names=TRUE)
#exclude the 0/1 land cover class rasters bc not informative at this scale
c <- c(2,3,4,6,7,9,10,12,14,18,19,22,24)
rastlist <- rastlist[c]

#get rid of developed_focal and elk
rastlist <- rastlist[-c(4,7)]

#select those that need scaling (i.e., remove allotment and disthwy)
rastlist_forstd <- rastlist[-c(3,4)]

#select those that don't need scaling
rastlist_notforstd <- rastlist[c(3)]

#disthwy as its own
disthwy <- rastlist[c(4)]

#create stack
stack_forstd <- stack(rastlist_forstd)

#getting mean and sd of data columns, and using those columns to standardize rasters
mean_and_sd <- function(x, column, rast) {
  mean <- mean(x[,column])
  sd <- sd(x[,column])
  rast_std <- (rast-mean) / sd
  return(rast_std)
}

head(data)

#mean(data[,"elev_focal"])
agric_focal_std <- mean_and_sd(data, "agriculture_focal", stack_forstd[[1]])
alldeer_focal_std <- mean_and_sd(data, "alldeer_final", stack_forstd[[2]])
elev_focal_std <-mean_and_sd(data, "elev_focal", stack_forstd[[3]])
forest_focal_std <-mean_and_sd(data, "forest_focal", stack_forstd[[4]])
grassland_focal_std <-mean_and_sd(data, "grassland_focal", stack_forstd[[5]])
popdens_focal_std <-mean_and_sd(data, "popdens_focal_truncated", stack_forstd[[6]])
roaddens_std <-mean_and_sd(data, "roaddensity_final", stack_forstd[[7]])
rugged_focal_std <-mean_and_sd(data, "rugged_focal_truncated", stack_forstd[[8]])
shrubland_focal_std <-mean_and_sd(data, "shrubland_focal", stack_forstd[[9]])

#doing disthwy somewhat differently
disthwy_log <- log(raster(disthwy) + 1)

#stacking them and giving them names
stack_all_forpred <- stack(agric_focal_std, alldeer_focal_std, elev_focal_std, forest_focal_std,
                           grassland_focal_std, popdens_focal_std, roaddens_std, rugged_focal_std, shrubland_focal_std,
                           disthwy_log, raster(rastlist_notforstd))
names(stack_all_forpred) <- c("agriculture_focal", "alldeer_focal", "elev_focal", "forest_focal", 
                              "grassland_focal", "popdens_focal", "roaddens", "rugged_focal", "shrubland_focal",
                              "disthwy", "allotment_01")
writeRaster(stack_all_forpred, filename=file.path("C:/My_Directory/a_Final_Rasters_Scaled", names(stack_all_forpred)), bylayer=TRUE, format='GTiff', overwrite=TRUE)

######------ READING IN SCALED RASTERS AND PREDICTING RSF ------######

wolf.rsf <- readRDS("C:/My_Directory/wolf_rsf_wallotment_wtruncatedpopdens.rds")
summary(wolf.rsf)
rastlist <- list.files(path = "C:/My_Directory/a_Final_Rasters_Scaled", pattern='.tif$', all.files=TRUE, full.names=TRUE)

#create stack
stack_std <- stack(rastlist)
names(stack_std)

#manually predicting model elk(0/1) removed and with TRUNCATED ruggedness
predict <- exp((-4.866e-01 * stack_std[[1]]) + (3.198e-01  * stack_std[[2]]) + (5.169e-01 * stack_std[[3]]) + (3.904e-01 * stack_std[[4]])+
                 (2.686e-01 * stack_std[[5]]) + (-4.741e-01 * (stack_std[[5]]^2)) + (9.496e-02 * stack_std[[6]])+ (-9.990e-02 * stack_std[[7]]) + 
                 (-1.157 * stack_std[[8]]) + (-1.238e-01 * stack_std[[9]]) + (-4.594e-01 * stack_std[[10]]) + (5.077e-02 * (stack_std[[10]]^2)) +
                 (1.119e-01 * stack_std[[11]]))
library(raster)
plot(predict)
#basically just dividing by max bc min is so small
max <- cellStats(predict, max) 

quantile(predict)

predict5 <- predict/max

#plot and write raster
plot(predict5)
writeRaster(predict5, "C:/My_Directory/RSF_noelk_trunc_rugged_newdeerinterp_newsample2021_trunpopdens.tif",overwrite=TRUE)

#making inverse version for resistance surface
orig_agg <- aggregate(predict5,15)
invert <- 1-orig_agg
plot(invert)
writeRaster(invert, "C:/My_Directory/RSF_aggregated_invert.tif",overwrite=TRUE)