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
#devtools::install_github("bbolker/broom")
library(nimble)

setwd(here())

######------ RUNNING AND PLOTTING RSF ------######

#testing for correlations among input

data <- as.data.frame(fread("Data/RSF/RSF_points_final.csv"))
dim(data) #start w 376047

#get rid of any cells that have NAs
data <- data[complete.cases(data), ]
dim(data) #now w 375993 pts (loss of 54 points, likely due to extent of WA and extent of some rasters being diff)
head(data)

hist(data$popdens_focal)
max(data$popdens_focal)
quantile(data$popdens_focal, probs = seq(0, 1, 0.01))

#this is at the 
data$popdens_focal_truncated <- data$popdens_focal
data$popdens_focal_truncated[data$popdens_focal_truncated>0.09065171] <- 0.09065171
#quantile(data$popdens_focal_truncated, probs = seq(0, 1, 0.01))

data$logdisthwy <- log(data$disthwy+1)

continuous_var <- data %>% dplyr::select(agriculture_focal, alldeer_final, developed_focal, logdisthwy, 
                                        elev_focal, forest_focal, grassland_focal, popdens_focal_truncated,
                                        roaddensity_final, rugged_focal_truncated, shrubland_focal)
dim(continuous_var)
head(continuous_var)
#seeing what the maxes are in the data
scaled <- as.data.frame(scale(continuous_var))
apply(scaled, 2, max) #ok so new popdens_focal_truncated reduces max from 64.278 to 7.88
apply(scaled, 2, mean)

#corrplot
M <- cor(scaled)
corrplot(M, method="number", type = "upper")

#and developed_focal and popdens_focal cannot go together (0.66)

#ok so let's test ruggedness relationship post facto
min(scaled$rugged_focal_truncated)
max(scaled$rugged_focal_truncated)
hist(scaled$rugged_focal_truncated)
rugged <- seq(from = -1.179554, to = 4.287925, length=1000)
output <- expit(-0.46*rugged + 0.05*(rugged^2))
plot(output) #selection is at low ruggedness, not at extremes

#get pack info
data$pack <- gsub("\\_.*","", data$packyear )
length(unique(data$pack))

covs <-cbind(scaled, data$allotment_01, data$elk_01)
head(covs)
covs <- covs %>% rename(allotment_01='data$allotment_01', elk_01='data$elk_01')

#get dataset for glmmTMB
data_glmmTMB <- cbind(data[,c("used", "pack")], covs)
head(data_glmmTMB)

#jags stuff
# pack <- as.numeric(as.factor(data$pack))
# y <- data$used

#which covariate to use considering they are correlated
#keep popdens_focal then; this took forever btw
# K      AIC Delta_AIC AICWt        LL
# Mod1 3 139826.3      0.00     1 -69910.16
# Mod2 3 143436.5   3610.13     0 -71715.23
# competing_mod <- list()
# competing_mod[[1]] <-  glmmTMB(used ~ scale(popdens_focal)  + (1|pack), family = "binomial", data = data)
# competing_mod[[2]] <- glmmTMB(used ~ scale(developed_focal) + (1|pack), family = "binomial", data = data)
# aictab(competing_mod, sort=F, second.ord=F)

#remove developed_focal from data_glmmTMB
data_glmmTMB <- data_glmmTMB %>% dplyr::select(-c(developed_focal))
head(data_glmmTMB)

######---- LETS TRY GLMMTMB LIKE MUFF ET AL 2020 ----#####
#using script from here based on muff et al. (2020)
#https://conservancy.umn.edu/bitstream/handle/11299/204737/fisher_rsf_and_ssf.R?sequence=26&isAllowed=y

head(data_glmmTMB)
dim(data_glmmTMB) #should be 375993

#assigns weights of 1000 to avail, 1 to used
data_glmmTMB$weight <- 1000^(1 - data_glmmTMB$used)
data_glmmTMB$weight2 <- 100^(1 - data_glmmTMB$used)

head(data_glmmTMB)
tail(data_glmmTMB)

#' As explained in the manuscript (Section 3.4), we recommend to manually fix the variance of the random intercept at a large value. 
#' This can be done in glmmTMB() by first setting up the model, but do not yet fit it:
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
#str(wolf.tmp)

#' Then fit the model and look at the results:
#+ echo=TRUE, message=FALSE, cache=TRUE 
wolf.rsf <- glmmTMB:::fitTMB(wolf.tmp)
summary(wolf.rsf)

saveRDS(wolf.rsf, "Outputs/RSF/wolf_rsf_wallotment_wtruncatedpopdens.rds")

wolf.rsf1 <- readRDS("Outputs/RSF/wolf_rsf_wallotment.rds")
summary(wolf.rsf1)

library(glmmTMB)
library(nloptr)
wolf.rsf <- readRDS("Outputs/RSF/wolf_rsf_wallotment_wtruncatedpopdens.rds")
summary(wolf.rsf)

######---- JAGS CODE STARTS FOR RSF ----#####

#UNSURE AB THIS. ITS NICE BC IT HAS THE PENALIZED COMPLEXITY PRIORS BUT
#ITS ALSO SLOW, AND I HAVE TO ADD WEIGHTS FOR WEIGHTED LIKELIHOOD

RSFmod <- nimbleCode({
  #use of penalized complexity priors with only two covariates to start
  
  #mean.RSF ~ dnorm(0, sd=5) # prior for mean RSF

  #forest and disthwy are out bc conflict w ruggedness/agric and elev, respectively
  beta.forest ~ dnorm(0,sd=sigma.beta.forest)
  sigma.beta.forest ~ dexp(1)
  beta.rugged ~ dnorm(0,sd=sigma.beta.rugged)
  sigma.beta.rugged ~ dexp(1)
  beta.rugged2 ~ dnorm(0,sd=sigma.beta.rugged2)
  sigma.beta.rugged2 ~ dexp(1)
  beta.agric ~ dnorm(0,sd=sigma.beta.agric)
  sigma.beta.agric ~ dexp(1)
  beta.deer ~ dnorm(0,sd=sigma.beta.deer)
  sigma.beta.deer ~ dexp(1)
  beta.allotment ~ dnorm(0,sd=sigma.beta.allotment)
  sigma.beta.allotment ~ dexp(1)
  beta.disthwy ~ dnorm(0,sd=sigma.beta.disthwy)
  sigma.beta.disthwy ~ dexp(1)
  beta.elev1 ~ dnorm(0,sd=sigma.beta.elev1)
  sigma.beta.elev1 ~ dexp(1)
  beta.elev2 ~ dnorm(0,sd=sigma.beta.elev2)
  sigma.beta.elev2 ~ dexp(1)
  beta.grassland ~ dnorm(0,sd=sigma.beta.grassland)
  sigma.beta.grassland ~ dexp(1)
  beta.popdens ~ dnorm(0,sd=sigma.beta.popdens)
  sigma.beta.popdens ~ dexp(1)
  beta.roaddens ~ dnorm(0,sd=sigma.beta.roaddens)
  sigma.beta.roaddens ~ dexp(1)
  beta.shrubland ~ dnorm(0,sd=sigma.beta.shrubland)
  sigma.beta.shrubland ~ dexp(1)
  beta.elk ~ dnorm(0,sd=sigma.beta.elk)
  sigma.beta.elk ~ dexp(1)
  
  # beta.year ~ dnorm(0,sd=sigma.beta.year)
  # sigma.beta.year ~ dexp(1)
  
  #random effect of cam
  for(i in 1:npacks){
    #there are 1343 cameras across 1303 camera clusters
    eps.pack[i] ~ dnorm(0, sd = sigma.pack) #need to specify that it is sd
  }
  sigma.pack ~ dexp(1)
  
  # for(i in 1:nyear){
  #   #there are 1343 cameras across 1303 camera clusters
  #   eps.year[i] ~ dnorm(0, sd = sigma.year) #need to specify that it is sd
  # }
  # sigma.year ~ dexp(1)

  for(i in 1:N){ #p is number of packs, 193
    logit(RSF[i])<- beta.forest*forest.focal[i] + 
      beta.rugged*rugged.focal.truncated[i] + (beta.rugged2*(rugged.focal.truncated[i]^2))+
      beta.agric*agriculture_focal[i] + beta.deer*alldeer_final[i] +
      beta.allotment*allotment_01[i] + beta.disthwy*disthwy[i] +
      beta.roaddens*roaddensity_final[i] +
      beta.elev1*elev_focal[i] + (beta.elev2*(elev_focal[i]^2)) + 
      beta.grassland*grassland_focal[i] +
      beta.elk*elk_01[i] +
      beta.popdens*popdens_focal[i] + beta.shrubland*shrubland_focal[i] +
      eps.pack[pack[i]]
      y[i] ~ dbern(RSF[i])
    }
})

#put the data together
dat1<-list(y=y, 
           forest.focal=covs$forest_focal,
           rugged.focal.truncated=covs$rugged_focal_truncated,
           agriculture_focal=covs$agriculture_focal,
           alldeer_final=covs$alldeer_final,
           allotment_01=covs$allotment_01,
           roaddensity_final=covs$roaddensity_final,
           elev_focal=covs$roaddensity_final,
           grassland_focal=covs$grassland_focal,
           elk_01=covs$elk_01,
           popdens_focal=covs$popdens_focal,
           shrubland_focal=covs$shrubland_focal,
           disthwy=covs$disthwy)
           #year=c(-2,-1,0,1,2), #year will be a vector with five values
           #pack=pack) #season will be a vector from 1:4 (=spring, summer, fall, winter))
#dist=array(scale(dist[1:N, 1:T, 1:J1]), dim=c(N,T,J1)))

rn<-sample(seq(1,10000, by=1))
set.seed(rn[2])
inits<-list(beta.forest=rexp(1),
            beta.rugged=rexp(1),
            beta.rugged2=rexp(1),
            beta.agric=rexp(1),
            beta.deer=rexp(1),
            beta.allotment=rexp(1),
            beta.disthwy=rexp(1),
            beta.elev1=rexp(1),
            beta.elev2=rexp(1),
            beta.grassland=rexp(1),
            beta.popdens=rexp(1),
            beta.shrubland=rexp(1),
            beta.elk=rexp(1),
            beta.roaddens=rexp(1)
            #beta.year=rexp(1)
            #beta.focal=rexp(1)
            )

nb <- 500 #100
ni <- 2000 # 50000+nb
nt <- 1
nc <- 2
# nb <- 150000 #100
# ni <- 300000 # 50000+nb
# nt <- 5
# nc <- 3
adaptInterval <- 100

constants<-list(pack=pack, N=length(pack), npacks=max(pack)) #J1=J1,Jm=Jmvec,ind=dataMod$ind[1:N,1:T,1:J1])
Rmodel<-nimbleModel(code=RSFmod, constants=constants, data=dat1, inits=inits)
conf<-configureMCMC(Rmodel, thin=nt, control=list(adaptInterval = adaptInterval))

conf$addMonitors(c("beta.forest", "beta.rugged", "beta.rugged2", "beta.agric", "beta.deer", "beta.allotment",
                   "beta.disthwy", "beta.elev1", "beta.elev2", "beta.grassland", "beta.popdens", "beta.roaddens", 
                   "beta.shrubland", "beta.elk",
                   "RSF", "eps.pack"))

Rmcmc<-buildMCMC(conf)
Cmodel<-compileNimble(Rmodel, showCompilerOutput = TRUE)

Cmodel$setData(dat1)
Cmodel$setInits(inits)

Cmcmc<-compileNimble(Rmcmc,project=Cmodel, showCompilerOutput = TRUE)
#Cmcmc$run(thin=10, reset=T, niter=500, nburnin=20) Hannah
#Cmcmc$run(thin=10, reset=T, niter=1500000, nburnin=700000)


# x1<-(as.data.frame(as.matrix(Cmcmc$mvSamples)))
# write.csv(x1, "test.csv")

t.start <- Sys.time()
out1 <- runMCMC(Cmcmc, niter = ni , nburnin = nb , nchains = nc, inits=inits,
                setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)  
t.end <- Sys.time()
runTime <- t.end-t.start

#saveRDS(out1, file = "Outputs/Occupancy/out1.rds")








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
ggsave("G:/My Drive/GitHub/Wolves/Outputs/a_composite_model/Scenarios/Figures/paper/Chapter1_AppendixS4_FigureS2.eps")
ggsave("G:/My Drive/GitHub/Wolves/Outputs/a_composite_model/Scenarios/Figures/paper/Chapter1_AppendixS4_FigureS2.jpg")

#plotting for glm
# require(GGally)
# ggcoef(tidy(multivariate_mod[[1]], conf.int = TRUE) %>% subset(std.error < 200), 
#        sort = "ascending", exclude_intercept = TRUE) +
#   ggtitle("Coefficients: Simple logistic model RSF")

#to make PREDICTIONS, will need to standardize the underlying rasters using the means and sd of the original covariates
#used to fit the models, not of the rasters themselves
#https://cran.r-project.org/web/packages/unmarked/vignettes/spp-dist.pdf


#####----- MAKING PREDICTIONS STARTS HERE -----#####

#read in rasters
rastlist <- list.files(path = "G:/My Drive/Data/Data/GIS/a_Final_Rasters", pattern='.tif$', all.files=TRUE, full.names=TRUE)
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
writeRaster(stack_all_forpred, filename=file.path("G:/My Drive/Data/Data/GIS/a_Final_Rasters_Scaled", names(stack_all_forpred)), bylayer=TRUE, format='GTiff', overwrite=TRUE)

######------ READING IN RASTERS AND PREDICTING RSF ------######

wolf.rsf <- readRDS("Outputs/RSF/wolf_rsf_wallotment_wtruncatedpopdens_weight100.rds")
summary(wolf.rsf)
rastlist <- list.files(path = "G:/My Drive/Data/Data/GIS/a_Final_Rasters_Scaled", pattern='.tif$', all.files=TRUE, full.names=TRUE)
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
#basically just dividing my max bc min is so small
#mean <- cellStats(predict, mean) 
max <- cellStats(predict, max) 

quantile(predict)

predict5 <- predict/max
#plot and write raster
plot(predict5)
writeRaster(predict5, "Outputs/RSF/RSF_noelk_trunc_rugged_newdeerinterp_newsample2021_trunpopdens.tif",overwrite=TRUE)

orig_agg <- aggregate(predict5,15)
invert <- 1-orig_agg
plot(invert)
writeRaster(invert, "Outputs/RSF/RSF_aggregated_invert.tif",overwrite=TRUE)