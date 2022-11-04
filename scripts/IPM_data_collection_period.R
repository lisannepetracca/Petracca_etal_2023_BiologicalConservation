#---
#"Wolf IPM for Washington State"
#"Lisanne Petracca"
#"Nov 2022"
#---

library(jagsUI)

#This is an IPM model for wolves in WA State that incorporates three different components: 
# (1) abundance (end-of-year counts -- minimum counts from winter aerial surveys)
# (2) reproduction (end-of-year pup counts -- pre-2015 single counts from winter surveys)
# (3) survival (from GPS-collared individuals)

#We are working with three ages at time of survey (6, 18, 30+ months)

#age class 1 in survival model is 6-11.99 mos
#age class 2 in survival model is 12-23.99 mos
#age class 3 in survival model is 24+ mos

#Abundance model has abundance following a log-normal distribution based on counts
#Reproduction is a multinomial process and is conditional on having 2+ reproductive-aged individuals (you will see use of step function)
#Survival model is known fate and multistate. Has random effect of six-month period.

######------ IPM MODEL ------######

library(jagsUI)
library(here)
library(stringr)
library(coda)
library(data.table)
library(tidyverse)

#load data for model
load("GitHub/Petracca_et_al_EcoApps/data/Data_for_IPM.Rdata")
load("GitHub/Petracca_et_al_EcoApps/data/Spatial_Information.RData")

#this stable age distribution comes from survival estimates, and median observed f and pack size
stableAge <- c(0.123, 0.114, 0.763)

#10/33 wolves dispersed at <=18 mo
probImmig <- c(0.3, 0.7) #there are only two entries here because pup class (6-11.99 mos) does not immigrate

#constants for IPM
G <- 3 #age groups
periods <- 2 #1 = December, 2 = June
sites <- 34 #number of sites
years <- 12 #number of years (2009-2020)

firstyear
surv.history <- as.matrix(surv.history)

sink("wolf_jags.txt")
cat("
    model {
  
  ##### ABUNDANCE STATE PROCESS ##### 
  
  # total abundance per pack-year needs to start somewhere
  pop_init ~ dunif(0,20) #this is at the pack-year level
  lambda.immig ~ dunif(0,5) #prior on immigration component

  #in general, lets set to 0 those nodes for which we don't have data
  #will make calcs of lambda and Nglobal easier such that there are no sums with NAs
  
      for (s in 1:length(terr_missing)){ #32 territories that do not have full data for 12 years; so only 2/34 packs have full data
        for(t in 1:occasions_missing[s]){ #number of years per pack with missing data
                    Ntot[yr_missing[s,t],terr_missing[s]] <-0 
                    N[1,1,yr_missing[s,t],terr_missing[s]] <- 0 
                    N[1,2,yr_missing[s,t],terr_missing[s]] <- 0 
                    N[2,1,yr_missing[s,t],terr_missing[s]] <- 0 
                    N[2,2,yr_missing[s,t],terr_missing[s]] <- 0 
                    N[3,1,yr_missing[s,t],terr_missing[s]] <- 0 
                    N[3,2,yr_missing[s,t],terr_missing[s]] <- 0 
           for(p in 1:periods){#there are two periods in each year
                    Tot.immig[p,yr_missing[s,t],terr_missing[s]] <- 0
           for(j in 1:G){ #age classes
                   N.immig[j,p,yr_missing[s,t],terr_missing[s]] <- 0
           }}}}
    
    #there are no 2-3 yo movers in very first period of first year
    #here, first year is 1 (rather than first year at a site) bc movement is latent in the data period
    #this is all movers across all sites in that period
     for(j in 2:3){
     N.movers[j,1,1] <- 0
     N.movers.newmove[j,1,1] <- 0
     N.movers.oldmove[j,1,1] <- 0
     }
    
    for(t in 1:years){
    N.movers[1,1,t] <- 0 #there are no 6-mo movers
    }
  
  #need to initialize Ntot only when count data first appear
  #this will be indexed by period (December is 1, June is 2), sampling year[1-12], and site s
  
  for (s in 1:sites){ #running through 34 territories w data
    
    Ntot[firstyear[s],s] ~ dpois(pop_init) 

    N[1,1,firstyear[s],s] ~ dbin(stableAge[1], Ntot[firstyear[s],s]) 
 
    #moving from December to June of first year here (2009 is December 2009 and June 2010)
    #the phis and eps are described below in survival model but BASICALLY
    #mu.phiA[1:3] is survival of resident in each of three age classes [1-3]
    #mu.phiB[1:3] is survival of mover in each of three age classes [1-3]
    #epsA[1:2] is probability of initiating movement for latter two age classes (as pup class can't move)
    #epsB[1:2] is probability of continuing movement for latter two age classes (as pup class can't move)
    #alpha is probability of staying in state while moving
    
    #phi is 2t-1 for params that have a period of 2
    #added the 1-epsA bc technically that age can move now
    N[1,2,firstyear[s],s] ~ dbin(phiA[2,(2*firstyear[s]-1)]*(1-epsA[1]),  N[1,1,firstyear[s],s])
    
    #now we are working with 18 mo old individuals
    # group 2:(G-1) prob of group g (2) given not in any previous group
    N[2,1,firstyear[s],s] ~ dbin(stableAge[2]/(1-stableAge[1]), 
                                   Ntot[firstyear[s],s] - N[1,1,firstyear[s],s])

    #stayers - intermediate class - the first period's 18-mo old residents survive and don't initiate movement
    N.stayers[2,2,firstyear[s],s] ~ dbin(phiA[2,(2*firstyear[s]-1)]*(1-epsA[1]), N[2,1,firstyear[s],s])
    
    #N is sum of stayers and immigrants
    N[2,2,firstyear[s],s] <- N.stayers[2,2,firstyear[s],s] + N.immig[2,2,firstyear[s],s] 
    
    #final age class here for first year
    N[3,1,firstyear[s],s] <- Ntot[firstyear[s],s] - sum(N[1:2,1,firstyear[s],s]) 

    #stayers - intermediate class - last period's 30-mo+ old residents survive and don't initiate movement  
    N.stayers[3,2,firstyear[s],s] ~ dbin(phiA[3,(2*firstyear[s]-1)]*(1-epsA[2]), N[3,1,firstyear[s],s])
    
    #N is sum of stayers and immigrants
    N[3,2,firstyear[s],s] <- N.stayers[3,2,firstyear[s],s] + N.immig[3,2,firstyear[s],s] 

    #no immigrants in first period of first year because pop initialized w stable age distribution
    Tot.immig[1,firstyear[s],s] <- 0
    for(g in 1:G){
      N.immig[g,1,firstyear[s],s] <- 0
    }
    
    #add immigration to year 1, period 2 
    #this is an immigration component that can allow some variation on estimated pack count
    
    Tot.immig[2,firstyear[s],s] ~ dpois(lambda.immig)
    
    #there are no 12-mo immigrants
    N.immig[1,2,firstyear[s],s] <- 0 
  
    #18 mo old class is the first class that can immigrate
    N.immig[2,2,firstyear[s],s] ~ dbin(probImmig[1], Tot.immig[2,firstyear[s],s])
    
    # group G: deterministic
    N.immig[3,2,firstyear[s],s] <- Tot.immig[2,firstyear[s],s] - N.immig[2,2,firstyear[s],s]
  }

  #random prior on settle, the probability of a mover settling in a given pack
  #this parameter had to change yearly based on number of sites available to be occupied, hence the indexing
  for(t in 1:years){
    settle[t, 1:sites_per_year[t]] ~ ddirch(alpha_dirichlet[t,1:sites_per_year[t]])
  } 
  
  #here are all the mover wolves
  #these wolves, importantly, do not belong to any pack but can then occupy a new pack in the next time step
  #these loops get N.movers in a given year [1:12] and period [1:2] across all sites
  #there are two loops here: one from 1-11 for 2-->1 transitions, and 1-12 for 1-->2 transitions
    
    for(t in 1:11){     
      #18mo movers
      #last periods 12-mo olds survive either keep moving (old move) or start moving (new move)
      N.movers.newmove[2,1,t+1] ~ dbin(phiA[2,(2*t)]*epsA[1]*alpha, sum(N[1,2,t,])) 
      N.movers.oldmove[2,1,t+1] ~ dbin(phiB[2,(2*t)]*epsB[1]*alpha, N.movers[1,2,t])
      N.movers[2,1,t+1] <- N.movers.newmove[2,1,t+1] + N.movers.oldmove[2,1,t+1]

      #30+mo movers
      #new and old movers - intermediate class - last period's 24-mo old and 36-mo+ old residents survive and start moving but stay in state AND last period's 24-mo old and 36-mo+ old movers continue moving but stay in state
      N.movers.newmove[3,1,t+1] ~ dbin(phiA[3,(2*t)]*epsA[2]*alpha, sum(N[2,2,t,]) + sum(N[3,2,t,])) 
      N.movers.oldmove[3,1,t+1] ~ dbin(phiB[3,(2*t)]*epsB[2]*alpha, N.movers[2,2,t]+ N.movers[3,2,t])
      N.movers[3,1,t+1] <- N.movers.newmove[3,1,t+1] + N.movers.oldmove[3,1,t+1]
    }
      
    for(t in 1:12){ #back to 12 bc we need estimates for june 2021
      #12 mo movers
      N.movers[1,2,t] ~ dbin(phiA[2,(2*t-1)]*epsA[1]*alpha, sum(N[1,1,t,]))

      #24mo movers
      #NOTE THAT T INDICES STAY T BECAUSE WITHIN SAME DEC-FEB SAMPLING PERIOD
      #new and old movers - intermediate class - last period's 18-mo old residents survive and start moving but stay in state AND last period's 18-mo old movers continue moving but stay in state 
      N.movers.newmove[2,2,t] ~ dbin(phiA[2,(2*t-1)]*epsA[1]*alpha, sum(N[2,1,t,]))
      N.movers.oldmove[2,2,t] ~ dbin(phiB[2,(2*t-1)]*epsB[1]*alpha, N.movers[2,1,t]) 
      N.movers[2,2,t] <- N.movers.newmove[2,2,t] + N.movers.oldmove[2,2,t]
             
      #36+mo movers
      #new movers - intermediate class - last period's 30-mo+ old residents survive and start moving but stay in state AND last period's 30-mo+ old movers continue moving but stay in state 
      N.movers.newmove[3,2,t] ~ dbin(phiA[3,(2*t-1)]*epsA[2]*alpha, sum(N[3,1,t,]))
      N.movers.oldmove[3,2,t] ~ dbin(phiB[3,(2*t-1)]*epsB[2]*alpha, N.movers[3,1,t])
      N.movers[3,2,t] <- N.movers.newmove[3,2,t] + N.movers.oldmove[3,2,t]
    }
  
  
  #assigning immigration and N to territories in second year and beyond
  #there are three territories that only have a single year of data (and it's the final year), so need to skip over those
      for (s in 1:length(morethanone)){ #running through 31 territories with >1 yr of data
      for (t in firstyear_morethanone[s]:11){ #t is first year to start

      #----------------------------- 
      
      #Dec 6 mo olds
      #had to use a max(0,N) thing to keep the model from having parent node issues
      N[1,1,t+1,morethanone[s]] <- max(0, lambda.pups[t+1,morethanone[s]] - removals[1,t+1,morethanone[s]])
      
      #----------------------------- 
      
      #Jun 12 mo 
      #last period's 6-mo olds survive
      #stayers - intermediate class - last period's 6-mo old residents survive and don't initiate movement  
      N[1,2,t+1,morethanone[s]] ~ dbin(phiA[2,(2*t-1)]*(1-epsA[1]), N[1,1,t+1,morethanone[s]]) #ADDED PLUS ONES HERE
      
      #----------------------------- 
      
      #Dec 18 mo 
      
      #last period's 12-mo olds survive and don't start moving
      N.stayers[2,1,t+1,morethanone[s]] ~ dbin(phiA[2,2*t]*(1-epsA[1]), N[1,2,t,morethanone[s]])
      
      #settlers - intermediate class - last periods 12-mo old new movers survive, stop moving, and settle at s  
      N.settlers[2,1,t+1,morethanone[s]] ~ dbin(phiB[2,2*t]*(1-epsB[1])*settle[t+1,settle_reps[t+1,s]], N.movers[1,2,t]) 
      
      N[2,1,t+1,morethanone[s]] <- max(0,N.stayers[2,1,t+1,morethanone[s]] + N.settlers[2,1,t+1,morethanone[s]] +
                                   N.immig[2,1,t+1,morethanone[s]] - removals[2,t+1,morethanone[s]])
                                   
      #----------------------------- 
      
      #Jun 24 mo 
      
      #stayers - intermediate class - last period's 18-mo old residents survive and don't initiate movement  
      N.stayers[2,2,t+1,morethanone[s]] ~ dbin(phiA[2,(2*t-1)]*(1-epsA[1]), N[2,1,t+1,morethanone[s]]) #ADDED PLUS ONES HERE
      
      #settlers - intermediate class - last periods 18-mo old new movers survive, stop moving, and settle at s  
      N.settlers[2,2,t+1,morethanone[s]] ~ dbin(phiB[2,(2*t-1)]*(1-epsB[1])*settle[t+1,settle_reps[t+1,s]], N.movers[2,1,t+1]) #ADDED PLUS ONES HERE
      
      #final age class - sum this period's stayers and settlers 
      #add immigrants here 
      N[2,2,t+1,morethanone[s]] <- N.stayers[2,2,t+1,morethanone[s]] + N.settlers[2,2,t+1,morethanone[s]] + N.immig[2,2,t+1,morethanone[s]] 
      
      #-----------------------------
      
      #Dec 30 mo+ (i.e., 30 mo, 42 mo, 54 mo, 66 mo, 78 mo, 90 mo...)
      #key here is that this group is 30 PLUS months; has 30, 42, 54 mos
      
      #stayers - intermediate class - last period's 24-mo old residents and 36-mo+ old residents survive and don't initiate movement 
      N.stayers[3,1,t+1,morethanone[s]] ~ dbin(phiA[3,2*t]*(1-epsA[2]), N[2,2,t,morethanone[s]]+N[3,2,t,morethanone[s]]) 
      
      #settlers - intermediate class - last periods 24-mo old movers and 36-mo+ old movers settle at s  
      N.settlers[3,1,t+1,morethanone[s]] ~ dbin(phiB[3,2*t]*(1-epsB[2])*settle[t+1,settle_reps[t+1,s]], N.movers[2,2,t]+N.movers[3,2,t])
      
      #final age class - sum this period's stayers and settlers  
      #add immigrants here - we shouldn't get earlier immigrants because we are assuming they don't move until 18 months 
      N[3,1,t+1,morethanone[s]] <- max(0,N.stayers[3,1,t+1,morethanone[s]] + N.settlers[3,1,t+1,morethanone[s]] + N.immig[3,1,t+1,morethanone[s]] - removals[3,t+1,morethanone[s]])
      
      #----------------------------- 
      
      #Jun 36 mo+ (i.e., 36 mo, 48 mo, 60 mo, 72 mo, 84 mo, 96 mo...) 
      
      #stayers - intermediate class - last period's 30-mo+ old residents survive and don't initiate movement  
      N.stayers[3,2,t+1,morethanone[s]] ~ dbin(phiA[3,(2*t-1)]*(1-epsA[2]), N[3,1,t+1,morethanone[s]]) #ADDED PLUS ONES HERE
      
      #settlers - intermediate class - last period's 30-mo+ old movers settle at s  
      N.settlers[3,2,t+1,morethanone[s]] ~ dbin(phiB[3,(2*t-1)]*(1-epsB[2])*settle[t+1,settle_reps[t+1,s]], N.movers[3,1,t+1]) #ADDED PLUS ONES HERE
      
      #final age class - sum this period's stayers and settlers 
      #add immigrants here 
      N[3,2,t+1,morethanone[s]] <- N.stayers[3,2,t+1,morethanone[s]] + N.settlers[3,2,t+1,morethanone[s]] + N.immig[3,2,t+1,morethanone[s]]
      
      #----------------------------- 
      
      #let's deal w immigrants now
    
      for(p in 1:2){
        Tot.immig[p,t+1,morethanone[s]] ~ dpois(lambda.immig)
        
        #there are no immigrant 6-11.99 mo olds
        N.immig[1,p,t+1,morethanone[s]] <- 0 
        
        #12-23.99 mo old class is the first class that can immigrate
        N.immig[2,p,t+1,morethanone[s]] ~ dbin(probImmig[1], Tot.immig[p,t+1,morethanone[s]])
        
        # group G: deterministic
        N.immig[3,p,t+1,morethanone[s]] <- Tot.immig[p,t+1,morethanone[s]] - N.immig[2,p,t+1,morethanone[s]]
      } 
      
      #Ntot is a sum of all 3 age classes within each pack-year
      Ntot[t+1,morethanone[s]] <-         
        N[1,1,t+1,morethanone[s]] + N[2,1,t+1,morethanone[s]] + N[3,1,t+1,morethanone[s]]
    }}
  

  ##### DERIVED TOTAL NUMBER OF WOLVES PER PACK-YEAR, AND POP GROWTH ESTIMATION #####
  
  #Nglobal refers to total number of wolves in WA per year
  for(t in 1:12){
    Nglobal[t] <- sum(Ntot[t,])  
    Nglobal_wmove[t] <- Nglobal[t] + sum(N.movers[,1,t])
  }
  
  #popgrowth is essentially lambda
  for(t in 1:(years-1)){ 
    popgrowth[t] <- Nglobal[t+1]/(Nglobal[t]+ 0.0001)
    popgrowth_wmove[t] <- Nglobal_wmove[t+1]/(Nglobal_wmove[t]+ 0.0001)
  }
  
    popgrowth_geomean <- prod(popgrowth_wmove[1:11]) ^ (1/11)

  ##### OBSERVATION MODEL ON ABUNDANCE ##### 
  
  precision.pack ~ dgamma(1, 1)

  #here we are modeling abundance as log-normal
  for(s in 1:sites){ 
    for (t in firstyear[s]:12){ #then over number of years with data
        N_singlesonly_plusone_log[s,t] ~ dnorm(log(Ntot[t,s]+1), precision.pack) #correct on log scale
      }}

  ##### OBSERVATION MODEL ON REPRODUCTION #####
  
  # This model uses WDFW end of year counts only
  
  alpha.pups <- c(1,1,1,1,1,1,1)
  probs.pup[1:7] ~ ddirch(alpha.pups[1:7])
  
  for (s in 1:sites){ #looping over all packs
    for (t in firstyear[s]:12){ #looping over all years
      
      #here is a step function that basically says that number of 6-mo old pups is lambda.pups when there are >=2 reproducing individuals,
      #and is set to 0 otherwise
      #if/else statement; 1 if e >= 0; 0 otherwise
      
      lambda.pups[t,s] <-  (pup_count_EOY_plusone[s,t]-1) * step(N[3,1,t,s]-2)

      pup_count_EOY_plusone[s,t] ~ dcat(probs.pup[1:7])
      
      
    }}
    
  ####### SURVIVAL MODEL ########
    
  alpha ~ dbeta(1,1) #probability of staying in state
  sigma.period ~ dunif(0, 10) # random effect of six-month period
  tau.period <- 1/(sigma.period*sigma.period)
  
  for(i in 1:2){ #only two age classes can move
    epsA[i] ~ dbeta(1,1) #epsA is prob of initiating movement
    epsB[i] ~ dbeta(36,84) #epsB is prob of continuing movement; prior taken from Jimenez paper
    int.surv1[i] ~ dnorm(0, 0.01)
    int.surv2[i] ~ dnorm(0, 0.01)}

  for (t in 1:surv.periods) {
    eps.period[t] ~ dnorm(0, tau.period)}
  
  for (t in 1:surv.periods) {
    phiA[1,t] <- 0 #setting these to 0 bc we are assigning age class 1 to phiA/B 2
    phiB[1,t] <- 0
    for(a in 1:2){
     logit(phiA[a+1,t]) <- int.surv1[a] + eps.period[t] #phiA is survival of non-mover, and has RE of six-month period
     logit(phiB[a+1,t]) <- int.surv2[a] + eps.period[t] #phiB is survival of mover; phiB[1,t] set to 0
  }}
    
#there are T-1 transitions, so t only goes to 22
for(t in 1:(surv.periods)){
  
  gamma[1,1,1,t] <- phiA[2,t] * (1 - epsA[1]) 
  gamma[1,1,2,t] <- 0
  gamma[1,1,3,t] <- 0 
  gamma[1,1,4,t] <- 1 - phiA[2,t]
  
  gamma[2,1,1,t] <- phiA[2,t] * (1 - epsA[1]) 
  gamma[2,1,2,t] <- phiA[2,t] * epsA[1] * alpha
  gamma[2,1,3,t] <- phiA[2,t] * epsA[1] * (1-alpha)
  gamma[2,1,4,t] <- 1 - phiA[2,t]
  
  gamma[3,1,1,t] <- phiA[3,t] * (1 - epsA[2]) 
  gamma[3,1,2,t] <- phiA[3,t] * epsA[2] * alpha
  gamma[3,1,3,t] <- phiA[3,t] * epsA[2] * (1-alpha)
  gamma[3,1,4,t] <- 1 - phiA[3,t]
  
  gamma[1,2,1,t] <- 0
  gamma[1,2,2,t] <- 0
  gamma[1,2,3,t] <- 0
  gamma[1,2,4,t] <- 0
  
  gamma[2,2,1,t] <- phiB[2,t] * (1-epsB[1])
  gamma[2,2,2,t] <- phiB[2,t] * epsB[1] * alpha
  gamma[2,2,3,t] <- phiB[2,t] * epsB[1] * (1-alpha)
  gamma[2,2,4,t] <- 1 - phiB[2,t]
  
  gamma[3,2,1,t] <- phiB[3,t] * (1-epsB[2])
  gamma[3,2,2,t] <- phiB[3,t] * epsB[2] * alpha
  gamma[3,2,3,t] <- phiB[3,t] * epsB[2] * (1-alpha)
  gamma[3,2,4,t] <- 1 - phiB[3,t]
  
  gamma[1,3,1,t] <- 0
  gamma[1,3,2,t] <- 0
  gamma[1,3,3,t] <- 0
  gamma[1,3,4,t] <- 1
  
  gamma[2,3,1,t] <- 0
  gamma[2,3,2,t] <- 0
  gamma[2,3,3,t] <- 0
  gamma[2,3,4,t] <- 1
  
  gamma[3,3,1,t] <- 0
  gamma[3,3,2,t] <- 0
  gamma[3,3,3,t] <- 0
  gamma[3,3,4,t] <- 1
  
  gamma[1,4,1,t] <- 0
  gamma[1,4,2,t] <- 0
  gamma[1,4,3,t] <- 0
  gamma[1,4,4,t] <- 1
  
  gamma[2,4,1,t] <- 0
  gamma[2,4,2,t] <- 0
  gamma[2,4,3,t] <- 0
  gamma[2,4,4,t] <- 1
  
  gamma[3,4,1,t] <- 0
  gamma[3,4,2,t] <- 0
  gamma[3,4,3,t] <- 0
  gamma[3,4,4,t] <- 1
}
  
  # likelihood 
  for (i in 1:nhist){
    for (t in (first[i]+1):last[i]){ 
      y[i,t] ~ dcat(gamma[age[i,t-1],y[i,t-1],1:4,t-1]) #multistate survival model
    }
  }
  
  #derived params
  for(a in 1:2){
  mu.phiA[a+1] <- exp(int.surv1[a])/(1+exp(int.surv1[a])) #getting avg phiA
  mu.phiB[a+1] <- exp(int.surv2[a])/(1+exp(int.surv2[a])) #getting avg phiB
}


  }",fill=TRUE)
sink()


win.data <- list(sites=sites, years=years, periods=periods, G=G,
                  firstyear=firstyear, 
                  morethanone=morethanone, removals=removal_array,
                  terr_missing=terr_missing, yr_missing=yr_missing, occasions_missing=occasions_missing,
                  first = first, last = last, 
                  age=age, firstyear_morethanone=firstyear_morethanone,
                  y = surv.history, nhist = nrow(surv.history),
                  surv.periods=23, #we are technically adding a period here, to estimate survival from dec to june in final year
                  N_singlesonly_plusone_log=N_singlesonly_plusone_log, pup_count_EOY_plusone=pup_count_EOY_plusone,
                  stableAge=stableAge, probImmig=probImmig,
                  alpha_dirichlet=alpha_dirichlet, sites_per_year=sites_per_year,
                  settle_reps=settle_reps)


#provide initial values
inits <- function(){list(
  lambda.immig=3,
  epsA = runif(2,0,1),
  epsB = runif(2,0.05,0.1),
  alpha = runif(1,0,1),
  int.surv1=runif(2,-5,5), 
  int.surv2=runif(2,-5,5),
  sigma.period = runif(1), 
  eps.period = runif(23,0,1))}

nb <- 120000 #100000 #100 #30000 #1000 #2000 #15000 #change iterations to low, burn in low, spit out
ni <- 240000 #300000 #1000 #150000 #100 #2000 #10000 #40000  
nt <- 6 #10 #3
nc <- 3

# Parameters to estimate
params <- c(  "Nglobal", "Nglobal_wmove",
              "Tot.immig", "N.immig", "N.movers", 
              "N.movers.oldmove", "N.movers.newmove",
              "N.stayers", "N.settlers",
              "Ntot", "popgrowth", "popgrowth_wmove", "popgrowth_geomean", "N",
              "precision.pack", #"precision.pups",
              "lambda.immig",
              "pup_count_EOY_plusone", "lambda.pups", "probs.pup", #"f",
              "mu.phiA", "mu.phiB", "epsA", "epsB", "alpha",
              "int.surv1", "int.surv2", "eps.period", "sigma.period", "phiA", "phiB")

out <- jags(win.data, inits, params, "wolf_jags.txt", n.adapt=100, n.chains=nc, n.iter=ni, n.burn = nb, n.thin=nt, parallel=T)
print(out,n=3)

write.csv(as.mcmc(out$summary),"Outputs/a_composite_model/jags_removals_ddirch_Oct25_2022_epsB.csv")
saveRDS(out, file = "Outputs/a_composite_model/jags_removals_ddirch_Oct25_2022_epsB.rds")