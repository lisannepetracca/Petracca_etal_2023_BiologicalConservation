get.move <- function(n.settlers.for.fxn,n.res,site_check,array_probs,array_siteID){

#as of January 31 2022, immigrants are no longer in this function at all
#we have also removed counting of solo vs. non-solo settlers
  #those former rules:
  #if occupancy draw is positive:
  #if less than two individuals (residents + immigrants + settlers), the settlers stay but are considered solo wolves
  #if there are >=2 individuals AND at least one is a settler, they stay and are counted as new wolves
  
#in this function, we are reassigning potential settlers to a new territory
#settlers pull a dispersal distance and a least cost path
#if a success, the settlers are assigned the new territory; otherwise settlers sent home

Nsamples <- dim(n.settlers.for.fxn)[1]

#these are the counts at i level
n.new.2.count <- n.new.3.count <- #n.new.2.count.wrejects <- n.new.3.count.wrejects <- 
  matrix(NA, nrow=Nsamples, ncol=224) #n.new.2.solo <- n.new.3.solo <-
   
 #now we are moving through the MCMC samples
 for(i in 1:Nsamples){    
  
  #make a vector that is the number of individuals by site
  #number of new settlers age 2 
  #number of new settlers age 3 
help.vector.2 <- help.vector.3 <- vector()

  #this is essentially when this becomes an IBM
  for(s in 1:224){
    help.vector.2 <- c(help.vector.2,rep(s,n.settlers.for.fxn[i,1,s])) #repeat site number for number of age 2 settlers
    help.vector.3 <- c(help.vector.3,rep(s,n.settlers.for.fxn[i,2,s])) #repeat site number for number of age 3 settlers
}
    
  #if there are no candidate movers, end the function and output all 0s 
  if((length(help.vector.2) + length(help.vector.3))==0) {
    n.new.2.temp <- n.new.3.temp <- rep(0,224) #n.new.2.temp.solo <- n.new.3.temp.solo #n.new.2.temp.wrejects <- n.new.3.temp.wrejects 
  } else {
    
  #if there are candidate movers, set up settlers matrix
  #settlers matrix includes the original site,
  #an indc for settler (1) vs immigrant (2), age 2 or 3, and 

   set.mat <- data.frame(orig.site = c(help.vector.2,help.vector.3),
                                  age = c(rep(2,length(help.vector.2)),rep(3,length(help.vector.3))))
    
   #have to account for possibility of getting an NA in that matrix for large distances
   for(ind in 1:nrow(set.mat)){
      success <- FALSE
       while (success==FALSE) {
         # do something
         temp_dist <- min(round(rgamma(n=1, shape=0.814, rate=0.005))+1, 633)
         test <- site_check[set.mat$orig.site[ind],temp_dist]
         x <- is.na(test)
         # check for success
         success <- x == FALSE 
       }
      
      #now we need to do multinomial thing
      probability_draw <- rmultinom(1, 1, 
                              array_probs[set.mat$orig.site[ind], temp_dist,
                                   1:length(which(!is.na(array_probs[set.mat$orig.site[ind],temp_dist,])))])
      site_choice <- which(as.vector(probability_draw==1))
      set.mat$new.site[ind] <- array_siteID[set.mat$orig.site[ind], temp_dist, site_choice]
      }
   
   #set up vectors to get the total numbers and total numbers of potential new guys by site 
   tots.news <- rep(NA,224)
   tots.site <- rep(NA,224)
   
   #total new guys by site  
   for(s in 1:224){
    tots.news[s] <- length(which(set.mat$new.site == s))
   }
  
   #tots.site is the total of everyone at the site
   #n.res is a summary of all stayers from the projection model
   tots.site <- tots.news + n.res[i,] 

   #do they stay or do they go?
   #going through 224 sites
   n.new.2.temp <- n.new.3.temp <- rep(0,224) #n.new.2.temp.solo <- n.new.3.temp.solon.new.2.temp.wrejects <- n.new.3.temp.wrejects 
   
   #rejected.2.settlers <- rejected.3.settlers <- vector() #solo.2.settlers <- solo.3.settlers 
   
   # ind.less2 is a 1 if there are less than 2 at the site, 0 otherwise
   # ind.less2 <- rep(0,224)
   # ind.less2[which(tots.site < 2)] <- 1 
   
   #get the occupancy probability for staying at that site
   # stay <- vector()
   # for(s in 1:224){
   #  stay[s] <- rbinom(1, 1, psi[i,s]) #sample a bernoulli outcome given psi for the site 
   # }
   
   for(s in 1:224){
    #identify site numbers that have single indiv, positive bernoulli draw, a new site, by age and indicator
    #THIS IS WHERE WE LET THEM STAY, BUT DONT COUNT THEM UNTIL THERE ARE 2
    # solo.2.settlers <- c(solo.2.settlers,which(ind.less2[set.mat$new.site] == 1 & stay[set.mat$new.site] == 1 &
    #                                              set.mat$new.site == s & set.mat$age == 2 & set.mat$indc == 1))
    # solo.3.settlers <- c(solo.3.settlers,which(ind.less2[set.mat$new.site] == 1 & stay[set.mat$new.site] == 1 &
    #                                              set.mat$new.site == s & set.mat$age == 3 & set.mat$indc == 1))
    
    #identify # settlers that have new site s and there are 2 inds there and it is occ = 1, by age 
    # n.new.2.temp[s] <- length(which(ind.less2[set.mat$new.site] == 0 & stay[set.mat$new.site] == 1 & set.mat$new.site == s & set.mat$age == 2))
    # n.new.3.temp[s] <- length(which(ind.less2[set.mat$new.site] == 0 & stay[set.mat$new.site] == 1 & set.mat$new.site == s & set.mat$age == 3))
     n.new.2.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 2)) #stay[set.mat$new.site] == 1 & 
     n.new.3.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 3)) #stay[set.mat$new.site] == 1 & 
    
    #identify row numbers that, regardless of whether 2 or not, were rejected
    #by age and indicator   
    # rejected.2.settlers <- c(rejected.2.settlers,which(stay[set.mat$new.site] == 0 & set.mat$new.site == s & set.mat$age == 2))
    # rejected.3.settlers <- c(rejected.3.settlers,which(stay[set.mat$new.site] == 0 & set.mat$new.site == s & set.mat$age == 3))
     } #close s loop 
   
   #send the rejected settlers back for the official count 
   #add in the lonely settlers in the population, but don't count them (solo)
   # for(s in 1:224){
   #  n.new.2.temp[s] <- n.new.2.temp[s] 
   #  n.new.3.temp[s] <- n.new.3.temp[s] 
   #  # n.new.2.temp.wrejects[s] <- n.new.2.temp[s] + length(which(set.mat$orig.site[rejected.2.settlers]==s))
   #  # n.new.3.temp.wrejects[s] <- n.new.3.temp[s] + length(which(set.mat$orig.site[rejected.3.settlers]==s))
   #  # n.new.2.temp.solo[s] <- n.new.2.temp[s]
   #  # n.new.3.temp.solo[s] <- n.new.3.temp[s]
   #  # n.new.2.temp.solo[s] <- n.new.2.temp.solo[s] + length(which(set.mat$orig.site[solo.2.settlers]==s)) 
   #  # n.new.3.temp.solo[s] <- n.new.3.temp.solo[s] + length(which(set.mat$orig.site[solo.3.settlers]==s)) 
   # }
   }
  n.new.2.count[i,] <- n.new.2.temp #these are the new settlers and rejected settlers from same site
  n.new.3.count[i,] <- n.new.3.temp #these are the new settlers and rejected settlers from same site
  # n.new.2.count.wrejects[i,] <- n.new.2.temp.wrejects #these are the new settlers and rejected settlers from same site
  # n.new.3.count.wrejects[i,] <- n.new.3.temp.wrejects #these are the new settlers and rejected settlers from same site
  #n.new.2.solo[i,] <- n.new.2.temp.solo #these are the counted new ones, plus the solo settlers
  #n.new.3.solo[i,] <- n.new.3.temp.solo #these are the counted new ones, plus the solo settlers
  
  }#close i loop

#to bring back to model    
return(list(n.new.2.count, n.new.3.count))    #n.new.2.solo, n.new.3.solon.new.2.count.wrejects, n.new.3.count.wrejects

}

# n.settlers.for.fxn <- array(rpois(22400,1), dim=c(50,2,224))
# n.res <- array(rpois(22400,3), dim=c(50,224))
