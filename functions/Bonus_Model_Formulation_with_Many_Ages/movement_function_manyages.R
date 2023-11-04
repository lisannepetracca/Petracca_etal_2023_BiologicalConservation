get.move <- function(N.settlers.for.fxn,n.res,site_check,array_probs,array_siteID){

#in this function, we are reassigning potential settlers to a new territory
#settlers pull a dispersal distance and a least cost path
#if a success, the settlers are assigned the new territory; otherwise settlers sent home

Nsamples <- dim(N.settlers.for.fxn)[1]

#these are the counts at i level
n.new.2.count <- n.new.3.count <- n.new.4.count <- 
  n.new.5.count <- n.new.6.count <- n.new.7.count <- 
  n.new.8.count <- n.new.9.count <- n.new.10.count <- 
  n.new.11.count <- n.new.12.count <- n.new.13.count <- 
  n.new.14.count <- n.new.15.count <- 
  matrix(NA, nrow=Nsamples, ncol=224) 
   
 #now we are moving through the MCMC samples
 for(i in 1:Nsamples){    
  
 #make a vector that is the number of individuals by site
help.vector.2 <- help.vector.3 <- help.vector.4 <- 
help.vector.5 <- help.vector.6 <- help.vector.7 <- 
help.vector.8 <- help.vector.9 <- help.vector.10 <- 
help.vector.11 <- help.vector.12 <- help.vector.13 <- 
help.vector.14 <- help.vector.15 <- vector()

  for(s in 1:224){
    help.vector.2 <- c(help.vector.2,rep(s,N.settlers.for.fxn[i,1,s])) #repeat site number for number of age 2 settlers
    help.vector.3 <- c(help.vector.3,rep(s,N.settlers.for.fxn[i,2,s])) #repeat site number for number of age 3 settlers
    help.vector.4 <- c(help.vector.4,rep(s,N.settlers.for.fxn[i,3,s])) #repeat site number for number of age 4 settlers
    help.vector.5 <- c(help.vector.5,rep(s,N.settlers.for.fxn[i,4,s])) #repeat site number for number of age 5 settlers
    help.vector.6 <- c(help.vector.6,rep(s,N.settlers.for.fxn[i,5,s])) #repeat site number for number of age 6 settlers
    help.vector.7 <- c(help.vector.7,rep(s,N.settlers.for.fxn[i,6,s])) #repeat site number for number of age 7 settlers
    help.vector.8 <- c(help.vector.8,rep(s,N.settlers.for.fxn[i,7,s])) #repeat site number for number of age 8 settlers
    help.vector.9 <- c(help.vector.9,rep(s,N.settlers.for.fxn[i,8,s])) #repeat site number for number of age 9 settlers
    help.vector.10 <- c(help.vector.10,rep(s,N.settlers.for.fxn[i,9,s])) #repeat site number for number of age 10 settlers
    help.vector.11 <- c(help.vector.11,rep(s,N.settlers.for.fxn[i,10,s])) #repeat site number for number of age 11 settlers
    help.vector.12 <- c(help.vector.12,rep(s,N.settlers.for.fxn[i,11,s])) #repeat site number for number of age 12 settlers
    help.vector.13 <- c(help.vector.13,rep(s,N.settlers.for.fxn[i,12,s])) #repeat site number for number of age 13 settlers
    help.vector.14 <- c(help.vector.14,rep(s,N.settlers.for.fxn[i,13,s])) #repeat site number for number of age 14 settlers
    help.vector.15 <- c(help.vector.15,rep(s,N.settlers.for.fxn[i,14,s])) #repeat site number for number of age 15 settlers
  }
    
  #if there are no candidate movers, end the function and output all 0s 
  if((length(help.vector.2) + length(help.vector.3)+length(help.vector.4)+
      length(help.vector.5) + length(help.vector.6)+length(help.vector.7)+
      length(help.vector.8) + length(help.vector.9)+length(help.vector.10)+
      length(help.vector.11) + length(help.vector.12)+length(help.vector.13)+
      length(help.vector.14) + length(help.vector.15))==0) {
      n.new.2.temp <- n.new.3.temp <- n.new.4.temp <- 
      n.new.5.temp <- n.new.6.temp <- n.new.7.temp <- 
      n.new.8.temp <- n.new.9.temp <- n.new.10.temp <- 
      n.new.11.temp <- n.new.12.temp <- n.new.13.temp <- 
      n.new.14.temp <- n.new.15.temp <- rep(0,224) 
  } else {
    
  #if there are candidate movers, set up settlers matrix
  #settlers matrix includes the original site and age

   set.mat <- data.frame(orig.site = c(help.vector.2, help.vector.3, help.vector.4, 
                                       help.vector.5, help.vector.6, help.vector.7,
                                       help.vector.8, help.vector.9, help.vector.10,
                                       help.vector.11,help.vector.12,help.vector.13,
                                       help.vector.14,help.vector.15),
                                  age = c(rep(2,length(help.vector.2)),rep(3,length(help.vector.3)),
                                          rep(4,length(help.vector.4)),rep(5,length(help.vector.5)),
                                          rep(6,length(help.vector.6)),rep(7,length(help.vector.7)),
                                          rep(8,length(help.vector.8)),rep(9,length(help.vector.9)),
                                          rep(10,length(help.vector.10)),rep(11,length(help.vector.11)),
                                          rep(12,length(help.vector.12)),rep(13,length(help.vector.13)),
                                          rep(14,length(help.vector.14)),rep(15,length(help.vector.15))))
   
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
   n.new.2.temp <- n.new.3.temp <- n.new.4.temp <- n.new.5.temp <- 
     n.new.6.temp <- n.new.7.temp <- n.new.8.temp <- n.new.9.temp <- 
     n.new.10.temp <- n.new.11.temp <- n.new.12.temp <- n.new.13.temp <- 
     n.new.14.temp <- n.new.15.temp <- rep(0,224) 
   
   for(s in 1:224){
     n.new.2.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 2)) #stay[set.mat$new.site] == 1 & 
     n.new.3.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 3)) #stay[set.mat$new.site] == 1 & 
     n.new.4.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 4)) #stay[set.mat$new.site] == 1 & 
     n.new.5.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 5)) #stay[set.mat$new.site] == 1 & 
     n.new.6.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 6)) #stay[set.mat$new.site] == 1 & 
     n.new.7.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 7)) #stay[set.mat$new.site] == 1 & 
     n.new.8.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 8)) #stay[set.mat$new.site] == 1 & 
     n.new.9.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 9)) #stay[set.mat$new.site] == 1 & 
     n.new.10.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 10)) #stay[set.mat$new.site] == 1 & 
     n.new.11.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 11)) #stay[set.mat$new.site] == 1 & 
     n.new.12.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 12)) #stay[set.mat$new.site] == 1 & 
     n.new.13.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 13)) #stay[set.mat$new.site] == 1 & 
     n.new.14.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 14)) #stay[set.mat$new.site] == 1 & 
     n.new.15.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 15)) #stay[set.mat$new.site] == 1 & 
     
     } #close s loop 
  }

  n.new.2.count[i,] <- n.new.2.temp #these are the new settlers and rejected settlers from same site
  n.new.3.count[i,] <- n.new.3.temp #these are the new settlers and rejected settlers from same site
  n.new.4.count[i,] <- n.new.4.temp #these are the new settlers and rejected settlers from same site
  n.new.5.count[i,] <- n.new.5.temp #these are the new settlers and rejected settlers from same site
  n.new.6.count[i,] <- n.new.6.temp #these are the new settlers and rejected settlers from same site
  n.new.7.count[i,] <- n.new.7.temp #these are the new settlers and rejected settlers from same site
  n.new.8.count[i,] <- n.new.8.temp #these are the new settlers and rejected settlers from same site
  n.new.9.count[i,] <- n.new.9.temp #these are the new settlers and rejected settlers from same site
  n.new.10.count[i,] <- n.new.10.temp #these are the new settlers and rejected settlers from same site
  n.new.11.count[i,] <- n.new.11.temp #these are the new settlers and rejected settlers from same site
  n.new.12.count[i,] <- n.new.12.temp #these are the new settlers and rejected settlers from same site
  n.new.13.count[i,] <- n.new.13.temp #these are the new settlers and rejected settlers from same site
  n.new.14.count[i,] <- n.new.14.temp #these are the new settlers and rejected settlers from same site
  n.new.15.count[i,] <- n.new.15.temp #these are the new settlers and rejected settlers from same site
  
  }#close i loop

#to bring back to model    
return(list(n.new.2.count, n.new.3.count, n.new.4.count, 
            n.new.5.count, n.new.6.count, n.new.7.count, 
            n.new.8.count, n.new.9.count, n.new.10.count, 
            n.new.11.count, n.new.12.count, n.new.13.count, 
            n.new.14.count, n.new.15.count))    

}

