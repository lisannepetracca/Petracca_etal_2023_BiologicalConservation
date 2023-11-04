get.move <- function(n.settlers.for.fxn,n.res,site_check,mat_singleLCP){
  
  #in this function, we are reassigning potential settlers to a new territory
  #settlers pull a dispersal distance and a least cost path
  #if a success, the settlers are assigned the new territory; otherwise settlers sent home
  
  Nsamples <- dim(n.settlers.for.fxn)[1]
  
  #these are the counts at i level
  n.new.2.count <- n.new.3.count <- matrix(NA, nrow=Nsamples, ncol=224) 
  
  #now we are moving through the MCMC samples
  for(i in 1:Nsamples){    
    
    #make a vector that is the number of individuals by site
    #number of new settlers age 2 
    #number of new settlers age 3 
    help.vector.2 <- help.vector.3 <- vector()
    
    for(s in 1:224){
      help.vector.2 <- c(help.vector.2,rep(s,n.settlers.for.fxn[i,1,s])) #repeat site number for number of age 2 settlers
      help.vector.3 <- c(help.vector.3,rep(s,n.settlers.for.fxn[i,2,s])) #repeat site number for number of age 3 settlers
    }
    
    #if there are no candidate movers, end the function and output all 0s 
    if((length(help.vector.2) + length(help.vector.3))==0) {
      n.new.2.temp <- n.new.3.temp <- rep(0,224) 
    } else {
      
      #if there are candidate movers, set up settlers matrix
      #settlers matrix includes the original site and age
      
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
        
        #now we need to pull the single LCP, see if it stays or pulls a new dist
        lcp.site <- mat_singleLCP[set.mat$orig.site[ind],temp_dist]
        if(rbinom(1, 1, 0.1471)==1){
          set.mat$new.site[ind] <- set.mat$orig.site[ind]} else{set.mat$new.site[ind] <- lcp.site}
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
      n.new.2.temp <- n.new.3.temp <- rep(0,224) 
      
      for(s in 1:224){
        n.new.2.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 2))  
        n.new.3.temp[s] <- length(which(set.mat$new.site == s & set.mat$age == 3))  
         }
    }
    
    n.new.2.count[i,] <- n.new.2.temp #these are the new settlers and rejected settlers from same site
    n.new.3.count[i,] <- n.new.3.temp #these are the new settlers and rejected settlers from same site
    
  }#close i loop
  
  #to bring back to model    
  return(list(n.new.2.count, n.new.3.count))    
}