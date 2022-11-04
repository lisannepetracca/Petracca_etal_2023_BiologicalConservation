get.translocations <- function(n.wolves.transl.fxn, EWash, site) {
  
  dim(n.wolves.transl.fxn)
  
  if(site=="St Helens"){
    packs_translocate <- c(14,10) #these are the packs that qualify based on the RSF [I personally did not trust occupancy]
  }
  if(site=="Olympic"){
    packs_translocate <- c(145,125)
  }
  
  #in this function, we are removing wolves by pack in E Washington based on a removal rate
  Nsamples <- dim(n.wolves.transl.fxn)[1]
  
  pup.translocated.all <- ad.translocated.all <- matrix(0, nrow=Nsamples, ncol=224)
  
  for(i in 1:Nsamples){    
    
    pup.translocated <- rep(0,224)
    ad.translocated <- rep(0,224)
    
      #THIS IS WHERE TRANSLOCATION HAPPENS, AT YEAR 3 {ARBITRARILY SET}
          #goal is to take 2 adults in age class 3 and 2 adults in pup class and put them in territories 14 and 23
          #let's make sure we'd leave at least two adults and one pup behind after the translocation
          packs_qualify_intermed <- which(n.wolves.transl.fxn[i,3,EWash]>=4 & n.wolves.transl.fxn[i,1,EWash]>=3)
          packs_qualify <- EWash[packs_qualify_intermed]  
          
          if(length(packs_qualify)<=1){
              next} else{
            packs_source <- sample(packs_qualify,2) #select two packs that quality
              pup.translocated[packs_source] <- -2
              ad.translocated[packs_source] <- -2
              pup.translocated[packs_translocate] <- 2
              ad.translocated[packs_translocate] <- 2
            
            pup.translocated.all[i,] <- pup.translocated #these sad immigs are those that have been rejected based on bernoulli draw
            ad.translocated.all[i,] <- ad.translocated 

              }#close else
} #close i loop
  #to bring back to model    
  return(list(pup.translocated.all, ad.translocated.all))    
  
}

# set.seed(101)
# n.wolves.transl.fxn <- array(rpois(33600,2), dim=c(50,3,224))
# sum(n.wolves.transl.fxn[5,,]) #1308
# sum(n.wolves.transl.fxn[5,,c(14,23)]) #12, should go up to 20