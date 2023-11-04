get.solos <- function(n.wolves.solo.fxn, neighbor_list){
  
  #in this function, we are assigning two "solo" wolves to the same territory if those territories are adjacent
  
  Nsamples <- dim(n.wolves.solo.fxn)[1]
  
  for(i in 1:Nsamples){
    
    #get the site numbers of sites w only one wolf
    packs_qualify <- which(n.wolves.solo.fxn[i,1,] + n.wolves.solo.fxn[i,2,] + n.wolves.solo.fxn[i,3,]==1)
    
    #if there are not 2+ sites with only one wolf, get out of loop 
     if(length(packs_qualify) <=1){
      next} else {#if there are 2+ sites w solo wolves
        while(length(packs_qualify) >1){ 
          
          pack_selected <- sample(packs_qualify,1) #select random pack
          packs_qualify <- packs_qualify[packs_qualify!= pack_selected] #remove selected pack from packs_qualify
          
          #find all solo neighbors -- i.e., it is in packs_qualify and also a neighbor
          solo_neighbors <- intersect(packs_qualify, neighbor_list[[pack_selected]]) 
          #if there are none that qualify, move on
          if(length(solo_neighbors)==0){
            next}
          else{ #so long as there is a solo neighbor
            #get the age class of the solo individual
            if (length(solo_neighbors) >1) { 
              chosen_neighbor <- sample(solo_neighbors,1) } else {chosen_neighbor <- solo_neighbors}
            packs_qualify <- packs_qualify[packs_qualify!= chosen_neighbor]
            age.vector <- which(n.wolves.solo.fxn[i,,chosen_neighbor] != 0, arr.ind = T)
            #add it to the selected pack and remove from neighbor pack
            n.wolves.solo.fxn[i,age.vector,pack_selected] <- n.wolves.solo.fxn[i,age.vector,pack_selected] + 1
            n.wolves.solo.fxn[i,age.vector,chosen_neighbor] <- n.wolves.solo.fxn[i,age.vector,chosen_neighbor] - 1
          } #closes else
        } #closes while
      } #closes loop of searching for solos
  } #closes i loop

  #to bring back to model    
  return (n.wolves.solo.fxn)    
  
} #closes function