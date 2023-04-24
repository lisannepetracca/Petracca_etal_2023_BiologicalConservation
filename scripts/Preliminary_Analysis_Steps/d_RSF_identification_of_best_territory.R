#---
#"Choosing 'best' territory at all distances from other territories, using output from RSF"
#"Also calculates multinomial probabilities of suitability for all territories at a given distance from another"
#"Lisanne Petracca"
#"Nov 2022"
#---

library(raster)
library(fasterize)
library(sf)

#read in shapefile of pack boundaries and assign coordinate system
packs <- st_read("C:/My_Directory/territory_polygons.shp")
st_crs(packs) <- 32610

#read in WA boundary
WA <- st_read("C:/My_Directory/washington_UTM_mainland.shp")

#rasterize extent of WA
r <- raster(extent(WA), resolution=1000)
ras_WA <- fasterize(WA, r)
crs(ras_WA) <- crs(packs)
plot(ras_WA)

output <- list()

#rasterize pack boundaries
for(i in 1:224){
r_pack <- raster(packs[i,], resolution=1000)
ras_pack <- fasterize(packs[i,], r_pack)
ras_pack[ras_pack==1] <- 2

#merge the two
merge <- raster::merge(ras_pack, ras_WA, tolerance = 0.5)
plot(merge)

#doing grid distance from each potential pack to all others
dist <- gridDistance(merge, 2, omit=NA)

dist_cells <- extract(dist, packs)

output[[i]] <- dist_cells
}

#make copy bc this step look a while
output_copy <- output

#gets distances in km
output_km <- list()
  for(i in 1:224){
  output_km[[i]] <- lapply(output[[i]], function(x) (round(x/1000)))
  output_km[[i]] <- lapply(output_km[[i]], function(x) (unique(x)))
  }

#this is collapsing list of lists. getting there.

#61 is the max number of distances (in km) between sites
output_array_site <- array(NA, dim=c(224,224,61))
nondup_res_df <- data.frame(Reduce(rbind, output_km))

#these are distances at which the sites are connected
for(i in 1:224){
  for(j in 1:224){
    for(k in 1:length(nondup_res_df[i,j][[1]])){
          output_array_site[i, j, k] <- as.numeric(nondup_res_df[i,j][[1]][k])
    }}}
head(output_array_site)

#ok, now let's get this array into 2d matrix like original approach
new.mat <- as.data.frame(matrix(NA, nrow=1500000,ncol=3))
colnames(new.mat) <- c("orig", "radius", "incl")

#this forms a "to" and "from" matrix
data_frame <- as.data.frame(array(aperm(output_array_site,c(1,2,3)), dim=c(50176,61)))
data_frame$orig <- rep(1:224, each=224)
data_frame$incl <- rep(1:224, times=224)

library(tidyverse)
data_frame <- data_frame %>% pivot_longer(
  cols = starts_with("V"),
  names_to = "site",
  values_to = "radius",
  values_drop_na = TRUE
)
data_frame <- data_frame[,c(1,4,2)]
head(data_frame)

#get rid of radius of 0, and distances over 632 km
data_frame_edit <- data_frame %>% dplyr::filter(radius>0 & radius <=632) 

#get rid of where incl==orig
data_frame_edit <- data_frame_edit[data_frame_edit$orig!=data_frame_edit$incl,]

#add column
data_frame_edit$median_RSF <- NA

#let's read in the median RSF values
median_RSF <- read.csv("C:/My_Directory/median_RSF_values.csv")
head(median_RSF)
median_RSF <- median_RSF[,c(2,3)]
colnames(median_RSF) <- c("ID", "median")
median <- median_RSF$median

#let's assume mean_RSF (or median) is a vector of 224 values
for(i in 1:nrow(data_frame_edit)){
  data_frame_edit$median_RSF[i] <-median[data_frame_edit$incl[i]]
}

#then add in possibility of returning to orig
#let's use out_new to find all pairs of operational distances from source
out_pairs <- data_frame_edit %>% group_by(orig, radius) %>% slice(1L)
out_pairs <- out_pairs[,c(1,2)]
out_pairs$incl <- out_pairs$orig
for(i in 1:nrow(out_pairs)){
  out_pairs$median_RSF[i] <-median[out_pairs$incl[i]]
}

tail(out_pairs)
subset <- subset(out_pairs,radius==632)

#now we need to merge out_new and out_pairs
out_all <- rbind(data_frame_edit, out_pairs)

#group them by origin and distance, get number in sequence as seq
out_seq <- out_all %>% group_by(orig, radius) %>% mutate(id = seq_len(n())) %>% arrange(orig, radius, id)
max <- max(out_seq$id)

#max number is 50 neighbors, and max distance is 632 km
array_RSF <- array(NA, dim=c(224, 632, max))
array_siteID <- array(NA, dim=c(224, 632, max))

#for origin site, distance, and number of sites, will produce mean rsf and site number
for (i in 1:nrow(out_seq)){
  array_RSF[out_seq$orig[i], out_seq$radius[i], out_seq$id[i]] <- as.numeric(out_seq$median_RSF[i])
  array_siteID[out_seq$orig[i], out_seq$radius[i], out_seq$id[i]] <- as.numeric(out_seq$incl[i])
}

#also come up w basic matrix re: whether there are territories at that distance 
site_check <- array_RSF[,,1]
dim(site_check)
site_check[site_check>=0] <- 1
unique(as.vector(as.matrix(site_check)))
site_check_add <- as.matrix(rep(NA, 224))
site_check <- cbind(site_check, site_check_add)

#calculate probabilities of choosing those RSF values
#do this across rows in the third dimension
get_probs <- function(array_row) {
  vector <- array_row/sum(array_row, na.rm=T)
  return(vector)
}

test<- c(1,2,3,4, NA, NA, NA)
get_probs(test)

#this will get probabilities of each site and distance
array_probs <- array(data=NA, dim=c(224,632,max))
for(i in 1:224){
  for(j in 1:632){
    array_probs[i,j,] <- get_probs(array_RSF[i,j,])
  }}

save(array_probs, array_siteID, site_check, file= "C:/My_Directory/territory_choice.RData")