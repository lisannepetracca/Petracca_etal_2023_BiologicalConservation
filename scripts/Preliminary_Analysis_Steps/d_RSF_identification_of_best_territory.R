#trying to see something with sf_distance
library(raster)
library(fasterize)
library(sf)

packs <- st_read("Outputs/spatial_model/territory_polygons.shp")
st_crs(packs) <- 32610

WA <- st_read("G:/My Drive/Data/Data/GIS/washington_UTM_mainland.shp")
r <- raster(extent(WA), resolution=1000)
ras_WA <- fasterize(WA, r)
crs(ras_WA) <- crs(packs)
plot(ras_WA)

output <- list()

for(i in 1:224){
r_pack <- raster(packs[i,], resolution=1000)
ras_pack <- fasterize(packs[i,], r_pack)
ras_pack[ras_pack==1] <- 2

merge <- raster::merge(ras_pack, ras_WA, tolerance = 0.5)
plot(merge)

dist <- gridDistance(merge, 2, omit=NA)

dist_cells <- extract(dist, packs)

output[[i]] <- dist_cells
}

output_copy <- output

output_km <- list()
  for(i in 1:224){
  output_km[[i]] <- lapply(output[[i]], function(x) (round(x/1000)))
  output_km[[i]] <- lapply(output_km[[i]], function(x) (unique(x)))
  }

#this is collapsing list of lists. getting there.

#61 is the max number of distances (in km) between sites
output_array_site <- array(NA, dim=c(224,224,61))
nondup_res_df <- data.frame(Reduce(rbind, output_km))
dim(nondup_res_df)
length(nondup_res_df[24,25][[1]])
nondup_res_df[24,25][[1]][5]

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
head(new.mat)

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
max(data_frame_edit$radius)
min(data_frame_edit$radius)

#get rid of where incl==orig
data_frame_edit <- data_frame_edit[data_frame_edit$orig!=data_frame_edit$incl,]

#add column
data_frame_edit$median_RSF <- NA

#let's read in the median RSF values
median_RSF <- read.csv("Outputs/RSF/median_RSF_values.csv")
head(median_RSF)
median_RSF <- median_RSF[,c(2,3)]
colnames(median_RSF) <- c("ID", "median")
median <- median_RSF$median

#let's assume mean_RSF (or median) is a vector of 224 values
for(i in 1:nrow(data_frame_edit)){
  data_frame_edit$median_RSF[i] <-median[data_frame_edit$incl[i]]
}

tail(data_frame_edit)

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
tail(out_all)
tail(out_pairs)

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

dim(array_probs)
sum(array_probs[1,385,], na.rm=T)

save(array_probs, array_siteID, site_check, file= "Outputs/spatial_model/territory_choice.RData")

#####---- GETTING RSF VALUES OF OCCUPIED CELLS ----####

library(tidyverse)

###figuring out min RSF median of occupied cells
#first, let's getting 1-34 pack ids of territories w wolves
packs_info <- read.csv("Data/spatial_model/packid_to_space.csv", header=T)
packs_info <- packs_info %>% dplyr::select(pack, packid, FID) %>% rename(Pack=pack) %>% 
  arrange(FID, desc=T) %>% mutate(FID_noNA = as.numeric(as.factor(FID)))
packs_info <- packs_info %>% group_by(FID_noNA) %>% slice_head()
#ok, so there are 34 territories that had wolves at any point
dim(packs_info)
head(packs_info)

#let's read in the median RSF values & produce vector called "median"
median_RSF <- read.csv("Outputs/RSF/median_RSF_values.csv")
head(median_RSF)
median_RSF <- median_RSF[,c(2,3)]
colnames(median_RSF) <- c("ID", "median")
median <- median_RSF$median

table(median[packs_info$FID])

#this is getting minimum median RSF of occupied cells (0.0006028379); this is low, and
#there are only 17 territories with lower median RSF
min(median[packs_info$FID])
length(which(median<median[19])) #there are 17 territories with lower median RSF than 19
median[19]

#now let's look at 5th percentile RSF of occupied cell
quantile(median[packs_info$FID], probs = seq(0, 1, 0.05))
length(which(median<0.0028739148 )) #there are 38 territories with lower RSF than 5th percentile

#let's get ids of cells with greater than that
terr_greater_5thpercent <- which(median>=0.0028739148) #there are 38 territories with lower RSF than 5th percentile
length(terr_greater_5thpercent)

#let's make map of cells that cannot be occupied

packs <- st_read("Outputs/spatial_model/territory_polygons.shp")
st_crs(packs) <- 32610

packs_bad <- packs[-terr_greater_5thpercent,]

ggplot() + 
  geom_sf(data = packs, size = 1, color = "black", fill = NA) + 
  geom_sf(data = packs_bad, size = 1, color = "black", fill = "orange") + 
  ggtitle("Sites < 5th percentile of median RSF")
ggsave("G:/My Drive/Data/Data/Figures/rejected_sites_5thpercentile.jpg")

save(terr_greater_5thpercent, file= "Outputs/spatial_model/terr_greater_5thpercent.RData")



#####---- PROBABILITY OF STAYING FOR EACH SOURCE AND DISTANCE BASED ON RSF ----####

load("Outputs/spatial_model/territory_choice.RData")
dim(array_probs)
dim(array_siteID)

array_siteID[10,400,]

#we are wanting to get the last non-NA for each source and distance (as that is source site)
get.last <-  function(x) max(which(!is.na(x)))

#creating empty matrices for storing the last non-NA, and the associated probability
test <- matrix(NA, nrow=224, ncol=632)
prob_stay <- matrix(NA, nrow=224, ncol=632)

#looping through and getting probability of staying
for(i in 1:224){
  for(j in 1:632){
test[i,j] <- get.last(array_siteID[i,j,])
prob_stay[i,j] <- array_probs[i,j,test[i,j]]
  }}

#checking outputs (yay!)
array_siteID[15,20,]
array_probs[15,20,]
prob_stay[15,20]

library(reshape2)
probstay <- as.data.frame(prob_stay)

#making figure of prob of staying for bad cell (106), medium cell (151) and good cell (173)
test <- cut(median, breaks = 3)

probstay$group <- test
head(probstay)
interest <- probstay %>% 
  melt(value.name = 'prob') %>% group_by(group, variable) %>% filter(!(is.na(prob))) %>% 
  summarise(mean_prob = mean(prob))

#seeing where these points are in space
ggplot(interest, aes(as.numeric(variable), mean_prob, color=as.factor(group))) +
  geom_line(lwd=1)+
  ylab("Probability of staying in territory") + xlab("Distance traveled (km)")+
  scale_color_discrete(name = "Quality of home site", labels = c("Lower Third", 
                                                                 "Middle Third", 
                                                                 "Upper Third"))
ggsave("G:/My Drive/Data/Data/Figures/Prob_Staying_DiffQualitySites.jpg")


#saving prob_stay
save(prob_stay, file= "Outputs/spatial_model/prob_stay_multinom.RData")

#let's summarize prob_stay a bit
load("Outputs/spatial_model/prob_stay_multinom.RData")
dim(prob_stay)

mean(prob_stay[,1:200],na.rm=T)
mean(prob_stay,na.rm=T)

#checking out an occupied cell
mean(prob_stay[197,],na.rm=T) #0.38; higher than .13 (collar)


#####----- OLD VERSION WHERE WE WERE USING INCREASING CIRCULAR BUFFERS -----#####

# library(sf)
# library(here)
# library(tidyverse)
# 
# centroids <- st_read("Outputs/spatial_model/territory_centroids.shp")
# packs <- st_read("Outputs/spatial_model/territory_polygons.shp")
# 
# st_crs(packs) <- st_crs(centroids)
# 
# packs$id <- 1:224
# centroids$id <- 1:224
# centroid_buff <- st_buffer(centroids,15000)
# centroid_buff$id <- 1:224
# 
# #making sure they read in ok
# plot(st_geometry(packs), col = sf.colors(12, categorical = TRUE), border = 'grey', 
#      axes = TRUE)
# plot(st_geometry(centroid_buff), #pch=16,
#      border = 'black', add = TRUE)
# 
# #another option
# #data    <- data.frame(id = paste0("ID_", 1:10), lon = runif(10, -10, 10), lat = runif(10, -10, 10), pop_size = runif(10, 0, 2000))
# #data_sf <- st_as_sf(data,coords = c("lon", "lat"), crs = 4326) %>% st_transform(CRS("+proj=laea"))
# 
# # plot(data_sf$geometry)
# 
# bufferR <- floor(seq(1000, 632000, length.out = 632)) # sequence of radii
# 
# out <- do.call("rbind", lapply(1:length(bufferR), function(y) {
#   bfr <- centroid_buff %>% st_buffer(bufferR[y]) %>% st_cast("LINESTRING") ## create Buffer
#   inters <- bfr %>% st_intersects(packs) 
#   do.call("rbind", lapply(which(sapply(inters, length)>0), 
#                           function(z) data.frame(orig = centroid_buff[z,]$id, radius = bufferR[y],
#                                                  incl = packs[unlist(inters[z]),]$id)))
# }))
# dim(out)
# head(out)
# 
# #let's read in the resistance values from the pairwise nodes
# median_RSF <- read.csv("Outputs/RSF/median_RSF_values.csv")
# head(median_RSF)
# median_RSF <- median_RSF[,c(2,3)]
# colnames(median_RSF) <- c("ID", "median")
# median <- median_RSF$median
# 
# #let's assume mean_RSF (or median) is a vector of 224 values
# for(i in 1:nrow(out)){
#   out$median_RSF[i] <-median[out$incl[i]]
# }
# 
# head(out)
# tail(out)
# 
# max(out$radius) #ok so 602 km is the maximum radius possible
# 
# #for now, let's get rid of rows where orig==incl
# out_new <- out[out$orig!=out$incl,]
# out_test <- subset(out_new, orig==1)
# dim(out)
# dim(out_new)
# 
# #let's use out_new to find all pairs of operational distances from source
# out_pairs <- out_new %>% group_by(orig, radius) %>% slice(1L)
# out_pairs <- out_pairs[,c(1,2)]
# out_pairs$incl <- out_pairs$orig
# for(i in 1:nrow(out_pairs)){
#   out_pairs$median_RSF[i] <-median[out_pairs$incl[i]]
# }
# 
# #now we need to merge out_new and out_pairs
# out_all <- rbind(out_new, out_pairs)
# tail(out_all)
# tail(out_pairs)
# 
# #group them by origin and distance, get number in sequence as seq
# out_seq <- out_all %>% group_by(orig, radius) %>% mutate(id = seq_len(n()))
# max <- max(out_seq$id)
# 
# #divide radii by 100
# max(out_seq$radius)
# out_seq$radius_km <- out_seq$radius/1000
# 
# #max number is 55 neighbors, and max distance is 602 km
# array_RSF <- array(NA, dim=c(224, 602, max))
# array_siteID <- array(NA, dim=c(224, 602, max))
# 
# #for origin site, distance, and number of sites, will produce mean rsf and site number
# for (i in 1:nrow(out_seq)){
#   array_RSF[out_seq$orig[i], out_seq$radius_km[i], out_seq$id[i]] <- as.numeric(out_seq$median_RSF[i])
#   array_siteID[out_seq$orig[i], out_seq$radius_km[i], out_seq$id[i]] <- as.numeric(out_seq$incl[i])
# }
# 
# #also come up w basic matrix re: whether there are territories at that distance 
# site_check <- array_RSF[,,1]
# dim(site_check)
# site_check[site_check>=0] <- 1
# unique(as.vector(as.matrix(site_check)))
# 
# #calculate probabilities of choosing those RSF values
# #do this across rows in the third dimension
# 
# get_probs <- function(array_row) {
#   vector <- array_row/sum(array_row, na.rm=T)
#   return(vector)
# }
# 
# test<- c(1,2,3,4, NA, NA, NA)
# get_probs(test)
# 
# #this will get probabilities of each site and distance
# array_probs <- array(data=NA, dim=c(224,602,55))
# for(i in 1:224){
#   for(j in 1:602){
#     array_probs[i,j,] <- get_probs(array_RSF[i,j,])
#   }}
# 
# dim(array_probs)
# sum(array_probs[1,385,], na.rm=T)
# 
# #this is matching and in site ID
# array_siteID
# 
# save(array_probs, array_siteID, site_check, file= "Outputs/spatial_model/territory_choice.RData")
