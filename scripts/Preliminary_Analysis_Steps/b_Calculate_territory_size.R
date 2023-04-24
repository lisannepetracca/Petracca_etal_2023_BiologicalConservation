library(adehabitatLT)
library(ggplot2)
library(adehabitatHR)
library(tidyverse)
library(sf)
library(units)
library(here)
library(lme4)
library(stringr)
library(lubridate)


#### ---- GETTING GPS DATA INTO PROPER FORMAT TO CONVERT TO MCPs ---- ####

#set working directory
here()

#read in GPS data post-segmentation
#this .csv came from "identify_dispersal_movements_standardized.R"
data <- read.csv("Data/spatial_model/GPS_data_forHR.csv", header=T)
dim(data) #23965
head(data)
tail(data)

#convert time column to POSIXct object & then get year/month
data$date <- as.POSIXct(data$date, tz="America/Los_Angeles")
data$year <- year(data$date)
data$month <- month(data$date)

#get rid of anything after the "." if there are more than one trajectory for a wolf
#as we are interested in all data by each wolf, rather than split by movement trajectory
unique(data$id)
data$id <- gsub("\\..*","", data$id)

#now let's ensure that we are using biological year to calculate HRs (May 2009-April 2010 = 2009)
for (i in 1:nrow(data)){
  if (data$month[i]<=4) {
    data$bioyear[i] <- data$year[i] - 1
  } else data$bioyear[i] <- data$year[i] 
}

#create grouping variable "id" with id and bioyear
data$id <- str_c(data$id,"_",data$bioyear)

#now let's subset to wolves with 60+ points
data_min60 <- data %>% group_by(id) %>% filter(n() > 60)

#and get number of points for each id
data_min60_grouped <- data_min60 %>% group_by(id) %>% dplyr::summarise(n=n()) %>% mutate(bioyear=word(id, 2, sep = "_"))
print(data_min60_grouped,n=100)
#let's transfer to an sf object so that we can write to a shapefile
data_sf <- st_as_sf(data_min60, coords = 
                         c("x", "y"), crs = 32610) 
saveRDS(data_min60, "Outputs/spatial_model/points_for_RSF_consideration.RDS")

#creating shapefiles for each individual wolf-year
# data_sf %>%
#   mutate(group = as.character(id)) %>%
#   group_by(group) %>%
#   nest() %>%
#   #move walk2 inside mutate as data and group were invisible for walk2
#   #also no need for data you can run walk2 directly
#   mutate(#data = map(data, ~st_as_sf(.x)),
#     txt = walk2(.x = data, .y = group, ~st_write(obj = .x, dsn = paste0(.y, ".shp"))))

####---- CALCULATING 99% MCP ----####

#IMPORTANT THAT WE REMOVE 19f (Diobsud)
data_min60 <- data_min60 %>% filter(id!="019f_2013" && id!="019f_2014")
length(unique(data_min60$id))

#converting data to spatial points data frame (package sp)
data_forspdf <- data_min60[,-1]
data_forspdf <- data_forspdf %>% dplyr::select(x,y,id)
spdf <- SpatialPointsDataFrame(coords = data_forspdf[,1:2], data = data_forspdf,
                               proj4string = CRS("+proj=utm +zone=10 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

#calculate 99% MCPs
mcp_99 <- mcp(spdf[,3], percent=99, unin = c("m"), unout = c("km2")) #ignore warning

####---- FILTERING THESE MCPs TO THOSE THAT HAVE >50% OF AREA IN STATE ----####

#bringing in source to add packid to this output
#packid was determined by WDFW monitoring team
packid <- read.csv("Data/spatial_model/packid.csv", header=T)
packid <- packid %>% dplyr::select(burst, pack) %>% dplyr::rename(id = burst)

#reading in WA to filter those that are out of state
WA <- st_read("G:/My Drive/Data/Data/GIS/washington_UTM_mainland.shp")

#let's write a function to project each MCP to UTM, add number of points, bring in wolf pack info, and calculate proportion
#of area overlapping WA state
#we then select those MCPs with >50% area in WA

prepare_HRs <- function(HR) {
  proj4string <- CRS("+proj=utm +zone=10 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")
  HR@proj4string <- proj4string #add coordinate system (UTM Zone 10N)
  HR.sf <- st_as_sf(HR) #convert to sf object
  HR.sf <- left_join(HR.sf, packid, by="id") #add pack identity
  HR.sf <- left_join(HR.sf, data_min60_grouped, by="id") #add number of points & bioyear
  HR.sf <- HR.sf[WA,] #remove those that don't overlap WA at all
  HR.sf$area <- st_area(HR.sf) #get area of MCPs
  HR.sf$area_km2 <- as.numeric(set_units(HR.sf$area, km^2))
  temp <- st_intersection(HR.sf, WA) #intersect those MCPs with WA
  temp$area <- as.numeric(st_area(temp)) #get that area of overlap
  HR.sf$prop <- as.numeric(temp$area/HR.sf$area * 100) #calculate proportion of area in state
  HR.sf <- HR.sf %>% dplyr::filter(prop>=50) %>% group_by(bioyear, pack) %>% #filter by >50% in state
    dplyr::mutate(packshare = cur_group_id()) #assign packshare variable to MCPs sharing pack and year
  return(HR.sf)
}

mcp.99.sf <- prepare_HRs(mcp_99)  #ignore warning

print(mcp.99.sf,n=116)

#let's plot the locations to see where they are
ggplot() +
  geom_sf(data=WA)+
  geom_sf(data = mcp.99.sf, color="blue")

ggplot() +
  geom_sf(data=WA)+
  geom_sf(data = mcp.99.sf_180min, color="blue")

####---- NOW LET'S FIND A MINIMUM NUMBER OF POINTS WHERE THERE IS NO ASSOC BTW AREA AND NUMBER OF POINTS ----####

summary(linear <- lm(area_km2~n, mcp.99.sf))
plot(mcp.99.sf$n, mcp.99.sf$area_km2)
abline(lm(area_km2~n, data=mcp.99.sf))

#removing HRs with <100 points
mcp.99.sf_100min <- mcp.99.sf[mcp.99.sf$n > 100,] 
summary(lm(area_km2~n, mcp.99.sf_100min))
plot(mcp.99.sf_100min$n, mcp.99.sf_100min$area_km2)
abline(lm(area_km2~n, data=mcp.99.sf_100min))

#removing HRs with <180 points seems to do the trick at the individual level
mcp.99.sf_180min <- mcp.99.sf[mcp.99.sf$n > 180,] 
summary(lm(area_km2~n, mcp.99.sf_180min))
plot(mcp.99.sf_180min$n, mcp.99.sf_180min$area_km2) #p-value 0.4115
abline(lm(area_km2~n, data=mcp.99.sf_180min))


####---- LET'S CALCULATE MEAN INDIV HOME RANGE USING INDIVIDUAL AS RE ----####

library(lme4)
#now we can calculate mean HR
fit.MCP.indiv <- lmer(area_km2 ~ (1 | pack), data = mcp.99.sf_180min)
summary(fit.MCP.indiv) #mean of 746.99 with SE of 58.26

#let's write these to shapefile
st_write(mcp.99.sf_180min,
          "Outputs/spatial_model/IndivHRs/MCP_99_territories.shp", driver = "ESRI Shapefile", append=F) #IGNORE WARNINGS

####---- NOW LET'S CALCULATE AREA OF SHARED TERRITORIES----####

#that is, if there is more than one wolf in a given territory in a given bioyear, let's combine those MCPs

#let's write a function that will summarize pack area by combining areas of shared territory
calculate_shared_territory <- function(SF) {
  SF.union <- SF %>%   group_by(packshare) %>% summarise(count = sum(n))
  #annoyingly have to just select first match using plyr
  SF.union.winfo <- plyr::join(as.data.frame(SF.union), as.data.frame(SF), by="packshare", type="left", match="first")
  SF.union[,4:6] <- SF.union.winfo[,c(5,6,8)]
  SF.union <- SF.union[WA,] 
  SF.union$area <- st_area(SF.union)
  SF.union$area_km2 <- as.numeric(set_units(SF.union$area, km^2))
  temp <- st_intersection(SF.union, WA)
  temp$area <- as.numeric(st_area(temp))
  SF.union$prop <- as.numeric(temp$area/SF.union$area * 100)
  SF.union <- SF.union %>% filter(prop>=50)
  return(SF.union)
}

#apply function
MCP_shared <- calculate_shared_territory(mcp.99.sf) #ignore warning

#see what it looks like
ggplot() +
  geom_sf(data=WA)+
  geom_sf(data = MCP_shared, color="blue")

####---- NOW LET'S FIND A MINIMUM NUMBER OF POINTS WHERE THERE IS NO ASSOC BTW AREA AND NUMBER OF POINTS ----####

summary(linear <- lm(area_km2~count, MCP_shared))
plot(MCP_shared$count, MCP_shared$area_km2)
abline(lm(area_km2~count, data=MCP_shared))

#removing HRs with <100 points
MCP_shared_100min <- MCP_shared[MCP_shared$count > 100,] 
summary(lm(area_km2~count, MCP_shared_100min))
plot(MCP_shared_100min$count, MCP_shared_100min$area_km2)
abline(lm(area_km2~count, data=MCP_shared_100min))

#removing HRs with <180 points seems to do the trick
MCP_shared_180min <- MCP_shared[MCP_shared$count > 180,] 
summary(lm(area_km2~count, MCP_shared_180min))
plot(MCP_shared_180min$count, MCP_shared_180min$area_km2) #p-value 0.393
abline(lm(area_km2~count, data=MCP_shared_180min))

####---- LET'S CALCULATE MEAN HOME RANGE USING PACK AS RE ----####

fit.MCP.shared <- lmer(area_km2 ~ (1 | pack), data = MCP_shared_180min)
summary(fit.MCP.shared) #mean of 760.03 with SE of 57.12
length(unique(MCP_shared_180min$pack))
length(unique(MCP_shared_180min$bioyear))

#write these MCPs to shapefile
st_write(MCP_shared_180min,
         "Outputs/spatial_model/SharedHRs/MCP_99_territories_shared.shp", driver = "ESRI Shapefile", append=F) #IGNORE WARNING

saveRDS(MCP_shared_180min, "Outputs/spatial_model/MCP_180min.rds")




















####---- LEGACY CODE ----####

#calculating kernel estimates

# #NOW LET'S DO KERNEL
# memory.limit(size=800000)
# kernel_99_bigger <- kernelUD(spdf[,3], h="href", grid=1000, same4all=FALSE, extent=1.5)
# kernel_99_vertices <- getverticeshr(kernel_99_bigger,99,unin = c("m"),unout = c("km2"))
# kernel_95_vertices <- getverticeshr(kernel_99_bigger,95,unin = c("m"),unout = c("km2"))
# 
# warnings()
# 
# plot(kernel_99_vertices)
# plot(kernel_95_vertices)

# kernel.99.sf <- prepare_HRs(kernel_99_vertices)
# kernel.95.sf <- prepare_HRs(kernel_95_vertices)
# 
# ggplot() +
#   geom_sf(data=WA)+
#   geom_sf(data = kernel.99.sf, color="blue")
# 
# ggplot() +
#   geom_sf(data=WA)+
#   geom_sf(data = kernel.95.sf, color="blue")

#let's look at assoc btw number of points and area for 99 percent kernel
#looks good
# summary(linear <- lm(numpts~area, kernel.99.sf))
# 
# #removing HRs at tails of dist re: numpts (10% was 80, 90% was 322)
# kernel.99.sf_tail <- kernel.99.sf[kernel.99.sf$numpts>80 & kernel.99.sf$numpts< 322,] 
# summary(lm(area~numpts, kernel.99.sf_tail))
# plot(kernel.99.sf_tail$numpts, kernel.99.sf_tail$area)
# abline(lm(area~numpts, data=kernel.99.sf_tail))
# 
# f <- ggplot(kernel.99.sf, aes(numpts, area)) + geom_smooth(method="lm") + geom_point()
# f
# 
# #removing HRs at tails of dist re: numpts (10% was 80, 90% was 322)
# kernel.95.sf_tail <- kernel.95.sf[kernel.95.sf$numpts>80 & kernel.95.sf$numpts< 322,] 
# summary(lm(area~numpts, kernel.95.sf_tail))
# plot(kernel.95.sf_tail$numpts, kernel.95.sf_tail$area)
# abline(lm(area~numpts, data=kernel.95.sf_tail))
# 
# f <- ggplot(kernel.95.sf, aes(numpts, area)) + geom_smooth(method="lm") + geom_point()
# f

#let's write these to shapefile
# st_write(kernel.99.sf_200min,
#          "Outputs/spatial_model/IndivHRs/kernel_99_territories.shp", driver = "ESRI Shapefile", append=F) #IGNORE WARNING
# st_write(kernel.95.sf_200min,
#          "Outputs/spatial_model/IndivHRs/kernel_95_territories.shp", driver = "ESRI Shapefile", append=F) #IGNORE WARNING
#
#apply function
# kernel_99_shared <- calculate_shared_territory(kernel.99.sf_tail)
# kernel_95_shared <- calculate_shared_territory(kernel.95.sf_tail)

# save(mcp.99.sf, kernel.99.sf, kernel.95.sf, 
#      MCP_shared, kernel_99_shared, kernel_95_shared, file = "Data/spatial_model/HR_outputs/data_dispremoved.RData")

#will write these to shapefile
# st_write(kernel_99_shared,
#          "Outputs/spatial_model/SharedHRs/kernel_99_territories_shared.shp", driver = "ESRI Shapefile", append=F) #IGNORE WARNING
# st_write(kernel_95_shared,
#          "Outputs/spatial_model/SharedHRs/kernel_95_territories_shared.shp", driver = "ESRI Shapefile", append=F) #IGNORE WARNING

######------ CALCULATING HRs FOR WHERE NO SEGMENTS REMOVED ------######

#goal was to compare HR size when you remove dispersals and when you don't

# #set working directory
# here()
# 
# #read in GPS data post-segmentation
# data <- read.csv("Data/spatial_model/GPS_data_forHR_noneremoved.csv", header=T)
# dim(data) #23722
# head(data)
# 
# #convert time column to POSIXlt object and then POSIXct object
# data$date <- as.POSIXct(data$date, tz="America/Los_Angeles")
# data$year <- year(data$date)
# 
# ####---- DOING 99% KERNEL AND MCP ----####
# #xy is spatial points
# 
# #LETS START W MCP
# data_forspdf <- data[,-1]
# xy <- data_forspdf[,c(1,2)]
# data_forspdf <- data_forspdf[,c(1,2,12)]
# spdf <- SpatialPointsDataFrame(coords = xy, data = data_forspdf,
#                                proj4string = CRS("+proj=utm +zone=10 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
# 
# mcp_99 <- mcp(spdf[,3], percent=99, unin = c("m"), unout = c("km2"))
# #plot(mcp_99)
# 
# #NOW LET'S DO KERNEL
# memory.limit(size=800000)
# kernel_99_bigger <- kernelUD(spdf[,3], h="href", grid=1000, same4all=FALSE, extent=1.5)
# kernel_99_vertices <- getverticeshr(kernel_99_bigger,99,unin = c("m"),unout = c("km2"))
# kernel_95_vertices <- getverticeshr(kernel_99_bigger,95,unin = c("m"),unout = c("km2"))
# 
# warnings()
# 
# plot(kernel_99_vertices)
# plot(kernel_95_vertices)
# plot(mcp_99)
# 
# #bringing in sources to add pack id and number of points to this output
# packid <- read.csv("Data/spatial_model/packid.csv", header=T)
# head(packid)
# packid <- packid %>% dplyr::select(burst, pack) %>% rename(id = burst)
# 
# numpts <- read.csv("Data/spatial_model/GPS_data_forHR_noneremoved_summary.csv")
# 
# #reading in WA to keep those that are in state
# WA <- st_read("G:/My Drive/Data/Data/GIS/washington_UTM_mainland.shp")
# 
# #let's write a function to project it to UTM, add number of points, and bring in wolf pack info
# prepare_HRs <- function(HR) {
#   proj4string <- CRS("+proj=utm +zone=10 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")
#   HR@proj4string <- proj4string
#   HR.sf <- st_as_sf(HR) 
#   HR.sf <- left_join(HR.sf, packid, by="id")
#   HR.sf$numpts <- numpts$nb.reloc
#   HR.sf$wolf <- word(HR.sf$id, 1, sep = fixed('_'))
#   HR.sf$year <- word(HR.sf$id, 2, sep = fixed('_'))
#   HR.sf <- HR.sf[WA,] 
#   HR.sf <- HR.sf %>% filter(id!="090m_2020") %>% #remove 90m from 2020 because 99% out of state 
#     group_by(year, pack) %>%
#     mutate(packshare = cur_group_id())
#   return(HR.sf)
# }
# 
# mcp.99.sf_allpts <- prepare_HRs(mcp_99) 
# kernel.99.sf_allpts <- prepare_HRs(kernel_99_vertices)
# kernel.95.sf_allpts <- prepare_HRs(kernel_95_vertices)
# 
# print(mcp.99.sf_allpts,n=116)
# print(kernel.99.sf_allpts,n=116)
# print(kernel.95.sf_allpts,n=116)
# 
# kernel.99.sf_allpts <- kernel.99.sf_allpts %>% filter(wolf!="OR49f") #remove bc largely out of state
# 
# print(kernel.99.sf_allpts,n=116)
# 
# #let's plot the locations to see where they are
# ggplot() +
#   geom_sf(data=WA)+
#   geom_sf(data = mcp.99.sf_allpts, color="blue")
# 
# ggplot() +
#   geom_sf(data=WA)+
#   geom_sf(data = kernel.99.sf_allpts, color="blue")
# 
# ggplot() +
#   geom_sf(data=WA)+
#   geom_sf(data = kernel.95.sf_allpts, color="blue")
# 
# #let's look at association btw number of points and area for MCP
# #raw data
# summary(linear <- lm(area~numpts, mcp.99.sf_allpts))
# plot(mcp.99.sf_allpts$numpts, mcp.99.sf_allpts$area)
# abline(lm(area~numpts, data=mcp.99.sf_allpts))
# 
# quantile(mcp.99.sf_allpts$numpts, probs = seq(0, 1, 0.05))
# 
# #removing HRs with <100 points
# mcp.99.sf_allpts_100min <- mcp.99.sf_allpts[mcp.99.sf_allpts$numpts > 100,] 
# summary(lm(area~numpts, mcp.99.sf_allpts_100min))
# plot(mcp.99.sf_allpts_100min$numpts, mcp.99.sf_allpts_100min$area)
# abline(lm(area~numpts, data=mcp.99.sf_allpts_100min))
# 
# #removing HRs with <150 points
# mcp.99.sf_allpts_150min <- mcp.99.sf_allpts[mcp.99.sf_allpts$numpts > 150,] 
# summary(lm(area~numpts, mcp.99.sf_allpts_150min))
# plot(mcp.99.sf_allpts_150min$numpts, mcp.99.sf_allpts_150min$area)
# abline(lm(area~numpts, data=mcp.99.sf_allpts_150min))
# 
# #removing HRs with <200 points
# mcp.99.sf_allpts_200min <- mcp.99.sf_allpts[mcp.99.sf_allpts$numpts > 200,] 
# summary(lm(area~numpts, mcp.99.sf_allpts_200min))
# plot(mcp.99.sf_allpts_200min$numpts, mcp.99.sf_allpts_200min$area)
# abline(lm(area~numpts, data=mcp.99.sf_allpts_200min))
# 
# #removing HRs at tails of dist re: numpts (10% was 80, 90% was 322)
# mcp.99.sf_allpts_tail <- mcp.99.sf_allpts[mcp.99.sf_allpts$numpts>96 & mcp.99.sf_allpts$numpts< 356,] 
# summary(lm(area~numpts, mcp.99.sf_allpts_tail))
# plot(mcp.99.sf_allpts_tail$numpts, mcp.99.sf_allpts_tail$area)
# abline(lm(area~numpts, data=mcp.99.sf_allpts_tail))
# 
# #getting unique number of wolves and packs
# length(unique(mcp.99.sf_allpts$wolf))
# length(unique(mcp.99.sf_allpts$pack))
# 
# f <- ggplot(mcp.99.sf_allpts, aes(numpts, area)) + geom_smooth(method="lm") + geom_point()
# f
# f <- ggplot(mcp.99.sf_allpts_tail, aes(numpts, area)) + geom_smooth(method="lm") + geom_point()
# f
# 
# #let's look at assoc btw number of points and area for 99 percent kernel
# #looks good
# summary(linear <- lm(numpts~area, kernel.99.sf_allpts))
# 
# #removing HRs at tails of dist re: numpts (10% was 80, 90% was 322)
# kernel.99.sf_allpts_tail <- kernel.99.sf_allpts[kernel.99.sf_allpts$numpts>80 & kernel.99.sf_allpts$numpts< 322,] 
# summary(lm(area~numpts, kernel.99.sf_allpts_tail))
# plot(kernel.99.sf_allpts_tail$numpts, kernel.99.sf_allpts_tail$area)
# abline(lm(area~numpts, data=kernel.99.sf_allpts_tail))
# 
# f <- ggplot(kernel.99.sf_allpts, aes(numpts, area)) + geom_smooth(method="lm") + geom_point()
# f
# 
# #removing HRs at tails of dist re: numpts (10% was 80, 90% was 322)
# kernel.95.sf_allpts_tail <- kernel.95.sf_allpts[kernel.95.sf_allpts$numpts>80 & kernel.95.sf_allpts$numpts< 322,] 
# summary(lm(area~numpts, kernel.95.sf_allpts_tail))
# plot(kernel.95.sf_allpts_tail$numpts, kernel.95.sf_allpts_tail$area)
# abline(lm(area~numpts, data=kernel.95.sf_allpts_tail))
# 
# f <- ggplot(kernel.95.sf_allpts, aes(numpts, area)) + geom_smooth(method="lm") + geom_point()
# f
# 
# #will write these to shapefile
# # st_write(mcp.99.sf_allpts_300min,
# #          "Outputs/spatial_model/IndivHRs/MCP_99_territories.shp", driver = "ESRI Shapefile", append=F) #IGNORE WARNING
# # st_write(kernel.99.sf_allpts_200min,
# #          "Outputs/spatial_model/IndivHRs/kernel_99_territories.shp", driver = "ESRI Shapefile", append=F) #IGNORE WARNING
# # st_write(kernel.95.sf_allpts_200min,
# #          "Outputs/spatial_model/IndivHRs/kernel_95_territories.shp", driver = "ESRI Shapefile", append=F) #IGNORE WARNING
# 
# 
# #let's write a function that will summarize pack area by combining areas of shared territory
# calculate_shared_territory <- function(SF) {
#   SF.union <- SF %>%   group_by(packshare) %>%   summarise()
#   #annoyingly have to just select first match using plyr
#   SF.union.winfo <- plyr::join(as.data.frame(SF.union), as.data.frame(SF), by="packshare", type="left", match="first")
#   #now can add to kernel.99.sf_allpts.200min.union
#   SF.union[,3:6] <- SF.union.winfo[,c(5,7:9)]
#   area_m2 <- st_area(SF.union)
#   SF.union$area_km2 <- as.numeric(set_units(area_m2, km^2))
#   return(SF.union)
# }
# 
# #apply function
# MCP_shared_allpts <- calculate_shared_territory(mcp.99.sf_allpts)
# kernel_99_shared_allpts <- calculate_shared_territory(kernel.99.sf_allpts)
# kernel_95_shared_allpts <- calculate_shared_territory(kernel.95.sf_allpts)
# 
# #just for fun
# MCP_shared_all <- calculate_shared_territory(mcp.99.sf_allpts)
# 
# save(mcp.99.sf_allpts, kernel.99.sf_allpts, kernel.95.sf_allpts, 
#      MCP_shared_allpts, kernel_99_shared_allpts, kernel_95_shared_allpts, file = "Data/spatial_model/HR_outputs/data_noneremoved.RData")
# load("Data/spatial_model/HR_outputs/data_noneremoved.RData")
# load("Data/spatial_model/HR_outputs/data_dispremoved.RData")
# 
# print(mcp.99.sf,n=100)
# print(mcp.99.sf_allpts,n=100)
# 
# mcps_disp <- as.data.frame(mcp.99.sf)
# 
# #replacing the .1 and the .2 with nothing
# mcps_disp$id<-gsub(".1_","_",mcps_disp$id)
# mcps_disp$id<-gsub(".2_","_",mcps_disp$id)
# 
# mcps_nondisp <- as.data.frame(mcp.99.sf_allpts)
# mcps_join <- left_join(mcps_disp, mcps_nondisp, by="id")
# mcps_join$diff <- mcps_join$area.y - mcps_join$area.x
# 
# mcps_join <- mcps_join %>% arrange(diff)
# mcps_join$row <- 1:116
# 
# p <- ggplot(mcps_join, aes(row, diff))
# p + geom_point() +  
#   scale_x_continuous(breaks=c(1:116),
#                    labels=unique(mcps_join$id))+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ylab("Difference in individual HR w and w/o dispersal (km2)")+
#   xlab("Wolf and Year")
# ggsave("G:/My Drive/Data/Data/Figures/HR_Diff_w_wo_Dispersal_detail.jpg")
# 
# quantile(mcps_join$diff, probs=seq(0,1,by=0.05))
# 
# #violin plot
# p <- ggplot(mcps_join, aes(1, diff))
# p + geom_violin() + geom_jitter(height = 0, width = 0.1)+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())+
#   ylab("Difference in individual HR w and w/o dispersal (km2)")
# ggsave("G:/My Drive/Data/Data/Figures/HR_Diff_w_wo_Dispersal.jpg")
# 
# #check to ensure non-dispersers same in both datasets
# disp.82f_2019 <- mcp.99.sf %>% filter(id=="082f_2019")
# nondisp.82f_2019 <- mcp.99.sf_allpts %>% filter(id=="082f_2019")
# 
# plot(disp.82f_2019)
# plot(nondisp.82f_2019)
# 
# ######------ CALCULATING MEANS WHEN USING INDIV OR TERRITORY AS RANDOM EFFECT ------######
# 
# # MCP_indiv <- st_read("Outputs/spatial_model/IndivHRs/MCP_99_territories.shp")
# # kernel_99_indiv <- st_read("Outputs/spatial_model/IndivHRs/kernel_99_territories.shp")
# # kernel_95_indiv <- st_read("Outputs/spatial_model/IndivHRs/kernel_95_territories.shp")
# # 
# # MCP_shared <- st_read("Outputs/spatial_model/SharedHRs/MCP_99_territories_shared.shp")
# # kernel_99_shared <- st_read("Outputs/spatial_model/SharedHRs/kernel_99_territories_shared.shp")
# # kernel_95_shared <- st_read("Outputs/spatial_model/SharedHRs/kernel_95_territories_shared.shp")
# 
# #for indiv calcs
# fit.MCP.indiv <- lmer(area ~ (1 | wolf), data = mcp.99.sf_tail)
# summary(fit.MCP.indiv)
# 
# fit.kernel99.indiv <- lmer(area ~ (1 | wolf), data = kernel.99.sf_tail)
# summary(fit.kernel99.indiv)
# 
# fit.kernel95.indiv <- lmer(area ~ (1 | wolf), data = kernel.95.sf_tail)
# summary(fit.kernel95.indiv)
# 
# #for shared calcs
# 
# fit.MCP.shared <- lmer(area_km2 ~ (1 | pack), data = MCP_shared)
# summary(fit.MCP.shared)
# 
# fit.kernel99.shared <- lmer(area_km2 ~ (1 | pack), data = kernel_99_shared)
# summary(fit.kernel99.shared)
# 
# fit.kernel95.shared <- lmer(area_km2 ~ (1 | pack), data = kernel_95_shared)
# summary(fit.kernel95.shared)
# 
# median(MCP_shared$area_km2)
# 
# #for fun
# fit.MCP.shared_all <- lmer(area_km2 ~ (1 | pack), data = MCP_shared_all)
# summary(fit.MCP.shared_all)
