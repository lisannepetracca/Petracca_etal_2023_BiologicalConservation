#---
#"Calculate territory size for wolves in Washington State"
#"Lisanne Petracca"
#"November 2023"
#---

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
#this .csv came from "a_Identify_dispersal_movements.R"
data <- read.csv("C:/Your_Directory/GPS_data_forHR.csv", header=T)

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
saveRDS(data_min60, "C:/Your_Directory/points_for_RSF_consideration.RDS")

####---- CALCULATING 99% MCP ----####

#IMPORTANT THAT WE REMOVE 19f BC ALL DISPERSAL (Diobsud)
data_min60 <- data_min60 %>% dplyr::filter(id!="019f_2013") 
data_min60 <- data_min60 %>% dplyr::filter(id!="019f_2014")
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
packid <- read.csv("C:/My_Directory/spatial_model/packid.csv", header=T)
packid <- packid %>% dplyr::select(burst, pack) %>% dplyr::rename(id = burst)

#reading in WA to filter those that are out of state
WA <- st_read("C:/My_Directory/washington_UTM_mainland.shp")

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
          "C:/Your_Directory/MCP_99_territories.shp", driver = "ESRI Shapefile", append=F) #IGNORE WARNINGS

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
         "C:/My_Directory/MCP_99_territories_shared.shp", driver = "ESRI Shapefile", append=F) #IGNORE WARNING

saveRDS(MCP_shared_180min, "C:/My_Directory/MCP_180min.rds")