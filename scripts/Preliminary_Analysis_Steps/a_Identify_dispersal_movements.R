
library(adehabitatLT)
library(ggplot2)
library(adehabitatHR)
library(tidyverse)
library(sf)
library(units)
library(raster)
library(cowplot)
library(here)
library(lubridate)
library(stringr)
library(raster)
library(data.table)

######------ LETS READ IN ALL PROCESSED DATA ------######

#set working directory
here()

#read in data, set data as POSIXct class, and generate columns for year, hour, month, and bioyear
data <- read.csv("Data/collar_data/wolf_gpsdata_all_processed.csv", header=T)
data$date <- as.POSIXct(data$date, tz="America/Los_Angeles")
data$year <- year(data$date)
data$hour <- hour(data$date)
data$month <- month(data$date)
data$bioyear <- NA

#coding bioyear (2009 is May 2009-April 2010)
data <- data %>%  mutate(bioyear = ifelse(month<=4, year-1, year))

head(data)
data[10000:10100,]

#create column for wolfyear (concatenating id and bioyear)
data$wolfyear <- str_c(data$id,"_",data$bioyear)

#format date
data$mdy <- format(data$date, "%Y-%m-%d")
dim(data)
head(data)
#start with 138487 rows

#let's subset to daily locs & then get locations closest to 7 am every day
data_daily <- data %>% group_by(id, mdy) %>% filter(abs(hour-7)==min(abs(hour-7))) %>% ungroup()
head(data_daily)
tail(data_daily)

#format date to just have m/d/y
data_daily$date <- as.POSIXct(format(as.POSIXct(data_daily$date,format='%Y-%m-%d" %H:%M:%S'),format='%Y-%m-%d'))

#remove duplicates for multiple fixes in same hour
#will keep the first entry for each day
data_daily <- data_daily %>% group_by(id,date) %>% slice(1) %>% ungroup()

#remove bursts with <60 days of data
data_daily <- data_daily %>% group_by(id) %>% filter(n()>=60)

#conversion to class ltraj
finaldata <- as.ltraj(xy = data_daily[,c("x","y")], date = data_daily$date, id=data_daily$id, burst=data_daily$id)
sum(summary(finaldata)$nb.reloc)
#now 31613 points #reduced to 31607 when the trajectory with 6 pts is removed

#convert to dataframe
dataframe <- ld(finaldata)

#let's write out this .csv for use in the survival model detection histories
write.csv(dataframe, "Data/survival/GPS_daily_data_for_survival.csv")

#now 31613 points, was 137631 points


#let's get max NSD by wolf here
#we end up using 50th percentile for max NSD to subset to potential dispersers

nsd <- list()
for(i in 1:74){
  nsd[i] <- max(finaldata[[i]]$R2n)
}
nsd <- as.vector(Reduce(rbind, nsd))
plot(nsd)
quantile(nsd, probs = seq(0,1,0.05))

length(which(nsd>11974379217, arr.ind=T))
length(which(nsd>5155952369, arr.ind=T)) #65th percentile #adds 21,36,37,38,43,44,69
length(which(nsd>3909502334, arr.ind=T)) #60th percentile #adds 18,30,31,40
length(which(nsd>3534604396, arr.ind=T)) #55th percentile #adds 8,17,27
length(which(nsd>2616581604, arr.ind=T)) #50th percentile #adds 16,22,42,64
length(which(nsd>1966890528, arr.ind=T)) #45th percentile #adds 10,13,34,55 (no bueno)

sqrt(2616581604) #equates to ~50 km

#dispersers
dispersers <- which(nsd>2616581604, arr.ind=T)
non_dispersers <- which(nsd<2616581604, arr.ind=T)
disp_data <- finaldata[dispersers]
nondisp_data <- finaldata[non_dispersers]

#remove those that entered from out of state
#disp_data <- disp_data[-c(35:37)]
#summary(disp_data)

sum(summary(disp_data)$nb.reloc) #19537 points from dispersers

#break down trajectories when there is gap of 14+ days
#60 * 60 * 24 * 14
foo2 <- function(dt) 
{
  return(dt> (2*1209600))
}

#removal of ONE point (yay!); go from 34 trajectories to 40 #Dec 6 I see 43
#total of 138487 points

dispdata_nogaps <- cutltraj(disp_data, "foo2(dt)", nextr = TRUE)
summary(dispdata_nogaps) 

nondispdata_nogaps <- cutltraj(nondisp_data, "foo2(dt)", nextr = TRUE)
summary(nondispdata_nogaps) 

alldata_nogaps <- rbind(ld(dispdata_nogaps), ld(nondispdata_nogaps))
alldata_nogaps$id <- alldata_nogaps$burst

alldata_nogaps <- dl(alldata_nogaps)
summary(alldata_nogaps)
sum(summary(alldata_nogaps)$nb.reloc) #31613 points from dispersers

# write.csv(dataframe, "Data/spatial_model/all_points_for_segmentation.csv")
# write.csv(summary, "Data/spatial_model/all_points_for_segmentation_summary.csv")

######------ DOING ALL DATA TOGETHER ------######

#important to have regularized trajectories for FPT, so let's do that with the dispersers
#now we can fill in missing days with setNA
refda <- all_data_regular <- list()
for(i in 1:84){
  refda[[i]] <-summary(alldata_nogaps)$date.begin[i]
  all_data_regular[i] <- setNA(alldata_nogaps[i], date.ref=refda[[i]][1],1,units="day")
  all_data_regular[[i]]$id <- summary(alldata_nogaps)$id[i]
  all_data_regular[[i]]$burst <- summary(alldata_nogaps)$burst[i]
}

all_data_regular <- as.data.frame(Reduce(rbind, all_data_regular))
all_data_regular$date <- as.POSIXct(round(all_data_regular$date, unit="day"))

all_data_regular <- dl(all_data_regular)

##***NEW LINE TO MAKE SHAPEFILE
id(all_data_regular) <- burst(all_data_regular)
summary <- summary(all_data_regular)

#remove burst with 6 relocs (87m)
all_data_regular <- all_data_regular[-which(summary$burst=="087m.1")]
summary(all_data_regular)
sum(summary(all_data_regular)$nb.reloc) 

#check that data separated by 1 day
#seems there is a v minor issue with use of daylight savings time but let's ignore
df <- ld(all_data_regular)
hist(df$dt)
#hist(disp_data_regular[1:10],"dt",freq=TRUE)

######------ WRITING DATA FROM INDIVIDUAL WOLVES & WOLF-YEARS TO SHAPEFILE  ------######

# df <- df[!is.na(df$x),]
# 
# disp_wolflocs_sf <- st_as_sf(df, coords =
#                                c("x", "y"), crs = 32610)
# head(disp_wolflocs_sf)

#writing as individual wolf
# disp_wolflocs_sf %>%
#   mutate(group = as.character(id)) %>%
#   group_by(group) %>%
#   nest() %>%
#   #move walk2 inside mutate as data and group were invisible for walk2
#   #also no need for data you can run walk2 directly
#   mutate(#data = map(data, ~st_as_sf(.x)),
#     txt = walk2(.x = data, .y = group, ~st_write(obj = .x, dsn = paste0(.y, ".shp"))))
# 
# #writing as individual wolf-year
# wolflocs_sf %>%
#   mutate(group = as.character(burst)) %>%
#   group_by(group) %>%
#   nest() %>%
#   #move walk2 inside mutate as data and group were invisible for walk2
#   #also no need for data you can run walk2 directly
#   mutate(#data = map(data, ~st_as_sf(.x)),
#     txt = walk2(.x = data, .y = group, ~st_write(obj = .x, dsn = paste0(.y, ".shp"))))

######------ START OF DATA SEGMENTATION USING LAVIELLE METHOD  ------######

fpt <- variance <- scale <- df <- ltraj <- list()

for (i in 1:83) {
  fpt[[i]] <- fpt(all_data_regular[i], seq(1000, 100000, by=1000), units="days")
  variance[[i]] <- varlogfpt(fpt[[i]])
  scale[[i]] <- as.numeric(apply(variance[[i]],1,which.max))
}

#getting median scale at which wolves are interacting with their environment
scale <- as.data.frame(Reduce(rbind, scale))
median(scale$V1, na.rm=T)
#ends up being 18

for (i in 1:83){
  df[[i]] <- ld(all_data_regular[i])
  df[[i]] <- df[[i]][!is.na(df[[i]]$x),]
  df[[i]] <- df[[i]][!is.na(df[[i]]$y),]
  df[[i]]$fpt_18 <- fpt[[i]][[1]]$r18
}
df <- as.data.frame(Reduce(rbind, df))

ltraj <- dl(df)
id(ltraj) <- burst(ltraj)
summary(scale$V1)

plot(ltraj[75])
#29 is only problem w fpt_14
#10, 29, 31, 64, 75 are problems w fpt_18

ltraj_qualify <- ltraj[-c(10,29,31,64,75)]
ltraj_notqualify <- ltraj[c(10,29,31,64,75)]

lav_wolf <- list()
for(i in 1:78){
  lav_wolf[[i]] <- lavielle(ltraj_qualify[i], Lmin=2, Kmax=8, type="mean", which="fpt_18")
}

#rule for getting K
#example where my method works better is lav_wolf[[2]]
#Laveille paper even says that 0.75 not great w observations <500
choice_K_Laveille <- function(lav) {
  temp <- chooseseg(lav)
  K <- max(which(temp$D>.75))
  return(K)
}

choice_K_Lisanne <- function(lav) {
  progress <- list()
  temp <- chooseseg(lav)
  range <- temp$Jk[1] - temp$Jk[7]
  for(i in 1:7){
    progress[i] <- (temp$Jk[1]-temp$Jk[i])/range
  }
  vec <- as.vector(Reduce(rbind, progress))
  K <- min(which(vec>0.8))
  return(K)
}

#applying function to list
K1 <- as.vector(Reduce(rbind, lapply(lav_wolf, choice_K_Laveille)))
K2 <- as.vector(Reduce(rbind, lapply(lav_wolf, choice_K_Lisanne)))
table(K2 - K1)

#compared to Laveille rule of thumb, my version
#agrees 30/78 times 
#overestimates by 1 14 times
#overestimates by 2 21 times
#overestimates by 3 6 times
#overestimates by 4 3 times
#overestimates by 5 2 time
#underestimates by 1 2 times

######------ CAN SKIP THIS SECTION, THIS WAS FOR INDIV CHECKS ------######

summary(ltraj_qualify)
plot(ltraj_qualify[25])
 #let's start
seg_4370 <- findpath(lav_wolf[[57]], K2[[57]])

plot(seg_4370)
plot(seg_4370[1])
plot(seg_4370[2])
plot(seg_4370[3])
plot(seg_4370[4])
plot(seg_4370[5])
plot(seg_4370[6])
plot(seg_4370[7])
#identifying segments to keep 
K2
######------ SKIPPED SECTION ENDS HERE ------######

#####----- SEGMENTATION OF PATHS AND CALCULATING NSD FROM FIRST SEGMENT ----####

path <- df <- list()
for (i in 1:78){
  path[[i]] <- findpath(lav_wolf[[i]], K2[[i]])
  df[[i]] <- ld(path[[i]])
}
df <- as.data.frame(Reduce(rbind, df))
head(df)

#first we need to get NSD by id rather than by segment
df <- df %>% group_by(id) %>% 
  mutate(NSD = (first(x) - x)^2 +(first(y) - y)^2) %>% ungroup()

#ensure that this calculation of NSD is correct when compared to R2n
df$R2n[1:3]
df$NSD[1:3]

#"segment" retrieves only the segment number from the "burst" string
df <- df %>% mutate(segment = as.numeric(gsub("Segment.", "", df$burst)))
head(df)

#calculate median fpt and max nsd along each segment
df <- df %>% group_by(id, segment) %>% 
  summarize(median_fpt = median(fpt_18, na.rm = TRUE),
            max_nsd = max(NSD, na.rm = TRUE))
print(df,n=100)

#looking at quantiles
#the quantiles that best separate are 50% FST, 75% NSD
quantile(df$max_nsd, seq(0,1,0.05)) #75th percentile is 2157521728 
quantile(df$median_fpt, seq(0,1,0.05)) #50th percentile is 26.574252  

#selecting those segments that meet our criteria
select <- df %>% filter(max_nsd>2157521728   & median_fpt<26.574252)
print(select,n=70)

#these create visual plots of each wolf's movement trajectories by NSD and FPT
#these plots are then written to file
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)

head(df)
plot_list1 <- plot_list2 <- plot_list3 <- list()
for (i in unique(df$id)){
  #this creates the FPT plot
  plot_fpt<- ggplot(data = df[df$id==i,], aes(y = median_fpt, x=segment)) +
    geom_point(size=3,col="black") + geom_line() +
    facet_wrap(~id)+
    theme_bw() + geom_hline(yintercept=26.574252  ,linetype="dashed", #20th is 11.422751, this is 50th
                            color = "red")
  #this creates the NSD plot
  plot_nsd<- ggplot(data = df[df$id==i,], aes(y = max_nsd, x=segment)) +
    geom_point(size=3,col="black") + geom_line() +
    facet_wrap(~id)+
    theme_bw() + geom_hline(yintercept=2157521728          ,linetype="dashed",  #75th 
                            color = "red")
  plot_list1[[i]]=plot_fpt
  plot_list2[[i]]=plot_nsd
  
  #then these plots are displayed together as plot_list3
  plot_list3[[i]]=plot_grid(plot_list1[[i]], plot_list2[[i]], 
                            labels = c("FPT", "NSD"),
                            ncol = 1, nrow = 2)
}

for (i in 1:length(unique(df$id))) {
  file_name = paste(i,"_",unique(df$id),".png")
  png(paste("G:/My Drive/GitHub/Wolves/Data/spatial_model/FPT_NSD_plots_all_5075","/",file_name[i],sep=""))
  grid.draw(plot_list3[[i]])
  dev.off()
}


#####----- KEEPING TRAJECTORIES FOR HR ANALYSIS ------#####

summary(path)
paths_HRonly <- path[-9] #this entire wolf removed bc all dispersal
paths_occupied <- path

#segments to keep as nondispersal segments
keep <- list(c(1:5), c(1:3), c(1:4), c(1:5), c(1:6),
             c(1:2), c(1,3), c(1:2), c(1:3), c(1:3), #KEEPING IN FINAL ONE (19f) BUT WILL BE REMOVED LATER
             c(1:3), c(1:5), c(1:3), c(1:5), c(1:4),
             c(1:5), c(1:4), c(1:3), c(1:2), c(1,2,4), 
             c(1:4), c(1:3), c(1:4), c(1:6), c(2:5), 
             c(1,2,4), c(1:4), c(1:5), c(1:5), c(1,2,4,5,6),
             c(1:5), c(2:3), c(1:3), c(1:3), c(1:3),
             c(1,4), c(1:4), c(1:6), 1, c(1,3),
             c(2,3,4), c(1,3), c(1:5), c(1,2,3,5), 3,
             c(1,3,4), c(1:4), c(1,3), c(2,3), 1,
             c(1:4), c(1:4), c(1:6), c(1:4), c(1:5),
             c(1,2,4), c(1:4), c(1:3), c(1:3), c(1:6),
             2, c(2,4), c(1:3), c(1:3), c(1:2),
             c(1:4), c(1:5), c(1:5), c(1,2,4), c(1:5),
             c(1:5), c(1:6), c(2:4), c(1:6), c(1:2),
             c(1,2,4), c(2,3), c(2:4))

#create df of "HR-only" trajectories (still has 19f for now)
ltraj_occupied <- df_occupied <- list()
for (i in 1:78) {
  ltraj_occupied[[i]] <- paths_occupied[[i]][c(keep[[i]])]
  df_occupied[[i]] <- ld(ltraj_occupied[[i]])
}

#####----- FOR OCCUPANCY STUFF -----#####

#add in the five that were excluded via FPT and we are golden
#this step reduces list to df
df_occupied_all <- Reduce(rbind, df_occupied)
#make the df an ltraj
ltraj_occupied_all <- as.ltraj(xy = df_occupied_all[,c("x","y")], date = df_occupied_all$date, id=df_occupied_all$id, burst=df_occupied_all$id)
ltraj_occupied_all <- removeinfo(ltraj_occupied_all) #removes infolocs
ltraj_notqualify <- removeinfo(ltraj_notqualify) #removes infolocs
ltraj_occupied_all[79:83] <- ltraj_notqualify

#let's convert this ltraj to a df and write to .csv for HR analysis
df_occupied_all <- ld(ltraj_occupied_all)
write.csv(df_occupied_all, "Data/spatial_model/GPS_data_forHR.csv")
length(unique(df_occupied_all$id))

#####----- KEEPING TRAJECTORIES FOR MOVEMENT PERIODS -----#####

#in order to get the durations for which wolves are moving, we need to remove non-movers first

paths_movers <- paths_HRonly[-c(1:5, 10:12, 14, 15, 17, 18, 20, 23,
                          27, 28, 30, 32, 33, 34, 36, 37,
                          42, 46, 50:54, 57:59, 62, 63, 
                          65:67, 69:71, 73)]
#these are the movement segments for these wolves
keep <- list(3, c(2,4), 3, 4, 
              6,5,3,4,5,
              1, c(3,5), 5, 3, 1,
              c(2,3,5), 2, c(2,4), 1, c(2,4),
              c(4,6), c(1,2,4), 2, 2, c(1,4),
              c(2:5), 3, 5, c(1,3), c(1,3,5),
              3, c(3,5), c(1,5,6), 3, 3,
              1, c(1,5))
  
#this step creates a df of the segments w movement
ltraj_movers <- df_movers <- list()
for (i in 1:36) {
  ltraj_movers[[i]] <- paths_movers[[i]][c(keep[[i]])]
  df_movers[[i]] <- ld(ltraj_movers[[i]])
}
#and reduces df list to simple df
df_movers_all <- Reduce(rbind, df_movers)

#let's create an id column with wolf id and unique movement period
df_movers_all$burst <- sub('.*\\.', '', df_movers_all$burst)
df_movers_all$pathid <- str_c(df_movers_all$id, df_movers_all$burst, sep="_")

#create ltraj
ltraj_movers_all <- as.ltraj(xy = df_movers_all[,c("x","y")], date = df_movers_all$date, id=df_movers_all$id, burst=df_movers_all$pathid)
#adding 19f back in bc this wolf apparently in movement state the entire time
ltraj_movers_all <- removeinfo(ltraj_movers_all) #remove infolocs
ltraj_qualify <- removeinfo(ltraj_qualify) #remove infolocs
ltraj_movers_all[57] <- ltraj_qualify[9]
movers_df <- ld(ltraj_movers_all)

summary_movers <- summary(ltraj_movers_all)
write.csv(summary_movers, "Data/survival/movement_data_for_survival.csv")
write.csv(movers_df, "Data/survival/movement_data_for_distancecalc.csv")

######------ CALCULATING MOVEMENT DISTANCE ------######

library(tidyverse)
movements <- read.csv("Data/survival/movement_data_for_distancecalc.csv", header=T)
summary <- read.csv("Data/survival/movement_data_for_survival.csv", header=T)
dim(summary)
dim(movements)
length(unique(movements$burst))
length(unique(movements$id))

#first we need to remove the three OR wolves bc we aren't interested in those distances
summary <- summary %>% filter(!burst %in% c("OR15m.1_3", "OR35f.1_1",
                                            "OR49f.1_1", "OR49f.1_5"))
movements <- movements %>% filter(!burst %in% c("OR15m.1_3", "OR35f.1_1",
                                                "OR49f.1_1", "OR49f.1_5"))
unique(summary$id[16])
unique(summary$id[22])

#convert to sf and take first/last coordinate by wolf
library(sf)
move_sf <- st_as_sf(movements, coords =
                      c("x", "y"), crs = 32610)
head(summary)
first_pts <- move_sf %>% group_by(id) %>% slice_head() %>% st_coordinates()
last_pts <- move_sf %>% group_by(id) %>% slice_tail() %>% st_coordinates()#whether individual wolves moved in or out

#getting distance
library(raster)
distances <- pointDistance(first_pts, last_pts, lonlat=F, allpairs=F)
distances_km = as.data.frame(distances/1000)

distances_km$location <- c("out", "out", "out", "in", "in",
                      "in", "in", "out","out","in", 
                      "in", "out", "in", "out", "in",
                      "out", "in", "in", "in", "out",
                      "in", "out","in", "out", "in",
                      "out","out", "out","out","out",
                      "out","in","in","out")

distances_km_instate <- distances_km %>% filter(location=="in")
distances_km_all <- distances_km
dim(distances_km_instate)
hist(distances_km_instate$`distances/1000`)
length(distances_km_instate$`distances/1000`)
summary(distances_km_instate$`distances/1000`)
max(distances_km_instate$`distances/1000`)
distances <- distances_km_instate$`distances/1000`
write.csv(distances, "Outputs/spatial_model/instate_distances.csv")

hist(distances_km_all$`distances/1000`)
length(distances_km_all$`distances/1000`)
summary(distances_km_all$`distances/1000`)
max(distances_km_all$`distances/1000`)
distances <- distances_km_all$`distances/1000`
write.csv(distances, "Outputs/spatial_model/all_distances.csv")

#random thing getting ids of dispersers
disperser_id <- move_sf %>% group_by(id) %>% slice(1L) %>% dplyr::select(date, id)
write.csv(disperser_id, "Outputs/spatial_model/disperser_id.csv")

#####----- OLD CODE STARTS HERE ----####

######------ PRODUCING FIGURE SHOWING DURATION OF MOVEMENT ------######
#THIS DIMENSIONALITY NO LONGER WORKS WITH DIMENSIONS OF THAT LTRAJ SUBSET BY BURST

library(RColorBrewer)
#for this figure, we are not interested in movements into the state by the three OR wolves
summary_movers <- summary(ltraj_movers_all[-c(34:36)])
summary_movers$location <- c("out", "out", "out", "in",
                           "in", "in", "out","out","in", 
                           "in", "out", "in", "in", "in",
                           "out", "in", "in", "in", "out",
                           "in", "out","in", "out", "in",
                           "out","out", "out","out","out",
                           "out","in","in","out", 
                           "in")
summary_movers$duration <- summary_movers$date.end - summary_movers$date.begin
c <- ggplot(summary_movers, aes(x=duration, color=location))
c + geom_freqpoly(size=1) + xlab("days of dispersal")+
  scale_color_brewer(palette = "Set1") + facet_grid(location~1)+
  theme(text = element_text(size=20))
ggsave("G:/My Drive/Data/Data/Figures/Dispersal_Times_byLocation.jpg")

######------ DOING JUST DISPERSERS ------#######

#important to have regularized trajectories for FPT, so let's do that with the dispersers
#now we can fill in missing days with setNA
refda <- disp_data_regular <- list()
for(i in 1:43){
  refda[[i]] <-summary(dispdata_nogaps)$date.begin[i]
  disp_data_regular[i] <- setNA(dispdata_nogaps[i], date.ref=refda[[i]][1],1,units="day")
  disp_data_regular[[i]]$id <- summary(dispdata_nogaps)$id[i]
  disp_data_regular[[i]]$burst <- summary(dispdata_nogaps)$burst[i]
}

disp_data_regular <- as.data.frame(Reduce(rbind, disp_data_regular))
disp_data_regular$date <- as.POSIXct(round(disp_data_regular$date, unit="day"))

disp_data_regular <- dl(disp_data_regular)

##***NEW LINE TO MAKE SHAPEFILE
id(disp_data_regular) <- burst(disp_data_regular)
summary <- summary(disp_data_regular)

#remove burst with 6 relocs (87m)
disp_data_regular <- disp_data_regular[-which(summary$burst=="087m.1")]
summary(disp_data_regular)

#check that data separated by 1 day
#seems there is a v minor issue with use of daylight savings time but let's ignore
df <- ld(disp_data_regular)
hist(df$dt)
#hist(disp_data_regular[1:10],"dt",freq=TRUE)

######------ WRITING DATA FROM INDIVIDUAL WOLVES & WOLF-YEARS TO SHAPEFILE  ------######

df <- df[!is.na(df$x),]

disp_wolflocs_sf <- st_as_sf(df, coords =
                          c("x", "y"), crs = 32610)
head(disp_wolflocs_sf)

#writing as individual wolf
disp_wolflocs_sf %>%
  mutate(group = as.character(id)) %>%
  group_by(group) %>%
  nest() %>%
  #move walk2 inside mutate as data and group were invisible for walk2
  #also no need for data you can run walk2 directly
  mutate(#data = map(data, ~st_as_sf(.x)),
    txt = walk2(.x = data, .y = group, ~st_write(obj = .x, dsn = paste0(.y, ".shp"))))
# 
# #writing as individual wolf-year
# wolflocs_sf %>%
#   mutate(group = as.character(burst)) %>%
#   group_by(group) %>%
#   nest() %>%
#   #move walk2 inside mutate as data and group were invisible for walk2
#   #also no need for data you can run walk2 directly
#   mutate(#data = map(data, ~st_as_sf(.x)),
#     txt = walk2(.x = data, .y = group, ~st_write(obj = .x, dsn = paste0(.y, ".shp"))))

######------ START OF DATA SEGMENTATION USING LAVIELLE METHOD  ------######

fpt <- variance <- scale <- df <- ltraj <- list()

for (i in 1:42) {
  fpt[[i]] <- fpt(disp_data_regular[i], seq(1000, 100000, by=1000), units="days")
  variance[[i]] <- varlogfpt(fpt[[i]])
  scale[[i]] <- as.numeric(apply(variance[[i]],1,which.max))
}

#getting median scale at which wolves are interacting with their environment
scale <- as.data.frame(Reduce(rbind, scale))
median(scale$V1, na.rm=T)
#ends up being 18.5 km

for (i in 1:42){
  df[[i]] <- ld(disp_data_regular[i])
  df[[i]] <- df[[i]][!is.na(df[[i]]$x),]
  df[[i]] <- df[[i]][!is.na(df[[i]]$y),]
  df[[i]]$fpt_18 <- fpt[[i]][[1]]$r18
}
df <- as.data.frame(Reduce(rbind, df))

ltraj <- dl(df)
id(ltraj) <- burst(ltraj)

lav_wolf <- list()
for(i in 1:42){
lav_wolf[[i]] <- lavielle(ltraj[i], Lmin=5, Kmax=10, type="mean", which="fpt_18")
}

#rule for getting K
#example where my method works better is lav_wolf[[2]]
#Laveille paper even says that 0.75 not great w observations <500
choice_K_Laveille <- function(lav) {
  temp <- chooseseg(lav)
  K <- max(which(temp$D>.75))
  return(K)
}

choice_K_Lisanne <- function(lav) {
  progress <- list()
  temp <- chooseseg(lav)
  range <- temp$Jk[1] - temp$Jk[9]
  for(i in 1:9){
    progress[i] <- (temp$Jk[1]-temp$Jk[i])/range
  }
  vec <- as.vector(Reduce(rbind, progress))
  K <- min(which(vec>0.8))
  return(K)
}


#applying function to list
K1 <- as.vector(Reduce(rbind, lapply(lav_wolf, choice_K_Laveille)))
K2 <- as.vector(Reduce(rbind, lapply(lav_wolf, choice_K_Lisanne)))

#compared to Laveille rule of thumb, my version
#agrees 18/42 times
#overestimates by 1 9 times
#overestimates by 2 8 times
#overestimates by 3 4 times
#overestimates by 5 1 time
#underestimates by 1 2 times

table(K2 - K1)
summary(ltraj)
seg_4370 <- findpath(lav_wolf[[42]], K2[[42]])

test <- ld(seg_4370)
meanfpt(seg_4370, graph = FALSE)
plot(seg_4370)
plot(seg_4370[1])
plot(seg_4370[2])
plot(seg_4370[3])
plot(seg_4370[4])
plot(seg_4370[5])
plot(seg_4370[6])
plot(seg_4370[7])
#identifying segments to keep 
K2

keep <- list(c(1:2), 1, c(1,2), c(1:3), c(1:4), #4 is resident (large HR)
         c(1,2), 1, c(1:3,5), c(1:2), c(1:4), 
         c(1:7), c(1:5), c(1:5), c(1:4), c(2:5), #11/12 is resident (long HR)
         c(1,2,4:6), c(1:3), c(1:3), c(2:3), c(1:3), #17 resident (long HR), #18 burst for that animal resident
         1, c(1:4), c(1:4), 1, 1, #22/23 is resident (large HR)
         1, c(1,3:5), c(1:4), 1, c(2:5), #28 burst for that animal resident
         1, c(1,2,4), c(1:4), c(1:2), c(2:5), 
         c(1:2), c(1:6), c(1:4), c(1:2), c(1,2,4),#37 is resident (large HR)
         c(2,3),c(1:4)) 

#this is creating a single df with fpt_18 as an infoloc across all wolves

summary(ltraj)
path <- df <- list()
for (i in 1:42){
  path[[i]] <- findpath(lav_wolf[[i]], K2[[i]])
  df[[i]] <- ld(path[[i]])
}
df <- as.data.frame(Reduce(rbind, df))
dim(df)
head(df)

df <- df %>% mutate(segment = as.numeric(gsub("Segment.", "", df$burst)))
head(df)

df <- df %>% group_by(id, segment) %>% 
  summarize(median_fpt = median(fpt_18, na.rm = TRUE),
            median_nsd = median(R2n, na.rm = TRUE))

max(df$median_nsd)
median(df$median_nsd)
head(df)
plot_list1 <- plot_list2 <- plot_list3 <- list()
for (i in unique(df$id)){
  plot_fpt<- ggplot(data = df[df$id==i,], aes(y = median_fpt, x=segment)) +
    geom_point(size=3,col="black") + geom_line() +
    facet_wrap(~id)+
    theme_bw() + geom_hline(yintercept=10,linetype="dashed", 
                            color = "red")
  plot_nsd<- ggplot(data = df[df$id==i,], aes(y = median_nsd, x=segment)) +
    geom_point(size=3,col="black") + geom_line() +
    facet_wrap(~id)+
    theme_bw() + geom_hline(yintercept=109451609,linetype="dashed", 
                            color = "red")
  plot_list1[[i]]=plot_fpt
  plot_list2[[i]]=plot_nsd
  plot_list3[[i]]=plot_grid(plot_list1[[i]], plot_list2[[i]], 
            labels = c("FPT", "NSD"),
            ncol = 1, nrow = 2)
}
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
for (i in 1:length(unique(df$id))) {
  file_name = paste(i,"_",unique(df$id),".png")
  png(paste("G:/My Drive/GitHub/Wolves/Data/spatial_model/FPT_NSD_plots","/",file_name[i],sep=""))
  grid.draw(plot_list3[[i]])
  dev.off()
}

######------ PERFORM SEGMENTATION AND CREATE NEW LTRAJ OBJECT  ------######

#create ltraj of edited wolfyears
lav_temp <- ltraj_edit <- df_temp <- list()
for (i in 1:42) {
  lav_temp[[i]] <- findpath(lav_wolf[[i]], K2[i])
  ltraj_edit[[i]] <- lav_temp[[i]][c(keep[[i]])]
  df_temp[[i]] <- ld(ltraj_edit[[i]])
}

#test <- lav_temp[[2]][c(keep[[2]])]

#get disperser data to data frame
disp_data_resident = do.call(rbind, df_temp)
disp_data_resident$date <- as.POSIXct(disp_data_resident$date, tz="America/Los_Angeles")
head(disp_data_resident)
disp_data_resident <- disp_data_resident %>% dplyr::select(-c(fpt_18, burst))

#plot(ltraj_edit[[40]])

#bringing the non-dispersers back in and making daily
#now we can fill in missing days with setNA

refda <- nondisp_data_regular <- list()
for(i in 1:41){
  refda[[i]] <-summary(nondispdata_nogaps)$date.begin[i]
  nondisp_data_regular[i] <- setNA(nondispdata_nogaps[i], date.ref=refda[[i]][1],1,units="day")
  nondisp_data_regular[[i]]$id <- summary(nondispdata_nogaps)$id[i]
  nondisp_data_regular[[i]]$burst <- summary(nondispdata_nogaps)$burst[i]
}

nondisp_data_regular <- as.data.frame(Reduce(rbind, nondisp_data_regular))
nondisp_data_regular$date <- as.POSIXct(round(nondisp_data_regular$date, unit="day"))

nondisp_data_regular_ltraj <- dl(nondisp_data_regular)
summary(nondisp_data_regular_ltraj)

#check that data separated by 1 day
#seems there is a v minor issue with use of daylight savings time but let's ignore
nondisp_data_regular <- ld(nondisp_data_regular_ltraj)
hist(nondisp_data_regular$dt)

resident_data <- nondisp_data_regular %>% filter(!is.na(x)) %>% dplyr::select(-burst)

#creating single data frame
HR_data <- rbind(disp_data_resident, resident_data)
HR_data$year <- year(HR_data$date)
HR_data$wolfyear <- str_c(HR_data$id,"_",HR_data$year)
HR_data <- as.ltraj(xy = HR_data[,c("x","y")], date = HR_data$date, id=HR_data$id, burst=HR_data$wolfyear)
summary <- summary(HR_data)

#need to remove 68m because not enough points in either HR in that year to support it
#should be 139 rows
HR_data <- HR_data[-which(summary$id=="068m.1")]

#get rid of those with <60 locs
summary <- summary(HR_data)
lower60 <- which(summary$nb.reloc<60)
HR_data <- HR_data[-lower60]
HR_data_df <- ld(HR_data)

write.csv(HR_data_df, "Data/spatial_model/GPS_data_forHR.csv")

summary <- summary(HR_data)
write.csv(summary, "Data/spatial_model/GPS_data_forHR_summary.csv")
#should be 120 rows

######------ VERSION WHERE NO SEGMENTS REMOVED ------######

#creating single data frame
resident_data <- nondisp_data_regular %>% filter(!is.na(x)) %>% dplyr::select(-burst)
disp_data <- ld(disp_data_regular)
disp_data <- disp_data %>% filter(!is.na(x)) %>% dplyr::select(-burst)

HR_data <- rbind(disp_data, resident_data)
HR_data$year <- year(HR_data$date)
HR_data$wolfyear <- str_c(HR_data$id,"_",HR_data$year)
HR_data <- as.ltraj(xy = HR_data[,c("x","y")], date = HR_data$date, id=HR_data$id, burst=HR_data$wolfyear)
id(HR_data) <- burst(HR_data)
summary <- summary(HR_data)

#need to remove 68m because not enough points in either HR in that year to support it
#should be 139 rows
HR_data <- HR_data[-which(summary$id=="068m")]

#get rid of those with <60 locs
summary <- summary(HR_data)
lower60 <- which(summary$nb.reloc<60)
HR_data <- HR_data[-lower60]
HR_data_df <- ld(HR_data)
write.csv(HR_data_df, "Data/spatial_model/GPS_data_forHR_noneremoved.csv")

summary <- summary(HR_data)
write.csv(summary, "Data/spatial_model/GPS_data_forHR_noneremoved_summary.csv")
#should be 120 rows

# #writing to shapefile
# wolflocs_postseg <- st_as_sf(big_data, coords =  c("x", "y"), crs = 32610)
# 
# #writing to shapefile by wolf
# wolflocs_postseg %>%
#   mutate(group = as.character(id)) %>%
#   group_by(group) %>%
#   nest() %>%
#   #move walk2 inside mutate as data and group were invisibale for walk2
#   #also no need for data you can ran walk2 directly
#   mutate(#data = map(data, ~st_as_sf(.x)),
#     txt = walk2(.x = data, .y = group, ~st_write(obj = .x, dsn = paste0(.y, ".shp"))))
# 
# #writing to shapefile by wolf-year
# wolflocs_postseg %>%
#   mutate(group = as.character(burst)) %>%
#   group_by(group) %>%
#   nest() %>%
#   #move walk2 inside mutate as data and group were invisibale for walk2
#   #also no need for data you can ran walk2 directly
#   mutate(#data = map(data, ~st_as_sf(.x)),
#     txt = walk2(.x = data, .y = group, ~st_write(obj = .x, dsn = paste0(.y, ".shp"), append=FALSE)))

######------ CHECKING NONDISPERSERS TO SEE IF THERE ARE ANY ABNORMALITIES ------######

#important to have regularized trajectories for FPT, so let's do that with the dispersers
#now we can fill in missing days with setNA
refda <- nondisp_data_regular <- list()
for(i in 1:41){
  refda[[i]] <-summary(nondispdata_nogaps)$date.begin[i]
  nondisp_data_regular[i] <- setNA(nondispdata_nogaps[i], date.ref=refda[[i]][1],1,units="day")
  nondisp_data_regular[[i]]$id <- summary(nondispdata_nogaps)$id[i]
  nondisp_data_regular[[i]]$burst <- summary(nondispdata_nogaps)$burst[i]
}

nondisp_data_regular <- as.data.frame(Reduce(rbind, nondisp_data_regular))
nondisp_data_regular$date <- as.POSIXct(round(nondisp_data_regular$date, unit="day"))

nondisp_data_regular <- dl(nondisp_data_regular)

# ##***NEW LINE TO MAKE SHAPEFILE
# id(disp_data_regular) <- burst(disp_data_regular)
# summary <- summary(disp_data_regular)

#check that data separated by 1 day
#seems there is a v minor issue with use of daylight savings time but let's ignore
df <- ld(nondisp_data_regular)
hist(df$dt)
#hist(disp_data_regular[1:10],"dt",freq=TRUE)

######------ WRITING DATA FROM INDIVIDUAL WOLVES & WOLF-YEARS TO SHAPEFILE  ------######
df <- df[!is.na(df$x),]

disp_wolflocs_sf <- st_as_sf(df, coords =
                               c("x", "y"), crs = 32610)
head(disp_wolflocs_sf)

#writing as individual wolf
disp_wolflocs_sf %>%
  mutate(group = as.character(id)) %>%
  group_by(group) %>%
  nest() %>%
  #move walk2 inside mutate as data and group were invisible for walk2
  #also no need for data you can run walk2 directly
  mutate(#data = map(data, ~st_as_sf(.x)),
    txt = walk2(.x = data, .y = group, ~st_write(obj = .x, dsn = paste0(.y, ".shp"))))
# 
# #writing as individual wolf-year
# wolflocs_sf %>%
#   mutate(group = as.character(burst)) %>%
#   group_by(group) %>%
#   nest() %>%
#   #move walk2 inside mutate as data and group were invisible for walk2
#   #also no need for data you can run walk2 directly
#   mutate(#data = map(data, ~st_as_sf(.x)),
#     txt = walk2(.x = data, .y = group, ~st_write(obj = .x, dsn = paste0(.y, ".shp"))))

######------ START OF DATA SEGMENTATION USING LAVIELLE METHOD  ------######

fpt <- variance <- scale <- df <- ltraj <- list()

for (i in 1:41) {
  fpt[[i]] <- fpt(nondisp_data_regular[i], seq(1000, 100000, by=1000), units="days")
  variance[[i]] <- varlogfpt(fpt[[i]])
  scale[[i]] <- as.numeric(apply(variance[[i]],1,which.max))
}

#getting median scale at which wolves are interacting with their environment
scale <- as.data.frame(Reduce(rbind, scale))
median(scale$V1, na.rm=T)
#ends up being 18.5 km

for (i in 1:41){
  df[[i]] <- ld(nondisp_data_regular[i])
  df[[i]] <- df[[i]][!is.na(df[[i]]$x),]
  df[[i]] <- df[[i]][!is.na(df[[i]]$y),]
  df[[i]]$fpt_14 <- fpt[[i]][[1]]$r14
}
df <- as.data.frame(Reduce(rbind, df))

ltraj <- dl(df)
id(ltraj) <- burst(ltraj)
summary <- summary(ltraj)
ltraj <- ltraj[-which(summary$id=="021m.1"|summary$id=="104f.1"|summary$id=="085m.1"|summary$id=="074m.2"|summary$id=="023f.1")]

plot(ltraj[38])
#6,7,15,16,17,24,31,38 is a problem
lav_wolf <- list()
for(i in 1:41){
  lav_wolf[[i]] <- lavielle(ltraj[i], Lmin=5, Kmax=10, type="mean", which="fpt_14")
}

test1 <- ltraj[which(summary$id=="021m.1")]
test <- ld(test)                  
lav_wolf <- lavielle(test1, Lmin=2, Kmax=10, type="mean", which="fpt_14")
seg_4370 <- findpath(lav_wolf, 2)
plot(seg_4370)


#rule for getting K
#example where my method works better is lav_wolf[[2]]
#Laveille paper even says that 0.75 not great w observations <500
choice_K_Laveille <- function(lav) {
  temp <- chooseseg(lav)
  K <- max(which(temp$D>.75))
  return(K)
}

choice_K_Lisanne <- function(lav) {
  progress <- list()
  temp <- chooseseg(lav)
  range <- temp$Jk[1] - temp$Jk[9]
  for(i in 1:9){
    progress[i] <- (temp$Jk[1]-temp$Jk[i])/range
  }
  vec <- as.vector(Reduce(rbind, progress))
  K <- min(which(vec>0.8))
  return(K)
}


#applying function to list
K1 <- as.vector(Reduce(rbind, lapply(lav_wolf, choice_K_Laveille)))
K2 <- as.vector(Reduce(rbind, lapply(lav_wolf, choice_K_Lisanne)))

#compared to Laveille rule of thumb, my version
#agrees 18/42 times
#overestimates by 1 9 times
#overestimates by 2 8 times
#overestimates by 3 4 times
#overestimates by 5 1 time
#underestimates by 1 2 times

table(K2 - K1)
summary(ltraj)
seg_4370 <- findpath(lav_wolf[[42]], K2[[42]])

test <- ld(seg_4370)
meanfpt(seg_4370, graph = FALSE)
plot(seg_4370)
plot(seg_4370[1])
plot(seg_4370[2])
plot(seg_4370[3])
plot(seg_4370[4])
plot(seg_4370[5])
plot(seg_4370[6])
plot(seg_4370[7])
#identifying segments to keep 
K2

keep <- list(c(1:2), 1, c(1,2), c(1:3), c(1:4), #4 is resident (large HR)
             c(1,2), 1, c(1:3,5), c(1:2), c(1:4), 
             c(1:7), c(1:5), c(1:5), c(1:4), c(2:5), #11/12 is resident (long HR)
             c(1,2,4:6), c(1:3), c(1:3), c(2:3), c(1:3), #17 resident (long HR), #18 burst for that animal resident
             1, c(1:4), c(1:4), 1, 1, #22/23 is resident (large HR)
             1, c(1,3:5), c(1:4), 1, c(2:5), #28 burst for that animal resident
             1, c(1,2,4), c(1:4), c(1:2), c(2:5), 
             c(1:2), c(1:6), c(1:4), c(1:2), c(1,2,4),#37 is resident (large HR)
             c(2,3),c(1:4)) 

######------ IDENTIFYING DISPERSAL AT WOLF LEVEL TO BUILD DISPERSAL DISTRIBUTION  ------######

#create new ltraj that comprising only the segments we wish to keep (generally first and last)
lav_temp <- list()
ltraj_edit <- list()
df_temp <- list()

for (i in 1:35) {
  lav_temp[[i]] <- findpath(lav[[cells[i]]], k[i])
  ltraj_edit[[i]] <- lav_temp[[i]][c(keep[[i]])]
  df_temp[[i]] <- ld(ltraj_edit[[i]])
}

big_data = do.call(rbind, df_temp)
dim(big_data) #38543 rows

edited_bursts <- as.ltraj(xy = big_data[,c("x","y")], date = big_data$date, id=big_data$id, burst=big_data$id)

write.csv(big_data, "Data/spatial_model/dispersal_points_pre_post.csv")

######------ USING DISPERSAL MOVEMENTS TO GENERATE LENGTH AND ANGLE DISTRIBUTIONS  ------######

#reading in points from beginning and end of dispersal movements
dispersal_points <- read.csv("Data/spatial_model/dispersal_points_pre_post.csv")
wolflocs_dispersal <- st_as_sf(dispersal_points, coords = 
                                 c("x", "y"), crs = 32610)

#writing these to shapefile by wolf id
wolflocs_dispersal %>%
  mutate(group = as.character(id)) %>%
  group_by(group) %>%
  nest() %>%
  #move walk2 inside mutate as data and group were invisibale for walk2
  #also no need for data you can ran walk2 directly
  mutate(#data = map(data, ~st_as_sf(.x)),
    txt = walk2(.x = data, .y = group, ~st_write(obj = .x, dsn = paste0(.y, ".shp"))))

#get centroids of each 50% kernel
head(dispersal_points)
dispersal_points$group <- paste(dispersal_points$id, dispersal_points$burst)
xy <- dispersal_points[,c(2,3)]
data <- dispersal_points[,c(2,3,15)]
spdf <- SpatialPointsDataFrame(coords = xy, data = data,
                               proj4string = CRS("+proj=utm +zone=10 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

memory.limit(size=800000)
kernel_UD <- kernelUD(spdf[,3], h="href", grid=1000, same4all=FALSE, extent=1.5)
kernel_50_vertices <- getverticeshr(kernel_UD,50,unin = c("m"),unout = c("km2"))
kernel_sf <- st_as_sf(kernel_50_vertices)
kernel_centroids <- st_centroid(kernel_sf)
plot(kernel_sf[1],reset = FALSE)
plot(kernel_centroids[1],add=T,pch=3,col="black",size=3)
plot(line[1], add=T)

#assigning wolf and segment to kernel centroids
kernel_centroids$wolf <- word(kernel_centroids$id,1,sep = "\\ ")
kernel_centroids$segment <- word(kernel_centroids$id,2,sep = "\\ ")

#writing kernel centroids to shapefile
st_write(kernel_centroids,
         "Outputs/spatial_model/dispersal_centroids.shp", driver = "ESRI Shapefile", append=F) #IGNORE WARNING

#reading it back in
kernel_centroids <- st_read("Outputs/spatial_model/dispersal_centroids.shp") #there should be 75 centroids
st_crs(kernel_centroids) <- 32610

#getting start and end centroids for each wolf
start <- kernel_centroids %>% filter(segment=="Segment.1")
end <- kernel_centroids %>% filter(segment!="Segment.1")

start$id
end$id
#need to add three more rows for the three occurrences (64f,79m,91f) where there are two directions taken
#also need to take those rows from "end" and tack them to the end, rather than where they are now
start[36:38,] <- start[c(18,24,29),]
substitute <- end[c(19,26,32),]
end <- end[-c(19,26,32),]
end[36:38,] <- substitute
tail(end)
tail(start)

#getting start and end coordinates
start_geo <- as.data.frame(st_coordinates(start))
end_geo <- as.data.frame(st_coordinates(end))

#calculating angle and creating rose diagram
dx <- c(end_geo$X - start_geo$X)
dy <- c(end_geo$Y - start_geo$Y)
abs.angle <- atan2(dy, dx)
degrees <- 180/pi * abs.angle
rose.diag(abs.angle, bins=24, prop=2)

#now working on distance

#one way to do distance (went w st_length in end, see later)
# dist <- st_distance(start, end, by_element=T)
# mean(dist)
# plot(dist)

#in order to make plotting the centroids work, need to create new wolf names for 64f,79m,91f
#basically, we are adding rows at the end that will union the start with first attempted destination of 64f, 79m, 91f
kernel_centroids$id
kernel_centroids$wolf
#first we are going to give the first attempted destination of those wolves a new name
kernel_centroids$wolf[c(37,50,61)] <- c("064f.2", "079m.2", "091f.2")
#now we will add rows at the end that duplicate the start of those three wolves
kernel_centroids[c(74:76),] <- kernel_centroids[c(35,48,59),]
#now we name those rows with new wolf names
kernel_centroids$wolf[c(74:76)] <- c("064f.2", "079m.2", "091f.2")
tail(kernel_centroids)

paths <- kernel_centroids %>% group_by(wolf) %>%  
  arrange(segment) %>%
  summarise(wolf = first(wolf), do_union = FALSE) %>%
  st_cast("LINESTRING")
plot(paths[1])

#assign group
#these group assignments are taken from .csv "Data/spatial_model/wolves_that_dispersed_classification.csv"
paths$group <- "assign"
paths$group[paths$wolf=="033m" | paths$wolf=="036f" | paths$wolf=="048f" | paths$wolf=="064f.2" |
              paths$wolf=="079m.2" | paths$wolf=="091f.2"] <- "attempt"
paths$group[paths$wolf=="015f" | paths$wolf=="018m" | paths$wolf=="039m" | paths$wolf=="040f" |
              paths$wolf=="054m" |paths$wolf=="064f" |paths$wolf=="079m"|
              paths$wolf=="085m" |paths$wolf=="091f" |paths$wolf=="93m"] <- "mortality"
paths$group[paths$wolf=="065m" | paths$wolf=="068m" | paths$wolf=="069m" | paths$wolf=="071f" |
              paths$wolf=="087m" | paths$wolf=="OR49f"] <- "unknown"
paths$group[paths$wolf=="017m" | paths$wolf=="041m" | paths$wolf=="043f" | paths$wolf=="044m" | paths$wolf=="047f" |
              paths$wolf=="051f" | paths$wolf=="052f" | paths$wolf=="061m" | paths$wolf=="062m" |
              paths$wolf=="070f" | paths$wolf=="088m" | paths$wolf=="090m" |
              paths$wolf=="102m"| paths$wolf=="107f"| paths$wolf=="OR15m" |
              paths$wolf=="OR35f"] <- "success"
unique(paths$group)

st_write(paths,
         "Outputs/spatial_model/dispersalpaths.shp", driver = "ESRI Shapefile", append=F) #IGNORE WARNING

#read in again as shapefile
library(sf)
library(units)
paths <- st_read("Outputs/spatial_model/dispersalpaths.shp")
st_crs(paths) <- 32610

#calculating length in km of dispersal paths
paths$length <- st_length(paths)
paths$len_km <- as.numeric(set_units(paths$length, km))

hist <- ggplot(paths, aes(x=len_km)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#E69F00")+
  ylab("Density") + xlab("Euclidean dispersal distance (km)")

median(paths$len_km) #77.05
mean(paths$len_km) #163.5925
sd(paths$len_km) #182.7631
min(paths$len_km) #13.89478
max(paths$len_km) #671.9307
