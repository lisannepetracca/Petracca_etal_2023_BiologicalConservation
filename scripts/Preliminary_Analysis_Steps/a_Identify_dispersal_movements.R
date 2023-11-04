#---
#"Identify dispersal movements from wolf GPS collar data"
#"Lisanne Petracca"
#"November 2023"
#---

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
data <- read.csv("C:/Your_Directory/wolf_gpsdata_all_processed.csv", header=T)
data$date <- as.POSIXct(data$date, tz="America/Los_Angeles")
data$year <- year(data$date)
data$hour <- hour(data$date)
data$month <- month(data$date)
data$bioyear <- NA

#coding bioyear (2009 is May 2009-April 2010)
data <- data %>%  mutate(bioyear = ifelse(month<=4, year-1, year))

#create column for wolfyear (concatenating id and bioyear)
data$wolfyear <- str_c(data$id,"_",data$bioyear)

#format date
data$mdy <- format(data$date, "%Y-%m-%d")

#let's subset to daily locs & then get locations closest to 7 am every day
data_daily <- data %>% group_by(id, mdy) %>% filter(abs(hour-7)==min(abs(hour-7))) %>% ungroup()

#format date to just have m/d/y
data_daily$date <- as.POSIXct(format(as.POSIXct(data_daily$date,format='%Y-%m-%d" %H:%M:%S'),format='%Y-%m-%d'))

#remove duplicates for multiple fixes in same hour
#will keep the first entry for each day
data_daily <- data_daily %>% group_by(id,date) %>% slice(1) %>% ungroup()

#remove bursts with <60 days of data, as two months of movement isn't necessarily representative of territory size
data_daily <- data_daily %>% group_by(id) %>% filter(n()>=60)

#conversion to class ltraj
finaldata <- as.ltraj(xy = data_daily[,c("x","y")], date = data_daily$date, id=data_daily$id, burst=data_daily$id)
sum(summary(finaldata)$nb.reloc)

#convert to dataframe
dataframe <- ld(finaldata)

#let's write out this .csv for use in the survival model detection histories
write.csv(dataframe, "C:/My_Directory/GPS_daily_data_for_survival.csv")

#let's get max NSD by wolf here
#we end up using 50th percentile for max NSD to subset to potential dispersers

nsd <- list()
for(i in 1:74){
  nsd[i] <- max(finaldata[[i]]$R2n)
}
nsd <- as.vector(Reduce(rbind, nsd))
plot(nsd)
quantile(nsd, probs = seq(0,1,0.05))

#seeing number of wolves falling into different bins of max NSD
length(which(nsd>11974379217, arr.ind=T))
length(which(nsd>5155952369, arr.ind=T)) #65th percentile #adds 21,36,37,38,43,44,69
length(which(nsd>3909502334, arr.ind=T)) #60th percentile #adds 18,30,31,40
length(which(nsd>3534604396, arr.ind=T)) #55th percentile #adds 8,17,27
length(which(nsd>2616581604, arr.ind=T)) #50th percentile #adds 16,22,42,64
length(which(nsd>1966890528, arr.ind=T)) #45th percentile #adds 10,13,34,55 

sqrt(2616581604) #equates to ~50 km

#dividing data into those with dispersal and those without
dispersers <- which(nsd>2616581604, arr.ind=T)
non_dispersers <- which(nsd<2616581604, arr.ind=T)
disp_data <- finaldata[dispersers]
nondisp_data <- finaldata[non_dispersers]

sum(summary(disp_data)$nb.reloc) #19537 points from dispersers

#break down trajectories when there is gap of 14+ days
#60 * 60 * 24 * 14
foo2 <- function(dt) 
{
  return(dt> (2*1209600))
}

#do this separately for those with and without dispersal movements
dispdata_nogaps <- cutltraj(disp_data, "foo2(dt)", nextr = TRUE)
summary(dispdata_nogaps) 

nondispdata_nogaps <- cutltraj(nondisp_data, "foo2(dt)", nextr = TRUE)
summary(nondispdata_nogaps) 

alldata_nogaps <- rbind(ld(dispdata_nogaps), ld(nondispdata_nogaps))
alldata_nogaps$id <- alldata_nogaps$burst

alldata_nogaps <- dl(alldata_nogaps)

######------ DOING ALL DATA TOGETHER ------######

#important to have regularized trajectories for FPT, so let's do that with all data
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

#line to make shapefile
id(all_data_regular) <- burst(all_data_regular)
summary <- summary(all_data_regular)

#remove burst with 6 relocs (87m)
all_data_regular <- all_data_regular[-which(summary$burst=="087m.1")]

#check that data separated by 1 day
#seems there is a v minor issue with use of daylight savings time but let's ignore
df <- ld(all_data_regular)
hist(df$dt)

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
#ends up being 18 km

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

#making ltraj of those trajectories that can and cannot be modeled at that scale
ltraj_qualify <- ltraj[-c(10,29,31,64,75)]
ltraj_notqualify <- ltraj[c(10,29,31,64,75)]

#perform segmentation on those trajectories
lav_wolf <- list()
for(i in 1:78){
  lav_wolf[[i]] <- lavielle(ltraj_qualify[i], Lmin=2, Kmax=8, type="mean", which="fpt_18")
}

#two different approaches for getting at number of segments (I went with choice_K_Lisanne)

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

#applying function to list, get at differences between two methods
K1 <- as.vector(Reduce(rbind, lapply(lav_wolf, choice_K_Laveille)))
K2 <- as.vector(Reduce(rbind, lapply(lav_wolf, choice_K_Lisanne)))
table(K2 - K1)

#####----- SEGMENTATION OF PATHS AND CALCULATING NSD FROM FIRST SEGMENT ----####

#performing the segmentation given specified number of segments
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
  png(paste("C:Your_Directory/FPT_NSD_plots_all_5075","/",file_name[i],sep=""))
  grid.draw(plot_list3[[i]])
  dev.off()
}


#####----- KEEPING TRAJECTORIES FOR HR ANALYSIS ------#####

summary(path)
paths_HRonly <- path[-9] #this entire wolf removed bc all dispersal
paths_occupied <- path

#segments to keep as nondispersal segments
keep <- list(c(1:5), c(1:3), c(1:4), c(1:5), c(1:6),
             c(1:2), c(1,3), c(1:2), c(1:3), c(1:3), 
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

#####----- FINAL DATASET FOR HR ANALYSIS -----#####

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
write.csv(df_occupied_all, "C:/Your_Directory/GPS_data_forHR.csv")
length(unique(df_occupied_all$id))

#####----- KEEPING TRAJECTORIES FOR MOVEMENT PERIODS -----#####

#this will help us estimate movement distances
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
write.csv(summary_movers, "C:/Your_Directory/movement_data_for_survival.csv")
write.csv(movers_df, "C:/Your_Directory/movement_data_for_distancecalc.csv")

######------ CALCULATING MOVEMENT DISTANCE ------######

library(tidyverse)
movements <- read.csv("C:/Your_Directory/movement_data_for_distancecalc.csv", header=T)
summary <- read.csv("C:/Your_Directory/movement_data_for_survival.csv", header=T)

#first we need to remove the three OR wolves bc we aren't interested in those distances
summary <- summary %>% filter(!burst %in% c("OR15m.1_3", "OR35f.1_1",
                                            "OR49f.1_1", "OR49f.1_5"))
movements <- movements %>% filter(!burst %in% c("OR15m.1_3", "OR35f.1_1",
                                                "OR49f.1_1", "OR49f.1_5"))

#convert to sf and take first/last coordinate by wolf
library(sf)
move_sf <- st_as_sf(movements, coords =
                      c("x", "y"), crs = 32610)
first_pts <- move_sf %>% group_by(id) %>% slice_head() %>% st_coordinates()
last_pts <- move_sf %>% group_by(id) %>% slice_tail() %>% st_coordinates()#whether individual wolves moved in or out

#getting distance btw first and last points as the crow flies
library(raster)
distances <- pointDistance(first_pts, last_pts, lonlat=F, allpairs=F)
distances_km = as.data.frame(distances/1000)

#say whether this distance is within or out of state
distances_km$location <- c("out", "out", "out", "in", "in",
                      "in", "in", "out","out","in", 
                      "in", "out", "in", "out", "in",
                      "out", "in", "in", "in", "out",
                      "in", "out","in", "out", "in",
                      "out","out", "out","out","out",
                      "out","in","in","out")

distances_km_instate <- distances_km %>% filter(location=="in")
distances_km_all <- distances_km
distances <- distances_km_instate$`distances/1000`
write.csv(distances, "C:/Your_Directory/instate_distances.csv")

distances <- distances_km_all$`distances/1000`
write.csv(distances, "C:/Your_Directory/spatial_model/all_distances.csv")

#random thing getting ids of dispersers
disperser_id <- move_sf %>% group_by(id) %>% slice(1L) %>% dplyr::select(date, id)
write.csv(disperser_id, "C:/Your_Directory/disperser_id.csv")
