library(wq)
library(plyr)
library(lubridate)
library(stringr)

source('//deqhq1/wqassessment/2012_WQAssessment/DO_Evaluation/R_scripts/00_DO_Evaluation_functions.R')

#### Read in the stations tables to recompile information necessary for segmentation ####
#First the usgs stations
usgs.DO.stations <- stations.w.data
rm(stations.w.data)
usgs.DO.stations$LAKE_LLID <- NA
usgs.DO.stations$LAKE_NAME <- NA
usgs.DO.stations <- usgs.DO.stations[,c("ID","Agency","Site_no","Site_name","Str_LLID","Str_name","Str_RM","LAKE_LLID","LAKE_NAME","DOcriteria")]

#i have to carry an id field for the usgs stations in order to relate them back to the gis layer
#this adds that column for the lasar stations but assigns an NA since it doesn't apply.
#This issue is due to the mismatch of site number to site name in the USGS layer. Since I use the name for USGS stations from
#the raw data table I didn't want to have to reprocess it to match the formatting from the layer table so I added an ID field to the 
#layer that I could use to reference back to.
lasar.DO.stations$ID <- NA
#clean up the lasar stations data frame to get it ready to combine with the usgs stations
lasar.DO.stations <- within(lasar.DO.stations, rm('spwnStart', 'spwnEnd', 'Site_elvft'))

#Combine the usgs and lasar stations information into a single data frame
DO.stations.info <- rbind(lasar.DO.stations, usgs.DO.stations)
#DO.stations.info$LLID_Stream_Lake <- ifelse(is.na(DO.stations.info$LAKE_LLID), DO.stations.info$Str_LLID, 
#                                    paste(DO.stations.info$Str_LLID, DO.stations.info$LAKE_LLID, sep = '/'))

#### Bring in the objects created from the 01b_USGS_DO.R and 01a_LASAR_DO.R files ####
ls.mgl <- rbind(lasar.ls.mgl, usgs.ls.mgl)
ls.temp <- rbind(lasar.ls.temp, usgs.ls.temp)

#### Calculate percent saturation ####
#Whittle down the temperature and concentration data frames to matches only
mgl.temp.match <- merge(ls.mgl, ls.temp[,c('id','Result_clean')], by = 'id')
#Converts site elevation from feet to meters
mgl.temp.match$Site_elvm <- mgl.temp.match$Site_elvft * .3048
#Field renaming for clearer calculation steps
mgl.temp.match <- rename(mgl.temp.match, c('Result_clean.x'='Result_DOmgl','Result_clean.y'='Result_tempC'))
#Field rounding to prevent future error
mgl.temp.match$Result_DOmgl <- round(as.numeric(mgl.temp.match$Result_DOmgl),2)
mgl.temp.match$Result_tempC <- round(as.numeric(mgl.temp.match$Result_tempC),2)
#Calculate theoretical dissolved oxygen concentration so we can determine percent saturation
#A salinity of 0 PPU is assumed since we are dealing with freshwater sites in the Willamette and the Umatilla basins
#The oxySol function is found in the wq package and was designed for sea level applications and therefore does not include
#an elevation variable. The elevation portion of the calculation was pulled from an equation used for TMDL calculations
mgl.temp.match$Result_theo_elvm <- oxySol(mgl.temp.match$Result_tempC, S = 0) * (1 - 0.0001148 * mgl.temp.match$Site_elvm)
#Calculate the percent saturation
mgl.temp.match$ps_calc <- (mgl.temp.match$Result_DOmgl/mgl.temp.match$Result_theo_elvm)*100
#This line applies logic in the standard that only applies to the calculation of daily means, which reads 'For the purpose of calculating the mean, concentrations in excess of 100 percent of saturation are valued at the saturation concentration.'
#We eventually decided not to calculate any daily means so it is not currently being run.
#mgl.temp.match$Result_use <- ifelse(mgl.temp.match$ps_calc > 100, mgl.temp.match$Result_theo_elvm, mgl.temp.match$Result_DOmgl)
mgl.temp.match$Result_use <- mgl.temp.match$Result_DOmgl
mgl.temp.match$Result_use <- round(mgl.temp.match$Result_use, 2)

#Subset mgl.temp.match to select only those columns we need to continue
data <- subset(mgl.temp.match, select = c('id', 'code', 'STATION', 'SAMPLE_DATE_TIME', 'day.POSIX', 'time',
                                          'DOcriteria', 'spwnStart', 'spwnEnd', 'Result_use', 'ps_calc'))
#data <- rbind(data, lasar.sub)

#rm(mgl.temp.match)
#### Split data sets into spawning and non-spawning seasons ####
#First, convert to POSIXlt to be able to access the yday attribute that gives the day of the year
#in order to determine if a date falls between given months
data$day <- as.POSIXlt(data$day)

#converting the month and day to a POSIX data type adds a year but we will only be using the
#yday attribute so it shouldn't matter
data$spwnStart <- strptime(data$spwnStart, '%B %d')
data$spwnEnd <- strptime(data$spwnEnd, '%B %d')

#First take out those sites which never have spawning
data.nospawn <- data[which(is.na(data$spwnStart)),]
data.spawn <- data[which(!is.na(data$spwnStart)),]

#rm(data)
#Then split up the sites with spawning time periods into spawning and non-spawning
#There are two sets of logic to apply because there are time periods that span within and across years
#so to extract the Jan - May/June criteria requires different logic than the rest
spawn.Jan <- data.spawn[yday(data.spawn$day) > yday(data.spawn$spwnStart) & 
                          yday(data.spawn$day) < yday(data.spawn$spwnEnd) & 
                          yday(data.spawn$spwnStart) == 1,]

spawn.Rest <- data.spawn[(yday(data.spawn$day) > yday(data.spawn$spwnStart) | 
                          yday(data.spawn$day) < yday(data.spawn$spwnEnd)) & 
                          yday(data.spawn$spwnStart) != 1,]

#Combine the spawning time periods into a single data frame
spawn <- rbind(spawn.Jan, spawn.Rest)
#Pulls out the non-spawning time periods for those stations that have both spawning and non-spawning time periods
nospawn <- data.spawn[!data.spawn$id %in% spawn$id,]

#Put all the non-spawning time period data together
nospawn <- rbind(nospawn, data.nospawn) 

#### Do the comparison to the standard ####
spawn$DOcriteria <- 11

#Treating every data point as a grab sample
#evluate function takes two fields. the first is the dataframe name and the second field 
#is a character string indicating the column you want to do the comparison on
nospawn$digress <- evaluate(nospawn, 'Result_use')
nospawn$ps_digress <- ifelse(nospawn$DOcriteria == 8, ifelse(nospawn$ps_calc < 90,1,0), NA)
spawn$digress <- evaluate(spawn, 'Result_use')
spawn$ps_digress <- ifelse(spawn$ps_calc < 95, 1, 0)

#final digression determination
nospawn$final_digress <- ifelse(nospawn$DOcriteria == 8, 
                                ifelse(nospawn$digress != 0 & nospawn$ps_digress != 0,1,0),
                                ifelse(nospawn$digress != 0,1,0))
spawn$final_digress <- ifelse(spawn$digress != 0 & spawn$ps_digress != 0,1,0)

#### Rolling up to station and determining 303d category ####
nospawn.totals <- assess.both(nospawn)
spawn.totals <- assess.both(spawn)

#### Bring them together and rename to get ready for segmentation ####
nospawn.crit <- DO.stations.info[,c('Site_no', 'DOcriteria')]
names(nospawn.crit) <- c('STATION', 'crit_val')
nospawn.totals <- merge(nospawn.totals, nospawn.crit, by = 'STATION', all.x = TRUE)
nospawn.totals$Season <- 'Year Around (Non-spawning)'

spawn.totals$crit_val <- 11
DO.Season <- ls.mgl[,c('STATION', 'spwnStart', 'spwnEnd')]
DO.Season <- DO.Season[!duplicated(DO.Season),]
DO.Season$Season <- ifelse(is.na(DO.Season$spwnEnd),
                           'Year Around', 
                           paste(as.character(DO.Season$spwnStart), as.character(DO.Season$spwnEnd), sep = ' - '))
DO.Season <- DO.Season[,c('STATION', 'Season')]
spawn.totals <- merge(spawn.totals, DO.Season, by = 'STATION', all.x = TRUE)

DO.final <- rbind(nospawn.totals, spawn.totals)

DO.final$Pollutant <- 'Dissolved Oxygen'

DO.final[DO.final$Season == 'Year Around (Non-spawning)','Season'] <- 'Year Round (Non-spawning)'

#merge(DO.final, DO.stations.info, by.x = 'STATION', by.y = 'Site_no')


rm(list=setdiff(ls(), c('DO.stations.info',"DO.final")))