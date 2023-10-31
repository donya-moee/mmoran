library(tidyverse)
library(dplyr)
library(lubridate)
library(wiqid)
library(unmarked)

# Prepare data
# read in data
COdata <- read.csv("CO9_07_23.csv")
# use MLE to make detection records
COdata$detect <- ifelse(COdata$MYOLUC.P>0.05,0,1)
# get julian day column
COdata$night_date <- ymd(COdata$night_date.x)
COdata$julian <- yday(COdata$night_date)
head(COdata)

# aggID is a character
COdata$aggID.x <- as.character(COdata$aggID.x)
class(COdata$aggID.x)

# select columns I'll actually be using
dh = COdata %>%
  select(aggID.x, julian, Year, open_water.x, developed_open.x, developed_low.x, developed_medium.x, developed_high.y, barren.x, deciduous_forest.x, evergreen_forest.x, mixed_forest.x, shrub.x, grassland.x, pasture.x, crops.x, woody_wetlands.x, herbaceous_wetlands.x, MeanElevation, MeanSlope, canopy.cover, detect) %>%
  unique()

# group certain variables together
dh$developed <- (dh$developed_open.x + dh$developed_low.x + dh$developed_medium.x + dh$developed_high.y)
dh$low_veg <- (dh$shrub.x + dh$grassland.x + dh$pasture.x + dh$crops.x)
dh$wetland <- (dh$woody_wetlands.x + dh$herbaceous_wetlands.x)
head(dh)

dh2 = dh %>%
  select(aggID.x, julian, Year, open_water.x, developed, low_veg, wetland, MeanElevation, MeanSlope, canopy.cover, detect)

# make pivot table
dhpivot <- pivot_wider(dh2, names_from = "julian", values_from = "detect")
View(dhpivot)

# write function to fix list in object cols
unlist_data <- function(x){
  output <- numeric(length(x))
  for(i in seq_along(x)){
    unlisted <- unlist(x[i])
    output[i] <- ifelse(1 %in% unlisted,1,0)
    if(is.null(unlisted)){output[i] <- NA}
  }
  return(output)
}

# fix detection history 
dhpivot_fixed1 <- dhpivot %>%
  mutate(across(10:ncol(dhpivot), unlist_data))

View(dhpivot_fixed1)

# randomly select one if site is used for multiple years
dhpivot_fixed <- dhpivot_fixed1 %>% 
  group_by(aggID.x) %>% 
  sample_n(1, replace = F)

# check to make sure the number of rows is correct
nrow(dhpivot_fixed)
View(dhpivot_fixed)

# use wiqid package to look at % of sites occupied 
wiqid::occSS(dhpivot_fixed[,10:ncol(dhpivot_fixed)])

# save an extra version of dhpivot_fixed for later
dhpivot_fixed2 <- dhpivot_fixed

# create detection matrix 
detection = dhpivot_fixed[,c(10:ncol(dhpivot_fixed))]
detectmatrix <- as.matrix(detection)
View(detectmatrix)
str(detectmatrix)
y <- detectmatrix

# quick and dirty bootstrap to see if anything comes up as significant
View(dhpivot_fixed2)

vip.data<- dhpivot_fixed2
vip.data$Presence<- ifelse(rowSums(vip.data[,c(10:147)], na.rm = T)>=1,1,0)
vip.data<- vip.data[,c(148,2:9)]
View(vip.data)

varimp.glm<- function(data){
  COVS<- data[,c(2:ncol(data))]
  ADJ_DEV<- data.frame(COVS = NA, AIC = NA, ADJ = NA)
  for (i in 1:ncol(COVS)){
    model <- glm(data[,1] ~ COVS[,i], na.action = na.omit, family = binomial)
    ADJ_DEV[i, 1]<- names(COVS)[i]
    ADJ_DEV[i, 2]<- round(wiqid::AICc(model), 3)
    ADJ_DEV[i, 3]<- round((1- (model$deviance/model$null.deviance)),6)
  }
  return(ADJ_DEV)
}

glm.covs<- data.frame(scale(vip.data[,2:ncol(vip.data)]))
VIP.DATA<- cbind(vip.data$Presence,glm.covs)
names(VIP.DATA)[1]<- 'Presence'
varimp.glm(data = VIP.DATA)

VIP.LIST<- list()
b<-100
for(i in 1:b){cat(round((i/b)*100,1),'% Complete','\n')
  boot<- dplyr::sample_n(VIP.DATA, nrow(VIP.DATA), replace = T)
  VIP.LIST[[i]]<- varimp.glm(data = boot)
}

VIP.LISTS<- bind_rows(VIP.LIST)
DEV.FIT<- VIP.LISTS %>% group_by(COVS) %>% summarise_all(mean)
DEV.SD<- VIP.LISTS %>% group_by(COVS) %>% summarise_all(sd)
DEV.FIT$SD<- DEV.SD$ADJ

library(ggplot2)

ggplot(DEV.FIT, aes(x=COVS, y = ADJ))+
  geom_point()+
  geom_errorbar(aes(ymin = ADJ + SD,
                    ymax = ADJ - SD), width = 0.5)

# make dataframe for occupancy covatiates (site covariates) 
site.cov <- vip.data[,3:ncol(vip.data)]
occ.covs <- site.cov
View(occ.covs)
# make matrix for detection covariates (covariates that depend on site and date)
# pivot longer back from dhpivot_fixed2
View(dhpivot_fixed2)
datepiv <- dhpivot_fixed2 %>%
  pivot_longer(
    cols = c(10:ncol(dhpivot_fixed2)),
    names_to = "julian",
    values_to = "detect"
  )
View(datepiv)

# make new year and julian columns so that NAs match up in matrix
datepiv$year_na <- ifelse(is.na(datepiv$detect), NA, datepiv$Year) 
datepiv$julian_na <- ifelse(is.na(datepiv$detect), NA, datepiv$julian)

juliandetect = datepiv %>%
  select(aggID.x, julian_na, julian) %>%
  unique()
juliandetect$julian_na <- as.numeric(juliandetect$julian_na)

julianpiv <-  pivot_wider(juliandetect, names_from = "julian", values_from = "julian_na")
View(julianpiv)

yeardetect = datepiv %>%
  select(aggID.x, year_na, julian) %>%
  unique()
View(yeardetect)
yeardetect$year_na <- as.numeric(yeardetect$year_na)

yearpiv <- pivot_wider(yeardetect, names_from = "julian", values_from = "year_na")
View(yearpiv)

julianmatrix <- as.matrix(julianpiv[c(2:ncol(julianpiv))])
View(julianmatrix)

yearmatrix <- as.matrix(yearpiv[c(2:ncol(julianpiv))])                                  
View(yearmatrix)

det.covs <- list(jday = julianmatrix, Year = yearmatrix)

# make unmarked object
umf <- unmarkedFrameOccu(y = y, siteCovs = occ.covs, obsCovs = det.covs)
summary(umf)

# scale data
umf@siteCovs$open_water.x <- scale(umf@siteCovs$open_water.x)
umf@siteCovs$developed <- scale(umf@siteCovs$developed)
umf@siteCovs$low_veg <- scale(umf@siteCovs$low_veg)
umf@siteCovs$wetland <- scale(umf@siteCovs$wetland)
umf@siteCovs$MeanElevation <- scale(umf@siteCovs$MeanElevation)
umf@siteCovs$MeanSlope <- scale(umf@siteCovs$MeanSlope)
umf@siteCovs$canopy.cover <- scale(umf@siteCovs$canopy.cover)

umf@obsCovs$jday <- scale(umf@obsCovs$jday)
umf@obsCovs$Year <- scale(umf@obsCovs$Year)

# fitting models
# most basic model
fm <- occu(formula = ~ 1 
           ~ 1,
           data = umf)

fm

# back transform to original scale
backTransform(fm, type = "state")
backTransform(fm, type = "det")

# add in site covariates 
fm1 <- occu(formula = ~ 1 
            ~ open_water.x + developed + low_veg + wetland + MeanElevation + MeanSlope + canopy.cover,
            data = umf)
fm1

# add in detection covariates 
fm2 <- occu(formula = ~ jday + Year 
            ~ open_water.x + developed + low_veg + wetland + MeanElevation + MeanSlope + canopy.cover,
            data = umf)

fm2

confint(fm2, type = "state")
confint(fm2, type = "det")

