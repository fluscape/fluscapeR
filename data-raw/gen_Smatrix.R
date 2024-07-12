#' ---
#' title: "File to prep S matrix for the gravity analysis of fluscape data"
#' author: Vivi Wang
#' date: 2024-07-04
#' output: html_document
#' ---

#' Clear memory, setup required libraries and set local data dir
rm(list=ls(all=TRUE))

library(devtools)
library(raster)
library(tidyverse)
library(raster)

load_all()

#' Assumes that the soure package fluscapeR is next to main private fluscape repo
local_data_dir <- "C:/Users/haowe/Desktop/fluscapeR" # Vivi
fluscape_top_dir <- "C:/Users/haowe/Desktop/fluscape/" # Vivi

source(paste0(fluscape_top_dir,"source/R/mob_utility_private.r"))
source(paste0(fluscape_top_dir,"source/R/fluscape_copy_stevensRfunctions.R"))
source(paste0(fluscape_top_dir,"source/R/GeneralUtility.r"))

#' Load the snapshot of landscan saved into the fluscape directory
x1 <- fsc.load.wide.raster(fluscapetopdir=fluscape_top_dir)

#' ## Make the contacts data base for visit 1 and a jittered version
#'
#' The non-jittered version is NOT used in the package, but is generated here so that the
#' analysis can be run by the fluscape team from exactly the same code that others would
#' use to test the methods.

#' Load the fluscape data sets, and trim to visit 1
participants <- load_particpant_data_long(topdir = fluscape_top_dir)
participants <- participants[participants$VISIT==1,]
participants$pid = paste(participants$LOC_ID, participants$HH_ID, participants$PARTICIPANT_ID, sep="_")
households <- load_hh_data_long(topdir = fluscape_top_dir)
# recreate locations information from households
locations <- unique( households[,c("LOC_ID","LOC_Lat","LOC_Long","URBAN","DIST_FRM_GZ")] )

#' Set area boundary
long.lim = c( 112.8, 114.2 )
lat.lim = c( 22.6, 24.0 )

#' Load the actual contact data
load(file=paste0(local_data_dir,"/data/contacts_fluscape_V1.rda"))

#' Generate actual potentially identifiable data and main jittered data
#' The approximate conversions are: Latitude: 1 deg = 110.574 km. Longitude: 1 deg = 111.320*cos(latitude) km. So at 23 degrees, 1 degree of logitude is 110*cos(23/180*pi) = 101 ~= 100
#' So a jitter error of 100 metres is a jitter of 0.001
contacts_jit <- mob_generate_contacts_data_from_scratch(
  locations,
  households,
  participants,
  jitterxy=TRUE,
  jitter_error=0.001,
  x1=x1,
  topdir=fluscape_top_dir
)

#' assess the jittering
hist(contacts_fluscape_V1$lat,breaks=c(0,seq(23,24,0.1),40))
hist(contacts_jit$lat-contacts_fluscape_V1$lat)
plot( contacts_fluscape_V1$HH_Long, contacts_fluscape_V1$HH_Lat, pch=16, cex=.3, col="blue", xlim=c(113.5,113.6), ylim=c(23.2,23.3))
points( contacts_jit$HH_Long, contacts_jit$HH_Lat, pch=16, cex=.3, col="red")

#' Rename and use 100m jittered data
contacts_V1_jittered_100m <- contacts_jit
usethis::use_data(contacts_V1_jittered_100m, overwrite = TRUE)

#' Prepare for making S matrix
contacts <- contacts_fluscape_V1
# contacts_tmp <- contacts[contacts$LOC_ID<=5 & contacts$LOC_ID!=2,]
contacts_tmp <- contacts

contacts_small <- contacts_tmp[!is.na(contacts_tmp$long) | !is.na(contacts_tmp$lat) | !is.na(contacts_tmp$HH_Long) | !is.na(contacts_tmp$HH_Lat),]
range(contacts_small$long)
range(contacts_small$lat)
originx <- contacts_small$HH_Long
originy <- contacts_small$HH_Lat
dim(contacts_small)

#' Find the spatial extent of destination long-lats, and crop landscan density raster.
ext <- extent(
  min(contacts_small$long),
  max(contacts_small$long),
  min(contacts_small$lat),
  max(contacts_small$lat)) # suitable for all locations
x2 <- crop(x1,ext)

#' Generating S matrix for radiation model will take very long time (~ a week).
#' For people who want to run the models by themselves, instead of generating full size
#' S matrix, we write a function to downsize the raw data and generate a sample S matrix.
#' Use agg_num to adjust the size - larger agg_num, smaller size
#' If agg_num = 0, it will calculate based on the full size of the data.

S.raiation.agg10 <- generate.agg.Smatrix(contacts_fluscape_V1, x2, 10)
S.raiation.agg5 <- generate.agg.Smatrix(contacts_fluscape_V1, x2, 5)

S.raiation.full <- generate.agg.Smatrix(contacts_fluscape_V1, x2, 0)

#' Check that the gravity model produces sensible results and try to calculate the likelihood of the
#' data

fnLog <- paste0(local_data_dir,"/gravity_log_debug.csv")
if (!file.exists(fnLog)) {
  dftmp <- make.data.df(nrow=0)
  write.csv(dftmp,fnLog,row.names=FALSE)
}

#' Fit radiation model to sample S matrix and actual S matrix 
fit.sampleS10 <- fit.mobility.model(
  contacts = contacts_fluscape_V1,
  popgrid = x2,
  Smat = S.raiation.agg,
  logfile = fnLog,
  optfun = fit.offset.radiation.optim,
  psToFit = c("offset"),
  psLB = c(1),
  psUB = c(20*1000),
  datasubset = "ALL",
  fdebug=TRUE,
  lognote = ""
)

fit.sampleS5 <- fit.mobility.model(
  contacts = contacts_fluscape_V1,
  popgrid = x2,
  Smat = S.raiation.agg,
  logfile = fnLog,
  optfun = fit.offset.radiation.optim,
  psToFit = c("offset"),
  psLB = c(1),
  psUB = c(20*1000),
  datasubset = "ALL",
  fdebug=TRUE,
  lognote = ""
)

fit.fullS <- fit.mobility.model(
  contacts = contacts_fluscape_V1,
  popgrid = x2,
  Smat = S.raiation.full,
  logfile = fnLog,
  optfun = fit.offset.radiation.optim,
  psToFit = c("offset"),
  psLB = c(1),
  psUB = c(20*1000),
  datasubset = "ALL",
  fdebug=TRUE,
  lognote = ""
)

###################################################
#' Generate a gravity model
# faster non-zero destination indices calculation:
tmp = getValues(x2)
cells = which( !is.na(tmp) & tmp>0 )
destindices <- data.frame(index=1:length(cells), cells=cells)

# faster origin indices calculation:
cells = unique( cellFromXY(x2, contacts_small[,c("HH_Long","HH_Lat")] ) )
originindices <- data.frame(index=1:length(cells), cells=cells)

# faster distances calculation:
a = xyFromCell( x2, destindices$cell )
b = xyFromCell( x2, originindices$cell )
distances = pointDistance( b, a, longlat=T )

# Next steps
# - generate a matrix of observations
# - generate a gravity model
# - compare the gravity models with the data

# number of origin and possible destination cells
noorig <- dim(originindices)[1] # number of origin cells
nodest <- dim(destindices)[1] # number of possible destination cells, large!

# make a matrix of observed contact origin-destinations
contactindices = cellFromXY( x2, contacts_small[,c("long","lat")] ) # in which cells did contact occur?
hhindices = cellFromXY( x2, contacts_small[,c("HH_Long","HH_Lat")] ) # which cells have participant households?
obs.tab = matrix( 0, nrow=noorig, ncol=nodest  ) # matrix of observations
pb = txtProgressBar(min = 0, max = noorig, initial = 0)
for (k in 1:noorig) {
  j = originindices$cell[k]
  i = contactindices[hhindices==j]
  i = i[!is.na(i)]
  z = table(i)
  i.indices = match(as.numeric(names(z)),destindices$cells)
  obs.tab[k,i.indices] = as.numeric(z)
  setTxtProgressBar(pb,k)
}
close(pb)
#' Faster gravity model generation
M <- 1
#' Make a matrix of predicted origin-destination weights for gravity model
#' with an offset of 500 and a power term M.
gravmodel <- matrix( nrow=noorig, ncol=nodest )
pb <- txtProgressBar(min = 0, max = noorig, initial = 0)
for (i in 1:noorig) {
  orcell 	<- originindices$cell[i]
  dscell 	<- destindices$cell
  n_o 	<- extract(x2,orcell)
  n_d 	<- extract(x2,dscell)
  dist 	<- distances[i, ]
  gravmodel[i, ] <- (n_o * n_d) / (500 + dist^M)
  setTxtProgressBar(pb,i)
}
close(pb)

pb <- txtProgressBar(min = 0, max = noorig, initial = 0)
for (i in 1:noorig) {
  gravmodel[i,] <- gravmodel[i,] / sum(gravmodel[i,])
  setTxtProgressBar(pb,i)
}
close(pb)



