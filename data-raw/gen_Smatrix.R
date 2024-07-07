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
load(file=paste0(local_data_dir,"/","contacts_fluscape_V1.rda"))

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

#' ## Prep for making the S matrix
lat.lim2 <- as.vector(c(22.11667,24.50833))
long.lim2 <- as.vector(c(112.2667,114.8000))
margin <- 0
ext <- extent(
  long.lim2[1]-margin,long.lim2[2]+margin,lat.lim2[1]-margin,lat.lim2[2]+margin
)
gz_pop_raster <- crop(x1,ext)

#' Brings in the previously generated (in 2017) S matrix for called pop_S_mat_fluscape
load("C:/Users/haowe/Desktop/fluscapeR/data/pop_S_mat_fluscape.rda") # vivi


