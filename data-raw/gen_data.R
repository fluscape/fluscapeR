#' ---
#' title: "File to prep data for the gravity analysis of fluscape data"
#' author: Steven Riley
#' date: 2024-01-04
#' output: html_document
#' ---

#' This file is designed to be run with the private fluscape repo sat next to this development 
#' version of the fluscapeR package. It will source some functions not available to fluscape 
#' and will make versions of th key data that are not available to the public. It will also create 
#' the jittered public versions of the data that are included in the package.

#' Clear memory, setup required libraries and set local data dir
rm(list=ls(all=TRUE))
library(devtools)
library(raster)
load_all()

#' Assumes that the soure package fluscapeR is next to main private fluscape repo
local_data_dir <- "~/tmp"
fluscape_top_dir <- "../fluscape/"

source("../fluscape/source/R/mob_utility_private.r")
source("../fluscape/source/R/fluscape_copy_stevensRfunctions.R")
source("../fluscape/source/R/GeneralUtility.r")

#' ## Load the population density matrix

#' Load the snapshot of landscan saved into the fluscape directory
x1 <- fsc.load.wide.raster(fluscapetopdir=fluscape_top_dir)

#' This is not yet saved as a package data object, until I see how it get used later

#' ## Make the contacts data base for visit 1 and a jittered version
#' 
#' The non-jittered version is NOT used in the package, but is generated here so that the 
#' analysis can be run by the fluscape team from exactly the same code that others would 
#' use to test the methods.

#' Load the standard visit 1 data sets
locations <- load.and.merge.locs.V1( topdir = fluscape_top_dir , make.corrections=TRUE)
households <- load.household.data.V1(topdir = fluscape_top_dir)
participants <- load.and.merge.part.V1(topdir = fluscape_top_dir)
participants$pid = paste(participants$LOC_ID, participants$HH_ID, participants$PARTICIPANT_ID, sep="_")

#' Set area boundary
long.lim = c( 112.8, 114.2 )
lat.lim = c( 22.6, 24.0 )

#' Keep old data line in just in case. Not run.
## contacts_0 <- mob_load_old_contact_data( locations, households, participants )

#' Generate actual potentially identifiable data and main jittered data
#' The approximate conversions are: Latitude: 1 deg = 110.574 km. Longitude: 1 deg = 111.320*cos(latitude) km. So at 23 degrees, 1 degree of logitude is 110*cos(23/180*pi) = 101 ~= 100
#' SO a jitter error of 100 metres is a jitter of 0.001
contacts_fluscape_V1 <- mob_generate_contacts_data_from_scratch(locations, households, participants,
                                                      jitterxy=FALSE, topdir=fluscape_top_dir)
contacts_jit <- mob_generate_contacts_data_from_scratch(locations, households, participants,
                                                      jitterxy=TRUE, jitter_error=0.001,x1=x1,
                                                      topdir=fluscape_top_dir )
                                                      
hist(contacts_fluscape_V1$lat,breaks=c(0,seq(23,24,0.1),40))
hist(contacts_jit$lat-contacts_fluscape_V1$lat)

#' rename and write the potentially personally identifiable version to the local temp directory
save(contacts_fluscape_V1,file=paste0(local_data_dir,"/","contacts_fluscape_V1.rda"))

#' Rename and use 100m jittered data
contacts_V1_jittered_100m <- contacts_jit
usethis::use_data(contacts_V1_jittered_100m, overwrite = TRUE)

#' ## Make the S matrix

#' When I get the point that I need it, I'll make a small version of the S matrix here and then rerun 
#' the code to make the larger version that will be stored on zenodo. 