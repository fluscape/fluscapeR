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
library(tidyverse)
library(raster)

load_all()

#' Assumes that the soure package fluscapeR is next to main private fluscape repo
#local_data_dir <- "~/tmp" # Steven
local_data_dir <- "D:/tmp" # jon
fluscape_top_dir <- "../fluscape/"

source(paste0(fluscape_top_dir,"source/R/mob_utility_private.r"))
source(paste0(fluscape_top_dir,"source/R/fluscape_copy_stevensRfunctions.R"))
source(paste0(fluscape_top_dir,"source/R/GeneralUtility.r"))

#' Load the snapshot of landscan saved into the fluscape directory
  x1 <- fsc.load.wide.raster(fluscapetopdir=fluscape_top_dir)

#' This is not yet saved as a package data object, until I see how it get used later

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

#' These are old versions of the data
part_all <- load_particpant_data_long()
cont_all <- read.csv("../fluscape/data/clean_datasets/ContactsAll.csv")
dim(part_all)
names(part_all)


#' Keep old data line in just in case. Not run.
## contacts_0 <- mob_load_old_contact_data( locations, households, participants )

#' Generate actual potentially identifiable data and main jittered data
#' The approximate conversions are: Latitude: 1 deg = 110.574 km. Longitude: 1 deg = 111.320*cos(latitude) km. So at 23 degrees, 1 degree of logitude is 110*cos(23/180*pi) = 101 ~= 100
#' SO a jitter error of 100 metres is a jitter of 0.001
  contacts_fluscape_V1 <- mob_generate_contacts_data_from_scratch(
    locations,
    households,
    participants,
    jitterxy=FALSE,
    topdir=fluscape_top_dir
  )
  contacts_jit <- mob_generate_contacts_data_from_scratch(
    locations,
    households,
    participants,
    jitterxy=TRUE,
    jitter_error=0.001,
    x1=x1,
    topdir=fluscape_top_dir
  )
# assess the jittering
  hist(contacts_fluscape_V1$lat,breaks=c(0,seq(23,24,0.1),40))
  hist(contacts_jit$lat-contacts_fluscape_V1$lat)

#' rename and write the potentially personally identifiable version to the local temp directory
  save(contacts_fluscape_V1,file=paste0(local_data_dir,"/","contacts_fluscape_V1.rda"))

#' Rename and use 100m jittered data
  contacts_V1_jittered_100m <- contacts_jit
  usethis::use_data(contacts_V1_jittered_100m, overwrite = TRUE)

#' ## Make the S matrix
  lat.lim2 <- as.vector(c(22.11667,24.50833))
  long.lim2 <- as.vector(c(112.2667,114.8000))
  margin <- 0
  ext <- extent(
    long.lim2[1]-margin,long.lim2[2]+margin,lat.lim2[1]-margin,lat.lim2[2]+margin
  )
  gz_pop_raster <- crop(x1,ext)

#' Brings in the S matrix for called pop_S_mat_fluscape
load("~/dbox/shares/me_jr_hm_dc_gravity/current/to_upload/pop_S_mat_fluscape.rda")

#' Firstly, we need to check that we can run this just for a small example. This test will also be
#' repeated in one of the package vingettes.
popsize_vector = values(popmatrix)
D=ncell(popsize_vector)
gz_pop_raster_agg <- aggregate(gz_pop_raster,popsize_vector,D,10)
pop_S_mat_fluscape_agg <- mob_calc_S_mat(gz_pop_raster_agg)
dim(pop_S_mat_fluscape_agg)

#' Now I need to load up the pre-calculated object and compare it with the size of this then
#' find the right size of matrix so I can check the radiation model runs.
pop_S_mat_fluscape <- readRDS( "~/dbox/projects/mobility/socio-spatial-behaviour/tmpdata/read_et_al_S_margin.rds")


#' ## Code below here from Jon and must be an alternative way of making the S matrix. I just need to figure ot out. I will only need a little bit of it

#Next - check that the gravity model produces sensible results and try to calculate the likelihood of the
# data

# Up to here XXXX needs to work below

fnLog <- paste0(local_data_dir,"/gravity_log_debug.csv")
if (!file.exists(fnLog)) {
    dftmp <- make.data.df(nrow=0)
    write.csv(dftmp,fnLog,row.names=FALSE)
}

tmp1 <- fit.mobility.model(
        contacts_fluscape_V1,
        gz_pop_raster,
        Smat = pop_S_mat_fluscape,
        logfile = fnLog,
        optfun = fit.offset.radiation.optim,
        psToFit = c("offset"),
        psLB = c(1),
        psUB = c(20*1000),
        datasubset = "ALL",
        fdebug=TRUE,
        lognote = ""
    )

comp <- fit.mobility.model(
  contacts_fluscape_V1,
  gz_pop_raster,
  logfile = "~/dbox/tmp/gravitylog.txt",
  optfun = fit.gravity.poppower.optim.nowithinregion,
  psToFit = c("Power", "DestPower"),
  psLB = c(0.1, 0.1),
  psUB = c(6, 4),
  datasubset = "ALL",
  pJustLike = c(2.646, 0.51375)
)


# Clear all other previous symbols and set a local working directory
rm(list=ls(all=TRUE))
setwd("/Users/sriley/Dropbox/svneclipse/fluscape/manuscripts/gravity")
# setwd("/home/sriley/Dropbox/svneclipse/fluscape/manuscripts/gravity")
# setwd("D:/Users/jon/Dropbox/work/Fluscape/manuscripts/gravity")

# source("../../../idsource/R/stevensRfunctions.R")
source("http://tinyurl.com/5t7gwnv")

datadir <- "../../data/"

# Load up the required study data
contacts <- read.csv(paste(datadir,"destination_data.csv",sep=""))
participants <- read.csv(paste(datadir,"Participants_V1.csv",sep=""))
households <- read.csv(paste(datadir,"HouseHolds_V1.csv",sep=""))

# Make the first population transect
require("raster") # help("raster-package")
require("rgdal")

#
contacts_tmp <- contacts[contacts$LOC_ID<=5 & contacts$LOC_ID!=2,]

contacts_small <- contacts_tmp[!is.na(contacts_tmp$long) | !is.na(contacts_tmp$lat) | !is.na(contacts_tmp$HH_Long) | !is.na(contacts_tmp$HH_Lat),]
range(contacts_small$long)
range(contacts_small$lat)
originx <- contacts_small$HH_Long
originy <- contacts_small$HH_Lat
dim(contacts_small)

x1 <- fsc.load.wide.raster()
 # JON: x1 <- fsc.load.wide.raster(fluscapetopdir="D:/Users/jon/svn folder/")
# ext <- extent(113.2,114.5,22.5,23.5) # suitable for LOC_ID==1
 ext <- extent(113.2,114.5,22.7,23.6) # suitable for LOC_ID<=5 & LOC_ID!=2
x2 <- crop(x1,ext)

#destindices <- fsc.gen.dest.index(x2)
#originindices <- fsc.gen.origin.index(originx,originy,x2)
#distances <- fsc.gen.dist.matrix(x2,originindices,destindices)

# faster non-zero destination indices calulation:
tmp = getValues(x2)
cells = which( !is.na(tmp) & tmp>0 )
destindices <- data.frame(index=1:length(cells), cells=cells)

# faster origin indices calulation:
cells = unique( cellFromXY(x2, contacts_small[,c("HH_Long","HH_Lat")] ) )
originindices <- data.frame(index=1:length(cells), cells=cells)

# faster distances calculation:
a = xyFromCell( x2, destindices$cell )
b = xyFromCell( x2, originindices$cell )
distances = pointDistance( b, a, longlat=T )



# save.image()
# rm(list=ls(all=TRUE))
# load(".RData")

# Next steps
# - generate a matrix of observations
# - generate a gravity model
# - compare the gravity models with the data

noorig 		<- dim(originindices)[1]
nodest 		<- dim(destindices)[1]

# observations
contactindices = cellFromXY( x2, contacts_small[,c("long","lat")] )
hhindices = cellFromXY( x2, contacts_small[,c("HH_Long","HH_Lat")] )
obs.tab = matrix( 0, nrow=noorig, ncol=nodest  )
for (k in 1:noorig) {
    j = originindices$cell[k]
	i = contactindices[hhindices==j]
	i = i[!is.na(i)]
	z = table(i)
	i.indices = match(as.numeric(names(z)),destindices$cells)
	obs.tab[k,i.indices] = as.numeric(z)
}


# faster gravity model generation
M <- 1
gravmodel <- matrix( nrow=noorig, ncol=nodest )
for (i in 1:noorig) {
		orcell 	<- originindices$cell[i]
		dscell 	<- destindices$cell
		n_o 	<- extract(x2,orcell)
		n_d 	<- extract(x2,dscell)
		dist 	<- distances[i, ]
		gravmodel[i, ] <- (n_o * n_d) / (500 + dist^M)
}
for (i in 1:noorig) {
	gravmodel[i,] <- gravmodel[i,] / sum(gravmodel[i,])
}


### WORK IN PROGRESS... NOT WORKING. IGNORE. START
#cells_in_range <- function( i, j, dist )
#{ destindices$cell[ which(dist<=dist[j]) ]
#}
#sum.n <- function( x2, dist, j )
#{
#  cells <- destindices$cell[ which(dist<=dist[j]) ]
#  sum(extract(x2,cells))
#}
#S <- matrix( nrow=noorig, ncol=nodest )
#i = 1
#	dist <- distances[i, ]
#	j = 5438
#	no_ij <- sum.n( x2, dist, j )
#
#	no_ij <- sum.n( x2, i, j )
### WORK IN PROGRESS... NOT WORKING. IGNORE. END



# find radiation density sums
S <- matrix( nrow=noorig, ncol=nodest )
pb = txtProgressBar(min = 0, max = noorig*nodest, style = 3); k = 0
for (i in 1:noorig) {
	dist <- distances[i, ]
	for (j in 1:nodest) {
	  k = k + 1; setTxtProgressBar( pb, k )
	  k = k + 1
		radcell <- destindices$cell[ which(dist<=dist[j]) ] # all dest cells within radius dist
		S[i,j] = sum( extract(x2,radcell) )
		setTxtProgressBar( pb, k )
	}
}
close(pb)
# write.csv(S, "S.csv")
# S <- as.matrix( read.csv("S.csv") ) # To read back in.



# fast radiation model, given S_ij
radmodel <- matrix( nrow=noorig, ncol=nodest)
for (i in 1:noorig) {
	orcell 	<- originindices$cell[i]
	dscell 	<- destindices$cell
	n_o 	<- extract(x2,orcell)
	n_d 	<- extract(x2,dscell)
	n_rad <- S[i,]
	radmodel[i,] <- (n_o * n_d) / (n_d + S[i,])*(n_d + n_o + S[i,])
}
for (i in 1:noorig) {
	radmodel[i, ] <- radmodel[i, ] / sum(radmodel[i, ])
}
# write.csv( radmodel, file="radmodel.csv", row.names=F, col.names=F )

par(mfcol=c(2,nrow(originindices)), mar=c(1,1,1,1) )
xy = xyFromCell(x2,destindices$cell)
for (i in 1:nrow(originindices))
{ col.seq = radmodel[i,]/max(radmodel[i,] )
  plot( xy, pch=15, cex=.6, col=rgb(col.seq,0,1-col.seq) )
  points( xyFromCell( x2, originindices$cell[i] ), pch=3, cex=1, col=3 )
  points( xy[obs.tab[i,]>0,], pch=1, cex=.75, col=3 )
  plot( distances[i,]/1000, radmodel[i,], pch=16, col=rgb(0,0,0,.1) )
 }


# errors/warnings:
# 1: In extract(x2, radcell) :
#  returning values at CELL NUMBERS (not coordinates) : 1919 and 1920
#2: In extract(x2, radcell) :
#  returning values at CELL NUMBERS (not coordinates) : 1921 and 1922
#3: In extract(x2, radcell) :
#  returning values at CELL NUMBERS (not coordinates) : 2075 and 2076
#4: In extract(x2, radcell) :
#  returning values at CELL NUMBERS (not coordinates) : 2077 and 2078




fsc.calc.model <- function()

# Why are there zeros above!
# The zeros above need to be fixed first

# !!!!! Starting code not run
# It should be possible to adapt the code below
# to calculate the log likelihood of the data given
# a gravity type model
dataMatrix <- gravmodel
dataMatrix[,] <- 1
veccurrentN <- colSums(dataMatrix)
veccurrentP <- vector(length=noorig)
veccurrentP[] <- 1
lnlike = 0
for (i in 1:(nodest-1)) {
	for (j in 1:noorig) {
		n <- dataMatrix[i,j]
		if (n > 0) {
			p <- gravmodel[i,j]
			lnlike <- lnlike + dbinom(n,veccurrentN[j],p/veccurrentP[j],log=TRUE)
			if (is.na(lnlike)) browser()
			veccurrentN[j] <- veccurrentN[j] - n
			veccurrentP[j] <- veccurrentP[j] - p
			if (veccurrentP[j] < -1e100) browser()
		}
	}
}
# !!!!! Ending code not run
