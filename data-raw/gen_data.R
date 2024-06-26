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
local_data_dir <- "~/tmp" # Steven
local_data_dir <- "D:/tmp" # Jon
local_data_dir <- "G:/OneDrive - Imperial College London/Postdoc/Fluscape/fluscapeR/fluscapeR/" # Vivi
fluscape_top_dir <- "../fluscape/"
fluscape_top_dir <- "G:/OneDrive - Imperial College London/Postdoc/Fluscape/fluscape/" # Vivi

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
# part_all <- load_particpant_data_long()
# cont_all <- read.csv("../fluscape/data/clean_datasets/ContactsAll.csv")
# dim(part_all)
# names(part_all)
#' Keep old data line in just in case. Not run.
## contacts_0 <- mob_load_old_contact_data( locations, households, participants )

#' If you have saved contacts_fluscape_V1 in your local,
#' you can load this contacts_fluscape_V1.rda directly and skip line 71-80,
#' otherwise you need to run line 71-80 and save it as a local file
#' In the future, we may intergrate it in the pacakge.
load(file=paste0(local_data_dir,"/","contacts_fluscape_V1.rda"))

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
  plot( contacts_fluscape_V1$HH_Long, contacts_fluscape_V1$HH_Lat, pch=16, cex=.3, col="blue", xlim=c(113.5,113.6), ylim=c(23.2,23.3))
  points( contacts_jit$HH_Long, contacts_jit$HH_Lat, pch=16, cex=.3, col="red")

#' rename and write the potentially personally identifiable version to the local temp directory
  save(contacts_fluscape_V1,file=paste0(local_data_dir,"/","contacts_fluscape_V1.rda"))

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

#' Brings in the previously generated S matrix for called pop_S_mat_fluscape
  load("~/dbox/shares/me_jr_hm_dc_gravity/current/to_upload/pop_S_mat_fluscape.rda") # steven
  load("G:/OneDrive - Imperial College London/Postdoc/Fluscape/fluscapeR/fluscapeR/pop_S_mat_fluscape.rda") # vivi
  load("D:/tmp/pop_S_mat_fluscape.rda") # jon
  # transfer this binary to jons machine
  str(pop_S_mat_fluscape)

#' Firstly, we need to check that we can run this just for a small example. This test will also be
#' repeated in one of the package vingettes.
  popsize_vector = values(popmatrix)
  D=ncell(popsize_vector)
  gz_pop_raster_agg <- aggregate(gz_pop_raster,popsize_vector,D,10)
  pop_S_mat_fluscape_agg <- mob_calc_S_mat(gz_pop_raster_agg)
  dim(pop_S_mat_fluscape_agg)

#' Now I need to load up the pre-calculated object and compare it with the size of this then
#' find the right size of matrix so I can check the radiation model runs.
#  pop_S_mat_fluscape <- readRDS( "~/dbox/projects/mobility/socio-spatial-behaviour/tmpdata/read_et_al_S_margin.rds")


#' ## Code below here from Jon and must be an alternative way of making the S matrix. I just need to figure ot out. I will only need a little bit of it

#Next - check that the gravity model produces sensible results and try to calculate the likelihood of the
# data

# Up to here XXXX needs to work below

fnLog <- paste0(local_data_dir,"/gravity_log_debug.csv")
if (!file.exists(fnLog)) {
    dftmp <- make.data.df(nrow=0)
    write.csv(dftmp,fnLog,row.names=FALSE)
}

# fit radiation model to previous existing instance of S matrix, pop_S_mat_fluscape
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


#------------------ works up to here.... ---------------------#


contacts <- contacts_fluscape_V1
  # contacts_tmp <- contacts[contacts$LOC_ID<=5 & contacts$LOC_ID!=2,]
  contacts_tmp <- contacts

  contacts_small <- contacts_tmp[!is.na(contacts_tmp$long) | !is.na(contacts_tmp$lat) | !is.na(contacts_tmp$HH_Long) | !is.na(contacts_tmp$HH_Lat),]
  range(contacts_small$long)
  range(contacts_small$lat)
  originx <- contacts_small$HH_Long
  originy <- contacts_small$HH_Lat
  dim(contacts_small)

  # Find the spatial extent of destination long-lats, and crop landscan density raster.
    # ext <- extent(113.2,114.5,22.5,23.5) # suitable for LOC_ID==1
    # ext <- extent(113.2,114.5,22.7,23.6) # suitable for LOC_ID<=5 & LOC_ID!=2
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
S.gravity <- generate.agg.Smatrix(contacts_fluscape_V1, x2, 0, "gravity")
S.raiation <- generate.agg.Smatrix(contacts_fluscape_V1, x2, 10, "radiation")
#' generate.agg.Smatrix includes codes from line 198 to line 299.
#' People can test if generate.agg.Smatrix() generates the same results as line 198-299 does
#' by using fundtion identical().

#destindices <- fsc.gen.dest.index(x2)
#originindices <- fsc.gen.origin.index(originx,originy,x2)
#distances <- fsc.gen.dist.matrix(x2,originindices,destindices)

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


# faster gravity model generation
M <- 1
# Make a matrix of predicted origin-destination weights for gravity model
#   with an offset of 500 and a power term M.
  gravmodel <- matrix( nrow=noorig, ncol=nodest )
  pb = txtProgressBar(min = 0, max = noorig, initial = 0)
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
  pb = txtProgressBar(min = 0, max = noorig, initial = 0)
  for (i in 1:noorig) {
  	gravmodel[i,] <- gravmodel[i,] / sum(gravmodel[i,])
  	setTxtProgressBar(pb,i)
  }
  close(pb)

  #' test if generate.agg.S.radiation() can generate the same result
  identical(S.gravity, gravmodel)

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

# Find radiation density sums, s_ij, for each orig-dest cell combination.
#   Finds the total population density within equivalent distance from the origin.
# Warning, this is slow and takes a long time.
#   Save to disk after generation for future use.
  S <- matrix( nrow=noorig, ncol=nodest )
  for (i in 1:noorig) {
    print(paste(i,"/",noorig))
  	dist <- distances[i, ]
  	pb = txtProgressBar(min = 0, max = nodest, style = 3)
  	for (j in 1:nodest) {
  	  within_d = which(dist<=dist[j])
  	  radcell <- destindices$cell[within_d] # all destination cells within radius dist
  		S[i,j] = sum( extract(x2,radcell) )
  		setTxtProgressBar( pb, j )
  	  close(pb)
  	}
  }
# write.csv(S, "S.csv", row.names=FALSE)
# S <- as.matrix( read.csv("S.csv") ) # To read back in.

# jon: see D:\Users\jon\Dropbox\work\Fluscape\source\R\gravy_functions.R
  gravy.fast.calc.S <- function( distances, x2, originindices, noorig, nodest,
                                 write.to.file=TRUE, filename="auxilliary/tmp_S.csv" ) {
    # Faster function to calculate S_ij, the sum of population densities within radius r_ij,
    # 	for each origin-destination pair.
    # Inputs:
    # 	distances
    #	x2
    #	originindices
    #	noorig
    #	nodest
    # Outputs:
    #	S -- matrix of S_ij, where rows correspond to originindices, and columns to destinations.
    max.r = max(distances) # maximum distance
    max.rij = which( distances==max.r, arr.ind=TRUE  ) # which pair is it?
    # ADD CHECK ***HERE*** THAT x2 IS BIG ENOUGH GIVEN ORIGINS AND max.rij
      print("Warning: no automatic check that x2 is large enough to generate correct S_ij.")
    allindices <- gravy.gen.dest.index( x2 )
    all.distances <- gravy.gen.dist.matrix( x2, originindices, allindices ) # distances from origins to all possible cells
    maxj = nrow(allindices)
    S = matrix( 0, nrow=noorig, ncol=maxj )
    for (i in 1:noorig) {
      png(filename="S_prog.png")
      plot( i, main=paste( "i =", i ) )
      dev.off()
      cell.i = originindices$cells[i]
      n_i = extract( x2, cell.i ) # density in cell i
      cell.j = destindices[ ,2]
      r.ij = all.distances[i, ]
      rank.r.ij <- rank( r.ij, ties.method="min" )
      d.ij <- data.frame( j=cell.j, rank.d=rank.r.ij, d=r.ij, n=0 )
      d.ij = d.ij[order(d.ij$d), ]
      k.seq = unique(d.ij$rank)
      # k: 1
      nk = length(k.seq)
      k1 = which( all.distances[i,]>=0 &
                    all.distances[i,]<=d.ij$d[1] )
      j = allindices$cells[ k1 ]
      n.j = extract( x2, j )
      d.ij$S[1] = sum( n.j )
      # k: 2 ... nk
      pb = txtProgressBar(min=0,max=nk,style=3)
      for (k in k.seq[2:nk]) {
        setTxtProgressBar(pb,k)
        k1 = which( all.distances[i,]>d.ij$d[k-1] &
                      all.distances[i,]<=d.ij$d[k] )
        j = allindices$cells[ k1 ]
        n.j = extract( x2, j )
        d.ij$n[k] = sum( n.j )
      }
      close(pb)
      d.ij$S = cumsum( d.ij$n )
      k2 = match( d.ij$j, allindices$cells )
      S[i,k2] <- d.ij$S

    }
    if (write.to.file==TRUE) {
      write.csv( S, filename, row.names=F )
    }
    return( S )
  }





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
