## Copyright Steven Riley (sr@stevenriley.net), Harriet L. Mills and
## Jonathan Read.
##
## This file is part of the library fluscapeR.
##
## fluscapeR is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This work is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this idsource.  If not, see <http://www.gnu.org/licenses/>.
## (at your option) any later version.

#' @import raster
fit.mobility.model <- function(contacts,
                               popgrid,
                               sd=928734924,
                               noRepeats = 1,
                               optfun = NULL,
                               psToFit = c("kernpower","offset"),
                               psLB = c(0.1,0.001),
                               psUB = c(6,100),
                               datasubset="ALL",
                               logfile="./gravity_logfile.csv",
                               Smat = NULL,
                               pJustLike = NULL,
                               lognote="",
                               fdebug=FALSE) {


    ## Check some preconditions for the function arguments
    if (is.null(optfun) && is.null(Smat)) {
        error("S matrix must be specified for the radiation model")
    }
    ## if (psToFit !%in% c("Power","Offset","DestPower")) {
    ##     error("Parameters must be one of: Power, Offset, DestPower")
    ## }

    ## Set he seed
    set.seed(sd)

    ## Check for the log file and make sure it has the correct headers
    if (!file.exists(logfile)) stop("Logfile must exist")
    testdf <- read.csv(logfile)
    testnames <- names(testdf)
    targetnames <- names(make.data.df())
    if (!all(testnames == targetnames)) {
        errormessage <- paste(
            "logfile anmes are not consistent with current",
            "version of the function")
        stop(errormessage)
    }

    ## Select type of contact
    if (datasubset=="ALL") {
        contacts_type <- contacts
    } else if (datasubset=="CHILDREN") {
        contacts_type <- contacts[contacts$CHILD==TRUE,]
    } else if (datasubset=="ADULTS") {
        contacts_type <- contacts[contacts$CHILD==FALSE,]
    } else if (datasubset=="URBAN") {
        contacts_type <- contacts[contacts$URBAN==TRUE,]
    } else if (datasubset=="RURAL") {
        contacts_type <- contacts[contacts$URBAN==FALSE,]
    } else {stop("datasubset not recognised")}

    ## Subset for contacts within the loaded gridsquare
    contacts_small <- gravy.subset.contact.data(
        contacts_type$pid, contacts_type, popgrid
    )

    ## Origin and destination indices calulation
    ## all_originindices <- gravy.gen.origin.index(popgrid, contacts_small)
    originx <- contacts_small$HH_Long
    originy <- contacts_small$HH_Lat
    destindices <- gravy.gen.dest.index(popgrid)
    originindices <- gravy.gen.origin.index(popgrid, contacts_small)
    all_originindices <- gravy.gen.origin.index(popgrid, contacts)

    ## Find the number of origin and destination cells
    noorig <- dim(originindices)[1]
    nodest <- dim(destindices)[1]

    ## Generate a table of the observations
    obs.tab <- gravy.gen.observations(
        popgrid, contacts_small, originindices, destindices, noorig, nodest
    )

    ## Calculation of distances could be cached
    ## But only a small proportion of typical run time
    distances <- gravy.gen.dist.matrix(popgrid, originindices, destindices)

    ## This seems to give the same results as the main table in the paper
    ## make this table be for all possible parameters and add the see and
    ## data subtype.
    ## Then append the table to the existing data
    ## Consider a switch to mcmc for this if the cis are still wrong
    nops <- length(psToFit)
    ## pvTab <- matrix(nrow=noRepeats,ncol=4*nops+1)
    ## colnames(pvTab) <- c(paste(
    ##    psToFit,
    ##    c(rep(".iv",nops),
    ##    rep(".pe",nops),
    ##    rep(".lb",nops),
    ##    rep(".ub",nops)),
    ##    sep=""),"lnlike")
    pvTab <- make.data.df(nrow=noRepeats)
    lnlike <- -999999999

    ## Insert an alternative table here, check that it works and then
    ## delete the other table only when done

    ## Run the radiation model with no repeats if no optimization
    ## function is specified
    if (nops < 1 && is.null(optfun)) {
        radiation.model = harriet.offset.radiation(
            popgrid, Smat, noorig, nodest, originindices, destindices,
            distances, obs.tab, all_originindices, offset=0
        )
        pvTab[,"lnlike"] <- rep(as.numeric(radiation.model),noRepeats)
        lnlike <- radiation.model

        ## If not pure ratiation, check to see if just likelihood needed for
        ## other model
    } else if (!is.null(pJustLike)) {

        ## Calculate just likelihood for the provided model
        ## output goes only into lnlike variable. Table ignored for these
        ## run options
        if (!fdebug) {
            lnlike <- optfun(
                pJustLike,
                psToFit=psToFit,
                originindices=originindices,
                all_originindices=all_originindices,
                destindices=destindices,
                distances=distances,
                noorig=noorig,
                nodest=nodest,
                obs.tab = obs.tab,
                x2=popgrid,
                S=Smat)
        } else {
            lnlike <- 0
        }

    } else {

        for (i in 1:noRepeats) {

            ## randomly chosen initial conditions
            psInitial = psLB + (psUB-psLB)*runif(length(psUB))

            ## Test if more than one parameter to choose optimization routine
            if (nops > 1) {

                ## Run multidimensional optimization
                if (fdebug) {

                    pests <- -psInitial[]
                    maxlike <-  0

                } else {

                    fit_gravity <- optimx(
                        psInitial,
                        optfun,
                        method="L-BFGS-B",
                        lower=psLB,
                        upper=psUB,
                        control=list(
                            trace=0,
                            fnscale=-1,
                            kkt=FALSE
                        ),
                        psToFit=psToFit,
                        originindices=originindices,
                        all_originindices=all_originindices,
                        destindices=destindices,
                        distances=distances,
                        noorig=noorig,
                        nodest=nodest,
                        obs.tab = obs.tab,
                        x2=popgrid,
                        S=Smat
                    )

                    ## Extract output for table
                    ## This can be changed for the alternative table
                    pests <- as.numeric(fit_gravity[1,1:nops])
                    maxlike <-  as.numeric(fit_gravity[1,"value"])

                }

                ## Need to use parameter names here to assign values, not
                ## based just on the location in the old pvTab
                ## There are about three other places I need to do this as well

                ## Below here maybe should be a function, but it seems like a waste to
                ## have it evaluated so often.
                pvTab[i,"lnlike"] <- maxlike
                for (j in 1:nops) {
                    colinit <- paste(psToFit[j],".iv",sep="")
                    colpe <- paste(psToFit[j],".pe",sep="")
                    pvTab[i,colinit] <- psInitial[j]
                    pvTab[i,colpe] <- pests[j]
                }

                ## pvTab[i,1:nops] <- psInitial[]
                ## pvTab[i,(nops+1):(2*nops)] <- pests
                ## pvTab[i,4*nops+1] <- maxlike

            } else {

                if (fdebug) {

                    ## Extract required variables from output
                    pests <- -psInitial[]
                    maxlike <-  0

                } else {

                    ## Run one dimensional optimisation
                    fit_gravity <- optimise(
                        optfun,
                        maximum = TRUE,
                        psToFit = psToFit,
                        originindices = originindices,
                        all_originindices = all_originindices,
                        destindices=destindices,
                        distances=distances,
                        noorig=noorig,
                        nodest=nodest,
                        obs.tab = obs.tab,
                        x2=popgrid,
                        S=Smat,
                        lower=psLB,
                        upper=psUB
                    )

                    ## Extract required variables from output
                    pests <- as.numeric(fit_gravity$maximum)
                    maxlike <-  as.numeric(fit_gravity$objective)

                }

                ## Below here maybe should be a function, but it seems like a waste to
                ## have it evaluated so often.
                pvTab[i,"lnlike"] <- maxlike
                for (j in 1:nops) {
                    colinit <- paste(psToFit[j],".iv",sep="")
                    colpe <- paste(psToFit[j],".pe",sep="")
                    pvTab[i,colinit] <- psInitial[j]
                    pvTab[i,colpe] <- pests[j]
                }

            }

            if (is.null(pJustLike) & !fdebug) {

                ## Add the univariate confidence bounds
                ## Defined to be enclosed in this function
                fCIs <- function(v,pind,maxl,offset=1.96) {
                    ps <- pests
                    ps[pind] <- v
                    val <- optfun(
                        ps,
                        psToFit=psToFit,
                        originindices=originindices,
                        all_originindices=all_originindices,
                        destindices=destindices,
                        distances=distances,
                        noorig=noorig,
                        nodest=nodest,
                        obs.tab = obs.tab,
                        x2=popgrid,
                        S=Smat
                    )

                    ## The value needs to be +offset if the value is equal
                    ## to max like. It needs to be 0 if the value is offset
                    ## below max like and it needs to be below zero if the
                    ## value is more than offset below max like. Remembering
                    ## that val is a log likelihood.
                    rtn <- maxl - val - offset
                    rtn

                }

                ## Find univariate confidence bounds
                ## Test for optimization bounds being inside CIs to
                ## Skecth this out in notebook to be sure it makes sense
                cioffset <- 1.96
                for (ip in 1:nops) {
                    colnamelb <- paste(psToFit[ip],".lb",sep="")
                    colnameub <- paste(psToFit[ip],".ub",sep="")
                    val_lb <- fCIs(psLB[ip],ip,maxlike,offset=cioffset)
                    val_ub <- fCIs(psUB[ip],ip,maxlike,offset=cioffset)
                    ## debug <- fCIs(pests[ip],ip,maxlike,offset=cioffset)
                    if (val_lb > -0.0000001) {
                        lb <- uniroot(
                            fCIs,
                            interval=c(psLB[ip],pests[ip]),
                            pind=ip,maxl=maxlike, offset=cioffset)
                        pvTab[i,colnamelb] <- lb$root
                    } else {
                        pvTab[i,colnamelb] <- psLB[ip]
                    }
                    if (val_ub > -0.0000001) {
                        ub <- uniroot(
                            fCIs,interval=c(pests[ip],psUB[ip]),
                            pind=ip,maxl=maxlike, offset=cioffset)
                        pvTab[i,colnameub] <- ub$root
                    } else {
                        pvTab[i,colnameub] <- psUB[ip]
                    }
                }
            }

            ## Close loop for number of repeats
        }
    }

    ## Record to log file if needed
    if (!is.null(logfile)) {

        ## If validation of column names needed
        ## can be done at the head of this function
        tmpinfo <- Sys.info()
        pvTab[,"day.time"] <- rep(as.character(Sys.time()),noRepeats)
        pvTab[,"node"] <- rep(as.character(tmpinfo["nodename"]),noRepeats)
        pvTab[,"datasubset"] <- rep(as.character(datasubset),noRepeats)
        pvTab[,"datatype"] <- rep(as.character(lognote),noRepeats)
        write.table(pvTab,
                    file=logfile,
                    row.names=FALSE,
                    col.names=FALSE,
                    append=TRUE,
                    sep=",")

        ## tmpTab <- make.data.df(nrow=1)
        ## tmpTab[1,"datasubset"] <- datasubset
        ## tmpTab[1,"likelihood"] <- as.numeric(1)

        ## write(datasubset,file=logfile,append=TRUE)
        ## write(psToFit,file=logfile,append=TRUE)
        ## write.table(pvTab,file=logfile,row.names=FALSE,
        ##            col.names=TRUE, append=TRUE)
        ## write(as.character(Sys.time()),file=logfile,append=TRUE)
        ## write(paste("Note: ",lognote,"\n"),file=logfile,append=TRUE)
        ## write("\n\n\n",file=logfile,append=TRUE)

    }

    ## Return function value
    list(tab=pvTab,like=lnlike)

}

make.data.df <- function(nrow=0) {
    colnames <- c(
        "node","day.time","datasubset","datatype","lnlike",
        "destpower.pe","destpower.lb","destpower.ub","destpower.iv",
        "kernpower.pe","kernpower.lb","kernpower.ub","kernpower.iv",
        "offset.pe","offset.lb","offset.ub","offset.iv"
    )
    ncols <- length(colnames)
    tmpmat <- matrix(nrow=nrow,ncol=ncols)
    rtn <- as.data.frame(tmpmat)
    names(rtn) <- colnames
    rtn
}

fit.gravity.optim.nowithinregion <- function(
  fitpars,
  psToFit ,
  originindices,
  all_originindices=NULL,
  destindices,
  distances,
  noorig,
  nodest,
  obs.tab,
  x2,
  S){
  # for use with optim, sorts out the parameters and calls "gravity" function
  #fitpars is vector of fitted parameter values
  #psToFit is vector of fitted parameter names in correct order c("Power", "Offset")
  #  contacts -- data.frame, contact data, can be subset to particular types of contact
  #  long.lim, lat.lim -- vectors of length 2, giving bounding limits for long/lat area
  #  x2 -- landscane density object
  #  pid -- IDs of participants to use in this fit.

  # USed to be an option
  gravitymodel="Harriet"

  #determine the parameters
  #power
  kernpower=fitpars[which(psToFit=="kernpower")]

  #offset
  if (gravitymodel == "Harriet"){ #correct model

    if (length(which(psToFit=="offset"))==0){ # gravity model
      OffsetAB=c(0, 1) # not fitting any offset
    } else { # offset gravity model
      OffsetAB = c(1, fitpars[which(psToFit=="offset")])
    }

  } else{ #Jon's model

    if (length(which(psToFit=="offset"))==0){ # gravity model
      OffsetAB=c(1, 1) # not fitting any offset
    } else {
      OffsetAB = c(fitpars[which(psToFit=="offset")], 1)
    }

  }
  #find the lnlikelihood
  lnlike = harriet.gravity.nowithinregion(
    originindices,
    destindices,
    distances,
    noorig,
    nodest,
    obs.tab,
    x2,
    kernpower = kernpower,
    OffsetAB= OffsetAB)

  cat(lnlike,kernpower,OffsetAB,"\n")

  return(lnlike)
}

fit.gravity.poppower.optim.nowithinregion <- function(
  fitpars,
  psToFit,
  originindices,
  all_originindices,
  destindices,
  distances,
  noorig,
  nodest,
  obs.tab,
  x2,
  S ){
  ## for use with optim, sorts out the parameters and calls "gravity" function
  ## fitpars is vector of fitted parameter values
  ## psToFit is vector of fitted parameter names in correct order c("Power", "Offset")
  ##  contacts -- data.frame, contact data, can be subset to particular types of contact
  ##  long.lim, lat.lim -- vectors of length 2, giving bounding limits for long/lat area
  ##  x2 -- landscane density object
  ##  pid -- IDs of participants to use in this fit.

  ## Used to be an option
  gravitymodel <- "Harriet"

    ## determine the parameters
    ## power

  kernpower=fitpars[which(psToFit=="kernpower")]

  if (length(which(psToFit=="OriginPower"))==0){ # origin powers
    OriginPower = 1 # not fitting any origin power
  } else { # fit origin powers
    OriginPower = fitpars[which(psToFit=="OriginPower")]
  }

    if (length(which(psToFit=="destpower"))==0){ # dest powers
        destpower = 1 # not fitting any dest power
  } else { # fit origin powers
    destpower = fitpars[which(psToFit=="destpower")]
  }

  #offset
  if (gravitymodel == "Harriet"){ #correct model

    if (length(which(psToFit=="offset"))==0){ # gravity model
      OffsetAB=c(0, 1) # not fitting any offset
    } else { # offset gravity model
      OffsetAB = c(1, fitpars[which(psToFit=="offset")])
    }

  } else{ #Jon's model

    if (length(which(psToFit=="offset"))==0){ # gravity model
      OffsetAB=c(1, 1) # not fitting any offset
    } else {
      OffsetAB = c(fitpars[which(psToFit=="offset")], 1)
    }

  }

                                        #find the lnlikelihood

    lnlike = harriet.gravity.poppower.nowithinregion( originindices,
                                                     destindices,
                                                     distances,
                                                     noorig,
                                                     nodest,
                                                     obs.tab,
                                                     x2,
                                                     kernpower = kernpower,
                                                     OffsetAB= OffsetAB,
                                                     OriginPower = OriginPower,
                                                     destpower = destpower)

    cat(lnlike,"kernpower",kernpower,"offset",OffsetAB,"Origin Power",
        OriginPower,"destpower",destpower,"\n")

  return(lnlike)

}

fit.offset.radiation.optim <- function(
  fitpars,
  psToFit,
  originindices,
  all_originindices,
  destindices,
  distances,
  noorig,
  nodest,
  obs.tab,
  x2,
  S){

  # for use with optim, sorts out the parameters and calls "gravity" function
  # fitpars is vector of fitted parameter values
  # psToFit is vector of fitted parameter names in correct order
  # c("vecPower", "vecOffset")
  # contacts -- data.frame, contact data, can be subset to particular types of contact
  # long.lim, lat.lim -- vectors of length 2, giving bounding limits for long/lat area
  # x2 -- landscane density object
  # pid -- IDs of participants to use in this fit.

  #determine the parameter
  offset=fitpars[which(psToFit=="offset")]

  #find the lnlikelihood
  lnlike = harriet.offset.radiation(x2,
                                    S,
                                    noorig,
                                    nodest,
                                    originindices,
                                    destindices,
                                    distances,
                                    obs.tab,
                                    all_originindices,
                                    offset)


  cat(lnlike," ",offset," ",fitpars,"\n")

  return(lnlike)
}


harriet.gravy.gravity.model.vcorrect.nowithinregion <- function(
  kernpower, OffsetAB, originindices, destindices,
  distances, noorig, nodest, poprast ) {

  # generate the probabilities for each destination cell, for each origin cell, based
  # on simple gravity model, with single parameter M.
  # NOTE, distances are expected to be in metres!

  # now calculate model
  OffsetA = OffsetAB[1]
  OffsetB = 10^(OffsetAB[2])

  gravmodel <- matrix( nrow=noorig, ncol=nodest )
  dscell <- destindices$cell

  n_d <- extract(poprast,dscell)
  for (i in 1:noorig) {
    orcell   <- originindices$cell[i]
    n_o 	<- extract(poprast, orcell)
    d_ij 	<- distances[i, ]/1000 # convert from metres to km
    gravmodel[i, ] <- (n_o * n_d) / (OffsetA+(d_ij/OffsetB)^kernpower) #our version
    # gravmodel[i, ] <- (n_o * n_d) / (OffsetA+d_ij/OffsetB)^Power #Truscott and Ferguson version
    # remove the within region mixing
    ind = dscell == orcell # gives indicator when same
    gravmodel[i, ind] = 0
  }

  # individual-level normalisation
  for (i in 1:noorig) {
    gravmodel[i,] <- gravmodel[i,] / sum(gravmodel[i,])
  }
  return( gravmodel )
}


mob_calc_S_mat_first <- function(popmatrix, A=1) {
  # Faster function to calculate S_ij, the sum of population densities
  # within radius r_ij,
  # for each origin-destination pair.
  # Inputs:
  # popmatrix is raster
  # D is number of districts
  # A is number of age groups
  # Outputs:
  #  S -- matrix of S_ij, where rows correspond
  # to originindices, and columns to destinations.

  require("raster")

  ## Re engineer arguments
  popsize_vector = values(popmatrix)
  D=ncell(popsize_vector)

  # within region mixing
  r_forfirstcell=raster::values(raster::distanceFromPoints(popmatrix, xyFromCell(popmatrix,1)))

  # distance between first and second
  xlength=r_forfirstcell[2]

  #distance between first and the second row
  ylength=r_forfirstcell[ncol(popmatrix)+1]
  dwr_r=distancewithinarectangle(xlength,ylength)

  S = matrix( 0, nrow=D*A, ncol=D*A)
  i=1
  while (i <= (D*A)) {
    #print(i)
    ai=i%%A; ai[ai==0]=A
    di=(i-ai)/A +1

    n_orig = popsize_vector[i] # density in cell i

    d.ij=raster::values(raster::distanceFromPoints(popmatrix, xyFromCell(popmatrix,di)))
    d.ij[di]=dwr_r #distance within that district
    popsize=popsize_vector
    destcells = 1:(D*A) #list of destination cells
    neworder=order(d.ij) #orders by the distance away
    d.ij = d.ij[neworder]
    popsize = popsize[neworder]
    destcells = destcells[neworder]

    S.ij = vector(len=D*A)
    rollingsum_ind = 1
    rollingsum = 0
    uni_dist = unique(d.ij)
    for (k in 1:length(uni_dist)){

      # finds cells which are the same distance away
      samedist_ind = which(d.ij==uni_dist[k])
      last_ind=samedist_ind[length(samedist_ind)]
      rollingsum = rollingsum + sum(popsize[rollingsum_ind:last_ind]) # the minimum sum up to this distance
      rollingsum_ind = last_ind + 1
      S.ij[samedist_ind] = rollingsum
    }
    S.ij = S.ij - n_orig - popsize # remove the origin and destination cell
    #NB there will be some negative cells
    S[i, ] = S.ij[order(destcells)] # reorder into correct order

    i=i+1
  }

  return( S )

}

harriet.offset.radiation <- function(
  x2, S, noorig, nodest, originindices,
  destindices, distances, obs.tab, all_originindices, offset) {

  # Inputs:

  # Outputs:
  #	rad.res -- data.frame, power|offset|loglikelihood.

  radmodel <- gravy.radiation.model.offset.harriet ( x2, S, noorig, nodest, originindices, destindices, distances, all_originindices, offset )

  lnlike <- gravy.calc.lnlike.nowithinregion.harriet( obs.tab, radmodel, noorig , originindices, destindices)

  # print("Finished.")

  return( lnlike )

}

gravy.subset.contact.data <- function(
  pid, contacts, x2
) {
  # Funciton to subset the contact data, according to pid and further subset
  #   ensure only participants
  # which participants to include in contact data?
  # long.lim = c( 113.1, 113.85 )  # default for 1st pass
  # lat.lim = c( 22.98, 23.62 )  # default for 1st pass
  print("Subsetting contact data...")
  long.lim <- c((extent(x2))@xmin,(extent(x2))@xmax)
  lat.lim <- c((extent(x2))@ymin,(extent(x2))@ymax)
  pid = gravy.pid.latlong.lim(pid,contacts,long.lim, lat.lim )
  contacts_small <- gravy.subset.contacts( contacts, pid );
  print( "dim(contacts_subset):")
  print( dim(contacts_small) )

  # return( list(pid, contacts_small)  )
  contacts_small

}


gravy.pid.latlong.lim <- function( pid, contacts, long.lim, lat.lim ) {
  # returns vector of participant IDs (pid) for which their contacts
  #  are within the limits defined by long.lim and lat.lim.
  tmp.c = contacts
  long.min = tapply( tmp.c$long, tmp.c$pid, FUN=min, na.rm=T )
  long.max = tapply( tmp.c$long, tmp.c$pid, FUN=max, na.rm=T )
  lat.min = tapply( tmp.c$lat, tmp.c$pid, FUN=min, na.rm=T )
  lat.max = tapply( tmp.c$lat, tmp.c$pid, FUN=max, na.rm=T )
  table( long.min>=long.lim[1], long.max<=long.lim[2], useNA="always" )
  table( lat.min>=lat.lim[1], lat.max<=lat.lim[2] , useNA="always" )
  i = which( long.min>=long.lim[1] & !is.na(long.min) &
               long.max<=long.lim[2] & !is.na(long.max) &
               lat.min>=lat.lim[1] & !is.na(lat.min) &
               lat.max<=lat.lim[2] & !is.na(lat.max) 		)
  x = names(long.min)
  tmp.pid = x[i]
  new.pid = intersect( pid, tmp.pid )
  return( new.pid )
}

gravy.subset.contacts <- function( contacts, pid ) {
  # subsets the contact data, according to pid, and ensures no NAs in lat/longs.
  contacts_small <- contacts[is.element(contacts$pid, pid), ] # subset to just these participants
  return( contacts_small )
}

gravy.gen.dest.index <- function( x2 ) {
  # create a dataframe of indices
  # for all non-zero density cells in the landscan 'world' of x2.
  tmp = getValues(x2)
  cells = which( !is.na(tmp) & tmp>0 ) #
  destindices <- data.frame(index=1:length(cells), cells=cells)
  return( destindices )
}

gravy.gen.origin.index <- function( x2, contacts_small ) {
  # create a dataframe of indices
  # of cells which contain origin long/lats.
  cells = unique( cellFromXY(x2, contacts_small[,c("HH_Long","HH_Lat")] ) )
  originindices <- data.frame(index=1:length(cells), cells=cells)
  return( originindices )
}

gravy.gen.observations <- function(
                                   x2,
                                   contacts_small,
                                   originindices,
                                   destindices,
                                   noorig,
                                   nodest) {

    ## Declare required libraries
    require(raster)

    ## make a matrix containing the observed counts of contact occuring
    ## in destination cell by participants coming from origin cell
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

    return(obs.tab)

}

gravy.gen.dist.matrix <- function( x2, originindices, destindices ) {
  # find the distance between all origin and destination cell combinations
  # Note, this uses the distance between cell centres, not the distance between
  # actual origin and contact destination long/lats.
  a = xyFromCell( x2, destindices$cell )
  b = xyFromCell( x2, originindices$cell )
  distances = pointDistance( b, a, longlat=T )
  return( distances )
}

harriet.gravity.nowithinregion <- function(
                                           originindices,
                                           destindices,
                                           distances,
                                           noorig,
                                           nodest,
                                           obs.tab,
                                           poprast,
                                           kernpower, # single value
                                           OffsetAB #2 values
                                           ) {

    ## finds gravity model with lowest loglikelihood,
    ## with base power and offset terms.
    ## Inputs:
    ## 	Power --  single value
    ## 	OffsetAB -- 2 values
    ## Outputs:
    ##	grav.res -- just the value of the loglikelihood is returned.

    gravmod <- harriet.gravy.gravity.model.vcorrect.nowithinregion(
        kernpower, OffsetAB, originindices, destindices, distances,
        noorig, nodest, poprast )

    lnlike = gravy.calc.lnlike.nowithinregion.harriet( obs.tab,
                                                      gravmod,
                                                      noorig ,
                                                      originindices,
                                                      destindices)

    ## cat(lnlike," ",Power," ",OffsetAB,"\n")

    return( lnlike )

}

harriet.gravity.poppower.nowithinregion <- function(
                                                    originindices,
                                                    destindices,
                                                    distances,
                                                    noorig,
                                                    nodest,
                                                    obs.tab,
                                                    poprast,
                                                    kernpower, # single value
                                                    OffsetAB, #2 values
                                                    OriginPower,
                                                    destpower)
{
    ## finds gravity model with lowest loglikelihood,
    ## with base power and offset terms.
    ## Inputs:
    ## Power --  single value
    ## OffsetAB -- 2 values
    ## Outputs:
    ##grav.res -- just the value of the loglikelihood is returned.

    ## Need to find out here why the likelihood is not changing
    ## So using the lines below in an interactive way

    gravmod <- poppower_model(kernpower,
                              OffsetAB,
                              OriginPower,
                              destpower,
                              originindices,
                              destindices,
                              distances,
                              noorig,
                              nodest,
                              poprast)

    lnlike = gravy.calc.lnlike.nowithinregion.harriet(obs.tab,
                                                      gravmod,
                                                      noorig,
                                                      originindices,
                                                      destindices)

    return(lnlike)

}

## Move this around the code to whereever working
version_check <- function() {cat("0.90001\n")}

poppower_model <- function(kernpower,
                           OffsetAB,
                           OriginPower,
                           destpower,
                           originindices,
                           destindices,
                           distances,
                           noorig,
                           nodest,
                           poprast) {
    ## generate the probabilities for each destination cell,
    ## for each origin cell, based
    ## on simple gravity model, with single parameter M.
    ## NOTE, distances are expected to be in metres!

    ## now calculate model

    OffsetA = OffsetAB[1]
    OffsetB = 10^(OffsetAB[2])

    gravmodel <- matrix( nrow=noorig, ncol=nodest )
    dscell <- destindices$cell
    n_d <- extract(poprast,dscell)
    for (i in 1:noorig) {

        orcell   <- originindices$cell[i]
        n_o 	<- extract(poprast, orcell)

        ## This is where the problem was
        ## Either convert to metres or do not
        ## d_ij 	<- distances[i, ]/1000
        d_ij 	<- distances[i, ] / 1000
        ## our version

        gravmodel[i, ] <- ((n_o^OriginPower) *
                           (n_d^destpower)) /
            (OffsetA+(d_ij/OffsetB)^kernpower)
        ## Truscott and Ferguson version
        ## gravmodel[i, ] <- ((n_o^OriginPower) * (n_d^DestPower)) /
        ## (OffsetA+d_ij/OffsetB)^Power

        ## Zero off within cell mixing
        ind <- (dscell == orcell)
        gravmodel[i, ind] = 0
    }

    ## individual-level normalisation
    for (i in 1:noorig) {
        gravmodel[i,] <- gravmodel[i,] / sum(gravmodel[i,])
    }

  return( gravmodel )

}

gravy.calc.lnlike.nowithinregion.harriet <- function( dataMatrix, model, noorig , originindices, destindices) {
    ## Inputs:
  #   dataMatrix (obs.tab)   matrix of observed events in each dest cell, row=origin, col=destination cells.
  # 	model - matrix of probabilities, with same dimensions as obs.tab.
  # 	noorig -  number of origin cells.
  # Outputs:
  # 	lnlike  log-likelihood value

  lnlike = 0
  dscell   <- destindices$cell
  for (i in 1:noorig) {

    orcell   <- originindices$cell[i]
    #check when origin == destination cell
    ind = dscell != orcell # gives indicator when different

    if (min(model[i, ind]) < 0) {
      #stop("Problem in gravy.calc.lnlike function: min(model[i,])<0" )
    }

    lnlike <- lnlike + dmultinom( dataMatrix[i, ind], prob=model[i, ind], log=TRUE )

  }

  return( lnlike )
}

mob_calc_S_mat <- function(popmatrix,popsize_vector,D,A=1) {
  # Faster function to calculate S_ij, the sum of population densities
  # within radius r_ij,
  # for each origin-destination pair.
  # Inputs:
  # popmatrix is raster
  # D is number of districts
  # A is number of age groups
  # Outputs:
  #  S -- matrix of S_ij, where rows correspond
  # to originindices, and columns to destinations.

  # within region mixing
  r_forfirstcell=raster::values(distanceFromPoints(popmatrix, xyFromCell(popmatrix,1)))

  # distance between first and second
  xlength=r_forfirstcell[2]

  #distance between first and the second row
  ylength=r_forfirstcell[ncol(popmatrix)+1]
  dwr_r=distancewithinarectangle(xlength,ylength)

  S = matrix( 0, nrow=D*A, ncol=D*A)
  i=1
  pb <- txtProgressBar(min=0,max=D*A,style=3)
  while (i <= (D*A)) {
    #print(i)
    ai=i%%A; ai[ai==0]=A
    di=(i-ai)/A +1

    n_orig = popsize_vector[i] # density in cell i

    d.ij=values(raster::distanceFromPoints(popmatrix, xyFromCell(popmatrix,di)))
    d.ij[di]=dwr_r #distance within that district
    popsize=popsize_vector
    destcells = 1:(D*A) #list of destination cells
    neworder=order(d.ij) #orders by the distance away
    d.ij = d.ij[neworder]
    popsize = popsize[neworder]
    destcells = destcells[neworder]

    S.ij = vector(len=D*A)
    rollingsum_ind = 1
    rollingsum = 0
    uni_dist = unique(d.ij)
    for (k in 1:length(uni_dist)){

      # finds cells which are the same distance away
      samedist_ind = which(d.ij==uni_dist[k])
      last_ind=samedist_ind[length(samedist_ind)]
      rollingsum = rollingsum + sum(popsize[rollingsum_ind:last_ind]) # the minimum sum up to this distance
      rollingsum_ind = last_ind + 1
      S.ij[samedist_ind] = rollingsum

    }
    S.ij = S.ij - n_orig - popsize # remove the origin and destination cell
    #NB there will be some negative cells
    S[i, ] = S.ij[order(destcells)] # reorder into correct order

    i=i+1

    setTxtProgressBar(pb,i)

  }

  return( S )

}

gravy.radiation.model.offset.harriet <- function(
  x2, S, noorig, nodest, originindices, destindices,
  distances, all_originindices, offset ) {
  # the offset radiation model allows within region mixing, by redefining the S matrix to be offset

  # offset S matrix
  offset.S=matrix(0, dim(S)[1], dim(S)[2])
  radmodel <- matrix(0, nrow=noorig, ncol=nodest)
  dscell   <- destindices$cell
  n_d   <- x2[dscell] 	# density at destination cell
  for (i in 1:noorig) {

    orcell   <- originindices$cell[i]

    n_o   <- x2[orcell] 	# density at origin cell

    # match i to correct row in S
    S_index = which(all_originindices[,2]==orcell)
    offset.S[S_index, ] = S[S_index, ]

    #redefine S matrix
    index_in_Offset = which(distances[i, ] <= offset)
    regions_in_Offset = destindices$cells[index_in_Offset] # indicates which regions are within the offset
    offset.S[S_index, index_in_Offset]=sum(x2[regions_in_Offset]) - n_o # all S entries less than Offset away from i are set to a maximum

    S_ij = offset.S[S_index, ]

    # None of the dimensions below here match
    radmodel[i, ] <- (n_o * n_d) / ( (n_d + S_ij)*(n_d + n_o + S_ij) )

    # consider only between region, set r to zero within region
    ind = dscell == orcell # gives indicator when same
    radmodel[i, ind] = 0

    if (min(radmodel[i,])<0) {
      # stop("model prediction < 0 in gravy.radiation.model function")
    }
  }
  # normalise probabilities
  for (i in 1:noorig) {
    radmodel[i, ] <- radmodel[i, ] / sum(radmodel[i, ])
  }
  return( radmodel )
}

distancewithinarectangle <- function(side1, side2){
  # the average distance between two points uniformly distributed in the rectange
  # with sides a, b with a>=b
  # http://www.math.uni-muenster.de/reine/u/burgstal/d18.pdf and Santalo, LA Integral
  # geometry and geometric probability pg 49 (book)
  if (side1>=side2){
    a=side1
    b=side2
  } else {
    a=side2
    b=side1
  }
  d=sqrt(a^2 + b^2)
    dist=(1/15)*(a^3/b^2 + b^3/a^2 + d*(3 - a^2/b^2 -
                                        b^2/a^2) + (5/2)*((b^2/a)*log((a+d)/b) +
                                                          (a^2/b)*log((b+d)/a)))
}
