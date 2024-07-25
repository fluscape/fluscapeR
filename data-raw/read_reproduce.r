#' The next items for this script are:
#' - run the full gravity model overnight and check that it completes without error
#' - run with at least two parameters oer night and check for completion
#' - if it runs OK - set many jobs going with the actual data
#' - once they are finished, set many jobs going with the jittered data
#' - think about uploading

#' The objective of this script is to load the first visit contact data from
#' t#he primary loading functions and to: recreate the old data set that was
#' used for the analysis initially, generate an updated set of analyses that
#' reflects the correction of some of the data between the start of this
#' proiject and now (Feb 2017), and to make a jittered version of that correct
#' data that can be given out with the manuscript as public data.

#' Set wd if needed
## setwd("~/Dropbox/shares/me_jr_hm_dc_gravity/elife_repro/act_data")
## setwd("~/Dropbox/shares/me_jr_hm_dc_gravity/elife_repro/act_data")
setwd("C:/Users/haowe/Desktop/fluscapeR/")

#' Clear memory
rm(list=ls(all=TRUE))

#' Install the required packages if needed
#' install.packages("R.utils")
library(optimx)
library(devtools)

## Messing around here
## install_github("jameshay218/lazymcmc")

## library("R.utils")

# install("~/dbox/git/lazymcmc",force=TRUE)
# devtools::install_github("jameshay218/lazymcmc")
library(lazymcmc)

#' Optional lines for installing developing packages. For batch mode,
#' the system needs to have the fluscapeR package installed already.
# if (FALSE) {
#     install_github("fluscape/fluscapeR",force=TRUE)
#     ## install("~/Dropbox/git/fluscapeR",force=TRUE)
#     detach("package:fluscapeR",unload=TRUE)
# }

load_all()
library(fluscapeR)

#' Temporary file debugging
fdb <- FALSE

#' Process the command line arguments for the log file
# args <- commandArgs(trailingOnly = TRUE)
jobType_opt <- c(1:6)
datasubset_opt <- c("ALL", "RURAL", "URBAN", "ADULTS", "CHILDREN")
CONTOPT_vec <- c("ACTUAL", "JITTERED")
args <- c("./gravity_log_debug.csv", jobType[2], datasubset[1], 20, CONTOPT[1])

noargs <- length(args)
if (noargs == 5) {
    fnLog <- trimws(args[1])
    jobType <- args[2]
    datasubset <- trimws(args[3])
    noreps <- as.numeric(args[4])
    conts_opt <- trimws(args[5])
    current_wd <- "./"
} else if (noargs == 0) {
    fnLog <- "./gravity_log_debug.csv"
    jobType <- 1
    datasubset="ALL"
    noreps <- 2
    conts_opt <- "JITTERED"
    current_wd <- "C:/Users/haowe/Desktop/fluscapeR/"
} else {
    stop("In script mode, must be either 0 or 3 arguments")
}
if (!file.exists(fnLog)) {
    dftmp <- make.data.df(nrow=0)
    write.csv(dftmp,fnLog,row.names=FALSE)
}

#' Change NAs for 0s in the population density and
#' load the actual contacts and the large population density
#' These lines shouldn't be in the main script
gz_pop_raster <- readRDS("./data/gz_pop_raster.rds")
zeromask <- is.na(as.matrix(gz_pop_raster[,,1]))
gz_pop_raster[zeromask] <- 0
load("./data/pop_S_mat_fluscape.rda") # S matrix
load("./data/contacts_fluscape_V1.rda") # actual contact 
load("./data/contacts_V1_jittered_100m.rda") # jittered contact

#' Select contacts option
if (conts_opt == "JITTERED") {
    contacts_used <- contacts_V1_jittered_100m
} else if (conts_opt == "ACTUAL") {
    cat("Actual contacts NOT available in public version\n")
    contacts_used <- contacts_fluscape_V1
} else {
    stop(paste("Invalid contacts option selected : ",conts_opt,conts_opt=="JITTERED"))
}

#' Some legacy code for S calculation that needs checking
#' This needs checking and perhaps putting into a separate script
#' that can be provided to help people make the S matrix
## gz_pop_raster_agg <- aggregate(gz_pop_raster,10)
## margin_coarse_vec <- values(gz_pop_raster_agg)
## pop_S_mat_fluscape_agg <- mob_calc_S_mat(
##   gz_pop_raster_agg,
##   margin_coarse_vec,
##   D=ncell(margin_coarse_vec),
##   A=1
## )

## Start the different job types
## This loop should become the top level function that has only few
## arguments to be changes routinely. One of them would be jobType
if (jobType == 0) {
    ## Get a specific likelihood value to check the univariate CIs
    fit.mobility.model(
        contacts_fluscape_V1,
        gz_pop_raster,
        Smat = pop_S_mat_fluscape,
        logfile=fnLog,
        noRepeats = noreps,
        optfun = fit.gravity.poppower.optim.nowithinregion,
        psToFit = c("destpower", "kernpower","offset"),
        ## pJustLike = c(0.5141599, 2.661739,  -1.68),
        pJustLike = NULL,
        psLB = c(0.1,0.1,-3),
        psUB = c(4,6,2),
        datasubset=datasubset,
        fdebug=fdb,
        lognote=conts_opt
    )
	## Get a specific likelihood value to check the univariate CIs
	for (val in seq(from=-3,to=-0.3,by=0.1)) {
		fit.mobility.model(
				contacts_used,
				gz_pop_raster,
				Smat = pop_S_mat_fluscape,
				logfile=fnLog,
				noRepeats = noreps,
				optfun = fit.gravity.poppower.optim.nowithinregion,
				psToFit = c("destpower", "kernpower","offset"),
				pJustLike = c(0.52592, 2.730,  val),
				## pJustLike = NULL,
				psLB = c(0.1,0.1,-3),
				psUB = c(4,6,2),
				datasubset=datasubset,
				fdebug=fdb,
				lognote=conts_opt
		)
	}
} else if (jobType == 1) {
    fit.mobility.model(
        contacts_used,
        gz_pop_raster,
        Smat = pop_S_mat_fluscape,
        logfile=fnLog,
        optfun = NULL,
        psToFit = NULL,
        psLB = NULL,
        psUB = NULL,
        datasubset=datasubset,
        fdebug=fdb,
        lognote=conts_opt
    )
} else if (jobType == 2) {
    fit.mobility.model(
        contacts_used,
        gz_pop_raster,
        Smat = pop_S_mat_fluscape,
        logfile = fnLog,
        optfun = fit.offset.radiation.optim,
        noRepeats = noreps,
        psToFit = c("offset"),
        psLB = c(1),
        psUB = c(20*1000),
        datasubset = datasubset,
        fdebug=fdb,
        lognote = conts_opt
    )
} else if (jobType == 3) {
    fit.mobility.model(
        contacts_used,
        gz_pop_raster,
        Smat = pop_S_mat_fluscape,
        logfile=fnLog,
        noRepeats = noreps,
        optfun = fit.gravity.poppower.optim.nowithinregion,
        psToFit = c("destpower", "kernpower","offset"),
        psLB = c(0.1,0.1,-3),
        psUB = c(4,6,2),
        datasubset=datasubset,
        fdebug=fdb,
        lognote= conts_opt
    )
} else if (jobType == 4) {
    fit.mobility.model(
        contacts_used,
        gz_pop_raster,
        Smat = pop_S_mat_fluscape,
        logfile=fnLog,
        noRepeats = noreps,
        optfun = fit.gravity.poppower.optim.nowithinregion,
        psToFit = c("kernpower","offset"),
        psLB = c(0.1,-3),
        psUB = c(6,2),
        datasubset=datasubset,
        fdebug=fdb,
        lognote= conts_opt
    )
} else if (jobType == 5) {
    fit.mobility.model(
        contacts_used,
        gz_pop_raster,
        Smat = pop_S_mat_fluscape,
        logfile=fnLog,
        noRepeats = noreps,
        optfun = fit.gravity.poppower.optim.nowithinregion,
        psToFit = c("kernpower","destpower"),
        psLB = c(0.1,0.1),
        psUB = c(6,4),
        datasubset=datasubset,
        fdebug=fdb,
        lognote= conts_opt
    )
} else if (jobType == 6) {
    fit.mobility.model(
        contacts_used,
        gz_pop_raster,
        Smat = pop_S_mat_fluscape,
        logfile=fnLog,
        noRepeats = noreps,
        optfun = fit.gravity.poppower.optim.nowithinregion,
        psToFit = c("kernpower"),
        psLB = c(0.1),
        psUB = c(6),
        datasubset=datasubset,
        fdebug=fdb,
        lognote=conts_opt
    )
} else {
    stop("job type not known")
}

## A few debug code scetions below here that might be of use if you are
## playing with the package, but are not yet refined or documented.

## A prototype for shared memory parallelization of the jobs
if (FALSE) {
    f <- function(d,s){
        rtn <- paste(d,s)
        list(d=d,s=s,ds=rtn)
    }
    wrapf <- function(vec) {
        f(vec[1],vec[2])
    }
    dfjobs <- data.frame(vecd=1:40,vecs=c("a","b","c","d"))
    wrapf(dfjobs[1,])
    testOutput1 <- apply(dfjobs,1,wrapf)
    cl <- makeCluster(4)
    clusterExport(cl,c("f"))
    testOutput2 <- parApply(cl,dfjobs,1,wrapf)
    stopCluster(cl)
    all(as.character(testOutput1)==as.character(testOutput2))
}

## Run a profile of likelihood
if (FALSE) {
    vecx <- c(500,200,1000,2000,5000,10000,20000,50000,60000,70000,80000,90000)
    vecy <- vector(mode="numeric",length=length(vecx))
    for (i in 1:length(vecx)){
        tmp <- fit.mobility.model(
            contacts_used,
            gz_pop_raster,
            Smat = pop_S_mat_fluscape,
            logfile=fnLog,
            optfun = fit.offset.radiation.optim,
            noRepeats = 5,
            psToFit = c("offset"),
            psLB = c(1),
            psUB = c(100*1000),
            datasubset=datasubset,
            lognote="XXXX",
            justLike = TRUE,
            pJustLike = c(vecx[i])
        )
        vecy[i] <- tmp$like
    }
    pdf()
    plot(vecx,vecy,log="x")
    dev.off()
}
