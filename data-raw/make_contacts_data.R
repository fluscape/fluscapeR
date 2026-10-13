#' ---
#' title: "Generate the actual and jittered contact datasets, and the GZ population raster"
#' ---
#'
#' This is the only script in data-raw/ that touches the private `fluscape`
#' repo. Its job is narrowly scoped: produce the two contact datasets
#' (`contacts_fluscape_V1`, real coordinates; `contacts_V1_jittered_100m`,
#' 100m-jittered) and the population-density raster (`gz_pop_raster`) that
#' the rest of the package (R/, vignettes/) consumes as package data.
#' Everything else that used to live in data-raw/ (S-matrix generation and
#' validation, the Table 1 model-fitting sweep) has moved to the vignette or
#' is exported R/ functions -- see vignettes/reproduce-mobility-ms.Rmd.
#'
#' Run from the fluscapeR package root, with the private fluscape repo
#' checked out as a sibling directory (../fluscape). Requires `devtools`,
#' `raster`, and `usethis`.

fluscape_top_dir <- "../fluscape/"  # override if fluscape is checked out elsewhere

devtools::load_all()
library(raster)
library(readr)  # load_participant_data_long()/load_hh_data_long() call read_csv() unqualified

source(file.path(fluscape_top_dir, "source/R/mob_utility_private.r"))
source(file.path(fluscape_top_dir, "source/R/fluscape_copy_stevensRfunctions.R"))
source(file.path(fluscape_top_dir, "source/R/GeneralUtility.r"))

## ---- Load participants/households/locations (visit 1) --------------------

participants <- load_participant_data_long(topdir = fluscape_top_dir)
participants <- participants[participants$VISIT == 1, ]
participants$pid <- paste(participants$LOC_ID, participants$HH_ID, participants$PARTICIPANT_ID, sep = "_")

households <- load_hh_data_long(topdir = fluscape_top_dir)
locations <- unique(households[, c("LOC_ID", "LOC_Lat", "LOC_Long", "URBAN", "DIST_FRM_GZ")])

x1 <- fsc.load.wide.raster(fluscapetopdir = fluscape_top_dir)

## ---- Generate the actual (real coordinates) and jittered contact data ----
##
## jitter_error = 0.001 degrees corresponds to ~100m at this latitude:
## 1 degree longitude ~= 111.32*cos(23 deg) ~= 102km, so 0.001 deg ~= 100m
## -- hence "_100m" in the dataset name.

contacts_fluscape_V1 <- mob_generate_contacts_data_from_scratch(
  locations, households, participants,
  jitterxy = FALSE,
  topdir = fluscape_top_dir
)

## The jitter seed lives in the private fluscape repo (source/R/gravity/jitter_seed.txt),
## not here: the jitter step (fluscape's jitter.xy(), called via
## mob_generate_contacts_data_from_scratch()) is a rejection-sampling loop over
## runif() draws, so fixing the seed makes it exactly reproducible -- but if
## that seed were public, anyone with this script, the published raster, and
## the published jittered coordinates could replay the RNG stream and recover
## every participant's real location exactly. fluscape is private, so it's
## safe to commit there; it must never be committed to this (public) repo.
set.seed(as.integer(readLines(file.path(fluscape_top_dir, "source/R/gravity/jitter_seed.txt"), n = 1)))

contacts_V1_jittered_100m <- mob_generate_contacts_data_from_scratch(
  locations, households, participants,
  jitterxy = TRUE,
  jitter_error = 0.001,
  x1 = x1,
  topdir = fluscape_top_dir
)

## ---- Crop the LandScan raster to the study region -------------------------
##
## This box (confirmed this session to match both the package's current
## inst/extdata/gz_pop_raster.rds and what the Table 1 fitting pipeline
## actually uses) is intentionally wider than the contact-filtering box used
## inside mob_generate_contacts_data_from_scratch() above (112.8-114.2 /
## 22.6-24.0): it needs to cover the full destination grid, not just the
## region contacts were filtered to.

gz_pop_raster <- crop(x1, extent(112.2667, 114.8000, 22.11667, 24.50833))

## ---- Save what ships with the package -------------------------------------

usethis::use_data(contacts_V1_jittered_100m, overwrite = TRUE)
saveRDS(gz_pop_raster, file = "inst/extdata/gz_pop_raster.rds")

## ---- contacts_fluscape_V1: deliberately NOT written to data/ here ---------
##
## data/contacts_fluscape_V1.rda already ships in this package (real,
## non-jittered household coordinates) -- a known issue, not introduced or
## fixed by this script (see AI_USAGE_LOG.md in the mobility-gz-ms repo for
## the session history behind it). This script does not perpetuate that: it
## deliberately does not call usethis::use_data(contacts_fluscape_V1),
## consistent with the convention that only public-safe objects get written
## to data/. The object is available in this R session after running this
## script for local, non-public analysis; if it needs to persist, save it
## outside the package (e.g. a gitignored local cache), never to data/.
