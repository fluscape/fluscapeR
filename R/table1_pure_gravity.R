#' Load the package's default population-density raster
#'
#' \code{gz_pop_raster.rds} ships in \code{inst/extdata/} (not \code{data/}),
#' per fluscapeR's existing convention, so it is not lazy-loaded
#' automatically -- this reads it in. Used as the default \code{popgrid}
#' argument throughout this package's figure/table regeneration functions.
#'
#' @return A \code{RasterLayer}.
#'
#' @examples
#' load_gz_pop_raster()
#'
#' @export
load_gz_pop_raster <- function() {
  readRDS(system.file("extdata", "gz_pop_raster.rds", package = "fluscapeR"))
}

#' Fit the pure-gravity mobility model for one data subset
#'
#' Thin wrapper around \code{fit.mobility.model()} for the "pure
#' gravity" model variant (distance-power kernel only, no offset, no
#' destination power) -- the simplest row of the manuscript's Table 1.
#' Matches the call made in \code{data-raw/read_reproduce.r}'s
#' \code{jobType == 6} branch.
#'
#' @param contacts A contacts data frame (e.g. \code{contacts_fluscape_V1} or
#'   \code{contacts_V1_jittered_100m}).
#' @param popgrid A population-density \code{RasterLayer} (e.g. from
#'   \code{\link{load_gz_pop_raster}}). \code{NA} cells are treated as zero
#'   population, matching the preprocessing done in
#'   \code{data-raw/read_reproduce.r}.
#' @param datasubset One of \code{"ALL"}, \code{"ADULTS"}, \code{"CHILDREN"},
#'   \code{"RURAL"}, \code{"URBAN"}.
#' @param noRepeats Number of random-start repeats passed to
#'   \code{fit.mobility.model()}.
#' @param seed Random seed passed to \code{fit.mobility.model()}.
#'
#' @return A one-row data frame with columns \code{datasubset}, \code{point},
#'   \code{lower}, \code{upper}, \code{lnlike}.
#'
#' @section Future work:
#' This covers only the pure-gravity row of Table 1. The other five model
#' rows (pure radiation, offset radiation, and the two gravity + offset /
#' destination-power combinations) are not yet implemented here.
#'
#' @examples
#' fit_pure_gravity(contacts_V1_jittered_100m, load_gz_pop_raster(), "ALL", noRepeats = 1)
#'
#' @export
fit_pure_gravity <- function(contacts, popgrid, datasubset = "ALL",
                              noRepeats = 2, seed = 928734924) {

  popgrid <- zero_na_cells(popgrid)

  logfile <- tempfile(fileext = ".csv")
  utils::write.csv(make.data.df(nrow = 0), logfile, row.names = FALSE)

  fit <- fit.mobility.model(
    contacts,
    popgrid,
    sd = seed,
    logfile = logfile,
    noRepeats = noRepeats,
    optfun = fit.gravity.poppower.optim.nowithinregion,
    psToFit = c("kernpower"),
    psLB = c(0.1),
    psUB = c(6),
    datasubset = datasubset,
    fdebug = FALSE,
    lognote = "pure_gravity"
  )

  tab <- fit$tab
  data.frame(
    datasubset = datasubset,
    point = tab$kernpower.pe[1],
    lower = tab$kernpower.lb[1],
    upper = tab$kernpower.ub[1],
    lnlike = tab$lnlike[1]
  )
}

#' Reproduce Table 1's pure-gravity row across all data subsets
#'
#' Calls \code{\link{fit_pure_gravity}} for each requested manuscript data
#' subset and returns one tidy table.
#'
#' @inheritParams fit_pure_gravity
#' @param datasubsets Character vector of subsets to fit.
#'
#' @return A data frame, one row per subset, with columns \code{datasubset},
#'   \code{point}, \code{lower}, \code{upper}, \code{lnlike}.
#'
#' @examples
#' table1_pure_gravity(contacts_V1_jittered_100m, load_gz_pop_raster(),
#'                      datasubsets = c("ALL", "ADULTS"), noRepeats = 1)
#'
#' @export
table1_pure_gravity <- function(contacts, popgrid,
                                 datasubsets = c("ALL", "ADULTS", "CHILDREN", "RURAL", "URBAN"),
                                 noRepeats = 2, seed = 928734924) {
  rows <- lapply(datasubsets, function(ds) {
    fit_pure_gravity(contacts, popgrid, datasubset = ds, noRepeats = noRepeats, seed = seed)
  })
  do.call(rbind, rows)
}

#' Replace NA cells in a population raster with zero
#'
#' \code{data-raw/read_reproduce.r} does this before any model fit; the same
#' preprocessing is applied here so callers can pass the package's raw
#' \code{gz_pop_raster} object directly.
#'
#' @param popgrid A \code{RasterLayer}.
#' @return The same raster with \code{NA} cells set to 0.
#' @keywords internal
zero_na_cells <- function(popgrid) {
  zeromask <- is.na(getValues(popgrid))
  if (any(zeromask)) popgrid[zeromask] <- 0
  popgrid
}
