#' Panel A: map of household and contact locations
#'
#' Reproduces panel A of the manuscript's "Study design" figure: a
#' population-density heatmap of the Guangzhou study region with household
#' and contact locations overlaid, Guangzhou/Dongguan/Foshan landmarks, a
#' distance scale bar, a population-density legend, and a China-location
#' inset map.
#'
#' Ported from the legacy \code{make.figure.1()} in \code{fluscape}'s
#' \code{source/R/gravy_functions.R}, rewritten to take only a contacts data
#' frame and a population raster -- no separate \code{locations}/
#' \code{participants} tables or raw-data loaders are needed, since
#' \code{contacts} already carries \code{HH_Long}/\code{HH_Lat} directly.
#' The legacy figure's numbered study-location markers are dropped, since
#' those need a \code{locations} table this package does not ship.
#'
#' @param contacts A contacts data frame (e.g. \code{contacts_fluscape_V1} or
#'   \code{contacts_V1_jittered_100m}) with columns \code{pid},
#'   \code{HH_Long}, \code{HH_Lat}, \code{long}, \code{lat}.
#' @param popgrid A population-density \code{RasterLayer} (e.g. from
#'   \code{\link{load_gz_pop_raster}}).
#' @param plot.landmarks Add Guangzhou/Dongguan/Foshan markers?
#' @param plot.country.map Add the China-location inset map?
#' @param fig The figure region to draw in, as \code{c(x1, x2, y1, y2)} (see
#'   \code{\link[graphics]{par}}'s \code{fig}). Defaults to the whole
#'   device, for standalone use; \code{\link{fig_study_design}} passes a
#'   smaller region when composing all three panels together. Sub-elements
#'   (the density legend, the country inset) are placed as fractions of this
#'   region, which is a deliberate generalisation of the legacy figure's
#'   fixed, absolute-device-fraction placement -- needed to make the panel
#'   composable at an arbitrary position.
#' @param add If \code{TRUE}, draw into the current plot region
#'   (\code{par(new = TRUE)}) instead of starting a fresh plot -- used when
#'   composing panels via \code{\link{fig_study_design}}. A device is opened
#'   automatically if none is currently open (e.g. for a standalone call in
#'   a fresh session), but an already-open device (including one opened via
#'   \code{\link[grDevices]{pdf}} to save to a file) is always respected and
#'   never overridden.
#'
#' @return \code{NULL}, invisibly. Called for its plotting side effect.
#'
#' @examples
#' fig_study_design_panel_map(contacts_V1_jittered_100m, load_gz_pop_raster())
#'
#' @importFrom graphics par box points rect axis mtext legend
#' @importFrom grDevices dev.new dev.cur colorRampPalette rgb
#' @export
fig_study_design_panel_map <- function(contacts, popgrid,
                                        plot.landmarks = TRUE,
                                        plot.country.map = TRUE,
                                        fig = c(0, 1, 0, 1),
                                        add = FALSE) {

  contacts_small <- gravy.subset.contact.data(contacts$pid, contacts, popgrid)

  long.lim <- c(extent(popgrid)@xmin, extent(popgrid)@xmax)
  lat.lim  <- c(extent(popgrid)@ymin, extent(popgrid)@ymax)

  hhxy <- unique(contacts_small[, c("HH_Long", "HH_Lat")])
  hhxy[, 1] <- hhxy[, 1] + stats::runif(nrow(hhxy), -0.001 * diff(long.lim), 0.001 * diff(long.lim))
  hhxy[, 2] <- hhxy[, 2] + stats::runif(nrow(hhxy), -0.001 * diff(lat.lim), 0.001 * diff(lat.lim))
  cxy <- contacts_small[, c("long", "lat")]
  cxy[, 1] <- cxy[, 1] + stats::runif(nrow(cxy), -0.001 * diff(long.lim), 0.001 * diff(long.lim))
  cxy[, 2] <- cxy[, 2] + stats::runif(nrow(cxy), -0.001 * diff(lat.lim), 0.001 * diff(lat.lim))

  if (dev.cur() == 1) dev.new(width = 14, height = 8)
  par(fig = fig, new = add)
  par(mar = c(0, 0, 0, 0), las = 1, lend = 1, cex.axis = .8, cex.lab = 1.2)

  pcol <- colorRampPalette(c("grey92", rgb(.65, .75, .65), rgb(.65, .60, .4)))(101)
  image(log10(popgrid), col = pcol,
        xlab = "", ylab = "", main = "", axes = FALSE,
        xlim = c(113.1, 113.9), ylim = c(22.98, 23.58),
        zlim = c(0.0001, 6))
  box()

  if (plot.landmarks) {
    pt.bg <- "gold"; pt.cex <- 1.5; txt.col <- rgb(0, 0, 0, .6); txt.cex <- 1.2
    points(113.25, 23.133, pch = 21, bg = pt.bg, cex = pt.cex)
    text(113.25, 23.133 + .015, "Guangzhou", pos = 2, cex = txt.cex, col = txt.col)
    points(113.7167, 23.0333, pch = 21, bg = pt.bg, cex = pt.cex)
    text(113.7167, 23.0333 - .01, "Dongguan", pos = 4, cex = txt.cex, col = txt.col)
    points(113.116667, 23.016667, pch = 21, bg = pt.bg, cex = pt.cex)
    text(113.116667, 23.016667 - .01, "Foshan", pos = 4, cex = txt.cex, col = txt.col)
  }

  points(cxy[, 1], cxy[, 2], pch = 16, cex = .7, col = rgb(0, 0, 1, .4))
  points(hhxy[, 1], hhxy[, 2], pch = 22, cex = .6, bg = "salmon", lwd = .5)

  xy1 <- c(113.4, 23.015); xy2 <- c(113.5, 23.015)
  d <- pointDistance(matrix(xy1, ncol = 2), matrix(xy2, ncol = 2), lonlat = TRUE) / 1000
  xy2[1] <- xy1[1] + ((10 / d) * 0.1)
  rect(xy1[1], xy1[2] - .0025, xy2[1], xy2[2] + .0025, col = "grey30", border = "grey30")
  x.seq <- seq(from = xy1[1] + (xy2[1] - xy1[1]) / 2, to = xy2[1], length.out = 6)
  rect(x.seq[1], xy1[2] - .0025, x.seq[2], xy2[2] + .0025, col = "grey70", border = "grey30")
  rect(x.seq[2], xy1[2] - .0025, x.seq[3], xy2[2] + .0025, col = "grey30", border = "grey30")
  rect(x.seq[3], xy1[2] - .0025, x.seq[4], xy2[2] + .0025, col = "grey70", border = "grey30")
  rect(x.seq[4], xy1[2] - .0025, x.seq[5], xy2[2] + .0025, col = "grey30", border = "grey30")
  rect(x.seq[5], xy1[2] - .0025, x.seq[6], xy2[2] + .0025, col = "grey70", border = "grey30")
  text(xy1[1] + ((5 / d) * 0.1), xy1[2], "10 km", pos = 1, col = rgb(0, 0, 0, .65))

  dens <- as.matrix(10^seq(0, 6, length.out = 101), nrow = 1)
  par(fig = c(fig[1], fig[1] + 0.3 * (fig[2] - fig[1]),
              fig[3] + 0.64 * (fig[4] - fig[3]), fig[3] + 0.75 * (fig[4] - fig[3])),
      new = TRUE)
  par(mar = rep(0, 4), xaxs = "i", yaxs = "i")
  plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", axes = FALSE)
  rect(0, 0, 1, 1, col = rgb(1, 1, 1, .25), border = NA)

  par(fig = c(fig[1], fig[1] + 0.3 * (fig[2] - fig[1]),
              fig[3] + 0.68 * (fig[4] - fig[3]), fig[3] + 0.74 * (fig[4] - fig[3])),
      new = TRUE)
  par(mar = c(1, 1, 1, 0.6))
  image(log10(dens), col = pcol, axes = FALSE); box()
  axis(side = 1, at = seq(0, 1, length.out = 4),
       labels = c(expression(10^0), expression(10^2), expression(10^4), expression(10^6)),
       cex.axis = 1)
  mtext(side = 3, expression(paste("population density, ", km^-1)), cex = 1)

  if (plot.country.map) {
    par(fig = c(fig[1], fig[1] + 0.3 * (fig[2] - fig[1]),
                fig[3] + 0.75 * (fig[4] - fig[3]), fig[4]),
        new = TRUE)
    par(mar = rep(0, 4), xaxs = "i", yaxs = "i")
    plot(c(107, 123), c(17, 33), type = "n", xlab = "", ylab = "", axes = FALSE)
    rect(105, 17, 123, 33, col = rgb(1, 1, 1, .85), border = "black")
    maps::map("world", "China", add = TRUE, fill = TRUE, col = rgb(.50, .45, .4, .5))
    txt.col <- rgb(0, 0, 0, .6)
    text(110.5, 29.5, "China", pos = 3, cex = 1.5, col = txt.col)
    points(121.473654, 31.229750, pch = 16, col = "grey30", cex = .8)
    text(121.473654 + .5, 31.229750, "Shanghai", pos = 2, cex = 1, col = txt.col)
    points(113.25, 23.133, pch = 16, col = "grey30", cex = 1)
    text(113.25 - 2.5, 23.133, "Guangzhou", pos = 3, cex = 1, col = txt.col)
    points(114.153291, 22.278629, pch = 16, col = "grey30", cex = .8)
    text(114.153291 + 2.5, 22.278629, "Hong Kong", pos = 1, cex = 1, col = txt.col)
    rect(113.1, 22.98, 113.9, 23.58, border = "red")
    box()
  }

  mtext(side = 2, "A", padj = -2.5, line = 1, cex = 2)
  invisible(NULL)
}

#' Panel B: example individual contact trajectories
#'
#' Reproduces panel B of the manuscript's "Study design" figure: a handful
#' of individuals' household-to-contact lines, drawn over a local
#' population-density crop.
#'
#' Selects individuals meeting the same eligibility criteria as the legacy
#' figure (more than \code{min_contacts} distinct contact locations, median
#' contact distance at least \code{min_median_km}), then samples \code{n} of
#' them at random. This is a general resampling, not the specific,
#' hand-curated set of individuals shown in the published figure.
#'
#' @inheritParams fig_study_design_panel_map
#' @param n Number of individuals to sample and display.
#' @param min_contacts Minimum number of distinct contact events required
#'   for an individual to be eligible for sampling.
#' @param min_median_km Minimum median contact distance (km) required for
#'   an individual to be eligible for sampling.
#' @param seed Random seed for sampling which individuals to display.
#'   \code{NULL} uses the current RNG state.
#'
#' @return \code{NULL}, invisibly.
#'
#' @examples
#' fig_study_design_panel_trajectories(contacts_V1_jittered_100m, load_gz_pop_raster(), seed = 1)
#'
#' @importFrom stats median
#' @export
fig_study_design_panel_trajectories <- function(contacts, popgrid, n = 10,
                                                 min_contacts = 5, min_median_km = 1.8,
                                                 seed = NULL,
                                                 fig = c(0, 1, 0, 1), add = FALSE) {

  contacts_small <- gravy.subset.contact.data(contacts$pid, contacts, popgrid)
  contacts_small$d <- pointDistance(
    contacts_small[, c("HH_Long", "HH_Lat")],
    contacts_small[, c("long", "lat")],
    lonlat = TRUE
  ) / 1000

  n_per_id <- tapply(contacts_small$d, contacts_small$pid, length)
  med_per_id <- tapply(contacts_small$d, contacts_small$pid, median)
  eligible <- intersect(
    names(n_per_id)[n_per_id > min_contacts],
    names(med_per_id)[med_per_id >= min_median_km]
  )
  if (length(eligible) == 0) stop("No individuals meet the trajectory-panel eligibility criteria")
  if (length(eligible) < n) {
    warning(sprintf("Only %d eligible individuals found (requested %d); showing all.",
                     length(eligible), n))
  }
  if (!is.null(seed)) set.seed(seed)
  ids <- sample(eligible, min(n, length(eligible)))

  tmp.c <- contacts_small[contacts_small$pid %in% ids, ]
  margin <- 0.04
  ext <- extent(min(c(tmp.c$HH_Long, tmp.c$long)) - margin,
                max(c(tmp.c$HH_Long, tmp.c$long)) + margin,
                min(c(tmp.c$HH_Lat, tmp.c$lat)) - margin,
                max(c(tmp.c$HH_Lat, tmp.c$lat)) + margin)
  x3 <- crop(popgrid, ext)

  if (dev.cur() == 1) dev.new(width = 14, height = 8)
  par(fig = fig, new = add)
  par(mar = c(0, 2.5, 0, 1), las = 1, lend = 1)

  pcol <- colorRampPalette(c("grey92", rgb(.65, .75, .65), rgb(.65, .65, .4)))(101)
  image(log10(x3), col = pcol, axes = FALSE, xlab = "", ylab = "", main = "")
  box()

  dcols <- rev(colorRampPalette(c("cyan", "blue", "yellow", "red"))(length(ids)))
  darkcols <- rev(colorRampPalette(c("blue", "darkblue", "gold", "darkred"))(length(ids)))
  tcols <- paste(dcols, "66", sep = "")

  for (j in seq_along(ids)) {
    tmp <- contacts_small[contacts_small$pid == ids[j], ]
    for (k in seq_len(nrow(tmp))) {
      lines(c(tmp$HH_Long[k], tmp$long[k]), c(tmp$HH_Lat[k], tmp$lat[k]),
            lwd = 3, col = tcols[j])
    }
    points(tmp[1, c("HH_Long", "HH_Lat")], pch = 22, bg = dcols[j], col = darkcols[j], cex = 1.5)
    points(tmp[, c("long", "lat")], pch = 16, col = dcols[j], cex = 1)
  }
  mtext(side = 2, "B", padj = -4.75, line = .75, cex = 2)
  invisible(NULL)
}

#' Panel C: reverse-CDF of contact distance from home
#'
#' Reproduces panel C of the manuscript's "Study design" figure: the
#' reverse-cumulative distribution of contact distance from home, split by
#' age group (child/adult) and urban/rural, plus the overall distribution.
#'
#' @param contacts A contacts data frame with columns \code{HH_Long},
#'   \code{HH_Lat}, \code{long}, \code{lat}, \code{CHILD} (logical),
#'   \code{URBAN} (0/1).
#' @param fig,add See \code{\link{fig_study_design_panel_map}}.
#'
#' @return \code{NULL}, invisibly.
#'
#' @examples
#' fig_study_design_panel_distcdf(contacts_V1_jittered_100m)
#'
#' @export
fig_study_design_panel_distcdf <- function(contacts, fig = c(0, 1, 0, 1), add = FALSE) {

  d <- pointDistance(
    contacts[, c("HH_Long", "HH_Lat")],
    contacts[, c("long", "lat")],
    lonlat = TRUE
  ) / 1000

  find_cum_prob <- function(x) {
    y <- table(x)
    z <- length(x) - cumsum(y) + 1
    out <- data.frame(x = as.numeric(names(z)), y = as.numeric(z))
    out$y <- out$y / out$y[1]
    out
  }

  p_all   <- find_cum_prob(d)
  p_child <- find_cum_prob(d[contacts$CHILD == TRUE])
  p_adult <- find_cum_prob(d[contacts$CHILD == FALSE])
  p_urban <- find_cum_prob(d[contacts$URBAN == 1])
  p_rural <- find_cum_prob(d[contacts$URBAN == 0])

  if (dev.cur() == 1) dev.new(width = 14, height = 8)
  par(fig = fig, new = add)
  par(mar = c(4.5, 5, 2, 1), las = 1, lend = 1)

  plot(p_all, log = "xy", type = "n",
       xlim = c(0.005, 100), ylim = c(0.0005, 1), axes = FALSE,
       xlab = "Distance from home (km), d", ylab = expression(p(D >= d)))
  box()
  axis(side = 1, at = 10^seq(-3, 2, by = 1),
       labels = c("0.001", "0.01", "0.1", "1", "10", "100"))
  axis(side = 2, at = c(0.001, 0.01, 0.1, 1.0),
       labels = c("0.001", "0.01", "0.1", "1.0"))

  pcols <- c(rgb(.5, .5, .5), rgb(1, 0, 0), rgb(0, 0, 1), rgb(.8, .6, 0), rgb(0, .4, 0))
  pt.pch <- c(16, 2, 6, 3, 4)
  points(p_urban, col = paste0(pcols[4], "99"), pch = pt.pch[4], cex = .5)
  points(p_rural, col = paste0(pcols[5], "99"), pch = pt.pch[5], cex = .5)
  points(p_child, col = paste0(pcols[2], "99"), pch = pt.pch[2], cex = .5)
  points(p_adult, col = paste0(pcols[3], "99"), pch = pt.pch[3], cex = .5)
  points(p_all,   col = paste0(pcols[1], "99"), pch = pt.pch[1], cex = .5)
  legend("bottomleft", bty = "n", cex = 1.3,
         legend = c("all", "children", "adults", "urban", "rural"),
         pt.cex = .75, pch = pt.pch, lty = 0, col = pcols)
  mtext(side = 2, "C", padj = -5.75, line = 2.75, cex = 2)
  invisible(NULL)
}

## contacts_V1_jittered_100m is a lazy-loaded package dataset, used here only
## as a default argument value -- R CMD check's static analysis can't see
## that, so it's declared here to silence the resulting NOTE.
utils::globalVariables("contacts_V1_jittered_100m")

#' Reproduce the manuscript's "Study design" figure (Figure 2)
#'
#' Composes \code{\link{fig_study_design_panel_map}} (A),
#' \code{\link{fig_study_design_panel_trajectories}} (B), and
#' \code{\link{fig_study_design_panel_distcdf}} (C) into the three-panel
#' figure. Defaults to the public jittered dataset so it runs for anyone
#' with just the installed package; pass \code{contacts_fluscape_V1} for the
#' figure as actually published.
#'
#' @inheritParams fig_study_design_panel_map
#' @param n_trajectories Passed to
#'   \code{\link{fig_study_design_panel_trajectories}} as \code{n}.
#' @param seed Passed to \code{\link{fig_study_design_panel_trajectories}}.
#'
#' @return \code{NULL}, invisibly.
#'
#' @examples
#' fig_study_design(contacts_V1_jittered_100m, load_gz_pop_raster(), seed = 1)
#'
#' @export
fig_study_design <- function(contacts = contacts_V1_jittered_100m,
                              popgrid = load_gz_pop_raster(),
                              n_trajectories = 10,
                              seed = NULL) {
  if (dev.cur() == 1) dev.new(width = 14, height = 8)
  fig_study_design_panel_map(contacts, popgrid, fig = c(0.05, 0.65, 0.02, 0.98), add = FALSE)
  fig_study_design_panel_trajectories(contacts, popgrid, n = n_trajectories, seed = seed,
                                       fig = c(0.65, 1, 0.6, 0.98), add = TRUE)
  fig_study_design_panel_distcdf(contacts, fig = c(0.65, 1.0, 0.0, 0.6), add = TRUE)
  invisible(NULL)
}
