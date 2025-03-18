## This file is for generating sample S matrix for radiation model by aggregating the original raster x2

generate.agg.Smatrix <- function(contacts,
                                 raster_object,
                                 agg_num
                                 ){
  # remove NA values of longitude and latitude
  contacts_small <- contacts[!is.na(contacts$long) | !is.na(contacts$lat) | !is.na(contacts$HH_Long) | !is.na(contacts$HH_Lat),]

  # aggregate large raster object by defining the number of aggregation times
  if (agg_num == 0){
    x <- raster_object
  }else{
    x <- aggregate(raster_object, agg_num)
  }

  # faster non-zero destination indices calculation:
  # get all values of aggregated raster object
  tmp <- getValues(x)
  # remove NA values
  cells_dest <- which(!is.na(tmp) & tmp>0 )
  # create a new data frame to store destination indexes
  destindices <- data.frame(index=1:length(cells_dest), cells=cells_dest)

  # faster origin indices calculation:
  cells_ori <- unique(cellFromXY(x, contacts_small[,c("HH_Long","HH_Lat")] ) )
  originindices <- data.frame(index=1:length(cells_ori), cells=cells_ori)

  # faster distances calculation:
  a <- xyFromCell(x, destindices$cell )
  b <- xyFromCell(x, originindices$cell )
  distances = pointDistance(b, a, longlat=T )

  # calculate number of origin and possible destination cells
  noorig <- dim(originindices)[1] # number of origin cells
  nodest <- dim(destindices)[1] # number of possible destination cells, will be very large if use original raster object!

  # make a matrix of observed contact origin-destinations
  contactindices <- cellFromXY(x, contacts[,c("long","lat")] ) # in which cells did contact occur?
  hhindices <- cellFromXY(x, contacts[,c("HH_Long","HH_Lat")] ) # which cells have participant households?
  obs.tab <- matrix(0, nrow=noorig, ncol=nodest  ) # matrix of observations
  pb <- txtProgressBar(min = 0, max = noorig, initial = 0)
  for (k in 1:noorig) {
    j <- originindices$cell[k]
    i <- contactindices[hhindices==j]
    i <- i[!is.na(i)]
    z <- table(i)
    i.indices <- match(as.numeric(names(z)),destindices$cells)
    obs.tab[k,i.indices] <- as.numeric(z)
    setTxtProgressBar(pb,k)
  }
  close(pb)

  # Find radiation density sums, s_ij, for each orig-dest cell combination.
  # Finds the total population density within equivalent distance from the origin.
  # Warning, this is slow and takes a long time.
  # Save to disk after generation for future use.
  S <- matrix(nrow=noorig, ncol=nodest )
  for (i in 1:noorig) {
    print(paste(i,"/",noorig))
    dist <- distances[i, ]
    pb <- txtProgressBar(min = 0, max = nodest, style = 3)
    for (j in 1:nodest) {
      within_d <- which(dist<=dist[j])
      radcell <- destindices$cell[within_d] # all destination cells within radius dist
      S[i,j] <- sum(raster::extract(x,radcell) )
      setTxtProgressBar(pb, j )
      close(pb)
    }
  }
   
  agg_objects <- list(S = S,
                      popgrid = x)
  return (agg_objects)
}
