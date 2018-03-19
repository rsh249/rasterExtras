#' @import raster
#' @import parallel
#' @import stats
NULL

#' Geographic Kernel Density Estimator using linear or Haversine distances
#' 
#' This function calculates a kernel density estimation for raster objects.
#' 
#' @param grid A raster object to match.
#' @param points A two column data frame in the form (lon,lat) or (x,y)
#' @param parallel TRUE or FALSE, should this code be executed in parallel.
#' @param nclus IF parallel==TRUE then how many cores in the cluster.
#' @param dist.method Which distance should we use? Haversine for lat/long projections,or Pythagorean for flat images and/or small areas.
#' @param maxram Maximum theoretical RAM usage. Will be divided by nclus for parallel jobs.
#' @param bw Bandwidth. Either 'nrd0' to find a default bandwidth (advised when points data are small), or a numeric value giving a suitable bandwidth (predetermined by you).
#'   
#' @export
#' @examples
#' require(raster)
#' grid = raster::raster(nrows=18, ncols=36, xmn=-180, xmx=180, ymn=-90, ymx=90, vals=NULL)
#' grid = raster::setValues(grid,values=(as.vector(seq(1:raster::ncell(grid)))))
#' points = cbind(seq(xmin(grid), xmax(grid), length.out=100), 
#'       seq(ymin(grid), ymax(grid), length.out=100))
#' plot(grid); points(points);
#' di = gkde(grid, points, parallel=TRUE, dist.method='Pythagorean');
#' plot(di)


gkde <-
  function(grid,
           points,
           parallel = TRUE,
           nclus = 4,
           dist.method = 'Haversine',
           maxram=4, bw = 'nrd0') {
    
    .gkde.core.h <- function(x) {
      require(rasterExtras)
      coords = latlonfromcell(as.vector(x),
                              as.vector(
                                c(
                                  raster::xmin(grid),
                                  raster::xmax(grid),
                                  raster::ymin(grid),
                                  raster::ymax(grid)
                                )
                              ),
                              nrow(grid),
                              ncol(grid))
      
      d = distance(as.matrix(coords[, 2:1]), as.matrix(points))
      
      di = vector()
      
      for (c in 1:raster::nrow(coords)) {
        di[c] = stats::density(
          d[c, ],
          n = 1,
          kernel = 'gaussian',
          from = 0,
          to = 0,
          bw = bw.gen,
          na.rm = TRUE
        )$y
        
      }
      return(di)
      
    }
    .gkde.core.p <- function(x) {
      coords = latlonfromcell(as.vector(x),
                              as.vector(
                                c(
                                  raster::xmin(grid),
                                  raster::xmax(grid),
                                  raster::ymin(grid),
                                  raster::ymax(grid)
                                )
                              ),
                              nrow(grid),
                              ncol(grid))
      
      d = pythagorean(as.matrix(coords[, 2:1]), as.matrix(points))
      
      di = vector()
      
      for (c in 1:raster::nrow(coords)) {
        di[c] = stats::density(
          d[c, ],
          n = 1,
          kernel = 'gaussian',
          from = 0,
          to = 0,
          bw = bw.gen,
          na.rm = TRUE
        )$y
        
      }
      return(di)
      
    }
    
    
    
    
    xx = seq(1:raster::ncell(grid))
    
    #one cell of a matrix should contain a double, so:
    dd = 32; ##Double double because Rcpp for some reason doubles ram on return.
    
    
    #16 bytes per cell. Estimates seem off using this value for memory of doubles in matrices
    vol = dd * (nrow(points) ^ 2) / 1024 / 1024 / 1024
    
    
    if(bw=='nrd0'){
      if (vol > maxram) {
        #if distance matrix will be > than ???
        ##Bootstrap bandwidth selection
        n = 1000
        
        bw = vector()
        
        for (i in 1:n) {
          sam = sample(c(1:nrow(points)), 100, replace = TRUE)
          p = points[sam, ]
          
          if (dist.method == "Pythagorean") {
            ps = as.vector(pythagorean(as.matrix(points[sam, ]), as.matrix(points[sam, ])))
            
          } else if (dist.method == "Haversine") {
            ps = as.vector(distance(as.matrix(points[sam, ]), as.matrix(points[sam, ])))
            
          }
          bw[i] = stats::bw.nrd(as.vector(ps))
          
        }
        bw.gen = stats::median(bw)
        
      } else {
        pbp = as.vector(pythagorean(as.matrix(points), as.matrix(points)))
        
        if (dist.method == "Pythagorean") {
          pbp = as.vector(pythagorean(as.matrix(points), as.matrix(points)))
          
        } else if (dist.method == "Haversine") {
          pbp = as.vector(distance(as.matrix(points), as.matrix(points)))
          pbp = stats::na.omit(pbp)
          
        }
        bw.gen = stats::bw.nrd(as.vector(pbp))
        
      }
    } else {
      bw.gen = bw;
    }
    
    ##Check grid x points matrix size.
    vol = (dd * (nrow(points) * raster::ncell(grid))) / 1024 / 1024 / 1024
    ramtarg= maxram/nclus;
    targ.n = ceiling((ramtarg / vol) * raster::ncell(grid))
    
    if (targ.n > raster::ncell(grid)) {
      splits = list(seq(1:raster::ncell(grid)))
      
    } else {
      n = targ.n
      f <- sort(rep(1:(trunc(length(
        xx
      ) / n) + 1), n))[1:length(xx)]
      splits = split(xx, f)
    }
    
    
    if (parallel == FALSE) {
      if (dist.method == "Pythagorean") {
        di = unlist(lapply(splits, .gkde.core.p))
        
      } else if (dist.method == "Haversine") {
        di = unlist(lapply(splits, .gkde.core.h))
        
      }
    } else {
      if (nclus > length(splits)) {
        n = length(xx) / nclus
        
        f <- sort(rep(1:(trunc(length(
          xx
        ) / n) + 1), n))[1:length(xx)]
        splits = split(xx, f)
      }
      
      ##Reporting bloc
      cat("BEGIN PARALLEL COMPUTATION\n");
      cat("Core count: ", nclus, "\n");
      cat("Cells/iteration: ", length(splits[[1]]), "of", raster::ncell(grid), "\n")
      cat("Points: ", nrow(points), "\n");
      cat("Maximum RAM per proc.: ", ramtarg/nclus, "\n");
      cat("Distance Method: ", dist.method, "\n\n");
      ###
      p = proc.time();
      cl = parallel::makeCluster(nclus, type = 'SOCK')
      
      parallel::clusterExport(cl,
                              c(
                                "grid",
                                "points",
                                "bw.gen",
                                "latlonfromcell",
                                "pythagorean"
                              ),
                              envir = environment())
      
      if (dist.method == "Pythagorean") {
        di = unlist(parallel::parLapply(cl, splits, .gkde.core.p))
        
      } else if (dist.method == "Haversine") {
        di = unlist(parallel::parLapply(cl, splits, .gkde.core.h))
        
      }
      parallel::stopCluster(cl)
      ep = proc.time() - p;
      cat("Time elapsed: ", ep[[3]], "\n\n");
      
    }
    r = raster::raster(
      nrows = nrow(grid),
      ncol = ncol(grid),
      crs = "+proj=longlat +datum=WGS84",
      ext = raster::extent(grid)
    )
    r = raster::setValues(r, values = di)
    
    r = raster::mask(r, grid)
    
    return(r)
    
  }









