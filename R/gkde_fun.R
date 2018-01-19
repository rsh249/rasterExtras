#' @import raster
#' @import parallel
#' @import stats
NULL

#' Geographic Kernel Density Estimator using linear or Haversine distances
#' 
#' This function requests data from a private service that
#' mirrors GBIF distribution data at cloud.diversityoflife.org. Account required.
#' 
#' @param grid A raster object to match.
#' @param points A two column data frame in the form (lon,lat) or (x,y)
#' @param parallel TRUE or FALSE, should this code be executed in parallel.
#' @param nclus IF parallel==TRUE then how many cores in the cluster.
#' @param dist.method Which distance should we use? Haversine for lat/long projections, 
#'   or Pythagorean for flat images and/or small areas.
#'   
#' @export
#' @examples
#' require(raster)
#' grid = raster::raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90, vals=NULL)
#' grid = raster::setValues(grid,values=(as.vector(seq(1:raster::ncell(grid)))))
#' points = cbind(seq(xmin(grid), xmax(grid), length.out=100), 
#'       seq(ymin(grid), ymax(grid), length.out=100))
#' plot(grid); points(points);
#' di = gkde(grid, points, parallel=FALSE, dist.method='Pythagorean');
#' plot(di)


gkde <- function(grid, points, parallel=TRUE, nclus = 4, dist.method = 'Haversine'){
	zz = seq(1:raster::ncell(grid)); ##Pass to gkde
	if(dist.method == "Haversine"){
  	pbp = distance(as.matrix(points), as.matrix(points));
  	bw.gen = stats::bw.nrd(as.vector(pbp));
  	x = seq(1:raster::ncell(grid))
  	n=25000;
  	f <- sort(rep(1:(trunc(length(x)/n)+1),n))[1:length(x)]
  	splits = split(x,f);
  	
  	##Hidden gkde core function
  	.gkde.core.h <- function(x){
  	  coords = latlonfromcell(as.vector(x), as.vector(raster::extent(grid)), nrow(grid), ncol(grid)); 
  	  d = distance(as.matrix(coords[,2:1]),as.matrix(points));
  	  di=vector();
  	  for(c in 1:nrow(coords)){
  	    di[c] = stats::density(d[c,], n = 1, kernel = 'gaussian', from = 0,  to = 0,  bw = bw.gen, na.rm = TRUE)$y;
  	  }
  	  return(di);
  	}

  	if(parallel == FALSE){
  	  #di = .gkde.core.h(zz);
  	  di = unlist(lapply(splits, .gkde.core.h));
  	} else {
  	  cl = parallel::makeCluster(nclus, type ='FORK');
  	  di = unlist(parallel::parLapply(cl, splits, .gkde.core.h));
  	  parallel::stopCluster(cl);
  	}
  	
	} else if(dist.method == "Pythagorean"){
	  pbp = pythagorean(as.matrix(points), as.matrix(points));
	  bw.gen = stats::bw.nrd(as.vector(pbp));
	  x = seq(1:raster::ncell(grid))
	  n=25000;
	  f <- sort(rep(1:(trunc(length(x)/n)+1),n))[1:length(x)]
	  splits = split(x,f);
	  
	  .gkde.core.p <- function(x){ 
	    coords = latlonfromcell(as.vector(x), as.vector(raster::extent(grid)), nrow(grid), ncol(grid)); 
	    d = pythagorean(as.matrix(coords[,2:1]),as.matrix(points));
	    di=vector();
	    for(c in 1:nrow(coords)){
	      di[c] = stats::density(d[c,], n = 1, kernel = 'gaussian', from = 0,  to = 0,  bw = bw.gen, na.rm = TRUE)$y;
	    }
	    return(di);
	  }
	  if(parallel == FALSE){
	    di = unlist(lapply(splits, .gkde.core.p));
	  } else {
	    cl = parallel::makeCluster(nclus, type ='FORK');
	    di = unlist(parallel::parLapply(cl, splits, .gkde.core.p));
	    parallel::stopCluster(cl);
	  }
	} else {
	  cat("ERR1: Bad Distance method declaration\n");
	  return(NULL);
	}

	r=raster::raster(nrows=nrow(grid),ncol=ncol(grid), crs="+proj=longlat +datum=WGS84",ext=raster::extent(grid))
	r=raster::setValues(r,values=di)
	r = raster::mask(r, grid)
	return(r);
}








