# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

distance <- function(start, end) {
    .Call('_rasterExtras_distance', PACKAGE = 'rasterExtras', start, end)
}

findcoord <- function(lon, lat, dist, brng) {
    .Call('_rasterExtras_findcoord', PACKAGE = 'rasterExtras', lon, lat, dist, brng)
}

latlonfromcell <- function(cells, extent, nrow, ncol) {
    .Call('_rasterExtras_latlonfromcell', PACKAGE = 'rasterExtras', cells, extent, nrow, ncol)
}

pythagorean <- function(start, end) {
    .Call('_rasterExtras_pythagorean', PACKAGE = 'rasterExtras', start, end)
}

