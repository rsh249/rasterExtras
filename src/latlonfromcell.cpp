#include<iostream>
#include<cmath>
#include <Rcpp.h>
using namespace std;

// [[Rcpp::export]]
Rcpp::NumericMatrix latlonfromcell(Rcpp::NumericVector cells, Rcpp::NumericVector extent, int nrow, int ncol)
{
    float uplon = extent[0];
    float uplat = extent[3];
    float xres = (extent[1] - extent[0])/ncol;
    float yres = (extent[3] - extent[2])/nrow;
    Rcpp::NumericMatrix m(cells.size(), 2);
    for(int i = 0; i < cells.size(); ++i) {
        int row = ceil(cells[i] / ncol);
        int col = cells[i] - ((row - 1) * ncol);
        float lat = uplat - (row * yres) + (0.5 * yres);
        float lon = uplon + (col * xres) - (0.5 * xres);
        m(i,0) = lat;
        m(i,1) = lon;
        
    }
    
    
    return m ;
}
