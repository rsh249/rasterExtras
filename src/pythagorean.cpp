#include<iostream>
#include<cmath>
#include <Rcpp.h>
using namespace std;
// [[Rcpp::export]]

Rcpp::NumericMatrix pythagorean(Rcpp::NumericMatrix start, Rcpp::NumericMatrix end)
{
    //Linear distance method
    Rcpp::NumericMatrix out((start.size()/2), (end.size()/2));
    int zz = 0;
    for(int i = 0; i < (start.size()/2); ++i) {
        for(int n = 0; n < (end.size()/2); ++n){
            
            float x1 = start(i,0);
            float x2 = end(n,0);
            float y1 = start(i,1);
            float y2 = end(n,1);
            float xdim = x1 - x2;
            float ydim = y1 - y2;
            double d = sqrt(pow(xdim,2) + pow(ydim, 2));
            
            out(i,n) = d;
            ++zz;
        }
    }
    
    return out ;
    

}






