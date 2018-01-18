#include<iostream>
#include<cmath>
#include <Rcpp.h>
using namespace std;
// [[Rcpp::export]]

Rcpp::NumericMatrix distance(Rcpp::NumericMatrix start, Rcpp::NumericMatrix end)
{
    //Haversine distance method
    float R = 6378.137;
    float toRad = 3.14159/180;
    Rcpp::NumericMatrix out((start.size()/2), (end.size()/2));
    int zz = 0;
    for(int i = 0; i < (start.size()/2); ++i) {
        for(int n = 0; n < (end.size()/2); ++n){
            
            float lon1 = start(i,0);
            float lon2 = end(n,0);
            float lat1 = start(i,1);
            float lat2 = end(n,1);
            
            lon1 = lon1 * toRad;
            lon2 = lon2 * toRad;
            lat1 = lat1 * toRad;
            lat2 = lat2 * toRad;
            float dlon = lon2 - lon1;
            float dlat = lat2 - lat1;
            
            double a = pow(sin(dlat / 2), 2) + (cos(lat1) * cos(lat2) * pow(sin(dlon / 2),2));
            double d = 2 * atan2(sqrt(a), sqrt(1 - a)) * R;
            out(i,n) = d;
            ++zz;
        }
    }
    
    return out ;
    

}






