#include<iostream>
#include<cmath>
#include <Rcpp.h>
using namespace std;
// [[Rcpp::export]]

Rcpp::NumericVector findcoord(double lon, double lat, double dist, double brng)
{
    float R = 6378.137;
    float pi = 3.14159;
    brng = brng * (pi / 180);
    lat = lat * pi / 180;
    lon = lon * pi / 180;
    
    float lat2 = asin(sin(lat) * cos(dist / R) + cos(lat) * sin(dist / R) * cos(brng));
    float lon2 = lon + atan2(sin(brng) * sin(dist / R) * cos(lat), cos(dist / R) - sin(lat) * sin(lat2));
    lat2 = lat2 / pi * 180;
    lon2 = lon2 / pi * 180;
    //  return lat2;
    
    Rcpp::NumericVector ll(2);
    ll[0] = lon2;
    ll[1] = lat2;
    return ll ;
    
}
