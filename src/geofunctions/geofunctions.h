#ifndef GEOFUNCTIONS_H
#define GEOFUNCTIONS_H

#include <assert.h>
#include <math.h>

namespace growing_balls
{
static const double DEG_TO_RAD = 0.017453292519943295769236907684886;
static const double EARTH_RADIUS_IN_CENTIMETERS = 637279756.0856;


/**
 * Compute the greater circle distance on the unit sphere from p1 to p2 given by
 * their lat and lon coordinate.
 *
 * Computing the distance with the haversine formula
 * Assume lat and lon are given in degrees:
 * - lat in range [-90, 90]
 * - lon is in range [-180, 180]
 */
double
haversine ( double p1_lat, double p1_lon, double p2_lat, double p2_lon )
{
    double lat_arc = ( p1_lat - p2_lat ) * DEG_TO_RAD;
    double lon_arc = ( p1_lon - p2_lon ) * DEG_TO_RAD;
    double lat_h = sin ( lat_arc * 0.5 );
    lat_h *= lat_h;
    double lon_h = sin ( lon_arc * 0.5 );
    lon_h *= lon_h;
    double tmp = cos ( p1_lat * DEG_TO_RAD ) * cos ( p2_lat * DEG_TO_RAD );
    return 2.0 * asin ( sqrt ( lat_h + tmp * lon_h ) );
}

/**
 * Compute the distance in centimeters for a spherical earth model between p1
 * and p2 given by their lat and lon coordinate.
 *
 * Computing the distance with the haversine formula
 * Assume lat and lon are given in degrees:
 * - lat in range [-90, 90]
 * - lon is in range [-180, 180]
 */
double
distance_in_centimeters ( double p1_lat, double p1_lon, double p2_lat, double p2_lon )
{
    double tmp = EARTH_RADIUS_IN_CENTIMETERS * haversine ( p1_lat, p1_lon, p2_lat, p2_lon );
    assert ( tmp >= 0 );
    return tmp;
}
}

#endif // GEOFUNCTIONS_H
