/**
 * @file sqms_geometry.h
 * @author Alptug Ulugol (a.ulugol@uu.nl)
 * @brief 
 * @version 0.1
 * @date 2021-09-07
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef SQMS_GEOMETRY_H
#define SQMS_GEOMETRY_H
#if defined(__cplusplus)
extern "C" {
#endif

double SqMS_distance_squared_rectangular_PBC(double *particle1, double *particle2 , double *Box, int num_dimensions);
//distance on rectangular box w/periodic boundaries
double SqMS_distance_rectangular_PBC(double *particle1, double *particle2 , double *Box, int num_dimensions);

#if defined(__cplusplus)
}
#endif
#endif /* SQMS_GEOMETRY_H */