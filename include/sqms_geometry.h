#ifndef SQMS_GEOMETRY_H
#define SQMS_GEOMETRY_H
#if defined(__cplusplus)
extern "C" {
#endif

//distance on rectangular box w/periodic boundaries
double SqMS_distance_squared_rectangular_PBC(double *particle1, double *particle2 , double *Box, int num_dimensions);
//distance on rectangular box w/periodic boundaries
double SqMS_distance_rectangular_PBC(double *particle1, double *particle2 , double *Box, int num_dimensions);

#if defined(__cplusplus)
}
#endif
#endif /* SQMS_GEOMETRY_H */