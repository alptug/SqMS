#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "../include/sqms_geometry.h"

double SqMS_distance_squared_rectangular_PBC(double *particle1, double *particle2 , double *Box, int num_dimensions)
{
    //dimension agnostic pbc distance
    double tmp = 0; 
    double distance_squared = 0;
    for (size_t i = 0; i < num_dimensions; i++)
    {
        tmp = fabs(particle1[i] - particle2[i]);

        if(tmp>0.5*Box[i])
        {
            tmp = Box[i]-tmp;
        }

        distance_squared += tmp*tmp;
    }
    
    return distance_squared;
}

double SqMS_distance_rectangular_PBC(double *particle1, double *particle2 , double *Box, int num_dimensions)
{
    return sqrt(SqMS_distance_squared_rectangular_PBC(particle1, particle2 , Box, num_dimensions));
}