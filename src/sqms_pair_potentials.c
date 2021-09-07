#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "../include/sqms_pair_potentials.h"


double SqMS_LJ_potential(double seperation_squared, double length_scale_squared, double energy_scale)
{
    double temp = (length_scale_squared * length_scale_squared * length_scale_squared) / (seperation_squared * seperation_squared * seperation_squared);
    return 4.0 * energy_scale * temp * (temp - 1.0);
}

double SqMS_truncated_LJ_potential_unsafe(double seperation_squared, double length_scale_squared, double energy_scale, double energy_cutoff)
{
    return SqMS_LJ_potential(seperation_squared, length_scale_squared, energy_scale) - energy_cutoff;
}

double SqMS_truncated_LJ_potential(double seperation_squared, double length_scale_squared, double energy_scale, double length_cutoff_squared, double energy_cutoff)
{
    double energy = 0;

    if (seperation_squared < length_cutoff_squared) 
        energy += SqMS_truncated_LJ_potential_unsafe(seperation_squared, length_scale_squared, energy_scale, energy_cutoff);

    return energy;
}

double SqMS_finite_well_potential(double seperation_squared, double well_depth, double well_width_squared)
{
    double energy = 0;

    if (seperation_squared < well_width_squared)
        energy -= well_depth;

    return energy;
}