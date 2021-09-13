#ifndef SQMS_PAIR_POTENTIALS_H
#define SQMS_PAIR_POTENTIALS_H
#if defined(__cplusplus)
extern "C" {
#endif

/**
 * @brief Calculates the value of Lennard-Jones potential
 * 
 * Using the notation in https://en.wikipedia.org/wiki/Lennard-Jones_potential
 * 
 * @param seperation_squared 
 * @param length_scale_squared 
 * @param energy_scale 
 * @return double 
 */
double SqMS_LJ_potential(double seperation_squared, double length_scale_squared, double energy_scale);
double SqMS_truncated_LJ_potential_unsafe(double seperation_squared, double length_scale_squared, double energy_scale, double energy_cutoff);
double SqMS_truncated_LJ_potential(double seperation_squared, double length_scale_squared, double energy_scale, double length_cutoff_squared, double energy_cutoff);

double SqMS_finite_well_potential(double seperation_squared, double well_depth, double well_width_squared);

#if defined(__cplusplus)
}
#endif
#endif /* SQMS_PAIR_POTENTIALS_H */