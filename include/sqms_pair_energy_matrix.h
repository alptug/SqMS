/**
 * @file sqms_pair_energy_matrix.h
 * @author Alptug Ulugol (a.ulugol@uu.nl)
 * @brief this file defines the functions to be used on an energy matrix
 * @version 0.1
 * @date 2021-09-07
 * 
 * @copyright Copyright (c) 2021
 * 
 * We define the energy matrix as a lower triangular matrix when i != j, 
 * a_{ij} is the energy due to the interaction between particle i and j.
 * When i==j, a_{ii} is the total energy of the particle pairs containing i.
 */

#ifndef SQMS_PAIR_ENERGY_MATRIX_H
#define SQMS_PAIR_ENERGY_MATRIX_H
#if defined(__cplusplus)
extern "C" {
#endif

/**
 * @brief this function initializes the energy matrix with zeros
 * 
 * @param N number of particles, also the dimension of the matrix (N*N)
 * @param energy_matrix a 1D array of size N(N+1)/2.
 */
void SqMS_init_energy_matrix(int N, double *energy_matrix);

/**
 * @brief Get the energy contribution of the particle pair (i,j)
 * if i==j returns total energy contribution of particle i
 * 
 * @param i index of a particle
 * @param j index of a particle
 * @param energy_matrix self explanatory
 * @return (double) energy contribution of the particle pair (i,j)
 */
double SqMS_get_energy_from_pair(int i, int j, double *energy_matrix);

/**
 * @brief Get the total energy of the system
 * 
 * @param N number of particles, also the dimension of the matrix (N*N)
 * @param energy_matrix a 1D array of size N(N+1)/2.
 * @return (double) total energy of the system
 */
double SqMS_get_total_energy(int N, double *energy_matrix);

/**
 * @brief Set the energy  contribution of the particle pair (i,j)
 * 
 * @param energy the energy value to be set
 * @param i index of a particle
 * @param j index of a particle
 * @param energy_matrix self explanatory
 */
void SqMS_set_energy_of_pair(double energy, int i, int j, double *energy_matrix);

/**
 * @brief Calculates 'i'th diagonal entry sets it in the matrix
 * 
 * The 'i'th diagonal entry is the total energy contribution of particle i
 * 
 * @param i index of a particle
 * @param N number of particles, also the dimension of the matrix (N*N)
 * @param energy_matrix self explanatory
 * @return (double) the total energy contribution of particle i
 */
double SqMS_calculate_total_energy_contribution_of_particle(int i, int N,
                                                        double *energy_matrix);

/**
 * @brief Calculates individual energy contributions of the particles
 * 
 * This function calculates individual energy contributions of the particles,
 * sets and checks the diagonal entry of the energy matrix
 * 
 * @param N number of particles, also the dimension of the matrix (N*N)
 * @param energy_matrix self explanatory
 * @return (double) total energy of the system 
 */
double SqMS_calculate_total_energy_contributions(int N, double *energy_matrix);

/**
 * @brief calculates the index of the (i,j) pair in the energy matrix
 * 
 * @param i index of a particle
 * @param j index of a particle
 * @return (int) index of the (i,j) pair in the energy matrix
 */
int SqMS_calculate_index_from_pair(int i, int j);


void SqMS_dump_all_interactions(size_t i, size_t N, double *energy_matrix);

#if defined(__cplusplus)
}
#endif
#endif /* SQMS_PAIR_ENERGY_MATRIX_H */