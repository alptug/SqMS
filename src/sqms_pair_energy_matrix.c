#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "../include/sqms_pair_energy_matrix.h"

void SqMS_init_energy_matrix(int N, double *energy_matrix)
{
    size_t length = (N*(N+1))/2;
    for (size_t i = 0; i < length; i++)
    {
        energy_matrix[i] = 0;
    }
    
}

int SqMS_calculate_index_from_pair(int i, int j)
{
    if (i>j)
    {
        return (i*(i+1))/2 + j;
    }
    else
    {
        return (j*(j+1))/2 + i;
    }
}

double SqMS_get_energy_from_pair(int i, int j, double *energy_matrix)
{
    return energy_matrix[SqMS_calculate_index_from_pair(i,j)];
}

double SqMS_get_total_energy(int N, double *energy_matrix)
{
    double twofold_energy = 0;
    for (size_t i = 0; i < N; i++)
    {
        twofold_energy += energy_matrix[SqMS_calculate_index_from_pair(i,i)];
    }
    
    return twofold_energy/2.;
}

void SqMS_set_energy_of_pair(double energy, int i, int j, double *energy_matrix)
{
    assert(i!=j); //this function is supposed to set the energy of a pair

    //get the old value
    double energy_old = SqMS_get_energy_from_pair(i, j, energy_matrix);
    //calculate the difference 
    double dE = energy - energy_old;
    
    //update the value
    energy_matrix[SqMS_calculate_index_from_pair(i,j)] = energy;
    //update the diagonals
    energy_matrix[SqMS_calculate_index_from_pair(i,i)] += dE;
    energy_matrix[SqMS_calculate_index_from_pair(j,j)] += dE;
    //energy_matrix[SqMS_calculate_index_from_pair(i,i)] -= energy_old;
    //energy_matrix[SqMS_calculate_index_from_pair(i,i)] += energy;
    //energy_matrix[SqMS_calculate_index_from_pair(j,j)] -= energy_old;
    //energy_matrix[SqMS_calculate_index_from_pair(j,j)] += energy;

}

double SqMS_calculate_total_energy_contribution_of_particle(int i, int N,
                                                        double *energy_matrix)
{
    double energy = 0;
    
    for (size_t j = 0; j < N; j++)
    {
        if (j==i) continue;
        energy += SqMS_get_energy_from_pair(i, j, energy_matrix);
    }
    
    if (energy != SqMS_get_energy_from_pair(i, i, energy_matrix))
    {
        printf("Warning: Calculated total energy of particle %d is not the "
        "same as the existing value.\nThe value is updated but you are strongly"
        " encouraged to debug your code!\n",i);
        energy_matrix[SqMS_calculate_index_from_pair(i,i)] = energy;
    }
    
    return energy;
}

double SqMS_calculate_total_energy_contributions(int N, double *energy_matrix)
{
    double twofold_total_energy = 0;
    for (size_t i = 0; i < N; i++)
    {
        twofold_total_energy += 
        SqMS_calculate_total_energy_contribution_of_particle(i,N,energy_matrix);
    }
    return twofold_total_energy/2.;
}

void SqMS_dump_all_interactions(size_t i, size_t N, double *energy_matrix)
{
    printf("Dumping all interactions for particle %ld\n", i);
    printf("----------------------------------------------\n");
    printf("Particle id \t pair energy contribution\n");

    double calculated_total = 0;
    double recorded_total = SqMS_get_energy_from_pair(i,i,energy_matrix);
    double tmp;
    for (size_t j=0; j < N; j++)
    {
        if (j == i) continue;
        tmp = SqMS_get_energy_from_pair(i,j,energy_matrix);
        if (tmp != 0)
        {
            printf("%ld\t\t%lf\n",j,tmp);
            calculated_total += tmp;
            tmp = 0;
        }

    }
    printf("----------------------------------------------\n");
    printf("Calculated total: %lf, Recorded total: %lf, Diff: %lf\n",
            calculated_total,
            recorded_total,
            fabs(calculated_total-recorded_total));
}