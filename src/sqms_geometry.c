#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "../include/sqms_geometry.h"

double SqMS_distance_squared_rectangular_PBC(double *particle1, double *particle2,
                                                 double *Box, int num_dimensions)
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

double SqMS_distance_rectangular_PBC(double *particle1, double *particle2,
                                             double *Box, int num_dimensions)
{
    return sqrt(SqMS_distance_squared_rectangular_PBC(particle1, particle2 , 
                                                        Box, num_dimensions));
}

void SqMS_init_right_cell_list(right_cell_list_t* cell_list, const double* Box, 
                                                const double cutoff_distance,
                                                const double particle_volume)
{
    //not tested yet
    double cell_volume = 1;
    double edge;
    cell_list->num_cell = 1;
    for (size_t i = 0; i < NDIM; i++)
    {
        //calculate how many cells per axis
        int size = Box[i]/cutoff_distance;
        size = (size > 3) ? size: 3;
        cell_list->list_shape[i] = size;
        cell_list->num_cell *= size;
        cell_list->limits[i] = calloc(sizeof(double),size+1);

        //calculate cell volume
        edge = Box[i]/size;
        cell_list->dimensions[i] = edge;
        cell_volume *= edge;
    }

    cell_list->max_population_per_cell = (int)(cell_volume/particle_volume) + 1;
    cell_list->max_population = cell_list->num_cell * cell_list->max_population_per_cell;

    cell_list->cell = calloc(sizeof(particle_t), cell_list->max_population);
    cell_list->cell_population = calloc(sizeof(int), cell_list->num_cell);
    cell_list->advised_insertion = calloc(sizeof(int), cell_list->num_cell);

    for (size_t i = 0; i < NDIM; i++)
    {
        for (size_t j = 0; j < cell_list->list_shape[i]; j++)
        {
            cell_list->limits[i][j] = j*cell_list->dimensions[i];
        }
        cell_list->limits[i][cell_list->list_shape[i]] = Box[i];
    }
    for (size_t i = 0; i < cell_list->num_cell; i++)
    {
        cell_list->cell_population[i] = 0;
        cell_list->advised_insertion[i] = 0;
    }
    for (size_t i = 0; i < cell_list->max_population; i++)
    {
        cell_list->cell[i].isplaceholder = 1;
    }

}

void SqMS_remove_particle_from_right_cell_list(right_cell_list_t* cell_list, 
                                    const int cell_id, const int particle_id)
{
    cell_list->cell[cell_id * cell_list->max_population_per_cell + particle_id].isplaceholder = 1;
    cell_list->advised_insertion[cell_id] = particle_id;
    cell_list->cell_population[cell_id]--;
}

int SqMS_add_particle_to_cell(right_cell_list_t* cell_list, 
                                const particle_t* particle, const int cell_id)
{
    int particle_id = cell_id*cell_list->max_population_per_cell+cell_list->cell_population[cell_id];
    int particle_id_advised = cell_id*cell_list->max_population_per_cell+cell_list->advised_insertion[cell_id];
    int last_stand=0, last_found=0;
    if (cell_list->cell[particle_id].isplaceholder == 1)
    {
        
    }
    else if (cell_list->cell[particle_id_advised].isplaceholder == 1)
    {
        particle_id = particle_id_advised;
    }
    else
    {
        last_stand = 1;
        int cellindex = cell_id*cell_list->max_population_per_cell;
        int index = 0;
        for (size_t i = 0; i < cell_list->max_population_per_cell; i++)
        {
            index = cellindex + i;
            if (cell_list->cell[index].isplaceholder == 1)
            {
                particle_id = index;
                last_found =1;
                break;
            }
            
        }
        
    }
    if (last_stand ==1) assert(last_found==1);

    cell_list->cell[particle_id].isplaceholder = 0;
    cell_list->cell[particle_id].radius = particle->radius;
    for (size_t i = 0; i < NDIM; i++)
    {
        cell_list->cell[particle_id].r[i]=particle->r[i];
    }
    cell_list->cell_population[cell_id]++;
    return particle_id;
}