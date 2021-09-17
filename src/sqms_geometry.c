#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "../include/sqms_geometry.h"
#include "../include/sqms_misc.h"

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

    cell_list->max_population_per_cell = (int)(4*cell_volume/particle_volume) + 1;
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
        cell_list->cell[i].uid = -1;
    }

}

void SqMS_remove_particle_from_right_cell_list(right_cell_list_t* cell_list,
                                    const int cell_id, const int particle_id)
{
    int index = cell_id * cell_list->max_population_per_cell + particle_id;
    int lastelement = cell_id * cell_list->max_population_per_cell + cell_list->cell_population[cell_id]-1;

    if (cell_list->cell[lastelement].isplaceholder == 0)
    {
        SqMS_copy_particle(&cell_list->cell[lastelement], &cell_list->cell[index]);
        cell_list->cell[lastelement].isplaceholder = 1;
        cell_list->advised_insertion[cell_id] = lastelement;
    }
    else
    {
        cell_list->cell[index].isplaceholder = 1;
        cell_list->advised_insertion[cell_id] = particle_id;
    }
    
    cell_list->cell_population[cell_id]--;
}

int SqMS_add_particle_to_cell(right_cell_list_t* cell_list,
                                const particle_t* particle, const int cell_id)
{
  
    assert(cell_list->cell_population[cell_id] < cell_list->max_population_per_cell);
    assert(cell_id < cell_list->num_cell);
    int particle_id = cell_id*cell_list->max_population_per_cell+cell_list->cell_population[cell_id];
    int particle_id_advised = cell_id*cell_list->max_population_per_cell+cell_list->advised_insertion[cell_id];
    int last_stand=0, last_found=0;
    if (cell_list->cell[particle_id].isplaceholder == 1)
    {
        
    }
    else if (particle_id_advised < particle_id && cell_list->cell[particle_id_advised].isplaceholder == 1)
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

    SqMS_copy_particle(particle, &cell_list->cell[particle_id]);
    assert(cell_list->cell[particle_id].isplaceholder == 0);
    
    cell_list->cell_population[cell_id]++;
    return particle_id;
}

int SqMS_find_where_particle_belongs(right_cell_list_t* cell_list, double *r)
{
    int cell_position[NDIM];
    for (size_t i = 0; i < NDIM; i++)
    {
        cell_position[i] = SqMS_biggest_lower_bound(cell_list->limits[i],cell_list->list_shape[i], r[i]);
    }
    int cell_id = SqMS_cell_position_to_id(cell_list, cell_position);
    
    return cell_id;
}

int SqMS_is_particle_in_cell(right_cell_list_t* cell_list,
                                const double *r, const int cell_id)
{
    int cell_position[NDIM];
    int tmp = cell_id;
    for (int i = NDIM-1; i>=0; i--)
    {
        cell_position[i] = tmp - (tmp/cell_list->list_shape[i]) * cell_list->list_shape[i];
        tmp /= cell_list->list_shape[i];
    }

    int in_cell = 1;
    for (size_t i = 0; i < NDIM; i++)
    {
        if (r[i] < cell_list->limits[i][cell_position[i]] || r[i] > cell_list->limits[i][cell_position[i]+1] )
        {
            in_cell = 0;
            break;
        }
    }
    return in_cell;
}


void SqMS_cell_id_to_position(right_cell_list_t* cell_list, const int cell_id, int cell_position[])
{
    int tmp = cell_id;
    for (int i = NDIM-1; i>=0; i--)
    {
        cell_position[i] = tmp - (tmp/cell_list->list_shape[i]) * cell_list->list_shape[i];
        tmp /= cell_list->list_shape[i];
    }
}

int SqMS_cell_position_to_id(right_cell_list_t* cell_list, int cell_position[])
{
    int cell_id = 0;
    for (size_t i = 0; i < NDIM; i++)
    { 
        cell_id *= cell_list->list_shape[i];
        cell_id += cell_position[i];
    }
    return cell_id;
}

void SqMS_copy_particle(const particle_t* source, particle_t* destination)
{
    destination->isplaceholder = source->isplaceholder;
    destination->radius = source->radius;
    destination->uid = source->uid;

    for (size_t i = 0; i < NDIM; i++)
    {
        destination->r[i] = source->r[i];
        //TODO: you can use memcpy there, it is faster but impact is too small
    }
    
}

void __SqMS_get_rigid_bounding_box_nestedloops(right_cell_list_t *cell_list, bounding_box_t* bbox,unsigned int N, int shifter[], int cell_position[], int *index)
{
    if (N>0) 
    {
        for (int i = -1; i < 2; i++)
        {
            shifter[N-1] = i; 
            __SqMS_get_rigid_bounding_box_nestedloops(cell_list, bbox, N-1, shifter,cell_position, index);
        }
    }
    else
    {
        int neighbour_cell[NDIM];
        for (size_t j = 0; j < NDIM; j++)
        {
            neighbour_cell[j] = shifter[j] + cell_position[j];
                if (neighbour_cell[j]<0)
                {
                    neighbour_cell[j] += cell_list->list_shape[j];
                }
                else if (neighbour_cell[j] >=cell_list->list_shape[j])
                {
                    neighbour_cell[j] -= cell_list->list_shape[j];
                }
        }

        bbox->cell_indices[*index] = SqMS_cell_position_to_id(cell_list,neighbour_cell);
        *index += 1; 
    }

}

void SqMS_get_rigid_bounding_box(right_cell_list_t *cell_list, bounding_box_t* bbox, const int cell_ind)
{
    bbox->num_cells = powi(3,NDIM);
    bbox->cell_indices = calloc(bbox->num_cells, sizeof(int));
    int cell_position[NDIM];
    int cell_position_shifter[NDIM];

    SqMS_cell_id_to_position(cell_list,cell_ind,cell_position);
    int index = 0;
    __SqMS_get_rigid_bounding_box_nestedloops(cell_list, bbox, NDIM,  cell_position_shifter, cell_position, &index);
    assert(index == bbox->num_cells);    
}

void SqMS_free_bounding_box(bounding_box_t* bbox)
{
    free(bbox->cell_indices);
}