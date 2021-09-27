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

void SqMS_init_right_cell_list(cell_list_t* cell_list, const double* Box,
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

void SqMS_remove_particle_from_cell_list(cell_list_t* cell_list,
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

int SqMS_add_particle_to_cell(cell_list_t* cell_list,
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

int SqMS_find_where_particle_belongs(cell_list_t* cell_list, double *r)
{
    int cell_position[NDIM];
    for (size_t i = 0; i < NDIM; i++)
    {
        cell_position[i] = SqMS_biggest_lower_bound(cell_list->limits[i],cell_list->list_shape[i], r[i]);
    }
    int cell_id = SqMS_cell_position_to_id(cell_list, cell_position);
    
    return cell_id;
}

int SqMS_is_particle_in_cell(cell_list_t* cell_list,
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


void SqMS_cell_id_to_position(cell_list_t* cell_list, const int cell_id, int cell_position[])
{
    int tmp = cell_id;
    for (int i = NDIM-1; i>=0; i--)
    {
        cell_position[i] = tmp - (tmp/cell_list->list_shape[i]) * cell_list->list_shape[i];
        tmp /= cell_list->list_shape[i];
    }
}

int SqMS_cell_position_to_id(cell_list_t* cell_list, int cell_position[])
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

void __SqMS_get_rigid_bounding_box_nestedloops(cell_list_t *cell_list, bounding_box_t* bbox,unsigned int N, int shifter[], int cell_position[], int *index)
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

void SqMS_get_rigid_bounding_box(cell_list_t *cell_list, bounding_box_t* bbox, const int cell_ind)
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


//==============================================================================
// IT IS TIME TO DO FLOPPY STUFF
//==============================================================================

double SqMS_distance_squared_floppy_PBC(double *particle1, double *particle2, 
                                        double *floppyBox, int num_dimensions)
{
    //dimension agnostic pbc distance
    double delta_r[num_dimensions];
    double distance_squared = 0, tmp = 0;
    for (size_t i = 0; i < num_dimensions; i++)
    {
        tmp = particle1[i] - particle2[i];
        tmp -= floor(tmp + 0.5);
        delta_r[i] = tmp;
    }
    for (size_t i = 0; i < num_dimensions; i++)
    {
        tmp = 0;
        for (size_t j = 0; j < num_dimensions; j++)
        {
            tmp += floppyBox[i*num_dimensions + j] * delta_r[j];
        }
        distance_squared += tmp * tmp;
    }
    return distance_squared;
}

double SqMS_distance_floppy_PBC(double *particle1, double *particle2,
                                        double *floppyBox, int num_dimensions)
{
    return sqrt(SqMS_distance_squared_floppy_PBC(particle1, particle2 ,
                                                    floppyBox, num_dimensions));
}

void SqMS_init_floppy_cell_list(cell_list_t* cell_list, const double* floppyBox,
                                                const double cutoff_distance,
                                                const double particle_volume)
{
    //not tested yet

    if (NDIM != 2)
    {
        printf("Only 2D floppy cell lists are supported");
        exit(-1);
    }
    
    double cell_volume = 1;
    double edge;
    cell_list->num_cell = 1;
    for (size_t i = 0; i < NDIM; i++)
    {
        //calculate how many cells per axis
        int size = floppyBox[i*(NDIM + 1)]/cutoff_distance; //initializer assumes we start with a rectangular box;
        size = (size > 3) ? size: 3;
        cell_list->list_shape[i] = size;
        cell_list->num_cell *= size;
        cell_list->limits[i] = calloc(sizeof(double),size+1);

        //calculate cell volume
        edge = 1./size; // fractional coord
        cell_list->dimensions[i] = edge;
        cell_volume *= edge;
    }
    cell_volume *= SqMS_determinant(floppyBox, NDIM);

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
        cell_list->limits[i][cell_list->list_shape[i]] = 1.; // fractional coord
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

void SqMS_get_ellipse_bounding_box(double* bounds, double cutoff_length, 
                                                            double* floppyBox)
{
    bounds[0] = cutoff_length * sqrt(floppyBox[1]*floppyBox[1] 
                + floppyBox[3]*floppyBox[3]) / fabs(floppyBox[0]*floppyBox[3]
                - floppyBox[1]*floppyBox[2]);

    bounds[1] = cutoff_length * sqrt(floppyBox[0]*floppyBox[0] 
                + floppyBox[2]*floppyBox[2]) / fabs(floppyBox[0]*floppyBox[3]
                - floppyBox[1]*floppyBox[2]);
}

void SqMS_get_ellipse_vertical_intersects(double* intersects, double x, 
                                        double cutoff_length, double* floppyBox)
{
    double b = x * x * (floppyBox[0]*floppyBox[1] + floppyBox[2]*floppyBox[3]);
    double a =  floppyBox[1]*floppyBox[1] + floppyBox[3]*floppyBox[3];
    double sqrtdiscriminant = x*(floppyBox[0]*floppyBox[3]
                                - floppyBox[1]*floppyBox[2]);
    sqrtdiscriminant *= -sqrtdiscriminant;
    sqrtdiscriminant += cutoff_length*cutoff_length*a;
    sqrtdiscriminant = sqrt(sqrtdiscriminant);

    intersects[0] = (-b + sqrtdiscriminant)/a;
    intersects[1] = (-b - sqrtdiscriminant)/a;

}

void SqMS_get_ellipse_horizontal_intersects(double* h_intersects, double y, 
                                        double cutoff_length, double* floppyBox)
{
    double b = y * y * (floppyBox[0]*floppyBox[1] + floppyBox[2]*floppyBox[3]);
    double a =  floppyBox[0]*floppyBox[0] + floppyBox[2]*floppyBox[2];
    double sqrtdiscriminant = y*(floppyBox[0]*floppyBox[3]
                                - floppyBox[1]*floppyBox[2]);
    sqrtdiscriminant *= -sqrtdiscriminant;
    sqrtdiscriminant += cutoff_length*cutoff_length*a;
    sqrtdiscriminant = sqrt(sqrtdiscriminant);

    h_intersects[0] = (-b + sqrtdiscriminant)/a;
    h_intersects[1] = (-b - sqrtdiscriminant)/a;

}

int SqMS_is_point_in_ellipse(double x, double y, double cutoff_length, 
                                                            double* floppyBox)
{
    double xx = floppyBox[0] * x + floppyBox[1] * y;
    double yy = floppyBox[2] * x + floppyBox[3] * y;
    return xx*xx + yy*yy < cutoff_length*cutoff_length;
}

void SqMS_get_interaction_shape(interaction_shape_t* bound, 
                cell_list_t *cell_list, double cutoff_length, double* floppyBox)
{
    /*
    let's talk what we are doing here, when the simulation volume is not a
    rectangular shape, we call it a 'floppy' box.

    at this point, our cells in the cell list are just scaled down copies of
    that shape, and working in the fractional coordinates are easier.

    there is a catch! the floppy box induces a metric! Due to that induced
    metric, circles (spheres) in Cartesian coordinates are deformed into
    ellipses (ellipsoids). 

    this has an important consequence! Now, it is NOT enough to check the
    nearest neighbours (NNs) of a given cell for interactions, the deformed 
    interaction circle (sphere), aka ellipse (ellipsoid), can span more or less
    than the NNs.

    In this function, We will learn how to get all the cells that contain a part
    of an ellipse located at the origin. Bear in mind that we (pragmatically) 
    assume the ellipse is just a sphere that is deformed by some linear operator
    */


    /*
    first get the top right corner of the rectangle that is bounding the
    ellipse, we can get the other corners by 2 reflections since the ellipse
    is centered at the origin.
    */
    double ellipse_bounds[2] = {0,0}; 
    SqMS_get_ellipse_bounding_box(ellipse_bounds, cutoff_length, floppyBox);
    size_t width = 2 * (size_t) ceil(ellipse_bounds[0] / cell_list->dimensions[0]);
    size_t height = 2 * (size_t) ceil(ellipse_bounds[0] / cell_list->dimensions[0]);
    bound->shape = (int*) calloc(width*height, sizeof(int));
    bound->height = height;
    bound->width = width;

    int contains = 0;
    for (size_t i = 0; i < height/2; i++)
    {
        for (size_t j = 0; j < width; j++)
        {
            
            double vline = cell_list->dimensions[0] * ((double)j - (double)width/2.);
            double hline = cell_list->dimensions[1] * (height/2 -1 - i);

            if (vline < 0)
            {
                contains = 
                    SqMS_is_point_in_ellipse(vline+cell_list->dimensions[0], 
                                            hline, cutoff_length, floppyBox);
            }
            else
            {
                contains = SqMS_is_point_in_ellipse(vline, hline, cutoff_length,
                                                                     floppyBox);
            }
            
            if (contains == 0)
            {
                double intersects[2];
                SqMS_get_ellipse_vertical_intersects(intersects, vline, 
                                                    cutoff_length, floppyBox);
                if (hline <= intersects[0] 
                            && intersects[0] < hline + cell_list->dimensions[1])
                {
                    contains = 1;
                }
                else if (hline <= intersects[1] 
                            && intersects[1] < hline + cell_list->dimensions[1])
                {
                    contains = 1;
                }
                else
                {
                    SqMS_get_ellipse_vertical_intersects(intersects, 
                    vline + cell_list->dimensions[0], cutoff_length, floppyBox);
                    if (hline <= intersects[0] 
                            && intersects[0] < hline + cell_list->dimensions[1])
                    {
                        contains = 1;
                    }
                    else if (hline <= intersects[1] 
                            && intersects[1] < hline + cell_list->dimensions[1])
                    {
                        contains = 1;
                    }
                    else
                    {
                        SqMS_get_ellipse_horizontal_intersects(intersects,hline, 
                                                    cutoff_length, floppyBox);
                        if (vline <= intersects[0] 
                            && intersects[0] < vline + cell_list->dimensions[0])
                        {
                            contains = 1;
                        }
                        else if (vline <= intersects[1] 
                                && intersects[1] < vline + cell_list->dimensions[0])
                        {
                            contains = 1;
                        }
                    }
                }
                
            }
            

            int index = i*width + j;
            bound->shape[index] = contains;
            bound->shape[width*height- 1 - index] = contains; //because of symmetry
        }
        
    }
}
void SqMS_free_interaction_shape(interaction_shape_t* x)
{
    free(x->shape);
}

void SqMS_interaction_shape_to_bounding_shape(interaction_shape_t* int_shape, int **bounding_shape)
{
    int bounding_shape_size = (int_shape->height+1)*(int_shape->width+1);
    *bounding_shape = (int*) calloc(bounding_shape_size, sizeof(int));

    for (size_t i = 0; i < bounding_shape_size; i++)
    {
        (*bounding_shape)[i] = 0;
    }

    for (size_t i = 0; i < int_shape->height; i++)
    {
        for (size_t j = 0; j < int_shape->width; j++)
        {
            int index = (i)*(int_shape->width+1) + (j);
            (*bounding_shape)[index] = ((*bounding_shape)[index] + int_shape->shape[(i)*(int_shape->width) + (j)])>0;
            index = (i)*(int_shape->width+1) + (j+1);
            (*bounding_shape)[index] = ((*bounding_shape)[index] + int_shape->shape[(i)*(int_shape->width) + (j)])>0;
            index = (i+1)*(int_shape->width+1) + (j);
            (*bounding_shape)[index] = ((*bounding_shape)[index] + int_shape->shape[(i)*(int_shape->width) + (j)])>0;
            index = (i+1)*(int_shape->width+1) + (j+1);
            (*bounding_shape)[index] = ((*bounding_shape)[index] + int_shape->shape[(i)*(int_shape->width) + (j)])>0;
        }
    }
}

int SqMS_bounding_shape_to_bounding_coordinates(int *bounding_shape, int height, int width, int* bounding_coordinates)
{
    int origin_i = height / 2;
    int origin_j = width / 2;
    int count = 0;
    for (size_t i = 0; i < height; i++)
    {
        for (size_t j = 0; j < width; j++)
        {
            if (bounding_shape[i*width+j]==1)
            {
                bounding_coordinates[2*count] = i - origin_i;;
                bounding_coordinates[2*count+1] = j - origin_j;
                count++;
            }
        }
    }
    return count;
}

void SqMS_get_floppy_bounding_shape(cell_list_t *cell_list, bounding_box_t* bbox, 
        const int cell_ind, int* bounding_coordinates, int bounding_cell_count)
{
    bbox->num_cells = bounding_cell_count;
    bbox->cell_indices = calloc(bounding_cell_count, sizeof(int));
    int cell_position[NDIM];
    SqMS_cell_id_to_position(cell_list, cell_ind, cell_position);

    int tmp[NDIM];
    for (size_t i = 0; i < bounding_cell_count; i++)
    {
        for (size_t j = 0; j < NDIM; j++)
        {
            tmp[j]= cell_position[j] + bounding_coordinates[NDIM*i+j];
            while (tmp[j]<0)
            {
                tmp[j] += cell_list->list_shape[j];
            }
            tmp[j] = tmp[j] % cell_list->list_shape[j];
        }
        bbox->cell_indices[i] = SqMS_cell_position_to_id(cell_list,tmp);
    }
    
}
