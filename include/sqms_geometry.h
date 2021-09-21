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
#ifndef NDIM
#define NDIM 2
#endif

#if defined(__cplusplus)
extern "C" {
#endif


/**
 * @brief type definition of particle
 *
 * @param r position of the particle
 * @param radius the radius of the particle (if applicable)
 * @param isplaceholder the indicator if the particle is real or not
 *
 */
typedef struct
{
    double r[NDIM];
    double radius;
    int isplaceholder;
    int uid;
} particle_t;

/**
 * @brief copies particle data
 * 
 * @param source the data to be copied from
 * @param destitanion the destination of the copied data
 */
void SqMS_copy_particle(const particle_t* source, particle_t* destitanion);

/**
 * @brief type definition of right angled cell list
 *
 * @param cell cells of the cell list that store the particles
 * @param dimensions dimensions of each cell
 * @param limits the limits of the cells that tesselate the simulation box
 * @param list_shape number of cells in each axis
 * @param cell_population the number of particles in each cell
 * @param max_population the maximum number of particles in the cell list
 * @param max_population_per_cell the maximum number of particles in a cell
 *
 */
typedef struct 
{
  int *cell_indices;
  int num_cells
} bounding_box_t;

typedef struct
{
    particle_t* cell;
    double dimensions[NDIM];
    double *limits[NDIM];
    int list_shape[NDIM];
    int *cell_population;
    int *advised_insertion;
    int max_population;
    int max_population_per_cell;
    int num_cell;
} right_cell_list_t;




/**
 * @brief The initializer of a right angled cell list
 *
 * @param cell_list the pointer to the cell list that is going to be initialized
 * @param Box the simulation box to be tesselated by the cells
 * @param cutoff_distance the cutoff distance of the pair potential
 */
void SqMS_init_right_cell_list(right_cell_list_t* cell_list, const double* Box,
                                                const double cutoff_distance,
                                                const double particle_volume);

void SqMS_free_right_cell_list(right_cell_list_t* cell_list);
/**
 * @brief remove particle from the cell
 *
 * @param cell_list
 * @param cell_id
 * @param particle_id
 */
void SqMS_remove_particle_from_right_cell_list(right_cell_list_t* cell_list,
                                    const int cell_id, const int particle_id);
 /**
  * @brief add particle to the cell
  *
  * @param cell_list
  * @param cell_id
  * @param particle
  * @return int particle id in the cell
  */
int SqMS_add_particle_to_cell(right_cell_list_t* cell_list,
                                const particle_t* particle, const int cell_id);
 /**
  * @brief find the suitable cell for the particle
  *
  * @param cell_list
  * @param r
  * @return int cell id of the particle
  */
int SqMS_find_where_particle_belongs(right_cell_list_t* cell_list, double *r);

/**
 * @brief is the particle supposed to be in that cell
 *
 * @param cell_list
 * @param particle
 * @param cell_id
 * @return int Yes/No
 */
int SqMS_is_particle_in_cell(right_cell_list_t* cell_list,
                                const double *r, const int cell_id);


void SqMS_cell_id_to_position(right_cell_list_t* cell_list, const int cell_id, int cell_position[]);
int SqMS_cell_position_to_id(right_cell_list_t* cell_list, int cell_position[]);
void SqMS_get_rigid_bounding_box(right_cell_list_t *cell_list, bounding_box_t* bbox, const int cell_ind);
void SqMS_free_bounding_box(bounding_box_t* bbox);
//==============================================================================
/**
 * @brief distance squared between two points on a flat torus of N dimensions
 *
 * @param particle1 the coordinats of the particle
 * @param particle2 the coordinates of the particle
 * @param Box the lengths of the dimensions
 * @param num_dimensions number of dimensions
 * @return (double) the distance squared between particles
 */
double SqMS_distance_squared_rectangular_PBC(double *particle1, double *particle2 , double *Box, int num_dimensions);

/**
 * @brief the distance between two points on a flat torus of N dimensions
 *
 * @param particle1 the coordinats of the particle
 * @param particle2 the coordinates of the particle
 * @param Box the lengths of the dimensions
 * @param num_dimensions number of dimensions
 * @return (double) the distance between particles
 */
double SqMS_distance_rectangular_PBC(double *particle1, double *particle2 , double *Box, int num_dimensions);




#if defined(__cplusplus)
}
#endif
#endif /* SQMS_GEOMETRY_H */