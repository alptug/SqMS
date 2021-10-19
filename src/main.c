//
//  main.c
//  NVT MCMC base file
//
//  Created by A. Ulugol on 02/09/2021.
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "../include/mt19937.h"
#include "../include/sqms_geometry.h"
#include "../include/sqms_pair_potentials.h"
#include "../include/sqms_pair_energy_matrix.h"
#include "../include/sqms_misc.h"
//#include <mt19937.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------  */
/* Initialization variables  */
/* ------------------------  */

#define NDIM 2 //number of dimensions
#define N 35
#define A_coef 2
#define B_coef 1
#define Box_coef 35
#define NMAX  ((A_coef + B_coef)*N)//max number of particles
#define ITER_MAX 100000000 //max number of iterations
#define DONT_WRITE 0

int n_particlesA = A_coef*N; //number of particles
int n_particlesB = B_coef*N; //number of particles
double pressure = 0.01; // NVT if negative
double well_depth = 5;
int n_particles; //number of particles
enum detection{YES, NO};

const int mc_steps = 500000; //it says steps but they are in fact sweeps
const int NAdjust = mc_steps / 50;
const int output_steps = 1000; //output a file every this many sweeps
const int SANITY_CHECK = 1;
const int SANITY_CHECK_PERIOD = 2000;

//const double diameter = 1.0; //particle diameter
double delta = 0.004; //MCMC move delta
double Vdelta = 0.01; //MCMC move delta
double beta = 1; //inverse temperature

enum detection collision_detection = YES; //detect hard cores
enum particle_type{A,B};
cell_list_t cell_list;

const double eps = 1e-1;
const double tolerance = 1e-6;

/* ---------------------  */
/*  Simulation variables  */
/* ---------------------  */

double radiusA = 1.426, radiusB = 1;
double particle_volumeA,particle_volumeB,max_length_cutoff;

double total_energy;
double *energy_matrix;
double floppyBox[NDIM*NDIM];
double volume;


int n_bounding_cells;
int *bounding_coordinates, *bounding_shape;
interaction_shape_t int_shape;

/* ---------------------------------  */
/*  FILE READER & WRITER SUBROUTINES  */
/* ---------------------------------  */

void read_parameters()
{
    FILE* fp = fopen("input.dat", "r"); //open  initial file
    if (fp == NULL) { //check if the  file exists
        printf("Cannot find parameters in directory, reverting to default values.\n");
    }
    else
    {
        //    first line is # of particles A
        fscanf(fp,"%d\n",&n_particlesA);
        printf("N_A = %d\n", n_particlesA);

        fscanf(fp,"%d\n",&n_particlesB);
        printf("N_B = %d\n", n_particlesB);

        double r_radii = 0;
        fscanf(fp,"%lf\n",&r_radii);
        if(r_radii > 1)
        {
            radiusA = 1;
            radiusB = 1/r_radii;
        }
        else
        {
            radiusA = r_radii;
            radiusB = 1;
        }

        printf("r_A = %lf, r_B = %lf\n", radiusA, radiusB);

        fscanf(fp,"%lf\n",&pressure);
        printf("P = %lf\n", pressure);

        fscanf(fp,"%lf\n",&well_depth);
        printf("E = %lf\n", well_depth);
        fclose(fp);
    }
}

/*void read_data(void){
    FILE* fp = fopen(init_filename, "r"); //open  initial file
    if (fp == NULL) { //check if the  file exists
        printf("Error: Cannot  open '%s'\n",init_filename);
        assert(fp != NULL); //if not raise an error and terminate!
    }
//    first line is # of particles
    fscanf(fp,"%d",&n_particles);
//    next 2/3 lines will give the box dimensions
    double tmp;
    int n,d;
    for(d = 0; d < NDIM; ++d){
        fscanf(fp, "%lf %lf\n",&tmp,&Box[d]);
    }
    //    the rest is n_particles *  NDIM matrix
    for(n = 0; n < n_particles; ++n){
        if (NDIM==3) {
            fscanf(fp, "%lf %lf %lf\n",&r[n][0],&r[n][1],&r[n][2]);
        } else {
            fscanf(fp, "%lf %lf",&r[n][0],&r[n][1]);
        }
    }
    fclose(fp);
}*/

void write_data_cell(int step){
    if (DONT_WRITE == 1)
    {
        return;
    }
    
    char buffer[64];
    sprintf(buffer, "data/coords_step%07d.dat", step);
    FILE* fp = fopen(buffer, "w");
    assert(fp != NULL);
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,Box_coef*max_length_cutoff);
    }
    if (NDIM == 2) fprintf(fp, "%lf %lf\n",0.0,0.0);

    for (size_t i = 0; i < cell_list.max_population; i++)
    {
        /* code */
        if (cell_list.cell[i].isplaceholder == 0)
        {
            for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", (floppyBox[d*NDIM] * 
                                cell_list.cell[i].r[0] + floppyBox[d*NDIM + 1] * 
                                                        cell_list.cell[i].r[1]));
            if (NDIM == 2) fprintf(fp, "%f\t", 0.);
            fprintf(fp, "%lf\n", 2.*cell_list.cell[i].radius);
            
        }
        
    }
    
    fclose(fp);
}

void write_data_cell_native(int step){
    if (DONT_WRITE == 1)
    {
        return;
    }
    
    char buffer[64];
    sprintf(buffer, "data_gl/coords_step%07d.dat", step);
    FILE* fp = fopen(buffer, "w");
    assert(fp != NULL);
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,Box_coef*max_length_cutoff);
    }
    if (NDIM == 2) fprintf(fp, "%lf %lf\n",0.0,0.0);

    for (size_t i = 0; i < cell_list.max_population; i++)
    {
        /* code */
        if (cell_list.cell[i].isplaceholder == 0)
        {
            for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", (floppyBox[d*NDIM] * 
                                cell_list.cell[i].r[0] + floppyBox[d*NDIM + 1] * 
                                                        cell_list.cell[i].r[1]));
            if (NDIM == 2) fprintf(fp, "%f\t", 0.);
            fprintf(fp, "%lf\n", 2.*cell_list.cell[i].radius);
            
        }
        
    }
    
    fclose(fp);
}

void write_data_cell_raw(int step){
    if (DONT_WRITE == 1)
    {
        return;
    }
    
    char buffer[64];
    sprintf(buffer, "data_raw/coords_step%07d.dat", step);
    FILE* fp = fopen(buffer, "w");
    assert(fp != NULL);
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,Box_coef*max_length_cutoff);
    }
    if (NDIM == 2) fprintf(fp, "%lf %lf\n",0.0,0.0);

    for (size_t i = 0; i < cell_list.max_population; i++)
    {
        /* code */
        if (cell_list.cell[i].isplaceholder == 0)
        {
            for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", (floppyBox[d*NDIM] * 
                                cell_list.cell[i].r[0] + floppyBox[d*NDIM + 1] * 
                                                        cell_list.cell[i].r[1]));
            if (NDIM == 2) fprintf(fp, "%f\t", 0.);
            fprintf(fp, "%lf\n", 2.*cell_list.cell[i].radius);
            
        }
        
    }
    
    fclose(fp);
}
/* ___________________________  */
/*  Euclidean metric with PBC   */
/* ---------------------------  */



double dist2(double *particle1, double *particle2)
{
    //dimension agnostic pbc distance
    return SqMS_distance_squared_floppy_PBC(particle1, particle2, 
                                                        floppyBox, NDIM);
}

double dist(double *particle1, double *particle2)
{
    //dimension agnostic pbc distance
    return sqrt(dist2(particle1, particle2));
}

/* ___________________________  */
/*       energy functions       */
/* ---------------------------  */

double pair_potential(double distance_square, double total_radius)
{
    return SqMS_finite_well_potential(distance_square, well_depth, 1.1*1.1*total_radius*total_radius);
}

void init_particle_energies_cell(void)
{
    SqMS_init_energy_matrix(n_particles,energy_matrix);
    total_energy = 0;
    double p_energy = 0;
    double p_dist2 = 0;

    for(size_t i = 0; i<cell_list.max_population; i++)
    {
        if (cell_list.cell[i].isplaceholder == 1)
        {
            continue;
        }
        
        for(size_t j = 0; j<i; j++)
        {
            if (cell_list.cell[j].isplaceholder == 1)
            {
                continue;
            }
            p_dist2 = dist2(cell_list.cell[i].r,cell_list.cell[j].r);
            p_energy = pair_potential(p_dist2, cell_list.cell[i].radius+cell_list.cell[j].radius);
            assert(p_energy <= 0);
            SqMS_set_energy_of_pair(p_energy,cell_list.cell[i].uid,cell_list.cell[j].uid,energy_matrix);
            total_energy += p_energy;
            
        }
    }
}

double particle_energy_cell_list(bounding_box_t *bbox, particle_t *particle)
{
    double p_energy = 0;
    double t_energy = 0;
    double p_dist2 = 0;

    for (size_t ii = 0; ii < bbox->num_cells; ii++)
    {
        for (size_t jj = 0; jj < cell_list.cell_population[bbox->cell_indices[ii]]; jj++)
        {
            int i = jj + cell_list.max_population_per_cell * bbox->cell_indices[ii];

            if (cell_list.cell[i].isplaceholder == 1 || cell_list.cell[i].uid == particle->uid) continue;

            p_dist2 = dist2(particle->r,cell_list.cell[i].r);
            p_energy = pair_potential(p_dist2, cell_list.cell[i].radius+particle->radius);
            assert(p_energy <= 0);
            SqMS_set_energy_of_pair(p_energy,cell_list.cell[i].uid,particle->uid,energy_matrix);
            t_energy += p_energy;
            
        }
        
    }

    double check_energy = SqMS_get_energy_from_pair(particle->uid,particle->uid,energy_matrix);
    double diff = fabs(t_energy - check_energy);
    if (diff > eps)
    {
        printf("fp error?\n");
        assert(diff < eps);
        
    }
    
    return t_energy;
}

double sanity_particle_energy(particle_t *particle, bounding_box_t *bbox)
{
    double p_energy = 0;
    double t_energy = 0;
    double p_dist2 = 0;
    for(size_t i = 0; i<cell_list.max_population; i++)
    {
        if (cell_list.cell[i].isplaceholder == 1) continue;
        if (particle->uid == cell_list.cell[i].uid) continue;

        p_dist2 = dist2(cell_list.cell[i].r,particle->r);
        if (p_dist2 < max_length_cutoff*max_length_cutoff)
        {
            p_energy = pair_potential(p_dist2, 
                    cell_list.cell[i].radius+particle->radius);

            if(collision_detection == YES && sqrt(p_dist2) < cell_list.cell[i].radius+particle->radius)
            {
                printf("Bug in collision detection\n");
                assert(sqrt(p_dist2) > cell_list.cell[i].radius+particle->radius);
            }

            if (fabs(p_energy) > tolerance)
            {
                int cell_ind = i/cell_list.max_population_per_cell;
                int in_box = 0;
                for (size_t j = 0; j < bbox->num_cells; j++)
                {
                    assert(bbox->cell_indices[j] > -1 && bbox->cell_indices[j] < cell_list.num_cell);
                    if (cell_ind == bbox->cell_indices[j])
                    {
                        in_box = 1;
                        
                        break;
                    }
                    
                }
                if(in_box==0) printf("BUG: Interaction outside of Interaction bound!\n");
                assert(in_box == 1);


            }
            t_energy += p_energy;
        }
    }
    
    return t_energy;
}

double dump_all_energies_related(particle_t *particle, bounding_box_t *bbox)
{
    double p_energy = 0;
    double t_energy = 0;
    double p_dist2 = 0;
    for(size_t i = 0; i<cell_list.max_population; i++)
    {
        if (cell_list.cell[i].isplaceholder == 1) continue;
        if (particle->uid == cell_list.cell[i].uid) continue;

        p_dist2 = dist2(cell_list.cell[i].r,particle->r);
        if (p_dist2 < max_length_cutoff*max_length_cutoff)
        {
            p_energy = pair_potential(p_dist2, 
                    cell_list.cell[i].radius+particle->radius);

            if(collision_detection == YES && sqrt(p_dist2) < cell_list.cell[i].radius+particle->radius)
            {
                printf("Bug in collision detection\n");
                assert(sqrt(p_dist2) > cell_list.cell[i].radius+particle->radius);
            }

            if (fabs(p_energy) > tolerance)
            {
                int cell_ind = i/cell_list.max_population_per_cell;
                int in_box = 0;
                printf("ID: %d, distance: %lf, cell: %d, calc cell: %d, 1.1*radii sum: %lf\n",
                        cell_list.cell[i].uid,
                        sqrt(p_dist2),
                        cell_ind,
                        SqMS_find_where_particle_belongs(&cell_list,cell_list.cell[i].r),
                        1.1*(cell_list.cell[i].radius +particle->radius)
                        );
                for (size_t j = 0; j < bbox->num_cells; j++)
                {
                    if (cell_ind == bbox->cell_indices[j])
                    {
                        in_box = 1;
                        break;
                    }
                    
                }
                if(in_box==0) printf("BUG: Interaction outside of Interaction bound!\n");
                assert(in_box == 1);

                bounding_box_t pbbox;
        
                SqMS_get_floppy_bounding_shape(&cell_list,&pbbox,cell_ind,
                                        bounding_coordinates,n_bounding_cells);
                
                int o_cell = SqMS_find_where_particle_belongs(&cell_list,particle->r);
                for (size_t j = 0; j < pbbox.num_cells; j++)
                {
                    if (o_cell == pbbox.cell_indices[j])
                    {
                        in_box = 1;
                        break;
                    }
                    
                }
                if(in_box==0) printf("BUG: Interaction outside of the other Interaction bound!\n");
                assert(in_box == 1);
                SqMS_free_bounding_box(&pbbox);

            }
            t_energy += p_energy;
        }
    }
    
    return t_energy;
}
/* ___________________________  */
/*  MOVE PARTICLE SUB_ROUTINE   */
/* ---------------------------  */

int move_particle(double del, int step)
{
    int particle_ind = (int)(dsfmt_genrand()*n_particles);
    int cell_ind;
    for (cell_ind = 0; cell_ind < cell_list.num_cell; cell_ind++)
    {
        if (particle_ind >= cell_list.cell_population[cell_ind])
        {
            particle_ind -= cell_list.cell_population[cell_ind];
        }
        else
        {
            break;
        }
        
    }
    
    int full_index = cell_ind * cell_list.max_population_per_cell + particle_ind;

    particle_t chosen;
    SqMS_copy_particle(&cell_list.cell[full_index],&chosen);
    assert(chosen.isplaceholder == 0);

    double old_energy = SqMS_get_energy_from_pair(chosen.uid,chosen.uid,energy_matrix);

    if (SANITY_CHECK == 1 && step%SANITY_CHECK_PERIOD==0){
        bounding_box_t sbbox;
        
        SqMS_get_floppy_bounding_shape(&cell_list,&sbbox,cell_ind,
                                        bounding_coordinates,n_bounding_cells);
        
        double sanity_old_energy = sanity_particle_energy(&chosen,&sbbox);

        if (fabs(old_energy - sanity_old_energy) > tolerance)
        {
            printf("old_energy: %lf, sanity_check: %lf, diff: %lf\n",
            old_energy, sanity_old_energy, fabs(old_energy - sanity_old_energy));

            SqMS_dump_all_interactions(chosen.uid,n_particles, energy_matrix);
            dump_all_energies_related(&chosen,&sbbox);

            for (size_t i = 0; i < sbbox.num_cells; i++)
            {
                printf("%d ",sbbox.cell_indices[i]);
            }
            printf("\n");
        }
        SqMS_free_bounding_box(&sbbox);
        assert(fabs(old_energy - sanity_old_energy) < tolerance);
    } 
    

    for (size_t i = 0; i < NDIM; i++) //propose a particle displacement
    {
        chosen.r[i] += (2*dsfmt_genrand()-1)*del;
    }

    //PBC
    for(int j=0;j<NDIM;j++)
    {
        if(chosen.r[j] > 1)
        {
            chosen.r[j] -= 1;
        }
        
        if(chosen.r[j] < 0)
        {
            chosen.r[j] += 1 ;
        }
    }
    int new_cell_ind;
    
    if (SqMS_is_particle_in_cell(&cell_list,chosen.r,cell_ind) == 0)
    {
        new_cell_ind = SqMS_find_where_particle_belongs(&cell_list,chosen.r);
    }
    else
    {
        new_cell_ind = cell_ind;
    }
    
    bounding_box_t bbox;
    SqMS_get_floppy_bounding_shape(&cell_list,&bbox,new_cell_ind,bounding_coordinates,n_bounding_cells);

    //collision detection
    if (collision_detection == YES)
    {
        for (size_t c = 0; c < bbox.num_cells; c++)
        {
            for (size_t pp = 0; pp < cell_list.cell_population[bbox.cell_indices[c]]; pp++)
            {
                int index = bbox.cell_indices[c] * cell_list.max_population_per_cell + pp;
                if (chosen.uid != cell_list.cell[index].uid && cell_list.cell[index].isplaceholder == 0)
                {
                    double seperation2 = dist2(chosen.r, cell_list.cell[index].r);
                    double min_seperation = chosen.radius + cell_list.cell[index].radius;
                    if (seperation2 < min_seperation*min_seperation)
                    {
                        /* code */
                        SqMS_free_bounding_box(&bbox);
                        return 0;
                    }
                    
                }
                
            }
            
        }
        
    }
    double energy_backup[n_particles];
    for (size_t i = 0; i < n_particles; i++)
    {
        energy_backup[i] = SqMS_get_energy_from_pair(chosen.uid,i,energy_matrix);        
    }

    if (new_cell_ind != cell_ind)
    {
        for (size_t i = 0; i < n_particles; i++)
        {
            if (i !=chosen.uid)
            {
                SqMS_set_energy_of_pair(0,i,chosen.uid,energy_matrix);
            }
        }
    }
    
    

    particle_energy_cell_list(&bbox, &chosen);
    double new_energy = SqMS_get_energy_from_pair(chosen.uid,chosen.uid,energy_matrix);
    double dE = new_energy - old_energy;
    
    if (SANITY_CHECK == 1 && step%SANITY_CHECK_PERIOD==0){
        double sanity_new_energy = sanity_particle_energy(&chosen, &bbox);
        if (fabs(new_energy - sanity_new_energy) >= tolerance)
        {
            printf("new_energy: %lf, sanity_check: %lf, diff: %lf\n",
            new_energy, sanity_new_energy, fabs(new_energy - sanity_new_energy));
            
        }
        assert(fabs(new_energy - sanity_new_energy) < tolerance);
    }

    double rnumber = dsfmt_genrand();

    if(dE < 0.0 || rnumber < exp(-beta * dE)){
        
        /*if (dE > 0)
        {
            
           printf("dE is %lf, rnumber is %lf, bfactor is %lf, total energy is %lf\n",dE,rnumber,exp(-beta * dE), total_energy);
        }*/
        total_energy += dE;
        if (cell_ind == new_cell_ind)
        {
            SqMS_copy_particle(&chosen, &cell_list.cell[full_index]);
        }
        else
        {
            SqMS_remove_particle_from_cell_list(&cell_list,cell_ind,particle_ind);
            SqMS_add_particle_to_cell(&cell_list, &chosen, new_cell_ind);
        }
        SqMS_free_bounding_box(&bbox);
        
        return 1;
    }
    else
    {
        //if (new_cell_ind != cell_ind) SqMS_get_rigid_bounding_box(&cell_list,&bbox,cell_ind);
        //particle_energy_cell_list(&bbox, &cell_list.cell[full_index]);
        //to update energy matrix THIS CAN BE OPTIMIZED!

        for (size_t i = 0; i < n_particles; i++)
        {
            if (i==chosen.uid)
            {
                continue;
            }
            
            SqMS_set_energy_of_pair(energy_backup[i],i,
                                    chosen.uid, energy_matrix);
        }
        assert(fabs(energy_backup[chosen.uid] - 
                SqMS_get_energy_from_pair(chosen.uid,chosen.uid,energy_matrix))
                                                                    <tolerance);

        SqMS_free_bounding_box(&bbox);
        
        return 0;
    }

        
}

int move_volume(double delV, int step)
{
    //backup
    double old_floppyBox[NDIM*NDIM];
    memcpy(old_floppyBox,floppyBox,NDIM*NDIM*sizeof(double));
    int old_bounding_coordinates[n_bounding_cells*NDIM];
    memcpy(old_bounding_coordinates,bounding_coordinates,
            n_bounding_cells*NDIM*sizeof(int));
    int old_n_bounding_cells = n_bounding_cells;

    if (SANITY_CHECK == 1 && step%SANITY_CHECK_PERIOD==0){
        for (size_t i = 0; i < cell_list.max_population; i++)
        {
            particle_t *chosen = &cell_list.cell[i];
            if (chosen->isplaceholder == 0)
            {
                bounding_box_t sbbox;
                int cellofchosen = SqMS_find_where_particle_belongs(&cell_list,chosen->r);
                assert (cellofchosen == i/cell_list.max_population_per_cell);

                SqMS_get_floppy_bounding_shape(&cell_list,&sbbox,cellofchosen,
                                            bounding_coordinates,
                                            n_bounding_cells);

                double old_penergy = SqMS_get_energy_from_pair(chosen->uid,
                                                        chosen->uid,
                                                        energy_matrix);
                double sanity_old_penergy = sanity_particle_energy(chosen, &sbbox);
                if (fabs(old_penergy - sanity_old_penergy) >= tolerance)
                {
                    printf("old_energy: %lf, sanity_check: %lf, diff: %lf\n",
                    old_penergy, sanity_old_penergy, fabs(old_penergy - sanity_old_penergy));
                    
                }
                assert(fabs(old_penergy - sanity_old_penergy) < tolerance);
                SqMS_free_bounding_box(&sbbox);
            }
            
        }
    }
    
    //propose a change
    double change = 2*(dsfmt_genrand()-0.5)*delV;

    int choice = (int) (NDIM*NDIM*dsfmt_genrand()) ;
    floppyBox[choice] += change;
    
    
    
    //update bounding coords
    SqMS_get_interaction_shape(&int_shape, 
                               &cell_list,
                               max_length_cutoff,
                               floppyBox);
    SqMS_interaction_shape_to_bounding_shape(&int_shape, &bounding_shape);
    
    free(bounding_coordinates);
    bounding_coordinates = calloc((int_shape.height+1)*(int_shape.width+1)*NDIM, 
                                  sizeof(int));
    n_bounding_cells = (
    SqMS_bounding_shape_to_bounding_coordinates(bounding_shape,
                                                int_shape.height+1,
                                                int_shape.width+1,
                                                bounding_coordinates));
    SqMS_free_interaction_shape(&int_shape);
    free(bounding_shape);


    bounding_box_t bbox;
    
    //collision detection
    if (collision_detection == YES)
    {
        for (size_t cell_i = 0; cell_i < cell_list.num_cell; cell_i++)
        {
            SqMS_get_floppy_bounding_shape(&cell_list,&bbox,cell_i,
                                            bounding_coordinates,
                                            n_bounding_cells);

            for (size_t particle_i = 0; 
                 particle_i < cell_list.cell_population[cell_i]; particle_i++)
            {
                int particle1_full_i = (cell_i*cell_list.max_population_per_cell
                                        +  particle_i);
                
                for (size_t c = 0; c < bbox.num_cells; c++)
                {
                    for (size_t pp = 0; 
                         pp < cell_list.cell_population[bbox.cell_indices[c]]; 
                         pp++)
                    {
                        int index = (bbox.cell_indices[c] * 
                                     cell_list.max_population_per_cell + pp);

                        if (cell_list.cell[particle1_full_i].uid > 
                            cell_list.cell[index].uid && 
                            cell_list.cell[index].isplaceholder == 0)
                        {
                            double seperation = dist(
                                cell_list.cell[particle1_full_i].r, 
                                cell_list.cell[index].r);
                            double min_seperation = (
                                cell_list.cell[particle1_full_i].radius 
                                + cell_list.cell[index].radius);
                            if (seperation < min_seperation)
                            {
                                /* code */
                                n_bounding_cells= old_n_bounding_cells;
                                memcpy(floppyBox,old_floppyBox,NDIM*NDIM*sizeof(double));

                                n_bounding_cells= old_n_bounding_cells;
                                free(bounding_coordinates);
                                bounding_coordinates = calloc(n_bounding_cells*NDIM,sizeof(int));
                                memcpy(bounding_coordinates,old_bounding_coordinates,
                                        n_bounding_cells*NDIM*sizeof(int));

                                SqMS_free_bounding_box(&bbox);
                                return 0;
                            }
                            
                        }
                        
                    }
                    
                }

            }
            SqMS_free_bounding_box(&bbox);
        }   
        
    }

    


    
    double old_volume = volume;
    volume = SqMS_determinant(floppyBox, NDIM);

    /*
    TODO: energy update
    */
    double old_energy = SqMS_get_total_energy(n_particles,energy_matrix);
    assert(fabs(old_energy-total_energy)<tolerance);

    double old_energy_matrix[ (n_particles * (n_particles + 1))/2 ];
    memcpy(old_energy_matrix,energy_matrix,
            sizeof(double)*(n_particles * (n_particles + 1))/2);
    
    for (size_t cell_i = 0; cell_i < cell_list.num_cell; cell_i++)
    {
        SqMS_get_floppy_bounding_shape(&cell_list,&bbox,cell_i,
                                        bounding_coordinates,
                                        n_bounding_cells);

        for (size_t particle_i = 0; 
                particle_i < cell_list.cell_population[cell_i]; particle_i++)
        {
            int particle1_full_i = (cell_i*cell_list.max_population_per_cell
                                    +  particle_i);

            assert(cell_list.cell[particle1_full_i].isplaceholder==0);
            
            for (size_t c = 0; c < bbox.num_cells; c++)
            {
                for (size_t pp = 0; 
                        pp < cell_list.cell_population[bbox.cell_indices[c]]; 
                        pp++)
                {
                    int index = (bbox.cell_indices[c] * 
                                    cell_list.max_population_per_cell + pp);

                    assert(cell_list.cell[index].isplaceholder==0);

                    if (cell_list.cell[particle1_full_i].uid > 
                        cell_list.cell[index].uid)
                    {
                        double p_dist2 = dist2(
                                        cell_list.cell[particle1_full_i].r,
                                        cell_list.cell[index].r);
                        double p_energy = pair_potential(p_dist2,
                            cell_list.cell[index].radius
                            +cell_list.cell[particle1_full_i].radius);
                        assert(p_energy <= 0);

                        SqMS_set_energy_of_pair(
                            p_energy,cell_list.cell[particle1_full_i].uid,
                            cell_list.cell[index].uid, energy_matrix);
                    }
                    
                }
                
            }

        }
        SqMS_free_bounding_box(&bbox);
        
    }
    double new_energy = SqMS_get_total_energy(n_particles,energy_matrix);

    if (SANITY_CHECK == 1 && step%SANITY_CHECK_PERIOD==0){
        for (size_t particle_i = 0; 
            particle_i < cell_list.max_population; 
            particle_i++)
        {
            particle_t *chosen = &cell_list.cell[particle_i];
            if (chosen->isplaceholder==1) continue;
            int chosen_cell = SqMS_find_where_particle_belongs(&cell_list,
                                                                chosen->r);
            double new_penergy = SqMS_get_energy_from_pair(chosen->uid,
                                                        chosen->uid,
                                                        energy_matrix);
            bounding_box_t sbbox;
            SqMS_get_floppy_bounding_shape(&cell_list,&sbbox,chosen_cell,
                                        bounding_coordinates,
                                        n_bounding_cells);
                                        
            double sanity_new_energy = sanity_particle_energy(chosen, &sbbox);
            if (fabs(new_penergy - sanity_new_energy) >= tolerance)
            {
                printf("vol new_energy: %lf, sanity_check: %lf, diff: %lf\n",
                new_penergy, sanity_new_energy, fabs(new_penergy - sanity_new_energy));
                
            }
            assert(fabs(new_penergy - sanity_new_energy) < tolerance);
            SqMS_free_bounding_box(&sbbox);

        }
    }

    double minus_beta_Delta_Free_energy = (-new_energy + old_energy 
                                - pressure * (volume - old_volume)
                                + n_particles * log(volume/old_volume));
   

   if (minus_beta_Delta_Free_energy>=0 || 
        dsfmt_genrand()< exp(minus_beta_Delta_Free_energy))
   {
       total_energy = new_energy;

       
       return 1;
   }
   else
   {
        volume = old_volume;
        n_bounding_cells= old_n_bounding_cells;
        memcpy(floppyBox,old_floppyBox,NDIM*NDIM*sizeof(double));

        n_bounding_cells= old_n_bounding_cells;
        free(bounding_coordinates);
        bounding_coordinates = calloc(n_bounding_cells*NDIM,sizeof(int));
        memcpy(bounding_coordinates,old_bounding_coordinates,
                n_bounding_cells*NDIM*sizeof(int));
        
        memcpy(energy_matrix,old_energy_matrix,
            sizeof(double)*(n_particles * (n_particles + 1))/2);

       return 0;
   }
   
}

/* --------------------------------------------------------------------- */


void init_positions_hot_cell_list(void)
{
    printf("Hot initialization has started.\n");
    double rr[NDIM];
    size_t placed_A = 0, placed_B=0;
    enum particle_type type;
    for (size_t i = 0; i < n_particles; i++)
    {
        /* particle loop */
        int is_collided = 0;
        double radius = 0;
        if (placed_A < n_particlesA)
        {
            radius = radiusA;
            type = A;
        }
        else if (placed_B < n_particlesB)
        {
            radius = radiusB;
            type = B ;
        }
        else
        {
            printf("something is wrong with particle numbers!\n");
            exit(-1);
        }
        
        

        size_t iter = 0;

        for (iter = 0; iter < ITER_MAX; iter++)
        {
            for (size_t d = 0; d < NDIM; d++)
            {
                /* dimension loop */
                rr[d]= dsfmt_genrand();
            }
            if (collision_detection == YES)
            {
                /* if particles have hard cores initialize accordingly */
                is_collided = 0;
                for (size_t j = 0; j < cell_list.max_population; j++)
                {
                    if(cell_list.cell[j].isplaceholder==0 && dist(rr, cell_list.cell[j].r) < radius + cell_list.cell[j].radius)
                    {
                        is_collided = 1;
                        break;
                    }
                }
            }
            if (is_collided == 0) break;
            
            
        }
        if (iter == ITER_MAX)
        {
            printf("Max iterations has been reached!\n");
            exit(-1);
        }

        int cell_index = SqMS_find_where_particle_belongs(&cell_list, rr);
        particle_t new_particle;
        new_particle.isplaceholder = 0;
        new_particle.radius = radius;
        new_particle.uid = i;
        for (size_t k = 0; k < NDIM; k++)
        {
            new_particle.r[k] = rr[k];
        }
        SqMS_add_particle_to_cell(&cell_list,&new_particle,cell_index);
        if (type == A)
        {
            placed_A++;
        }
        else if (type == B)
        {
            placed_B++;
        }
        else
        {
            printf("something is wrong with particle types!\n");
            exit(-1);
        }
        
        
    }
    printf("Hot initialization has ended.\n");
    if (SANITY_CHECK == 1)
    {
        int tot=0;
        for (size_t i = 0; i < cell_list.num_cell; i++)
        {
            tot += cell_list.cell_population[i];
        }
        assert(tot == n_particles);

        for (size_t i = 0; i < cell_list.num_cell; i++)
        {
            tot = 0;
            for (size_t j = 0; j < cell_list.cell_population[i]; j++)
            {
                if (cell_list.cell[i * cell_list.max_population_per_cell + j].isplaceholder == 0)
                {
                    tot++;
                    assert(SqMS_find_where_particle_belongs(&cell_list, cell_list.cell[i * cell_list.max_population_per_cell + j].r)==i);
                }
                
            }
            assert(tot == cell_list.cell_population[i]);
        }
        
    }
    
}

double DeltaAdjuster(double DD, double acceptance, double target)
{
    double newDelta = DD;
    double A = -log(3)/log(target);
    double scale = 1.5 * pow(acceptance, A) + 0.5;
    newDelta *= scale;
//    printf("Move Acceptence is %lf. Delta is rescaled by %lf\n",acceptance,scale);
    
    return newDelta;
}


void find_delta()
{
    double backup_floppyBox[NDIM*NDIM];
    memcpy(backup_floppyBox,floppyBox,NDIM*NDIM*sizeof(double));

    printf("Adjusting Delta:\n");
    double m_acceptance=0, m_mean_acc=0, m_mean2_acc=0;
    double v_acceptance=0, v_mean_acc=0, v_mean2_acc=0;
    init_positions_hot_cell_list();
    init_particle_energies_cell();
    
    int count = 5;
    double tol = 0.01;
    double target = 0.6;

    double move_arr_acc[count];
    double volume_arr_acc[count];
    int move_accepted = 0,volume_accepted = 0,move_trials = 0,volume_trials = 0, step;

    for(step = 0; step <= (count-1)*output_steps; ++step){
        for(int n = 0; n < n_particles; ++n){
            if (pressure > 0 && dsfmt_genrand() < 1./(double)(n_particles+1))
            {
                volume_accepted += move_volume(Vdelta, 1);
                volume_trials++;
            }
            else
            {
                move_accepted += move_particle(delta, 1);
                move_trials++;
            }
        }
        if(step % output_steps == 0 && step != 0)
        {
            m_acceptance = (double)move_accepted / (double) move_trials;
            move_arr_acc[step / output_steps] = m_acceptance;
            delta = DeltaAdjuster(delta, m_acceptance, target);
            move_accepted = 0;
            move_trials = 0;

            v_acceptance = (double)volume_accepted / (double) volume_trials;
            volume_arr_acc[step / output_steps] = v_acceptance;
            Vdelta = DeltaAdjuster(Vdelta, v_acceptance, target);
            volume_accepted = 0;
            volume_trials = 0;
        }
    }
    move_accepted = 0;
    move_trials = 0;
    volume_accepted = 0;
    volume_trials = 0;
    for(step = 0; step <= NAdjust; ++step){
        for(int n = 0; n < n_particles; ++n){
            if (pressure > 0 && dsfmt_genrand() < 1./(double)(n_particles+1))
            {
                volume_accepted += move_volume(Vdelta, 1);
                volume_trials++;
            }
            else
            {
                move_accepted += move_particle(delta, 1);
                move_trials++;
            }
        }
        if(step % output_steps == 0 && step != 0)
        {
            m_acceptance = (double)move_accepted / (double) move_trials;
            move_arr_acc[(step/output_steps)%count] = m_acceptance;

            v_acceptance = (double)volume_accepted / (double) volume_trials;
            volume_arr_acc[(step/output_steps)%count] = v_acceptance;
            
            move_accepted = 0;
            move_trials = 0;
            volume_accepted = 0;
            volume_trials = 0;
            
            printf("Trial Step %d, Move acceptance: %f, Volume acceptance: %f.\n",
                step / output_steps, m_acceptance,v_acceptance
            );
            m_mean_acc = 0;
            m_mean2_acc = 0;
            v_mean_acc = 0;
            v_mean2_acc = 0;
            for (size_t ii = 0; ii < count; ii++)
            {
                m_mean_acc += move_arr_acc[ii];
                m_mean2_acc += move_arr_acc[ii]*move_arr_acc[ii];
                v_mean_acc += volume_arr_acc[ii];
                v_mean2_acc += volume_arr_acc[ii]*volume_arr_acc[ii];
            }
            m_mean_acc /= count;
            m_mean2_acc /= count;
            v_mean_acc /= count;
            v_mean2_acc /= count;

            double m_variance = m_mean2_acc - m_mean_acc * m_mean_acc;
            double v_variance = v_mean2_acc - v_mean_acc * v_mean_acc;

            if (fabs(target - m_mean_acc)<tol 
                && sqrt(m_variance) < tol
                && fabs(target-v_mean_acc)<4*tol 
                && sqrt(v_variance) < 16*tol )
            {
                printf("Deltas found! Delta = %lf,mean acceptance: %lf, variance: %lf.\nVDelta = %lf,mean acceptance: %lf, variance: %lf.\n\nStarting MC Simulation.\n",
                delta, m_mean_acc, m_variance,
                Vdelta, v_mean_acc, v_variance);
                printf("----------------------------------------\n");
                break;
            }
            else
            {
                delta = DeltaAdjuster(delta, m_acceptance, target);
                delta = fmod(delta, 1.);
                Vdelta = DeltaAdjuster(Vdelta, v_acceptance, target);
            }
        }
    }
    if (step == NAdjust+1)
    {
        printf("Acceptence(Move: %lf, Volume: %lf) could not be adjusted to fall in the desired interval!\n Delta = %lf, VDelta = %lf.\n\nStarting MC Simulation.\n",m_acceptance,v_acceptance,delta,Vdelta);
        printf("----------------------------------------\n");
    }

    memcpy(floppyBox,backup_floppyBox,NDIM*NDIM*sizeof(double));
    volume = SqMS_determinant(floppyBox, NDIM);
}

/* _____________________________________________________________________ */
/* _____________________________________________________________________ */
/* _____________________________________________________________________ */
/* _____________________________________________________________________ */

int main(int argc, const char * argv[])
{
    read_parameters();
    
    assert(radiusA > 0.0);
    assert(radiusB > 0.0);
    assert(delta > 0.0);
    n_particles = n_particlesA + n_particlesB;
    energy_matrix = calloc((n_particles * (n_particles + 1))/2,sizeof(double));
    max_length_cutoff = 1.*2.2*fmax(radiusB,radiusA);

    for (size_t i = 0; i < NDIM*NDIM; i++)
    {
        floppyBox[i] = 0;
    }
    for (size_t i = 0; i < NDIM; i++)
    {
        floppyBox[NDIM*i+i] = Box_coef*max_length_cutoff;
    }

    if(NDIM == 3) 
    {
        particle_volumeA = 4. * M_PI * pow(radiusA, 3.0) / 3.0;
        particle_volumeB = 4. * M_PI * pow(radiusB, 3.0) / 3.0;
    }
    else if(NDIM == 2)
    {
        particle_volumeA = M_PI * pow(radiusA, 2.0);
        particle_volumeB = M_PI * pow(radiusB, 2.0);
    }
    else
    {
        particle_volumeA = pow(M_PI,(double) NDIM/2.) * pow(radiusA, NDIM) / tgamma(1.+(double) NDIM/2.);
        particle_volumeB = pow(M_PI,(double) NDIM/2.) * pow(radiusB, NDIM) / tgamma(1.+(double) NDIM/2.);
    }
    
    volume = SqMS_determinant(floppyBox, NDIM);
    assert(volume > particle_volumeA*n_particlesA + particle_volumeB*n_particlesB);
    
    
    
    //write_data(-1); //writes initial config
    // read_data(); //reads initial config (just a proof of concept)
    
    //dsfmt_seed ((int) time (NULL) ) ;
    dsfmt_seed (314159) ;
    
    
    SqMS_init_floppy_cell_list(&cell_list, floppyBox, max_length_cutoff, fmin(particle_volumeB,particle_volumeA));
    assert(cell_list.max_population > n_particles);
    
    SqMS_get_interaction_shape(&int_shape,&cell_list,max_length_cutoff,floppyBox);
    SqMS_interaction_shape_to_bounding_shape(&int_shape, &bounding_shape);
    bounding_coordinates = calloc((int_shape.height+1)*(int_shape.width+1)*NDIM, sizeof(int));
    n_bounding_cells = SqMS_bounding_shape_to_bounding_coordinates(bounding_shape,int_shape.height+1,int_shape.width+1,bounding_coordinates);
    SqMS_free_interaction_shape(&int_shape);
    free(bounding_shape);

    find_delta();
    SqMS_free_cell_list(&cell_list);
    total_energy = 0 ;
    SqMS_init_floppy_cell_list(&cell_list, floppyBox, max_length_cutoff, fmin(particle_volumeB,particle_volumeA));
    init_positions_hot_cell_list();
    init_particle_energies_cell();
    write_data_cell(0);
    



    
    int move_accepted = 0;
    int volume_accepted = 0;
    int move_trials = 0;
    int volume_trials = 0;
    
    int step, n;
    for(step = 1; step < mc_steps+1; ++step){
        if (SANITY_CHECK == 1 && step % SANITY_CHECK_PERIOD == 0)
        {
            printf("This sweep(%d) is a sanity check step!\n",step);
        }
        
        for(n = 0; n < n_particles+1; ++n){
            if (pressure > 0 && dsfmt_genrand() < 1./(double)(n_particles+1))
            {
                volume_accepted += move_volume(Vdelta, step);
                volume_trials++;
            }
            else
            {
                move_accepted += move_particle(delta, step);
                move_trials++;
            }
        }

        if(step % output_steps == 0){
            printf("Sweep %d. Move acc: %.2lf (%d/%d), Volume acc: %.2lf (%d/%d), total energy: %lf.\n",
            step, (double)move_accepted / (double) move_trials,
            move_accepted,move_trials,
            (double)volume_accepted / (double) volume_trials,
            volume_accepted, volume_trials,
            total_energy);
            move_accepted = 0;
            volume_accepted = 0;
            move_trials = 0;
            volume_trials = 0;
            write_data_cell(step);
        }
    }

    free(bounding_coordinates);
    SqMS_free_cell_list(&cell_list);
    return 0;
}

