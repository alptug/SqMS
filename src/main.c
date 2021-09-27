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

//#define NDIM 2 //number of dimensions
#define N 15
#define A_coef 1
#define B_coef 24
#define Box_coef 10
#define NMAX  ((A_coef + B_coef)*N)//max number of particles
#define ITER_MAX 100000000 //max number of iterations
#define SANITY_CHECK 1

const int n_particlesA = A_coef*N; //number of particles
const int n_particlesB = B_coef*N; //number of particles
int n_particles; //number of particles
enum detection{YES, NO};
const int mc_steps = 20000; //it says steps but they are in fact sweeps
const int NAdjust = mc_steps / 5;
const int output_steps = 100; //output a file every this many sweeps
//const double diameter = 1.0; //particle diameter
double delta = 0.004; //MCMC move delta
double beta = 5; //inverse temperature
const char* init_filename = "coords_step-000001.dat";
enum detection collision_detection = YES; //detect hard cores
enum particle_type{A,B};
cell_list_t cell_list;

const double eps = 1e-1;
const double tolerance = 1e-6;

/* ---------------------  */
/*  Simulation variables  */
/* ---------------------  */

const double radiusA = 1., radiusB = .345;
double particle_volumeA,particle_volumeB,max_length_cutoff;

double total_energy;
double energy_matrix[ (NMAX * (NMAX + 1))/2 ];
double floppyBox[NDIM*NDIM];

int n_bounding_cells;
int *bounding_coordinates, *bounding_shape;
interaction_shape_t int_shape;

/* ---------------------------------  */
/*  FILE READER & WRITER SUBROUTINES  */
/* ---------------------------------  */

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
    char buffer[128];
    sprintf(buffer, "data/coords_step%07d.dat", step);
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,Box_coef*max_length_cutoff);
    }
    if (NDIM == 2) fprintf(fp, "%lf %lf\n",0.0,Box_coef*max_length_cutoff);

    for (size_t i = 0; i < cell_list.max_population; i++)
    {
        /* code */
        if (cell_list.cell[i].isplaceholder == 0)
        {
            for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", cell_list.cell[i].r[d]*Box_coef*max_length_cutoff);
            if (NDIM == 2) fprintf(fp, "%f\t", 0.);
            fprintf(fp, "%lf\n", 2.*cell_list.cell[i].radius);
        }
        
    }
    
    fclose(fp);
}

/* ___________________________  */
/*  Euclidean metric with PBC   */
/* ---------------------------  */

double dist(double *particle1, double *particle2)
{
    //dimension agnostic pbc distance
    return SqMS_distance_floppy_PBC(particle1, particle2, floppyBox, NDIM);
}

double dist2(double *particle1, double *particle2)
{
    //dimension agnostic pbc distance
    return SqMS_distance_squared_floppy_PBC(particle1, particle2, floppyBox, NDIM);
}

/* ___________________________  */
/*       energy functions       */
/* ---------------------------  */

double pair_potential(double distance_square, double total_radius)
{
    return SqMS_finite_well_potential(distance_square, 1., 1.1*1.1*total_radius*total_radius);
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

    /*particle_t particle;
    SqMS_copy_particle(&cell_list.cell[cell_index * cell_list.max_population_per_cell + particle_ind], &particle);

    bounding_box_t bbox;
    SqMS_get_bounding_box(&cell_list,&bbox,cell_index);*/
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
    /* RESERVED FOR SANITY CHECK
    for(size_t i = 0; i<cell_list.max_population; i++)
    {
        if (cell_list.cell[i].isplaceholder == 1 || cell_list.cell[i].uid == particle.uid) continue;

        p_dist2 = dist2(particle.r,cell_list.cell[i].r);
        if (p_dist2 < LJ_length_cutoff2)
        {
            p_energy = pair_potential(p_dist2);
            SqMS_set_energy_of_pair(p_energy,cell_list.cell[i].uid,particle.uid,energy_matrix);
            t_energy += p_energy;
        }
    }*/
    double check_energy = SqMS_get_energy_from_pair(particle->uid,particle->uid,energy_matrix);
    double diff = fabs(t_energy - check_energy);
    if (diff > eps)
    {
        //assert(diff < eps);
        //printf("fp error?\n");
    }
    
    return t_energy;
}

/* ___________________________  */
/*  MOVE PARTICLE SUB_ROUTINE   */
/* ---------------------------  */

int move_particle(double del)
{
    
    int cell_ind = (int)(dsfmt_genrand()*cell_list.num_cell);
    int cell_population = cell_list.cell_population[cell_ind];
    size_t dummy = 0;
    for (dummy = 0; dummy < ITER_MAX; dummy++)
    {
        if (cell_population == 0)
        {
            cell_ind = (int)(dsfmt_genrand()*cell_list.num_cell);
            cell_population = cell_list.cell_population[cell_ind];
        }
        else
        {
            break;
        }
    }
    if (dummy == ITER_MAX)
    {
        printf("Max iterations has been reached in cell selection\n");
        exit(-2);
    }
    int particle_ind = (int)(dsfmt_genrand()* cell_population ); //choose a random particle

    int full_index = cell_ind * cell_list.max_population_per_cell + particle_ind;

    particle_t chosen;
    SqMS_copy_particle(&cell_list.cell[full_index],&chosen);
    assert(chosen.isplaceholder == 0);
    //double old_energy = SqMS_get_energy_from_pair(n,n,energy_matrix);
    //double old_energy = particle_energy(n);

    double old_energy = SqMS_get_energy_from_pair(chosen.uid,chosen.uid,energy_matrix);

    bounding_box_t sbbox;
    SqMS_get_floppy_bounding_shape(&cell_list,&sbbox,cell_ind,bounding_coordinates,n_bounding_cells);

    double sanity_old_energy = 0;
    if (SANITY_CHECK == 1) sanity_old_energy = particle_energy_cell_list(&sbbox, &chosen);
    

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
    SqMS_get_floppy_bounding_shape(&cell_list,&bbox,cell_ind,bounding_coordinates,n_bounding_cells);

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
                    double seperation = dist(chosen.r, cell_list.cell[index].r);
                    double min_seperation = chosen.radius + cell_list.cell[index].radius;
                    if (seperation < min_seperation)
                    {
                        /* code */
                        SqMS_free_bounding_box(&bbox);
                        SqMS_free_bounding_box(&sbbox);
                        return 0;
                    }
                    
                }
                
            }
            
        }
        
    }

    //pair potential
    //double new_energy = particle_energy(n);
    //particle_energy(n);

    double sanity_new_energy = particle_energy_cell_list(&bbox, &chosen);
    double new_energy = SqMS_get_energy_from_pair(chosen.uid,chosen.uid,energy_matrix);
    double dE = new_energy - old_energy;
    //double dE = sanity_new_energy - sanity_old_energy;
    double rnumber = dsfmt_genrand();

    if(dE < 0.0 || rnumber < exp(-beta * dE)){
        
        if (dE > 0)
        {
            
           //printf("dE is %lf, rnumber is %lf, bfactor is %lf, total energy is %lf\n",dE,rnumber,exp(-beta * dE), total_energy);
        }
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
        SqMS_free_bounding_box(&sbbox);
        return 1;
    }
    else
    {
        if (new_cell_ind != cell_ind) SqMS_get_rigid_bounding_box(&cell_list,&bbox,cell_ind);
        particle_energy_cell_list(&bbox, &cell_list.cell[full_index]);
        //to update energy matrix THIS CAN BE OPTIMIZED!
        SqMS_free_bounding_box(&bbox);
        SqMS_free_bounding_box(&sbbox);
        return 0;
    }

        
}


/* --------------------------------------------------------------------- */
/*      EXPERIMENT FUNCTIONS --> Run MCMC for particular objectives      */
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
    double minl = 0;
    newDelta = fmod(newDelta, 1.);
//    printf("Move Acceptence is %lf. Delta is rescaled by %lf\n",acceptance,scale);
    
    return newDelta;
}


void find_delta()
{
    printf("Adjusting Delta:\n");
    double acceptance=0, mean_acc=0, mean2_acc=0;
    init_positions_hot_cell_list();
    init_particle_energies_cell();
    
    int count = 5;
    double tol = 0.02;
    double target = 0.6;

    double arr_acc[count];
    int accepted = 0, step;

    for(step = 0; step <= (count-1)*output_steps; ++step){
        for(int n = 0; n < n_particles; ++n){
            accepted += move_particle(delta);
        }
        if(step % output_steps == 0)
        {
            acceptance = (double)accepted / (n_particles * output_steps);
            arr_acc[step / output_steps] = acceptance;
            delta = DeltaAdjuster(delta, acceptance, target);
            accepted = 0;
        }
    }

    for(step = 0; step <= NAdjust; ++step){
        for(int n = 0; n < n_particles; ++n){
            accepted += move_particle(delta);
        }
        if(step % output_steps == 0 && step != 0)
        {
            acceptance = (double)accepted / (n_particles * output_steps);
            arr_acc[(step/output_steps)%count] = acceptance;
            printf("Trial Step %d. Move acceptance: %f.\n",
                step / output_steps, acceptance
            );
            mean_acc = 0;
            mean2_acc = 0;
            for (size_t ii = 0; ii < count; ii++)
            {
                mean_acc += arr_acc[ii];
                mean2_acc += arr_acc[ii]*arr_acc[ii];
            }
            mean_acc /= count;
            mean2_acc /= count;
            double variance = mean2_acc - mean_acc * mean_acc;

            if (fabs(target-mean_acc)<tol && sqrt(variance) < tol)
            {
                printf("Delta found! Delta = %lf,mean acceptance: %lf, variance: %lf.\n\nStarting MC Simulation.\n",delta, mean_acc, variance);
                printf("----------------------------------------\n");
                break;
            }
            else
            {
                delta = DeltaAdjuster(delta, acceptance, target);
            }
            accepted = 0;
        }
    }
    if (step == NAdjust+1)
    {
        printf("Acceptence(%lf) could not be found in the desired interval!\n Delta = %lf.\n\nStarting MC Simulation.\n",acceptance,delta);
        printf("----------------------------------------\n");
    }
    
}

/* _____________________________________________________________________ */
/* _____________________________________________________________________ */
/* _____________________________________________________________________ */
/* _____________________________________________________________________ */

int main(int argc, const char * argv[])
{
    assert(radiusA > 0.0);
    assert(radiusB > 0.0);
    assert(delta > 0.0);
    n_particles = n_particlesA + n_particlesB;
    
    max_length_cutoff = 1.*2.2*fmax(radiusB,radiusA);

    for (size_t i = 0; i < NDIM*NDIM; i++)
    {
        floppyBox[NDIM*i+i] = 0;
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
    if (collision_detection == YES)
    {
        double box_volume = SqMS_determinant(floppyBox, NDIM);
        assert(box_volume > particle_volumeA*n_particlesA + particle_volumeB*n_particlesB);
    }
    
    
    //write_data(-1); //writes initial config
    // read_data(); //reads initial config (just a proof of concept)
    
    //dsfmt_seed ((int) time (NULL) ) ;
    dsfmt_seed (314159) ;
    
    
    SqMS_init_floppy_cell_list(&cell_list, floppyBox, max_length_cutoff, fmin(particle_volumeB,particle_volumeA));

    //init_positions_hot();
    assert(cell_list.max_population > n_particles);
    
    SqMS_get_interaction_shape(&int_shape,&cell_list,max_length_cutoff,floppyBox);
    SqMS_interaction_shape_to_bounding_shape(&int_shape, &bounding_shape);
    bounding_coordinates = calloc((int_shape.height+1)*(int_shape.width+1)*NDIM, sizeof(int));
    n_bounding_cells = SqMS_bounding_shape_to_bounding_coordinates(bounding_shape,int_shape.height+1,int_shape.width+1,bounding_coordinates);
    SqMS_free_interaction_shape(&int_shape);
    free(bounding_shape);

    find_delta();
    SqMS_init_floppy_cell_list(&cell_list, floppyBox, max_length_cutoff, fmin(particle_volumeB,particle_volumeA));
    init_positions_hot_cell_list();
    //init_particle_energies();
    init_particle_energies_cell();
    write_data_cell(0);

    //total_energy = 0 ;

    //write_data(0); //writes initial config

    
    int accepted = 0;
    int step, n;
    for(step = 1; step < mc_steps+1; ++step){
        for(n = 0; n < n_particles; ++n){
            accepted += move_particle(delta);
        }

        if(step % output_steps == 0){
            printf("Step %d. Move acceptance: %lf, total energy: %lf.\n",
            step, (double)accepted / (n_particles * output_steps), total_energy);
            accepted = 0;
            write_data_cell(step);
        }
    }


    return 0;
}

