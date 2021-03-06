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
//#include <mt19937.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------  */
/* Initialization variables  */
/* ------------------------  */

#define NDIM 3 //number of dimensions
#define NMAX 32 //max number of particles
#define ITER_MAX 100000000 //max number of iterations

int n_particles = NMAX; //number of particles
enum detection{YES, NO};
const int mc_steps = 100000; //it says steps but they are in fact sweeps
const int output_steps = 100; //output a file every this many sweeps
const double diameter = 1.0; //particle diameter
double delta = 0.02; //MCMC move delta
double beta = 100; //inverse temperature
const char* init_filename = "coords_step-000001.dat";
enum detection collision_detection = NO; //detect hard cores

const double eps = 1e-1;

/* ---------------------  */
/*  Simulation variables  */
/* ---------------------  */


double radius, LJ_length_cutoff2, LJ_energy_cutoff,LJ_length_scale_squared;
double particle_volume;
double r[NMAX][NDIM];
//double particle_energies[NMAX];
double total_energy;
double energy_matrix[ (NMAX * (NMAX + 1))/2 ];
double Box[NDIM];

/* ---------------------------------  */
/*  FILE READER & WRITER SUBROUTINES  */
/* ---------------------------------  */

void read_data(void){
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
}

void write_data(int step){
    char buffer[128];
    sprintf(buffer, "data/coords_step%07d.dat", step);
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,Box[d]);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", r[n][d]);
        fprintf(fp, "%lf\n", diameter);
    }
    fclose(fp);
}



/* ___________________________  */
/*  Euclidean metric with PBC   */
/* ---------------------------  */

double dist(double *particle1, double *particle2)
{
    //dimension agnostic pbc distance
    return SqMS_distance_rectangular_PBC(particle1, particle2, Box, NDIM);
}

double dist2(double *particle1, double *particle2)
{
    //dimension agnostic pbc distance
    return SqMS_distance_squared_rectangular_PBC(particle1, particle2, Box, NDIM);
}

/* ___________________________  */
/*       energy functions       */
/* ---------------------------  */

double pair_potential(double distance_square)
{
    return SqMS_truncated_LJ_potential_unsafe(distance_square, LJ_length_scale_squared, 1., LJ_energy_cutoff);
}

void init_particle_energies(void)
{
    SqMS_init_energy_matrix(n_particles,energy_matrix);
    total_energy = 0;
    double p_energy = 0;
    double p_dist2 = 0;

    for(size_t i = 0; i<n_particles; i++)
    {
        for(size_t j = 0; j<i; j++)
        {
            p_dist2 = dist2(r[i],r[j]);
            if (p_dist2 < LJ_length_cutoff2)
            {
                p_energy = pair_potential(p_dist2);
                SqMS_set_energy_of_pair(p_energy,i,j,energy_matrix);
                total_energy += p_energy;
            }
        }
    }
}

double particle_energy(int particle_id)
{
    double p_energy = 0;
    double t_energy = 0;
    double p_dist2 = 0;
    for(size_t i = 0; i<n_particles; i++)
    {
        if (particle_id == i) continue;

        p_dist2 = dist2(r[i],r[particle_id]);
        if (p_dist2 < LJ_length_cutoff2)
        {
            p_energy = pair_potential(p_dist2);
            SqMS_set_energy_of_pair(p_energy,i,particle_id,energy_matrix);
            t_energy += p_energy;
        }
    }
    double check_energy = SqMS_get_energy_from_pair(particle_id,particle_id,energy_matrix);
    double diff = fabs(t_energy - check_energy);
    if (diff > eps)
    {
        assert(diff < eps);
    }
    
    return t_energy;
}

/* ___________________________  */
/*  MOVE PARTICLE SUB_ROUTINE   */
/* ---------------------------  */

int move_particle(double del)
{
    int n = (int)(dsfmt_genrand()*n_particles); //choose a random particle
    
    double backup[NDIM]; //backup the position in case the proposed displacement gets rejected
    memcpy(backup, r[n], NDIM*sizeof(double));
    double old_energy = SqMS_get_energy_from_pair(n,n,energy_matrix);
    //double old_energy = particle_energy(n);

    for (size_t i = 0; i < NDIM; i++) //propose a particle displacement
    {
        r[n][i] += (2*dsfmt_genrand()-1)*del;
    }

    //PBC
    for(int j=0;j<NDIM;j++)
    {
        if(r[n][j] > Box[j])
        {
            r[n][j] -= Box[j] ;
        }
        
        if(r[n][j] < 0)
        {
            r[n][j] += Box[j] ;
        }
    }
    //collision detection
    if (collision_detection == YES)
    {
        
        for(int p=0; p<n_particles;p++)
        {
            if(n==p) continue ;
            
            if(dist(r[n], r[p]) < diameter)
            {
                memcpy(r[n], backup, NDIM*sizeof(double));
                return 0 ;
            }
        }
    }
    //pair potential
    //double new_energy = particle_energy(n);
    particle_energy(n);
    double new_energy = SqMS_get_energy_from_pair(n,n,energy_matrix);
    double dE = new_energy - old_energy;
    double rnumber = dsfmt_genrand();

    if(dE < 0.0 || rnumber < exp(-beta * dE)){
        
       /* if (dE > 0)
        {
            
            printf("dE is %lf, rnumber is %lf, bfactor is %lf, total energy is %lf\n",dE,rnumber,exp(-beta * dE), total_energy);
        }*/
        total_energy += dE;
        return 1;
    }
    else
    {
        memcpy(r[n], backup, NDIM*sizeof(double));
        particle_energy(n); //to update energy matrix THIS CAN BE OPTIMIZED!
        return 0;
    }

        
}


/* --------------------------------------------------------------------- */
/*      EXPERIMENT FUNCTIONS --> Run MCMC for particular objectives      */
/* --------------------------------------------------------------------- */

double MCsim(double del)
{
    int accepted = 0;
    int step, n;
    for(step = 0; step < mc_steps; ++step)
    {
        for(n = 0; n < n_particles; ++n)
        {
            accepted += move_particle(del);
        }
        if(step % output_steps == 0){
            printf("Step %d. Move acceptance: %lf, total energy: %lf.\n", 
            step, (double)accepted / (n_particles * output_steps), total_energy);
            accepted = 0;
            write_data(step);
        }
    }
    double acceptance = (double)accepted / (n_particles * mc_steps);
    printf("MC complete. Move acceptance: %lf.\n", acceptance);
    return acceptance;
}

void init_positions_hot()
{
    printf("Hot initialization has started.\n");
    for (size_t i = 0; i < n_particles; i++)
    {
        /* particle loop */
        int is_collided = 0;
        size_t iter = 0;

        for (iter = 0; iter < ITER_MAX; iter++)
        {
            for (size_t d = 0; d < NDIM; d++)
            {
                /* dimension loop */
                r[i][d]= Box[d] * dsfmt_genrand();                
            }
            if (collision_detection == YES)
            {
                /* if particles have hard cores initialize accordingly */
                is_collided = 0;
                for (size_t j = 0; j < i; j++)
                {
                    if(dist(r[i], r[j]) < diameter)
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
        
    }
    printf("Hot initialization has ended.\n");
}

/* _____________________________________________________________________ */
/* _____________________________________________________________________ */
/* _____________________________________________________________________ */
/* _____________________________________________________________________ */

int main(int argc, const char * argv[])
{
    assert(diameter > 0.0);
    assert(delta > 0.0);
    
    radius = 0.5 * diameter;
    LJ_length_scale_squared = diameter*diameter;
    LJ_length_cutoff2 = LJ_length_scale_squared*6.25;
    LJ_energy_cutoff = SqMS_LJ_potential(LJ_length_cutoff2, LJ_length_scale_squared, 1.);


    for (size_t i = 0; i < NDIM; i++)
    {
        Box[i] = 3.55;
    }
    

    if(NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
    else if(NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
    else{
        particle_volume = pow(M_PI,(double) NDIM/2.) * pow(radius, NDIM) / tgamma(1.+(double) NDIM/2.);
    }
    if (collision_detection == YES)
    {
        double box_volume = 1;
        for (size_t i = 0; i < NDIM; i++)
        {
            box_volume *= Box[i];
        }
        assert(box_volume > particle_volume*n_particles);
    }
    
    //write_data(-1); //writes initial config
    // read_data(); //reads initial config (just a proof of concept)
    
    //dsfmt_seed ((int) time (NULL) ) ;
    dsfmt_seed (314159) ;
    
    init_positions_hot();
    init_particle_energies();



    write_data(0); //writes initial config
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
            write_data(step);
        }
    }


    return 0;
}
