#ifndef LENNARD_JONES_CUDA
#define LENNARD_JONES_CUDA

__device__ inline void d_add_pair_force_lj( int i , int j , double cutoff2 ,
					    double *d_x , double *d_y , double *d_z ,
					    double *d_fx , double *d_fy , double *d_fz ) ;
__device__ inline void d_add_pair_force_lj( int i , int j , double cutoff2 ,
					    double *d_x , double *d_y ,
					    double *d_fx , double *d_fy ) ;
__device__ inline void d_add_pair_virial_lj( int i , int j , double cutoff2 ,
					     double *d_x , double *d_y , double *d_z , double *virial ) ;
__device__ inline void d_add_pair_virial_lj( int i , int j , double cutoff2 ,
					     double *d_x , double *d_y , double *virial ) ;
__device__ inline void d_add_pair_energy_lj( int i , int j , double cutoff2 , double shift ,
					     double *d_x , double *d_y , double *d_z , double *energy ) ;
__device__ inline void d_add_pair_energy_lj( int i , int j , double cutoff2 , double shift ,
					     double *d_x , double *d_y , double *energy ) ;

__global__ void d_compute_forces_lj( double *d_x , double *d_y , double *d_z ,
				     double *d_fx , double *d_fy , double *d_fz ,
				     int *cellStart , int *cellEnd , double cutoff2 ) ;
__global__ void d_compute_forces_lj( double *d_x , double *d_y ,
				     double *d_fx , double *d_fy ,
				     int *cellStart , int *cellEnd , double cutoff2 ) ;
__global__ void d_compute_energy_lj( double *d_x , double *d_y , double *d_z ,
				     double *d_energy , int *cellStart , int *cellEnd , double cutoff2 , double shift ) ;
__global__ void d_compute_energy_lj( double *d_x , double *d_y ,
				     double *d_energy , int *cellStart , int *cellEnd , double cutoff2 , double shift ) ;
__global__ void d_compute_virial_lj( double *d_x , double *d_y , double *d_z ,
				     double *d_virial , int *cellStart , int *cellEnd , double cutoff2 ) ;
__global__ void d_compute_virial_lj( double *d_x , double *d_y ,
				     double *d_virial , int *cellStart , int *cellEnd , double cutoff2 ) ;


#include "lennard_jones_kernel.cu"

#endif
