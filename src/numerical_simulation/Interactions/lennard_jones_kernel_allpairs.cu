/********************************************************************************************/
#include "assert.h"

__device__ inline void d_add_pair_force_lj( int i , int j , double cutoff2 ,
					    double *d_x , double *d_y , double *d_z ,
					    double *d_fx , double *d_fy , double *d_fz ) {
  double3 dr ;
  dr.x = d_x[i] - d_x[j] ;
  dr.y = d_y[i] - d_y[j] ;
  dr.z = d_z[i] - d_z[j] ;
  if( dr.x > d_midside[0] ) dr.x -= d_box_sides[0] ;
  else if( dr.x < -d_midside[0] ) dr.x += d_box_sides[0] ;
  if( dr.y > d_midside[1] ) dr.y -= d_box_sides[1] ;
  else if( dr.y < -d_midside[1] ) dr.y += d_box_sides[1] ;
  if( dr.z > d_midside[2] ) dr.z -= d_box_sides[2] ;
  else if( dr.z < -d_midside[2] ) dr.z += d_box_sides[2] ;

  double r2_inv = ( dr.x*dr.x + dr.y*dr.y + dr.z*dr.z ) ;
  if( r2_inv <= cutoff2 ) {
    r2_inv = 1.0 / r2_inv ;
    double multiplier = r2_inv * r2_inv * r2_inv ;
    multiplier = ( 48.0 * multiplier * multiplier - 24.0 * multiplier ) * r2_inv ;
    d_fx[i] += ( dr.x * multiplier ) ;
    d_fy[i] += ( dr.y * multiplier ) ;
    d_fz[i] += ( dr.z * multiplier ) ;
  }
};
/********************************************************************************************/


__device__ inline void d_add_pair_force_lj( int i , int j , double cutoff2 ,
					    double *d_x , double *d_y ,
					    double *d_fx , double *d_fy ) {
  double3 dr ;
  dr.x = d_x[i] - d_x[j] ;
  dr.y = d_y[i] - d_y[j] ;
  if( dr.x > d_midside[0] ) dr.x -= d_box_sides[0] ;
  else if( dr.x < -d_midside[0] ) dr.x += d_box_sides[0] ;
  if( dr.y > d_midside[1] ) dr.y -= d_box_sides[1] ;
  else if( dr.y < -d_midside[1] ) dr.y += d_box_sides[1] ;

  double r2_inv = ( dr.x*dr.x + dr.y*dr.y ) ;
  if( r2_inv <= cutoff2 ) {
    r2_inv = 1.0 / r2_inv ;
    double multiplier = r2_inv * r2_inv * r2_inv ;
    multiplier = ( 48.0 * multiplier * multiplier - 24.0 * multiplier ) * r2_inv ;
    d_fx[i] += dr.x * multiplier ;
    d_fy[i] += dr.y * multiplier ;
  }
};
/********************************************************************************************/


__device__ inline void d_add_pair_virial_lj( int i , int j , double cutoff2 ,
					     double *d_x , double *d_y , double *d_z , double *virial ) {
  double3 dr ;
  dr.x = d_x[i] - d_x[j] ;
  dr.y = d_y[i] - d_y[j] ;
  dr.z = d_z[i] - d_z[j] ;
  if( dr.x > d_midside[0] ) dr.x -= d_box_sides[0] ;
  else if( dr.x < -d_midside[0] ) dr.x += d_box_sides[0] ;
  if( dr.y > d_midside[1] ) dr.y -= d_box_sides[1] ;
  else if( dr.y < -d_midside[1] ) dr.y += d_box_sides[1] ;
  if( dr.z > d_midside[2] ) dr.z -= d_box_sides[2] ;
  else if( dr.z < -d_midside[2] ) dr.z += d_box_sides[2] ;

  double r2_inv = ( dr.x*dr.x + dr.y*dr.y + dr.z*dr.z ) ;
  if( r2_inv <= cutoff2 ) {
    r2_inv = 1.0 / r2_inv ;
    double multiplier = r2_inv * r2_inv * r2_inv ;
    multiplier = ( 48.0 * multiplier * multiplier - 24.0 * multiplier ) * r2_inv ;
    virial[i] += ( dr.x * dr.x + dr.y * dr.y + dr.z * dr.z ) * multiplier ;
  }
};
/********************************************************************************************/


__device__ inline void d_add_pair_virial_lj( int i , int j , double cutoff2 ,
					     double *d_x , double *d_y , double *virial ) {
  double3 dr ;
  dr.x = d_x[i] - d_x[j] ;
  dr.y = d_y[i] - d_y[j] ;
  if( dr.x > d_midside[0] ) dr.x -= d_box_sides[0] ;
  else if( dr.x < -d_midside[0] ) dr.x += d_box_sides[0] ;
  if( dr.y > d_midside[1] ) dr.y -= d_box_sides[1] ;
  else if( dr.y < -d_midside[1] ) dr.y += d_box_sides[1] ;

  double r2_inv = ( dr.x*dr.x + dr.y*dr.y ) ;
  if( r2_inv <= cutoff2 ) {
    r2_inv = 1.0 / r2_inv ;
    double multiplier = r2_inv * r2_inv * r2_inv ;
    multiplier = ( 48.0 * multiplier * multiplier - 24.0 * multiplier ) * r2_inv ;
    virial[i] += ( dr.x * dr.x + dr.y * dr.y ) * multiplier ;
  }
};
/********************************************************************************************/


__device__ inline void d_add_pair_energy_lj( int i , int j , double cutoff2 , double shift ,
					     double *d_x , double *d_y , double *d_z , double *energy ) {
  double3 dr ;
  dr.x = d_x[i] - d_x[j] ;
  dr.y = d_y[i] - d_y[j] ;
  dr.z = d_z[i] - d_z[j] ;
  if( dr.x > d_midside[0] ) dr.x -= d_box_sides[0] ;
  else if( dr.x < -d_midside[0] ) dr.x += d_box_sides[0] ;
  if( dr.y > d_midside[1] ) dr.y -= d_box_sides[1] ;
  else if( dr.y < -d_midside[1] ) dr.y += d_box_sides[1] ;
  if( dr.z > d_midside[2] ) dr.z -= d_box_sides[2] ;
  else if( dr.z < -d_midside[2] ) dr.z += d_box_sides[2] ;

  double r2_inv = ( dr.x*dr.x + dr.y*dr.y + dr.z*dr.z ) ;
  if( r2_inv <= cutoff2 ) {
    r2_inv = 1.0 / r2_inv ;
    double r6_inv = r2_inv * r2_inv * r2_inv ;
    energy[i] += ( 4.0 * ( r6_inv * r6_inv - r6_inv ) + shift ) ;
  }
};
/********************************************************************************************/


__device__ inline void d_add_pair_energy_lj( int i , int j , double cutoff2 , double shift ,
					     double *d_x , double *d_y , double *energy ) {
  double3 dr ;
  dr.x = d_x[i] - d_x[j] ;
  dr.y = d_y[i] - d_y[j] ;
  if( dr.x > d_midside[0] ) dr.x -= d_box_sides[0] ;
  else if( dr.x < -d_midside[0] ) dr.x += d_box_sides[0] ;
  if( dr.y > d_midside[1] ) dr.y -= d_box_sides[1] ;
  else if( dr.y < -d_midside[1] ) dr.y += d_box_sides[1] ;

  double r2_inv = ( dr.x*dr.x + dr.y*dr.y ) ;
  if( r2_inv <= cutoff2 ) {
    r2_inv = 1.0 / r2_inv ;
    double r6_inv = r2_inv * r2_inv * r2_inv ;
    energy[i] += ( 4.0 * ( r6_inv * r6_inv - r6_inv ) + shift ) ;
  }
};
/********************************************************************************************/


__global__ void d_compute_energy_lj( double *d_x , double *d_y ,
				     double *d_energy ,
				     int *cellStart , int *cellEnd , double cutoff2 , double shift ) {
  int index = blockIdx.x * blockDim.x + threadIdx.x ;

  if ( index < d_numberOfPart ) {
    for (int j=0; j<d_numberOfPart; j++) {
      if (j != index) {               // check not colliding with self
	d_add_pair_energy_lj( index , j , cutoff2 , shift , d_x , d_y , d_energy ) ;
      }
    }
  }
}
/********************************************************************************************/

__global__ void d_compute_energy_lj( double *d_x , double *d_y , double *d_z ,
				     double *d_energy ,
				     int *cellStart , int *cellEnd , double cutoff2 , double shift ) {
  int index = blockIdx.x * blockDim.x + threadIdx.x ;

  if ( index < d_numberOfPart ) {
    for (int j=0; j<d_numberOfPart; j++) {
      if (j != index) {               // check not colliding with self
	d_add_pair_energy_lj( index , j , cutoff2 , shift , d_x , d_y , d_z , d_energy ) ;
      }
    }
  }
}
/********************************************************************************************/


__global__ void d_compute_forces_lj( double *d_x , double *d_y ,
				     double *d_fx , double *d_fy ,
				     int *cellStart , int *cellEnd , double cutoff2 ) {
  int index = blockIdx.x * blockDim.x + threadIdx.x ;

  if ( index < d_numberOfPart ) {
    for (int j=0; j<d_numberOfPart; j++) {
      if (j != index) {               // check not colliding with self
	d_add_pair_force_lj( index , j , cutoff2 , d_x , d_y , d_fx , d_fy ) ;
      }
    }
  }
}
/********************************************************************************************/

__global__ void d_compute_forces_lj( double *d_x , double *d_y , double *d_z ,
				     double *d_fx , double *d_fy , double *d_fz ,
				     int *cellStart , int *cellEnd , double cutoff2 ) {
  int index = blockIdx.x * blockDim.x + threadIdx.x ;

  if ( index < d_numberOfPart ) {
    for (int j=0; j<d_numberOfPart; j++) {
      if (j != index) {               // check not colliding with self
	d_add_pair_force_lj( index , j , cutoff2 , d_x , d_y , d_z , d_fx , d_fy , d_fz ) ;
      }
    }
  }
}
/********************************************************************************************/


__global__ void d_compute_virial_lj( double *d_x , double *d_y , double *d_z ,
				     double *d_virial ,
				     int *cellStart , int *cellEnd , double cutoff2 ) {
  int index = blockIdx.x * blockDim.x + threadIdx.x ;

  if ( index < d_numberOfPart ) {
    for (int j=0; j<d_numberOfPart; j++) {
      if (j != index) {               // check not colliding with self
	d_add_pair_virial_lj( index , j , cutoff2 , d_x , d_y , d_z , d_virial ) ;
      }
    }
  }
}
/********************************************************************************************/

__global__ void d_compute_virial_lj( double *d_x , double *d_y ,
				     double *d_virial ,
				     int *cellStart , int *cellEnd , double cutoff2 ) {
  int index = blockIdx.x * blockDim.x + threadIdx.x ;

  if ( index < d_numberOfPart ) {
    for (int j=0; j<d_numberOfPart; j++) {
      if (j != index) {               // check not colliding with self
	d_add_pair_virial_lj( index , j , cutoff2 , d_x , d_y , d_virial ) ;
      }
    }
  }
}
/********************************************************************************************/
