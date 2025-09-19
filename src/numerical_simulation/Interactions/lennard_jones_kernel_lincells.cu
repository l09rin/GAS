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
    // get address in grid
    int2 cellPos ;
    cellPos.x = floor( d_x[index] / d_CellSide[0] ) ;
    cellPos.y = floor( d_y[index] / d_CellSide[1] ) ;

    for (int y=-1; y<=1; y++) {
      for (int x=-1; x<=1; x++) {
	int cellHash = ( ( cellPos.y + y + d_CellSize[1] )%d_CellSize[1] ) * d_CellSize[0] + \
	  ( ( cellPos.x + x + d_CellSize[0] )%d_CellSize[0] ) ;
	int startIndex = cellStart[cellHash] ;
	if( startIndex != 0xffffffff ) {         // cell is not empty
	  for (int j=startIndex; j<cellEnd[cellHash]; j++) {
	    if (j != index) {               // check not colliding with self
	      d_add_pair_energy_lj( index , j , cutoff2 , shift , d_x , d_y , d_energy ) ;
	    }
	  }
	}
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
    // get address in grid
    int3 cellPos ;
    cellPos.x = floor( d_x[index] / d_CellSide[0] ) ;
    cellPos.y = floor( d_y[index] / d_CellSide[1] ) ;
    cellPos.z = floor( d_z[index] / d_CellSide[2] ) ;
    int cellHash0 = cellPos.z * d_CellSize[1] * d_CellSize[0]  +  cellPos.y * d_CellSize[0]  +  cellPos.x ;
    int startIndex = cellStart[cellHash0] ;
    if( startIndex != 0xffffffff ) {         // cell is not empty
      for( int j=startIndex; j<cellEnd[cellHash0]; j++ ) {
	if( j != index ) {               // check not colliding with self
	  d_add_pair_energy_lj( index , j , cutoff2 , shift , d_x , d_y , d_z , d_energy ) ;
	}
      }
    }
    // QUIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII gattibbaff

    int cellHash ;
    int cborder_x = (cellPos.x==(d_CellSize[0]-1)) - (cellPos.x==0) ;
    int cborder_y = (cellPos.y==(d_CellSize[1]-1)) - (cellPos.y==0) ;
    int cborder_z = (cellPos.z==(d_CellSize[2]-1)) - (cellPos.z==0) ;
    for( int cn=0; cn<d_NneighCells; cn++ ) {
      cellHash = cellHash0 + d_neighDeltaHashes[cn] - cborder_x * (cborder_x==d_neighMask_x[cn]) * d_CellHashPBC[0] \
	- cborder_y * (cborder_y==d_neighMask_y[cn]) * d_CellHashPBC[1] \
	- cborder_z * (cborder_z==d_neighMask_z[cn]) * d_CellHashPBC[2] ;
      startIndex = cellStart[cellHash] ;
      if( startIndex != 0xffffffff ) {         // cell is not empty
	for( int j=startIndex; j<cellEnd[cellHash]; j++ ) {
	  d_add_pair_energy_lj( index , j , cutoff2 , shift , d_x , d_y , d_z , d_energy ) ;
	}
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
    // get address in grid
    int2 cellPos ;
    cellPos.x = floor( d_x[index] / d_CellSide[0] ) ;
    cellPos.y = floor( d_y[index] / d_CellSide[1] ) ;

    for (int y=-1; y<=1; y++) {
      for (int x=-1; x<=1; x++) {
	int cellHash = ( ( cellPos.y + y + d_CellSize[1] )%d_CellSize[1] ) * d_CellSize[0] + \
	  ( ( cellPos.x + x + d_CellSize[0] )%d_CellSize[0] ) ;
	int startIndex = cellStart[cellHash] ;
	if( startIndex != 0xffffffff ) {         // cell is not empty
	  for (int j=startIndex; j<cellEnd[cellHash]; j++) {
	    if (j != index) {               // check not colliding with self
	      d_add_pair_force_lj( index , j , cutoff2 , d_x , d_y , d_fx , d_fy ) ;
	    }
	  }
	}
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
    // get address in grid
    int3 cellPos ;
    cellPos.x = floor( d_x[index] / d_CellSide[0] ) ;
    cellPos.y = floor( d_y[index] / d_CellSide[1] ) ;
    cellPos.z = floor( d_z[index] / d_CellSide[2] ) ;

    for (int z=-1; z<=1; z++) {
      for (int y=-1; y<=1; y++) {
	for (int x=-1; x<=1; x++) {
	  int cellHash = ( ( cellPos.z + z + d_CellSize[2] )%d_CellSize[2] ) * d_CellSize[1] * d_CellSize[0] + \
	    ( ( cellPos.y + y + d_CellSize[1] )%d_CellSize[1] ) * d_CellSize[0] + \
	    ( ( cellPos.x + x + d_CellSize[0] )%d_CellSize[0] ) ;
	  int startIndex = cellStart[cellHash] ;
	  if( startIndex != 0xffffffff ) {         // cell is not empty
	    for (int j=startIndex; j<cellEnd[cellHash]; j++) {
	      if (j != index) {               // check not colliding with self
                d_add_pair_force_lj( index , j , cutoff2 , d_x , d_y , d_z , d_fx , d_fy , d_fz ) ;
	      }
	    }
	  }
	}
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
    // get address in grid
    int3 cellPos ;
    cellPos.x = floor( d_x[index] / d_CellSide[0] ) ;
    cellPos.y = floor( d_y[index] / d_CellSide[1] ) ;
    cellPos.z = floor( d_z[index] / d_CellSide[2] ) ;

    for (int z=-1; z<=1; z++) {
      for (int y=-1; y<=1; y++) {
	for (int x=-1; x<=1; x++) {
	  int cellHash = ( ( cellPos.z + z + d_CellSize[2] )%d_CellSize[2] ) * d_CellSize[1] * d_CellSize[0] + \
	    ( ( cellPos.y + y + d_CellSize[1] )%d_CellSize[1] ) * d_CellSize[0] + \
	    ( ( cellPos.x + x + d_CellSize[0] )%d_CellSize[0] ) ;
	  int startIndex = cellStart[cellHash] ;
	  if( startIndex != 0xffffffff ) {         // cell is not empty
	    for (int j=startIndex; j<cellEnd[cellHash]; j++) {
	      if (j != index) {               // check not colliding with self
                d_add_pair_virial_lj( index , j , cutoff2 , d_x , d_y , d_z , d_virial ) ;
	      }
	    }
	  }
	}
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
    // get address in grid
    int2 cellPos ;
    cellPos.x = floor( d_x[index] / d_CellSide[0] ) ;
    cellPos.y = floor( d_y[index] / d_CellSide[1] ) ;

    for (int y=-1; y<=1; y++) {
      for (int x=-1; x<=1; x++) {
	int cellHash = ( ( cellPos.y + y + d_CellSize[1] )%d_CellSize[1] ) * d_CellSize[0] + \
	  ( ( cellPos.x + x + d_CellSize[0] )%d_CellSize[0] ) ;
	int startIndex = cellStart[cellHash] ;
	if( startIndex != 0xffffffff ) {         // cell is not empty
	  for (int j=startIndex; j<cellEnd[cellHash]; j++) {
	    if (j != index) {               // check not colliding with self
	      d_add_pair_virial_lj( index , j , cutoff2 , d_x , d_y , d_virial ) ;
	    }
	  }
	}
      }
    }
  }
}
/********************************************************************************************/
