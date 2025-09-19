#include <cooperative_groups.h>
namespace cg = cooperative_groups ;

__global__ void calcCellHash( int *m_dCellParticleHash , int *m_dCellParticleIndex ,
			      double *d_x , double *d_y ,
			      int *d_ix , int *d_iy ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  int3 cell_idx ;

  if( i < d_numberOfPart ) {
    if( d_x[i] >= d_box_sides[0] ){
      d_x[i] -= d_box_sides[0] ;
      d_ix[i] ++ ;
    } else if( d_x[i] < 0 ){
      d_x[i] += d_box_sides[0] ;
      d_ix[i] -- ;
    }
    cell_idx.x = floor( d_x[i] / d_CellSide[0] ) ;

    if( d_y[i] >= d_box_sides[1] ){
      d_y[i] -= d_box_sides[1] ;
      d_iy[i] ++ ;
    } else if( d_y[i] < 0 ){
      d_y[i] += d_box_sides[1] ;
      d_iy[i] -- ;
    }
    cell_idx.y = floor( d_y[i] / d_CellSide[1] ) ;

    m_dCellParticleHash[i] = cell_idx.y * d_CellSize[0] + cell_idx.x ;
    m_dCellParticleIndex[i] = i ;
  }
}
/****************************************************************************************/

__global__ void calcCellHash( int *m_dCellParticleHash , int *m_dCellParticleIndex ,
			      double *d_x , double *d_y , double *d_z ,
			      int *d_ix , int *d_iy , int *d_iz ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  int3 cell_idx ;

  if( i < d_numberOfPart ) {
    if( d_x[i] >= d_box_sides[0] ){
      d_x[i] -= d_box_sides[0] ;
      d_ix[i] ++ ;
    } else if( d_x[i] < 0 ){
      d_x[i] += d_box_sides[0] ;
      d_ix[i] -- ;
    }
    cell_idx.x = floor( d_x[i] / d_CellSide[0] ) ;

    if( d_y[i] >= d_box_sides[1] ){
      d_y[i] -= d_box_sides[1] ;
      d_iy[i] ++ ;
    } else if( d_y[i] < 0 ){
      d_y[i] += d_box_sides[1] ;
      d_iy[i] -- ;
    }
    cell_idx.y = floor( d_y[i] / d_CellSide[1] ) ;

    if( d_z[i] >= d_box_sides[2] ){
      d_z[i] -= d_box_sides[2] ;
      d_iz[i] ++ ;
    } else if( d_z[i] < 0 ){
      d_z[i] += d_box_sides[2] ;
      d_iz[i] -- ;
    }
    cell_idx.z = floor( d_z[i] / d_CellSide[2] ) ;

    m_dCellParticleHash[i] = cell_idx.z * d_CellSize[1] * d_CellSize[0]  + cell_idx.y * d_CellSize[0] + cell_idx.x ;
    m_dCellParticleIndex[i] = i ;
  }
}
/****************************************************************************************/

// find the start of each cell
// in the sorted hash array
__global__ void findCellStartEnd( int *m_dCellParticleHash , int *m_dCellStart , int *m_dCellEnd ) {
  // Handle to thread block group
  cg::thread_block cta = cg::this_thread_block();
  extern __shared__ int sharedHash[] ;    // blockSize + 1 elements
  int i = blockIdx.x * blockDim.x + threadIdx.x ;

  int hash ;
  // handle case when no. of particles not multiple of block size
  if( i < d_numberOfPart ) {
    hash = m_dCellParticleHash[i];
    // Load hash data into shared memory so that we can look
    // at neighboring particle's hash value without loading
    // two hash values per thread
    sharedHash[threadIdx.x+1] = hash;

    if( i > 0 && threadIdx.x == 0 ) {
      // first thread in block must load neighbor particle hash
      sharedHash[0] = m_dCellParticleHash[i-1];
    }
  }
  cg::sync(cta);

  if( i < d_numberOfPart ) {
    // If this particle has a different cell i to the previous
    // particle then it must be the first particle in the cell,
    // so store the i of this particle in the cell.
    // As it isn't the first particle, it must also be the cell end of
    // the previous particle's cell
    if( i == 0 || hash != sharedHash[threadIdx.x] ) {
      m_dCellStart[hash] = i ;
      if( i > 0 ) m_dCellEnd[sharedHash[threadIdx.x]] = i ;
    }

    if( i == d_numberOfPart - 1 )  m_dCellEnd[hash] = i + 1;
  }
}
/****************************************************************************************/

// reorder the arrays according to the reordered hashes
__global__ void reorderDeviceArrays( double *arr , double *d_sorted_double , int *m_dCellParticleIndex ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;

  if( i < d_numberOfPart ) {
    int sidx = m_dCellParticleIndex[i] ;
    d_sorted_double[i] = arr[sidx] ;
  }
}
/****************************************************************************************/

// reorder the arrays according to the reordered hashes
__global__ void reorderDeviceArrays( int *arr , int *d_sorted_int , int *m_dCellParticleIndex ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;

  if( i < d_numberOfPart ) {
    int sidx = m_dCellParticleIndex[i] ;
    d_sorted_int[i] = arr[sidx] ;
  }
}
/****************************************************************************************/

__global__ void zero_vectors( double *d_fx ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  if( i < d_numberOfPart )  d_fx[i] = 0.0 ;
}
/****************************************************************************************/

__global__ void zero_vectors( double *d_fx , double *d_fy ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  if( i < d_numberOfPart ) {
    d_fx[i] = 0.0 ;
    d_fy[i] = 0.0 ;
  }
}
/****************************************************************************************/

__global__ void zero_vectors( double *d_fx , double *d_fy , double *d_fz ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  if( i < d_numberOfPart ) {
    d_fx[i] = 0.0 ;
    d_fy[i] = 0.0 ;
    d_fz[i] = 0.0 ;
  }
}
/****************************************************************************************/

__global__ void VelocityVerletH1( double *d_x , double *d_y ,
				  double *d_vx , double *d_vy ,
				  double *d_fx , double *d_fy ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;

  if( i < d_numberOfPart ) {
    d_vx[i] += ( d_fx[i] * d_half_time_step ) ;
    d_vy[i] += ( d_fy[i] * d_half_time_step ) ;
    // positions update
    d_x[i] += ( d_vx[i] * d_time_step ) ;
    d_y[i] += ( d_vy[i] * d_time_step ) ;
  }
}
/****************************************************************************************/

__global__ void VelocityVerletH1( double *d_x , double *d_y , double *d_z ,
				  double *d_vx , double *d_vy , double *d_vz ,
				  double *d_fx , double *d_fy , double *d_fz ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;

  if( i < d_numberOfPart ) {
    d_vx[i] += ( d_fx[i] * d_half_time_step ) ;
    d_vy[i] += ( d_fy[i] * d_half_time_step ) ;
    d_vz[i] += ( d_fz[i] * d_half_time_step ) ;
    // positions update
    d_x[i] += ( d_vx[i] * d_time_step ) ;
    d_y[i] += ( d_vy[i] * d_time_step ) ;
    d_z[i] += ( d_vz[i] * d_time_step ) ;
  }
}
/****************************************************************************************/

__global__ void VelocityVerletH2( double *d_vx , double *d_vy ,
				  double *d_fx , double *d_fy ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;

  if( i < d_numberOfPart ) {
    d_vx[i] += ( d_fx[i] * d_half_time_step ) ;
    d_vy[i] += ( d_fy[i] * d_half_time_step ) ;
  }
}
/****************************************************************************************/

__global__ void VelocityVerletH2( double *d_vx , double *d_vy , double *d_vz ,
				  double *d_fx , double *d_fy , double *d_fz ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;

  if( i < d_numberOfPart ) {
    d_vx[i] += ( d_fx[i] * d_half_time_step ) ;
    d_vy[i] += ( d_fy[i] * d_half_time_step ) ;
    d_vz[i] += ( d_fz[i] * d_half_time_step ) ;
  }
}
/****************************************************************************************/

__global__ void compute_v2( double *d_vx , double *d_vy , double *d_vz , double *d_compute_vec ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  if( i < d_numberOfPart )  d_compute_vec[i] = d_vx[i] * d_vx[i] + d_vy[i] * d_vy[i] + d_vz[i] * d_vz[i] ;
}
/****************************************************************************************/

__global__ void compute_v2( double *d_vx , double *d_vy , double *d_compute_vec ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  if( i < d_numberOfPart )  d_compute_vec[i] = d_vx[i] * d_vx[i] + d_vy[i] * d_vy[i] ;
}
/****************************************************************************************/

__global__ void fill_vectors( double *d_x , double val ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  if( i < d_numberOfPart )  d_x[i] = val ;
}
/****************************************************************************************/

__global__ void fill_vectors( double *d_x , double *d_y , double val ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  if( i < d_numberOfPart ) {
    d_x[i] = val ;
    d_y[i] = val ;
  }
}
/****************************************************************************************/

__global__ void fill_vectors( double *d_x , double *d_y , double *d_z , double val ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  if( i < d_numberOfPart ) {
    d_x[i] = val ;
    d_y[i] = val ;
    d_z[i] = val ;
  }
}
/****************************************************************************************/

__global__ void fill_vectors( int *d_x , int val ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  if( i < d_numberOfPart )  d_x[i] = val ;
}
/****************************************************************************************/

__global__ void fill_vectors( int *d_x , int *d_y , int val ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  if( i < d_numberOfPart ) {
    d_x[i] = val ;
    d_y[i] = val ;
  }
}
/****************************************************************************************/

__global__ void fill_vectors( int *d_x , int *d_y , int *d_z , int val ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  if( i < d_numberOfPart ) {
    d_x[i] = val ;
    d_y[i] = val ;
    d_z[i] = val ;
  }
}
/****************************************************************************************/

__global__ void fill_vectors_range( double *d_x , double start ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  if( i < d_numberOfPart )  d_x[i] = start + i ;
}
/****************************************************************************************/

__global__ void fill_vectors_range( double *d_x , double *d_y , double start ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  if( i < d_numberOfPart ) {
    d_x[i] = start + i ;
    d_y[i] = start + i ;
  }
}
/****************************************************************************************/

__global__ void fill_vectors_range( double *d_x , double *d_y , double *d_z , double start ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  if( i < d_numberOfPart ) {
    d_x[i] = start + i ;
    d_y[i] = start + i ;
    d_z[i] = start + i ;
  }
}
/****************************************************************************************/

__global__ void fill_vectors_range( int *d_x , int start ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
 if( i < d_numberOfPart )  d_x[i] = start + i ;
}
/****************************************************************************************/

__global__ void fill_vectors_range( int *d_x , int *d_y , int start ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  if( i < d_numberOfPart ) {
    d_x[i] = start + i ;
    d_y[i] = start + i ;
  }
}
/****************************************************************************************/

__global__ void fill_vectors_range( int *d_x , int *d_y , int *d_z , int start ) {
  int i = blockIdx.x * blockDim.x + threadIdx.x ;
  if( i < d_numberOfPart ) {
    d_x[i] = start + i ;
    d_y[i] = start + i ;
    d_z[i] = start + i ;
  }
}
/****************************************************************************************/
