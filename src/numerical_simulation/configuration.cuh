
// Device data arrays in case of GPU compilation:
struct _GPU_arrays {
  int *d_idx ;
  double *d_x , *d_y , *d_z ;     // positions
  int *d_ix , *d_iy , *d_iz ;     // box images
  double *d_vx , *d_vy , *d_vz ;  // velocities
  double *d_fx , *d_fy , *d_fz ;  // forces
  double *d_compute_vec ;
};

// Global quantities on GPU
__constant__ int d_numberOfPart ;
__constant__ double d_box_sides[3] , d_midside[3] ;
__constant__ double d_time_step , d_half_time_step ;

// Cell lists on GPU
struct _GPU_cells {
  int numCells , *m_dCellStart , *m_dCellEnd ;
  int3 CellSize ;
  double3 CellSide ;
  int *m_dCellParticleHash , *m_dCellParticleIndex ;
  int *d_sorted_int ;
  double *d_sorted_double ;
};
__constant__ int d_CellSize[3] , d_CellHashPBC[3] ;
__constant__ double d_CellSide[3] ;
__constant__ int d_NneighCells , d_neighDeltaHashes[26] , d_neighMask_x[26] , d_neighMask_y[26] , d_neighMask_z[26] ;

struct _GPU_params {
  dim3 block , grid ;
  dim3 cell_block , cell_grid ;
};



__global__ void calcCellHash( int *m_dCellParticleHash , int *m_dCellParticleIndex ,
			      double *d_x , double *d_y , double *d_z ,
			      int *d_ix , int *d_iy , int *d_iz ) ;
__global__ void calcCellHash( int *m_dCellParticleHash , int *m_dCellParticleIndex ,
			      double *d_x , double *d_y ,
			      int *d_ix , int *d_iy ) ;
__global__ void findCellStartEnd( int *m_dCellParticleHash , int *m_dCellStart , int *m_dCellEnd ) ;
__global__ void reorderDeviceArrays( double *arr , double *d_sorted_double , int *m_dCellParticleIndex ) ;
__global__ void reorderDeviceArrays( int *arr , int *d_sorted_int , int *m_dCellParticleIndex ) ;
__global__ void VelocityVerletH1( double *d_x , double *d_y , double *d_z ,
				  double *d_vx , double *d_vy , double *d_vz ,
				  double *d_fx , double *d_fy , double *d_fz ) ;
__global__ void VelocityVerletH1( double *d_x , double *d_y ,
				  double *d_vx , double *d_vy ,
				  double *d_fx , double *d_fy ) ;
__global__ void VelocityVerletH2( double *d_vx , double *d_vy , double *d_vz ,
				  double *d_fx , double *d_fy , double *d_fz ) ;
__global__ void VelocityVerletH2( double *d_vx , double *d_vy ,
				  double *d_fx , double *d_fy ) ;
__global__ void zero_vectors( double *d_fx , double *d_fy , double *d_fz ) ;
__global__ void zero_vectors( double *d_fx , double *d_fy ) ;
__global__ void zero_vectors( double *d_fx ) ;
__global__ void compute_v2( double *d_vx , double *d_vy , double *d_vz , double *d_compute_vec ) ;
__global__ void compute_v2( double *d_vx , double *d_vy , double *d_compute_vec ) ;
__global__ void fill_vectors( double *d_x , double val ) ;
__global__ void fill_vectors( double *d_x , double *d_y , double val ) ;
__global__ void fill_vectors( double *d_x , double *d_y , double *d_z , double val ) ;
__global__ void fill_vectors( int *d_x , int val ) ;
__global__ void fill_vectors( int *d_x , int *d_y , int val ) ;
__global__ void fill_vectors( int *d_x , int *d_y , int *d_z , int val ) ;
__global__ void fill_vectors_range( double *d_x , double start ) ;
__global__ void fill_vectors_range( double *d_x , double *d_y , double start ) ;
__global__ void fill_vectors_range( double *d_x , double *d_y , double *d_z , double start ) ;
__global__ void fill_vectors_range( int *d_x , int start ) ;
__global__ void fill_vectors_range( int *d_x , int *d_y , int start ) ;
__global__ void fill_vectors_range( int *d_x , int *d_y , int *d_z , int start ) ;
