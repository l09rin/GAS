/****************************************************************************************/
#include <vector>

// GPU initialization
template <typename particle>
void configuration<particle>::allocateDeviceBox( void ) {
  double box[3] ;
  box[0] = box_sides.position.x ;
  box[1] = box_sides.position.y ;
  box[2] = 0.0 ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_box_sides , box , sizeof(double) * 3 )) ;
  box[0] *= 0.5 ;
  box[1] *= 0.5 ;
  box[2] *= 0.5 ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_midside , box , sizeof(double) * 3 )) ;
}

template <>
void configuration<particle_3D>::allocateDeviceBox( void ) {
  double box[3] ;
  box[0] = box_sides.position.x ;
  box[1] = box_sides.position.y ;
  box[2] = box_sides.position.z ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_box_sides , box , sizeof(double) * 3 )) ;
  box[0] *= 0.5 ;
  box[1] *= 0.5 ;
  box[2] *= 0.5 ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_midside , box , sizeof(double) * 3 )) ;
}
/****************************************************************************************/
/****************************************************************************************/

template <typename particle>
void configuration<particle>::allocateDeviceNeighHashes( void ) {
  int NneighCells = 8 ;
  int neighDeltaHashes[8] , neighMask_x[8] , neighMask_y[8] , neighMask_z[8] ;

  neighDeltaHashes[0] = -1 - GPU_cells.CellSize.x ; // [ -1  -1 ]
  neighDeltaHashes[1] = -1 ;                        // [ -1   0 ]
  neighDeltaHashes[2] = -1 + GPU_cells.CellSize.x ; // [ -1  +1 ]
  neighDeltaHashes[3] = -1 * GPU_cells.CellSize.x ; // [  0  -1 ]
  neighDeltaHashes[4] =      GPU_cells.CellSize.x ; // [  0  +1 ]
  neighDeltaHashes[5] = 1 - GPU_cells.CellSize.x ;  // [ +1  -1 ]
  neighDeltaHashes[6] = 1 ;                         // [ +1   0 ]
  neighDeltaHashes[7] = 1 + GPU_cells.CellSize.x ;  // [ +1  +1 ]

  neighMask_x[0] = -1 , neighMask_y[0] = -1 , neighMask_z[0] = 0 ;
  neighMask_x[1] = -1 , neighMask_y[1] =  0 , neighMask_z[1] = 0 ;
  neighMask_x[2] = -1 , neighMask_y[2] =  1 , neighMask_z[2] = 0 ;
  neighMask_x[3] =  0 , neighMask_y[3] = -1 , neighMask_z[3] = 0 ;
  neighMask_x[4] =  0 , neighMask_y[4] =  1 , neighMask_z[4] = 0 ;
  neighMask_x[5] =  1 , neighMask_y[5] = -1 , neighMask_z[5] = 0 ;
  neighMask_x[6] =  1 , neighMask_y[6] =  0 , neighMask_z[6] = 0 ;
  neighMask_x[7] =  1 , neighMask_y[7] =  1 , neighMask_z[7] = 0 ;

  CUDA_CHECK(cudaMemcpyToSymbol( d_NneighCells , &NneighCells , sizeof(int) )) ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_neighDeltaHashes , neighDeltaHashes , NneighCells*sizeof(int) )) ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_neighMask_x , neighMask_x , NneighCells*sizeof(int) )) ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_neighMask_y , neighMask_y , NneighCells*sizeof(int) )) ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_neighMask_z , neighMask_z , NneighCells*sizeof(int) )) ;
}

template <>
void configuration<particle_3D>::allocateDeviceNeighHashes( void ) {
  int NneighCells = 26 ;
  int neighDeltaHashes[26] , neighMask_x[26] , neighMask_y[26] , neighMask_z[26] ;

  neighDeltaHashes[0] =  -1 - GPU_cells.CellSize.x - GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [ -1 -1 -1 ]
  neighDeltaHashes[1] =  -1 - GPU_cells.CellSize.x                                               ;  // [ -1 -1  0 ]
  neighDeltaHashes[2] =  -1 - GPU_cells.CellSize.x + GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [ -1 -1  1 ]
  neighDeltaHashes[3] =  -1                        - GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [ -1  0 -1 ]
  neighDeltaHashes[4] =  -1                                                                      ;  // [ -1  0  0 ]
  neighDeltaHashes[5] =  -1                        + GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [ -1  0  1 ]
  neighDeltaHashes[6] =  -1 + GPU_cells.CellSize.x - GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [ -1  1 -1 ]
  neighDeltaHashes[7] =  -1 + GPU_cells.CellSize.x                                               ;  // [ -1  1  0 ]
  neighDeltaHashes[8] =  -1 + GPU_cells.CellSize.x + GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [ -1  1  1 ]
  neighDeltaHashes[9] =   0 - GPU_cells.CellSize.x - GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [  0 -1 -1 ]
  neighDeltaHashes[10] =  0 - GPU_cells.CellSize.x                                               ;  // [  0 -1  0 ]
  neighDeltaHashes[11] =  0 - GPU_cells.CellSize.x + GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [  0 -1  1 ]
  neighDeltaHashes[12] =  0                        - GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [  0  0 -1 ]
  neighDeltaHashes[13] =  0                        + GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [  0  0  1 ]
  neighDeltaHashes[14] =  0 + GPU_cells.CellSize.x - GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [  0  1 -1 ]
  neighDeltaHashes[15] =  0 + GPU_cells.CellSize.x                                               ;  // [  0  1  0 ]
  neighDeltaHashes[16] =  0 + GPU_cells.CellSize.x + GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [  0  1  1 ]
  neighDeltaHashes[17] =  1 - GPU_cells.CellSize.x - GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [  1 -1 -1 ]
  neighDeltaHashes[18] =  1 - GPU_cells.CellSize.x                                               ;  // [  1 -1  0 ]
  neighDeltaHashes[19] =  1 - GPU_cells.CellSize.x + GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [  1 -1  1 ]
  neighDeltaHashes[20] =  1                        - GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [  1  0 -1 ]
  neighDeltaHashes[21] =  1                                                                      ;  // [  1  0  0 ]
  neighDeltaHashes[22] =  1                        + GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [  1  0  1 ]
  neighDeltaHashes[23] =  1 + GPU_cells.CellSize.x - GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [  1  1 -1 ]
  neighDeltaHashes[24] =  1 + GPU_cells.CellSize.x                                               ;  // [  1  1  0 ]
  neighDeltaHashes[25] =  1 + GPU_cells.CellSize.x + GPU_cells.CellSize.x * GPU_cells.CellSize.y ;  // [  1  1  1 ]

  neighMask_x[0] =  -1 , neighMask_y[0] =  -1 , neighMask_z[0] =  -1 ;
  neighMask_x[1] =  -1 , neighMask_y[1] =  -1 , neighMask_z[1] =   0 ;
  neighMask_x[2] =  -1 , neighMask_y[2] =  -1 , neighMask_z[2] =   1 ;
  neighMask_x[3] =  -1 , neighMask_y[3] =   0 , neighMask_z[3] =  -1 ;
  neighMask_x[4] =  -1 , neighMask_y[4] =   0 , neighMask_z[4] =   0 ;
  neighMask_x[5] =  -1 , neighMask_y[5] =   0 , neighMask_z[5] =   1 ;
  neighMask_x[6] =  -1 , neighMask_y[6] =   1 , neighMask_z[6] =  -1 ;
  neighMask_x[7] =  -1 , neighMask_y[7] =   1 , neighMask_z[7] =   0 ;
  neighMask_x[8] =  -1 , neighMask_y[8] =   1 , neighMask_z[8] =   1 ;
  neighMask_x[9] =   0 , neighMask_y[9] =  -1 , neighMask_z[9] =  -1 ;
  neighMask_x[10] =  0 , neighMask_y[10] = -1 , neighMask_z[10] =  0 ;
  neighMask_x[11] =  0 , neighMask_y[11] = -1 , neighMask_z[11] =  1 ;
  neighMask_x[12] =  0 , neighMask_y[12] =  0 , neighMask_z[12] = -1 ;
  neighMask_x[13] =  0 , neighMask_y[13] =  0 , neighMask_z[13] =  1 ;
  neighMask_x[14] =  0 , neighMask_y[14] =  1 , neighMask_z[14] = -1 ;
  neighMask_x[15] =  0 , neighMask_y[15] =  1 , neighMask_z[15] =  0 ;
  neighMask_x[16] =  0 , neighMask_y[16] =  1 , neighMask_z[16] =  1 ;
  neighMask_x[17] =  1 , neighMask_y[17] = -1 , neighMask_z[17] = -1 ;
  neighMask_x[18] =  1 , neighMask_y[18] = -1 , neighMask_z[18] =  0 ;
  neighMask_x[19] =  1 , neighMask_y[19] = -1 , neighMask_z[19] =  1 ;
  neighMask_x[20] =  1 , neighMask_y[20] =  0 , neighMask_z[20] = -1 ;
  neighMask_x[21] =  1 , neighMask_y[21] =  0 , neighMask_z[21] =  0 ;
  neighMask_x[22] =  1 , neighMask_y[22] =  0 , neighMask_z[22] =  1 ;
  neighMask_x[23] =  1 , neighMask_y[23] =  1 , neighMask_z[23] = -1 ;
  neighMask_x[24] =  1 , neighMask_y[24] =  1 , neighMask_z[24] =  0 ;
  neighMask_x[25] =  1 , neighMask_y[25] =  1 , neighMask_z[25] =  1 ;

  CUDA_CHECK(cudaMemcpyToSymbol( d_NneighCells , &NneighCells , sizeof(int) )) ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_neighDeltaHashes , neighDeltaHashes , NneighCells*sizeof(int) )) ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_neighMask_x , neighMask_x , NneighCells*sizeof(int) )) ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_neighMask_y , neighMask_y , NneighCells*sizeof(int) )) ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_neighMask_z , neighMask_z , NneighCells*sizeof(int) )) ;
  //  create a mask for the PBC (0 or abs pbc value to sum for x,y,z, then the computations in the kernel will be only (c.x==0 - c.x==d_CellSide.x) . . . )
  // Servono 3 vettori deltaHashPBC, per x (+o- Cx), y(+o- Cx*Cy), z (+o- Cx*Cy*Cz)... o forse 6, per + e per - ... e si calcolano solo grid.i==0oN-1... ma forse non 6, se si definiscono le mask come scritto sotto...
  //			  (-c.x==0 + c.x==d_CellSide.x) == mask (verde sopra)
}
/****************************************************************************************/

// GPU initialization
template <typename particle>
void configuration<particle>::initialize_GPU_data( void ) {
  double max_cutoff = 0 ;

  // Allocation of global parameters
  CUDA_CHECK(cudaMemcpyToSymbol( d_numberOfPart , &numberOfPart , sizeof(int) )) ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_time_step , &time_step , sizeof(double) )) ;
  double half_time_step = 0.5*time_step ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_half_time_step , &half_time_step , sizeof(double) )) ;
  allocateDeviceBox() ;
  allocateDeviceNeighHashes() ;

  // Allocation of particles arrays
  allocateDeviceArrays( 1 , 1 ) ;
  copyParticlesToDevice( 1 , 1 , 1 ) ;  

  // Allocation of Cell list structures
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
    if( ! INTERACTIONS.type[k]->ON_GPU ) {
      cout << "*** ERROR: Potential " << INTERACTIONS.type[k]->idx << ", " << INTERACTIONS.type[k]->name <<
	" cannot be computed on GPU!" << endl ;
      exit( EXIT_FAILURE ) ;
    }
    max_cutoff = (INTERACTIONS.type[k]->cutoff > max_cutoff)  ?  INTERACTIONS.type[k]->cutoff  :  max_cutoff ;
  }
  allocateDeviceCellList( max_cutoff ) ;

  /* Determination of the optimal block and grid sizes for GPU kernels execution */
  cout << "***********************************************************************" << endl ;
  cout << "=======================================================================" << endl ;
  cout << "_____ Find Cells" << endl ;
  findOptimalGrid( GPU_params.cell_grid , GPU_params.cell_block , findCellStartEnd , numberOfPart , 0 ) ;
  for(int i=0; i<3; i++) findOptimalGrid( GPU_params.cell_grid , GPU_params.cell_block , findCellStartEnd , numberOfPart , GPU_params.cell_block.x+1 ) ;
  cout << "_____ Compute Hash" << endl ;
  if( DIMENSION == 3 ) findOptimalGrid( GPU_params.grid , GPU_params.block ,
					static_cast<void (*)(int*, int*, double*, double*, double*, int*, int*, int*)>(calcCellHash) ,
					numberOfPart , 0 ) ;
  else if( DIMENSION == 2 ) findOptimalGrid( GPU_params.grid , GPU_params.block ,
					     static_cast<void (*)(int*, int*, double*, double*, int*, int*)>(calcCellHash) ,
					     numberOfPart , 0 ) ;
  cout << "_____ Reorder Arrays" << endl ;
  findOptimalGrid( GPU_params.grid , GPU_params.block ,
		   static_cast<void (*)(double*, double*, int*)>(reorderDeviceArrays) ,
		   numberOfPart , 0 ) ;
  cout << "=======================================================================" << endl ;
  cout << "***********************************************************************" << endl ;
  buildDeviceCellList() ;
}
/****************************************************************************************/

// Alloc device arrays
template <typename particle>
void configuration<particle>::allocateDeviceArrays( bool VEL , bool FORCE ) {
  size_t bytes = numberOfPart * sizeof(int) ;
  CUDA_CHECK(cudaMalloc( &(GPU_arrays.d_idx) ,  bytes )) ;
  CUDA_CHECK(cudaMalloc( &GPU_arrays.d_ix ,  bytes )) ;
  CUDA_CHECK(cudaMalloc( &GPU_arrays.d_iy ,  bytes )) ;
  bytes = numberOfPart * sizeof(double) ;
  CUDA_CHECK(cudaMalloc( &GPU_arrays.d_x ,  bytes )) ;
  CUDA_CHECK(cudaMalloc( &GPU_arrays.d_y ,  bytes )) ;
  if( VEL ){
    CUDA_CHECK(cudaMalloc( &GPU_arrays.d_vx , bytes )) ;
    CUDA_CHECK(cudaMalloc( &GPU_arrays.d_vy , bytes )) ;
  }
  if( FORCE ){
    CUDA_CHECK(cudaMalloc( &GPU_arrays.d_fx , bytes )) ;
    CUDA_CHECK(cudaMalloc( &GPU_arrays.d_fy , bytes )) ;
  }
  CUDA_CHECK(cudaMalloc( &GPU_arrays.d_compute_vec ,  bytes )) ;
}
/****************************************************************************************/

template <>
void configuration<particle_3D>::allocateDeviceArrays( bool VEL , bool FORCE ) {
  size_t bytes = numberOfPart * sizeof(int) ;
  CUDA_CHECK(cudaMalloc( &(GPU_arrays.d_idx) ,  bytes )) ;
  CUDA_CHECK(cudaMalloc( &GPU_arrays.d_ix ,  bytes )) ;
  CUDA_CHECK(cudaMalloc( &GPU_arrays.d_iy ,  bytes )) ;
  CUDA_CHECK(cudaMalloc( &GPU_arrays.d_iz ,  bytes )) ;
  bytes = numberOfPart * sizeof(double) ;
  CUDA_CHECK(cudaMalloc( &GPU_arrays.d_x ,  bytes )) ;
  CUDA_CHECK(cudaMalloc( &GPU_arrays.d_y ,  bytes )) ;
  CUDA_CHECK(cudaMalloc( &GPU_arrays.d_z ,  bytes )) ;
  if( VEL ){
    CUDA_CHECK(cudaMalloc( &GPU_arrays.d_vx , bytes )) ;
    CUDA_CHECK(cudaMalloc( &GPU_arrays.d_vy , bytes )) ;
    CUDA_CHECK(cudaMalloc( &GPU_arrays.d_vz , bytes )) ;
  }
  if( FORCE ){
    CUDA_CHECK(cudaMalloc( &GPU_arrays.d_fx , bytes )) ;
    CUDA_CHECK(cudaMalloc( &GPU_arrays.d_fy , bytes )) ;
    CUDA_CHECK(cudaMalloc( &GPU_arrays.d_fz , bytes )) ;
  }
  CUDA_CHECK(cudaMalloc( &GPU_arrays.d_compute_vec ,  bytes )) ;
}
/****************************************************************************************/

// Free device arrays
template <typename particle>
void configuration<particle>::freeDeviceArrays() {
  if( !GPU_arrays.d_idx ) CUDA_CHECK(cudaFree( GPU_arrays.d_idx )) ;
  if( !GPU_arrays.d_x ) CUDA_CHECK(cudaFree( GPU_arrays.d_x )) ;
  if( !GPU_arrays.d_y ) CUDA_CHECK(cudaFree( GPU_arrays.d_y )) ;
  if( !GPU_arrays.d_z ) CUDA_CHECK(cudaFree( GPU_arrays.d_z )) ;
  if( !GPU_arrays.d_ix ) CUDA_CHECK(cudaFree( GPU_arrays.d_ix )) ;
  if( !GPU_arrays.d_iy ) CUDA_CHECK(cudaFree( GPU_arrays.d_iy )) ;
  if( !GPU_arrays.d_iz ) CUDA_CHECK(cudaFree( GPU_arrays.d_iz )) ;
  if( ! GPU_arrays.d_vx ) CUDA_CHECK(cudaFree( GPU_arrays.d_vx )) ;
  if( ! GPU_arrays.d_vy ) CUDA_CHECK(cudaFree( GPU_arrays.d_vy )) ;
  if( ! GPU_arrays.d_vz ) CUDA_CHECK(cudaFree( GPU_arrays.d_vz )) ;
  if( ! GPU_arrays.d_fx ) CUDA_CHECK(cudaFree( GPU_arrays.d_fx )) ;
  if( ! GPU_arrays.d_fy ) CUDA_CHECK(cudaFree( GPU_arrays.d_fy )) ;
  if( ! GPU_arrays.d_fz ) CUDA_CHECK(cudaFree( GPU_arrays.d_fz )) ;
  if( !GPU_arrays.d_compute_vec ) CUDA_CHECK(cudaFree( GPU_arrays.d_compute_vec )) ;
}
/****************************************************************************************/

// Copy host -> device
template <typename particle>
void configuration<particle>::copyParticlesToDevice( bool POS , bool VEL , bool FORCE ) {
    std::vector<double> hx(numberOfPart) , hy(numberOfPart) ;
    std::vector<int> hix(numberOfPart) , hiy(numberOfPart) ;
    size_t bytes ;

    if( GPU_params.block.x*GPU_params.block.y*GPU_params.block.z * GPU_params.block.x*GPU_params.block.y*GPU_params.block.z >= static_cast<unsigned int>(numberOfPart) ) {
      fill_vectors_range<<<GPU_params.grid, GPU_params.block>>>( GPU_arrays.d_idx , 0 ) ;
      cudaDeviceSynchronize() ;
    } else {
      int block = 1024 ;
      int grid = (numberOfPart + block - 1) / block ;
      fill_vectors_range<<<grid, block>>>( GPU_arrays.d_idx , 0 ) ;
      cudaDeviceSynchronize() ;
    }

    bytes = numberOfPart * sizeof(double) ;
    if( POS ) {
      for(int i=0; i<numberOfPart; i++) {
        hx[i]  = particles[i]->position.x ;
        hy[i]  = particles[i]->position.y ;
        hix[i]  = particles[i]->periodic_box.x ;
        hiy[i]  = particles[i]->periodic_box.y ;
      }
      bytes = numberOfPart * sizeof(int) ;
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_ix , hix.data() , bytes , cudaMemcpyHostToDevice )) ;
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_iy , hiy.data() , bytes , cudaMemcpyHostToDevice )) ;
      bytes = numberOfPart * sizeof(double) ;
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_x , hx.data() , bytes , cudaMemcpyHostToDevice )) ;
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_y , hy.data() , bytes , cudaMemcpyHostToDevice )) ;
    }

    if( VEL ) {
      for(int i=0; i<numberOfPart; i++) {
        hx[i] = velocities[i]->position.x ;
        hy[i] = velocities[i]->position.y ;
      }
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_vx , hx.data() , bytes , cudaMemcpyHostToDevice )) ;
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_vy , hy.data() , bytes , cudaMemcpyHostToDevice )) ;
    }

    if( FORCE ) {
      for(int i=0; i<numberOfPart; i++) {
        hx[i] = forces_t2[i]->position.x ;
        hy[i] = forces_t2[i]->position.y ;
      }
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_fx , hx.data() , bytes , cudaMemcpyHostToDevice )) ;
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_fy , hy.data() , bytes , cudaMemcpyHostToDevice )) ;
    }
}
/****************************************************************************************/

template <>
void configuration<particle_3D>::copyParticlesToDevice( bool POS , bool VEL , bool FORCE ) {
    std::vector<double> hx(numberOfPart) , hy(numberOfPart) , hz(numberOfPart) ;
    std::vector<int> hix(numberOfPart) , hiy(numberOfPart) , hiz(numberOfPart) ;
    size_t bytes ;

    if( GPU_params.block.x*GPU_params.block.y*GPU_params.block.z * GPU_params.block.x*GPU_params.block.y*GPU_params.block.z >= static_cast<unsigned int>(numberOfPart) ) {
      fill_vectors_range<<<GPU_params.grid, GPU_params.block>>>( GPU_arrays.d_idx , 0 ) ;
      cudaDeviceSynchronize() ;
    } else {
      int block = 1024 ;
      int grid = (numberOfPart + block - 1) / block ;
      fill_vectors_range<<<grid, block>>>( GPU_arrays.d_idx , 0 ) ;
      cudaDeviceSynchronize() ;
    }

    bytes = numberOfPart * sizeof(double) ;
    if( POS ) {
      for( int i=0; i<numberOfPart; i++ ) {
        hx[i]  = particles[i]->position.x ;
        hy[i]  = particles[i]->position.y ;
        hz[i]  = particles[i]->position.z ;
        hix[i]  = particles[i]->periodic_box.x ;
        hiy[i]  = particles[i]->periodic_box.y ;
        hiz[i]  = particles[i]->periodic_box.z ;
      }
      bytes = numberOfPart * sizeof(int) ;
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_ix , hix.data() , bytes , cudaMemcpyHostToDevice )) ;
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_iy , hiy.data() , bytes , cudaMemcpyHostToDevice )) ;
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_iz , hiz.data() , bytes , cudaMemcpyHostToDevice )) ;
      bytes = numberOfPart * sizeof(double) ;
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_x , hx.data() , bytes , cudaMemcpyHostToDevice )) ;
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_y , hy.data() , bytes , cudaMemcpyHostToDevice )) ;
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_z , hz.data() , bytes , cudaMemcpyHostToDevice )) ;
    }

    if( VEL ) {
      for( int i=0; i<numberOfPart; i++ ) {
        hx[i] = velocities[i]->position.x ;
        hy[i] = velocities[i]->position.y ;
        hz[i] = velocities[i]->position.z ;
      }
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_vx , hx.data() , bytes , cudaMemcpyHostToDevice )) ;
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_vy , hy.data() , bytes , cudaMemcpyHostToDevice )) ;
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_vz , hz.data() , bytes , cudaMemcpyHostToDevice )) ;
    }

    if( FORCE ) {
      for( int i=0; i<numberOfPart; i++ ) {
        hx[i] = forces_t2[i]->position.x ;
        hy[i] = forces_t2[i]->position.y ;
        hz[i] = forces_t2[i]->position.z ;
      }
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_fx , hx.data() , bytes , cudaMemcpyHostToDevice )) ;
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_fy , hy.data() , bytes , cudaMemcpyHostToDevice )) ;
      CUDA_CHECK(cudaMemcpy( GPU_arrays.d_fz , hz.data() , bytes , cudaMemcpyHostToDevice )) ;
    }
}
/****************************************************************************************/


// Copy device -> host
template <typename particle>
void configuration<particle>::copyParticlesToHost( bool POS , bool VEL , bool FORCE , bool POS_IMG ) {
    std::vector<int> hidx(numberOfPart) ;
    std::vector<double> hx(numberOfPart) , hy(numberOfPart) ;
    std::vector<double> hvx(numberOfPart) , hvy(numberOfPart) ;
    std::vector<double> hfx(numberOfPart) , hfy(numberOfPart) ;
    size_t bytes ;
    int idx ;

    bytes = numberOfPart * sizeof(int) ;
    CUDA_CHECK(cudaMemcpy( hidx.data() , GPU_arrays.d_idx , bytes , cudaMemcpyDeviceToHost )) ;

    if( POS_IMG ) {
      std::vector<int> hix(numberOfPart) , hiy(numberOfPart) ;
      CUDA_CHECK(cudaMemcpy( hix.data() ,  GPU_arrays.d_ix ,  bytes , cudaMemcpyDeviceToHost )) ;
      CUDA_CHECK(cudaMemcpy( hiy.data() ,  GPU_arrays.d_iy ,  bytes , cudaMemcpyDeviceToHost )) ;
      for( int i=0; i<numberOfPart; i++ ){
	idx = hidx[i] ;
        particles[idx]->periodic_box.x = hix[i];
        particles[idx]->periodic_box.y = hiy[i];
      }
    }

    bytes = numberOfPart * sizeof(double) ;
    if( POS ) {
      CUDA_CHECK(cudaMemcpy( hx.data() ,  GPU_arrays.d_x ,  bytes , cudaMemcpyDeviceToHost )) ;
      CUDA_CHECK(cudaMemcpy( hy.data() ,  GPU_arrays.d_y ,  bytes , cudaMemcpyDeviceToHost )) ;
    }
    if( VEL ) {
      CUDA_CHECK(cudaMemcpy( hvx.data() , GPU_arrays.d_vx , bytes , cudaMemcpyDeviceToHost )) ;
      CUDA_CHECK(cudaMemcpy( hvy.data() , GPU_arrays.d_vy , bytes , cudaMemcpyDeviceToHost )) ;
    }
    if( FORCE ) {
      CUDA_CHECK(cudaMemcpy( hfx.data() , GPU_arrays.d_fx , bytes , cudaMemcpyDeviceToHost )) ;
      CUDA_CHECK(cudaMemcpy( hfy.data() , GPU_arrays.d_fy , bytes , cudaMemcpyDeviceToHost )) ;
    }

    // merge back into AoS
    if( POS && VEL && FORCE ) {
      for( int i=0; i<numberOfPart; i++ ){
	idx = hidx[i] ;
        particles[idx]->position.x = hx[i];
        particles[idx]->position.y = hy[i];

        velocities[idx]->position.x = hvx[i];
        velocities[idx]->position.y = hvy[i];

        forces_t1[idx]->position.x = hfx[i];
        forces_t1[idx]->position.y = hfy[i];
      }
    } else if( POS && VEL ) {
      for( int i=0; i<numberOfPart; i++ ){
	idx = hidx[i] ;
        particles[idx]->position.x = hx[i];
        particles[idx]->position.y = hy[i];

        velocities[idx]->position.x = hvx[i];
        velocities[idx]->position.y = hvy[i];
      }
    } else if( POS && FORCE ) {
      for( int i=0; i<numberOfPart; i++ ){
	idx = hidx[i] ;
        particles[idx]->position.x = hx[i];
        particles[idx]->position.y = hy[i];

        forces_t1[idx]->position.x = hfx[i];
        forces_t1[idx]->position.y = hfy[i];
      }
    } else if( VEL && FORCE ) {
      for( int i=0; i<numberOfPart; i++ ){
	idx = hidx[i] ;
        velocities[idx]->position.x = hvx[i];
        velocities[idx]->position.y = hvy[i];

        forces_t1[idx]->position.x = hfx[i];
        forces_t1[idx]->position.y = hfy[i];
      }
    } else if( POS ) {
      for( int i=0; i<numberOfPart; i++ ){
	idx = hidx[i] ;
        particles[idx]->position.x = hx[i];
        particles[idx]->position.y = hy[i];
      }
    } else if( VEL ) {
      for( int i=0; i<numberOfPart; i++ ){
	idx = hidx[i] ;
        velocities[idx]->position.x = hvx[i];
        velocities[idx]->position.y = hvy[i];
      }
    } else if( FORCE ) {
      for( int i=0; i<numberOfPart; i++ ){
	idx = hidx[i] ;
        forces_t1[idx]->position.x = hfx[i];
        forces_t1[idx]->position.y = hfy[i];
      }
    }
}
/****************************************************************************************/

template <>
void configuration<particle_3D>::copyParticlesToHost( bool POS , bool VEL , bool FORCE , bool POS_IMG ) {
    std::vector<int> hidx(numberOfPart) ;
    std::vector<double> hx(numberOfPart) , hy(numberOfPart) , hz(numberOfPart) ;
    std::vector<double> hvx(numberOfPart) , hvy(numberOfPart) , hvz(numberOfPart) ;
    std::vector<double> hfx(numberOfPart) , hfy(numberOfPart) , hfz(numberOfPart) ;
    size_t bytes ;
    int idx ;

    bytes = numberOfPart * sizeof(int) ;
    CUDA_CHECK(cudaMemcpy( hidx.data() , GPU_arrays.d_idx , bytes , cudaMemcpyDeviceToHost )) ;
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess)  std::cerr << "CUDA error after reorderDeviceArrays: " << cudaGetErrorString(error) << std::endl;

    if( POS_IMG ) {
      std::vector<int> hix(numberOfPart) , hiy(numberOfPart) , hiz(numberOfPart) ;
      CUDA_CHECK(cudaMemcpy( hix.data() ,  GPU_arrays.d_ix ,  bytes , cudaMemcpyDeviceToHost )) ;
      CUDA_CHECK(cudaMemcpy( hiy.data() ,  GPU_arrays.d_iy ,  bytes , cudaMemcpyDeviceToHost )) ;
      CUDA_CHECK(cudaMemcpy( hiz.data() ,  GPU_arrays.d_iz ,  bytes , cudaMemcpyDeviceToHost )) ;
      for( int i=0; i<numberOfPart; i++ ){
	idx = hidx[i] ;
        particles[idx]->periodic_box.x = hix[i];
        particles[idx]->periodic_box.y = hiy[i];
        particles[idx]->periodic_box.z = hiz[i];
      }
    }

    bytes = numberOfPart * sizeof(double) ;
    if( POS ) {
      CUDA_CHECK(cudaMemcpy( hx.data() ,  GPU_arrays.d_x ,  bytes , cudaMemcpyDeviceToHost )) ;
      CUDA_CHECK(cudaMemcpy( hy.data() ,  GPU_arrays.d_y ,  bytes , cudaMemcpyDeviceToHost )) ;
      CUDA_CHECK(cudaMemcpy( hz.data() ,  GPU_arrays.d_z ,  bytes , cudaMemcpyDeviceToHost )) ;
    }
    if( VEL ) {
      CUDA_CHECK(cudaMemcpy( hvx.data() , GPU_arrays.d_vx , bytes , cudaMemcpyDeviceToHost )) ;
      CUDA_CHECK(cudaMemcpy( hvy.data() , GPU_arrays.d_vy , bytes , cudaMemcpyDeviceToHost )) ;
      CUDA_CHECK(cudaMemcpy( hvz.data() , GPU_arrays.d_vz , bytes , cudaMemcpyDeviceToHost )) ;
    }
    if( FORCE ) {
      CUDA_CHECK(cudaMemcpy( hfx.data() , GPU_arrays.d_fx , bytes , cudaMemcpyDeviceToHost )) ;
      CUDA_CHECK(cudaMemcpy( hfy.data() , GPU_arrays.d_fy , bytes , cudaMemcpyDeviceToHost )) ;
      CUDA_CHECK(cudaMemcpy( hfz.data() , GPU_arrays.d_fz , bytes , cudaMemcpyDeviceToHost )) ;
    }

    // merge back into AoS
    if( POS && VEL && FORCE ) {
      for( int i=0; i<numberOfPart; i++ ){
	idx = hidx[i] ;
        particles[idx]->position.x = hx[i];
        particles[idx]->position.y = hy[i];
        particles[idx]->position.z = hz[i];

        velocities[idx]->position.x = hvx[i];
        velocities[idx]->position.y = hvy[i];
        velocities[idx]->position.z = hvz[i];

        forces_t1[idx]->position.x = hfx[i];
        forces_t1[idx]->position.y = hfy[i];
        forces_t1[idx]->position.z = hfz[i];
      }
    } else if( POS && VEL ) {
      for( int i=0; i<numberOfPart; i++ ){
	idx = hidx[i] ;
        particles[idx]->position.x = hx[i];
        particles[idx]->position.y = hy[i];
        particles[idx]->position.z = hz[i];

        velocities[idx]->position.x = hvx[i];
        velocities[idx]->position.y = hvy[i];
        velocities[idx]->position.z = hvz[i];
      }
    } else if( POS && FORCE ) {
      for( int i=0; i<numberOfPart; i++ ){
	idx = hidx[i] ;
        particles[idx]->position.x = hx[i];
        particles[idx]->position.y = hy[i];
        particles[idx]->position.z = hz[i];

        forces_t1[idx]->position.x = hfx[i];
        forces_t1[idx]->position.y = hfy[i];
        forces_t1[idx]->position.z = hfz[i];
      }
    } else if( VEL && FORCE ) {
      for( int i=0; i<numberOfPart; i++ ){
	idx = hidx[i] ;

        velocities[idx]->position.x = hvx[i];
        velocities[idx]->position.y = hvy[i];
        velocities[idx]->position.z = hvz[i];

        forces_t1[idx]->position.x = hfx[i];
        forces_t1[idx]->position.y = hfy[i];
        forces_t1[idx]->position.z = hfz[i];
      }
    } else if( POS ) {
      for( int i=0; i<numberOfPart; i++ ){
	idx = hidx[i] ;
        particles[idx]->position.x = hx[i];
        particles[idx]->position.y = hy[i];
        particles[idx]->position.z = hz[i];
      }
    } else if( VEL ) {
      for( int i=0; i<numberOfPart; i++ ){
	idx = hidx[i] ;
        velocities[idx]->position.x = hvx[i];
        velocities[idx]->position.y = hvy[i];
        velocities[idx]->position.z = hvz[i];
      }
    } else if( FORCE ) {
      for( int i=0; i<numberOfPart; i++ ){
	idx = hidx[i] ;
        forces_t1[idx]->position.x = hfx[i];
        forces_t1[idx]->position.y = hfy[i];
        forces_t1[idx]->position.z = hfz[i];
      }
    }
}
/****************************************************************************************/

// Alloc device cells
template <typename particle>
void configuration<particle>::allocateDeviceCellList( double cutoff ) {
  GPU_cells.CellSize.x = (int)ceil( box_sides.position.x / cutoff ) ;
  GPU_cells.CellSize.y = (int)ceil( box_sides.position.y / cutoff ) ;
  GPU_cells.CellSize.z = 1 ;
  GPU_cells.numCells = GPU_cells.CellSize.x * GPU_cells.CellSize.y * GPU_cells.CellSize.z ;
  GPU_cells.CellSide.x = box_sides.position.x / GPU_cells.CellSize.x ;
  GPU_cells.CellSide.y = box_sides.position.y / GPU_cells.CellSize.y ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_CellSide , &GPU_cells.CellSide , sizeof(double)*3 )) ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_CellSize , &GPU_cells.CellSize , sizeof(int)*3 )) ;
  int CellHashPBC[3] = { GPU_cells.CellSize.x , GPU_cells.CellSize.x * GPU_cells.CellSize.y , GPU_cells.CellSize.x * GPU_cells.CellSize.y * GPU_cells.CellSize.z } ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_CellHashPBC , &CellHashPBC , sizeof(int)*3 )) ;

  CUDA_CHECK(cudaMalloc( (void **)&GPU_cells.m_dCellStart , sizeof(int)*GPU_cells.numCells ));
  CUDA_CHECK(cudaMalloc( (void **)&GPU_cells.m_dCellEnd , sizeof(int)*GPU_cells.numCells ));
  CUDA_CHECK(cudaMalloc( (void **)&GPU_cells.m_dCellParticleIndex , sizeof(int)*numberOfPart ));
  CUDA_CHECK(cudaMalloc( (void **)&GPU_cells.m_dCellParticleHash , sizeof(int)*numberOfPart ));
  CUDA_CHECK(cudaMalloc( &GPU_cells.d_sorted_int ,  sizeof(int)*numberOfPart )) ;
  CUDA_CHECK(cudaMalloc( &GPU_cells.d_sorted_double ,  sizeof(double)*numberOfPart )) ;

  allocateDeviceNeighHashes() ;
}
/****************************************************************************************/

// Alloc device cells
template <>
void configuration<particle_3D>::allocateDeviceCellList( double cutoff ) {
  GPU_cells.CellSize.x = (int)ceil( box_sides.position.x / cutoff ) ;
  GPU_cells.CellSize.y = (int)ceil( box_sides.position.y / cutoff ) ;
  GPU_cells.CellSize.z = (int)ceil( box_sides.position.z / cutoff ) ;
  GPU_cells.numCells = GPU_cells.CellSize.x * GPU_cells.CellSize.y * GPU_cells.CellSize.z ;
  GPU_cells.CellSide.x = box_sides.position.x / GPU_cells.CellSize.x ;
  GPU_cells.CellSide.y = box_sides.position.y / GPU_cells.CellSize.y ;
  GPU_cells.CellSide.z = box_sides.position.z / GPU_cells.CellSize.z ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_CellSide , &GPU_cells.CellSide , sizeof(double)*3 )) ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_CellSize , &GPU_cells.CellSize , sizeof(int)*3 )) ;
  int CellHashPBC[3] = { GPU_cells.CellSize.x , GPU_cells.CellSize.x * GPU_cells.CellSize.y , GPU_cells.CellSize.x * GPU_cells.CellSize.y * GPU_cells.CellSize.z } ;
  CUDA_CHECK(cudaMemcpyToSymbol( d_CellHashPBC , &CellHashPBC , sizeof(int)*3 )) ;

  CUDA_CHECK(cudaMalloc( (void **)&GPU_cells.m_dCellStart , sizeof(int)*GPU_cells.numCells ));
  CUDA_CHECK(cudaMalloc( (void **)&GPU_cells.m_dCellEnd , sizeof(int)*GPU_cells.numCells ));
  CUDA_CHECK(cudaMalloc( (void **)&GPU_cells.m_dCellParticleIndex , sizeof(int)*numberOfPart ));
  CUDA_CHECK(cudaMalloc( (void **)&GPU_cells.m_dCellParticleHash , sizeof(int)*numberOfPart ));
  CUDA_CHECK(cudaMalloc( &GPU_cells.d_sorted_int ,  sizeof(int)*numberOfPart )) ;
  CUDA_CHECK(cudaMalloc( &GPU_cells.d_sorted_double ,  sizeof(double)*numberOfPart )) ;

  allocateDeviceNeighHashes() ;
}
/****************************************************************************************/

// Free device cells
template <typename particle>
void configuration<particle>::freeDeviceCellList( void ) {
  if( ! GPU_cells.m_dCellStart ) CUDA_CHECK(cudaFree( GPU_cells.m_dCellStart )) ;
  if( ! GPU_cells.m_dCellEnd ) CUDA_CHECK(cudaFree( GPU_cells.m_dCellEnd )) ;
  if( ! GPU_cells.m_dCellParticleIndex ) CUDA_CHECK(cudaFree( GPU_cells.m_dCellParticleIndex )) ;
  if( ! GPU_cells.m_dCellParticleHash ) CUDA_CHECK(cudaFree( GPU_cells.m_dCellParticleHash )) ;
  if( ! GPU_cells.d_sorted_int ) CUDA_CHECK(cudaFree( GPU_cells.d_sorted_int )) ;
  if( ! GPU_cells.d_sorted_double ) CUDA_CHECK(cudaFree( GPU_cells.d_sorted_double )) ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::buildDeviceCellList() {
  int *swapi ;
  double *swapd ;
  calcCellHash<<<GPU_params.grid, GPU_params.block>>>( GPU_cells.m_dCellParticleHash , GPU_cells.m_dCellParticleIndex ,
				 GPU_arrays.d_x , GPU_arrays.d_y , GPU_arrays.d_ix , GPU_arrays.d_iy ) ;
  cudaDeviceSynchronize() ;
  // sort particles array
  thrust::sort_by_key(thrust::device_ptr<int>(GPU_cells.m_dCellParticleHash),
		      thrust::device_ptr<int>(GPU_cells.m_dCellParticleHash) + numberOfPart,
		      thrust::device_ptr<int>(GPU_cells.m_dCellParticleIndex) );

  // set all cells to empty
  CUDA_CHECK(cudaMemset( GPU_cells.m_dCellStart , 0xffffffff , GPU_cells.numCells*sizeof(int) )) ;

  int smemSize = sizeof(int) * (GPU_params.cell_block.x+1) ;
  findCellStartEnd<<<GPU_params.cell_grid, GPU_params.cell_block, smemSize>>>( GPU_cells.m_dCellParticleHash , GPU_cells.m_dCellStart , GPU_cells.m_dCellEnd ) ;
  cudaDeviceSynchronize();
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess)  std::cerr << "CUDA error after reorderDeviceArrays: " << cudaGetErrorString(error) << std::endl;

  // reordering all the state vectors
  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_x , GPU_cells.d_sorted_double , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapd = GPU_arrays.d_x ;
  GPU_arrays.d_x = GPU_cells.d_sorted_double ;
  GPU_cells.d_sorted_double = swapd ;
  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_y , GPU_cells.d_sorted_double , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapd = GPU_arrays.d_y ;
  GPU_arrays.d_y = GPU_cells.d_sorted_double ;
  GPU_cells.d_sorted_double = swapd ;
  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_vx , GPU_cells.d_sorted_double , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapd = GPU_arrays.d_vx ;
  GPU_arrays.d_vx = GPU_cells.d_sorted_double ;
  GPU_cells.d_sorted_double = swapd ;
  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_vy , GPU_cells.d_sorted_double , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapd = GPU_arrays.d_vy ;
  GPU_arrays.d_vy = GPU_cells.d_sorted_double ;
  GPU_cells.d_sorted_double = swapd ;

  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_idx , GPU_cells.d_sorted_int , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapi = GPU_arrays.d_idx ;
  GPU_arrays.d_idx = GPU_cells.d_sorted_int ;
  GPU_cells.d_sorted_int = swapi ;
  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_ix , GPU_cells.d_sorted_int , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapi = GPU_arrays.d_ix ;
  GPU_arrays.d_ix = GPU_cells.d_sorted_int ;
  GPU_cells.d_sorted_int = swapi ;
  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_iy , GPU_cells.d_sorted_int , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapi = GPU_arrays.d_iy ;
  GPU_arrays.d_iy = GPU_cells.d_sorted_int ;
  GPU_cells.d_sorted_int = swapi ;
}
/****************************************************************************************/

template <>
void configuration<particle_3D>::buildDeviceCellList() {
  int *swapi ;
  double *swapd ;
  calcCellHash<<<GPU_params.grid, GPU_params.block>>>( GPU_cells.m_dCellParticleHash , GPU_cells.m_dCellParticleIndex ,
				 GPU_arrays.d_x , GPU_arrays.d_y , GPU_arrays.d_z , GPU_arrays.d_ix , GPU_arrays.d_iy , GPU_arrays.d_iz ) ;
  cudaDeviceSynchronize() ;
  // sort particles array
  thrust::sort_by_key(thrust::device_ptr<int>(GPU_cells.m_dCellParticleHash),
		      thrust::device_ptr<int>(GPU_cells.m_dCellParticleHash) + numberOfPart,
		      thrust::device_ptr<int>(GPU_cells.m_dCellParticleIndex) );

  // set all cells to empty
  CUDA_CHECK(cudaMemset( GPU_cells.m_dCellStart , 0xffffffff , GPU_cells.numCells*sizeof(int) )) ;

  int smemSize = sizeof(int) * (GPU_params.cell_block.x+1) ;
  findCellStartEnd<<<GPU_params.cell_grid, GPU_params.cell_block, smemSize>>>( GPU_cells.m_dCellParticleHash , GPU_cells.m_dCellStart , GPU_cells.m_dCellEnd ) ;
  cudaDeviceSynchronize() ;
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess)  std::cerr << "CUDA error after reorderDeviceArrays: " << cudaGetErrorString(error) << std::endl;

  // reordering all the state vectors
  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_x , GPU_cells.d_sorted_double , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapd = GPU_arrays.d_x ;
  GPU_arrays.d_x = GPU_cells.d_sorted_double ;
  GPU_cells.d_sorted_double = swapd ;
  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_y , GPU_cells.d_sorted_double , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapd = GPU_arrays.d_y ;
  GPU_arrays.d_y = GPU_cells.d_sorted_double ;
  GPU_cells.d_sorted_double = swapd ;
  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_z , GPU_cells.d_sorted_double , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapd = GPU_arrays.d_z ;
  GPU_arrays.d_z = GPU_cells.d_sorted_double ;
  GPU_cells.d_sorted_double = swapd ;
  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_vx , GPU_cells.d_sorted_double , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapd = GPU_arrays.d_vx ;
  GPU_arrays.d_vx = GPU_cells.d_sorted_double ;
  GPU_cells.d_sorted_double = swapd ;
  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_vy , GPU_cells.d_sorted_double , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapd = GPU_arrays.d_vy ;
  GPU_arrays.d_vy = GPU_cells.d_sorted_double ;
  GPU_cells.d_sorted_double = swapd ;
  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_vz , GPU_cells.d_sorted_double , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapd = GPU_arrays.d_vz ;
  GPU_arrays.d_vz = GPU_cells.d_sorted_double ;
  GPU_cells.d_sorted_double = swapd ;

  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_idx , GPU_cells.d_sorted_int , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapi = GPU_arrays.d_idx ;
  GPU_arrays.d_idx = GPU_cells.d_sorted_int ;
  GPU_cells.d_sorted_int = swapi ;
  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_ix , GPU_cells.d_sorted_int , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapi = GPU_arrays.d_ix ;
  GPU_arrays.d_ix = GPU_cells.d_sorted_int ;
  GPU_cells.d_sorted_int = swapi ;
  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_iy , GPU_cells.d_sorted_int , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapi = GPU_arrays.d_iy ;
  GPU_arrays.d_iy = GPU_cells.d_sorted_int ;
  GPU_cells.d_sorted_int = swapi ;
  reorderDeviceArrays<<<GPU_params.grid,GPU_params.block>>>( GPU_arrays.d_iz , GPU_cells.d_sorted_int , GPU_cells.m_dCellParticleIndex ) ;
  cudaDeviceSynchronize() ;
  swapi = GPU_arrays.d_iz ;
  GPU_arrays.d_iz = GPU_cells.d_sorted_int ;
  GPU_cells.d_sorted_int = swapi ;
}
/****************************************************************************************/
