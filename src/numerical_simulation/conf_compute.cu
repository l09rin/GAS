
template <typename particle>
void configuration<particle>::calculate_forces_gpu() {
  if( strcmp( simulation_mode , "MC" ) == 0 ) {
    cout << " *** ERROR : function calculate_forces_gpu() cannot be used in MC mode !" << endl ;
    exit( EXIT_FAILURE ) ;
  }
  buildDeviceCellList() ;
  if( DIMENSION == 2 ) zero_vectors<<< GPU_params.grid , GPU_params.block >>>( GPU_arrays.d_fx , GPU_arrays.d_fy ) ;
  else zero_vectors<<< GPU_params.grid , GPU_params.block >>>( GPU_arrays.d_fx , GPU_arrays.d_fy , GPU_arrays.d_fz ) ;
  cudaDeviceSynchronize() ;
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) INTERACTIONS.type[k]->compute_forces_contribution() ;
}
/*****************************************************************************************/
