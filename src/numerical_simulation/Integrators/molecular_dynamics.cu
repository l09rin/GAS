/* methods to integrate the system in time with molecular dynamics */
#include "../../../lib/system_comm.h"

template <typename particle>
void configuration<particle>::MD_GPU( void ) {
  int current_step = 1 ;
  int starting_time = time(0) ;

  if( strcmp( simulation_mode , "MD" ) != 0 ) {
    cout << "ERROR: The simulation is not in the Molecular Dynamics mode !" << endl ;
    exit( EXIT_FAILURE ) ;
  }
  set_saving_options( "EQUILIBRATION" ) ;
  save_simulation_conditions( "MD_simulation_conditions.dat" ) ;

  cout << "_____ Velocity Verlet Kernel" << endl ;
  if( DIMENSION == 3 ) findOptimalGrid( GPU_params.grid , GPU_params.block ,
					static_cast<void (*)(double*, double*, double*, double*, double*, double*, double*, double*, double*)>(VelocityVerletH1) ,
					numberOfPart , 0 ) ;
  else if( DIMENSION == 2 ) findOptimalGrid( GPU_params.grid , GPU_params.block ,
					     static_cast<void (*)(double*, double*, double*, double*, double*, double*)>(VelocityVerletH1) ,
					     numberOfPart , 0 ) ;

  /*************                  equilibration stage                  *************/
  if( EQUILIBRATION_STEPS != 0 ) {
    if( SAVE.backup_interval == 0 ) SAVE.backup_interval = EQUILIBRATION_STEPS + 10 ;
                 // generation of the sub-directory containing the equilibration data
    char ROOT[100] ;
    if( generate_equilibration_subdirectory( ROOT ) == NULL ) cout << " Impossible to open equilibration sub-directory\n";
    initialize_storage_file( -1 , seed ) ;

    printf( "\n  EQUILIBRATION STAGE . . .\n" ) ; fflush(0) ;
    int equilibration_steps = EQUILIBRATION_STEPS ;
    int total_step_div_10 = floor( equilibration_steps/10. ) ;
    if( total_step_div_10 == 0 ) total_step_div_10 = 1 ;
    check_action( 0 , starting_time , total_step_div_10 ) ;

    if( strcmp( simulation_ensemble , "NVE" ) == 0 ) {
      for( current_step=current_step ; current_step<=equilibration_steps ; current_step++ ) {                   // integration steps
	VelocityVerletIntegrationStep_GPU() ;
	check_action( current_step , starting_time , total_step_div_10 ) ;
      }
    }else{
      cout << "*** ERROR: This integration scheme is not available on GPU." << endl ;
      exit( EXIT_FAILURE ) ;
    }

    close_storage_file() ;
    if( strcpy( _DIRECTORY_NAME , ROOT ) == NULL ) cout << " Impossible to open simulation directory" ;
  }

  /*************                  production stage                  *************/
  starting_time = time(0) ;
  set_saving_options( "PRODUCTION" ) ;
  if( PRODUCTION_STEPS != 0 ) {
    int production_steps = PRODUCTION_STEPS ;
    int total_step_div_10 = floor( production_steps/10. ) ;
    if( total_step_div_10 == 0 ) total_step_div_10 = 1 ;

    printf( "\n  PRODUCTION STAGE . . .\n" ) ; fflush(0) ;
    initialize_storage_file( production_steps , seed ) ;     // open the file in which configurations are periodically saved
    if( SAVE.backup_interval == 0 ) SAVE.backup_interval = PRODUCTION_STEPS + 10 ;

    check_action( 0 , starting_time , total_step_div_10 ) ;
    if( strcmp( simulation_ensemble , "NVE" ) == 0 ) {
      for( current_step=1 ; current_step<=production_steps ; current_step++ ) {
	VelocityVerletIntegrationStep_GPU() ;
	check_action( current_step , starting_time , total_step_div_10 ) ;
      }
    }else{
      cout << "*** ERROR: This integration scheme is not available on GPU." << endl ;
      exit( EXIT_FAILURE ) ;
    }
  
    cout << endl ;
    fflush( NULL ) ;
    close_storage_file() ;
  }
}
/*****************************************************************************************/

template <typename particle>
inline void configuration<particle>::VelocityVerletIntegrationStep_GPU(void) {
  VelocityVerletH1<<< GPU_params.grid , GPU_params.block >>>( GPU_arrays.d_x , GPU_arrays.d_y , GPU_arrays.d_vx , GPU_arrays.d_vy , GPU_arrays.d_fx , GPU_arrays.d_fy ) ;
  cudaDeviceSynchronize() ;
  calculate_forces_gpu() ;
  VelocityVerletH2<<< GPU_params.grid , GPU_params.block >>>( GPU_arrays.d_vx , GPU_arrays.d_vy , GPU_arrays.d_fx , GPU_arrays.d_fy ) ;
  cudaDeviceSynchronize() ;
  global_time += time_step ;
}
/****************************************************************************************/

template <>
inline void configuration<particle_3D>::VelocityVerletIntegrationStep_GPU(void) {
  VelocityVerletH1<<< GPU_params.grid , GPU_params.block >>>( GPU_arrays.d_x , GPU_arrays.d_y , GPU_arrays.d_z , GPU_arrays.d_vx , GPU_arrays.d_vy , GPU_arrays.d_vz , GPU_arrays.d_fx , GPU_arrays.d_fy , GPU_arrays.d_fz ) ;
  cudaDeviceSynchronize() ;
  calculate_forces_gpu() ;
  VelocityVerletH2<<< GPU_params.grid , GPU_params.block >>>( GPU_arrays.d_vx , GPU_arrays.d_vy , GPU_arrays.d_vz , GPU_arrays.d_fx , GPU_arrays.d_fy , GPU_arrays.d_fz ) ;
  cudaDeviceSynchronize() ;
  global_time += time_step ;
}
/****************************************************************************************/
