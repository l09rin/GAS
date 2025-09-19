/* methods to integrate the system in time with molecular dynamics */
#include "../../../lib/system_comm.h"

/**************************   BLOCK 3 - MD integrators   ********************************/
/****************************************************************************************/

template <typename particle>
void configuration<particle>::MD( const char *questions ) {
  char control, command_output ;
  int current_step = 1 ;
  int starting_time = time(0) ;
  // Variables needed for NVT integrators
  double half_timestep = 0.5 * time_step , N3kT = 3.0 * ((double)numberOfPart) * k_T ;

  if( strcmp( simulation_mode , "MD" ) != 0 ) {
    cout << "ERROR: The simulation is not in the Molecular Dynamics mode !" << endl ;
    fflush( 0 ) ;
    exit( EXIT_FAILURE ) ;
  }

  set_saving_options( "EQUILIBRATION" ) ;
  save_simulation_conditions( "MD_simulation_conditions.dat" ) ;

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
    int input_int = 0 ;
    check_action( 0 , starting_time , total_step_div_10 ) ;
    do{
      lists_rebuilds = 0 ;  // number of events of recalculation of the verlet lists
      equilibration_steps += input_int ;              //  if system needs another equilibration run I can specify n

      if( strcmp( simulation_ensemble , "NVE" ) == 0 ) {
	for( current_step=current_step ; current_step<=equilibration_steps ; current_step++ ) {                   // integration steps
	  VelocityVerletIntegrationStep() ;
	  check_action( current_step , starting_time , total_step_div_10 ) ;
	}
      }

      cout << "\tNumber of events of verlet lists recalculation : " << lists_rebuilds << endl ;
      input_int = 0 ;
      if( strcmp( questions, "NO_STOP" ) != 0 ) {
	control = request( "\n Do you want to continue the equilibration stage?" , "y/n" , '/' ) ; // if system needs more time to equilibrate
	if( control == 'y' ) {
	  printf( "  Insert the number of further equilibration steps: _" ) ;
	  command_output = scanf( "%d%*c" , &input_int ) ;
	  if (command_output == 0) cout << "*** Warning: input steps: " << input_int << endl ;
	}
      }

    } while( input_int != 0 ) ;
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
    if( strcmp( questions, "NO_STOP" ) != 0 ) {
      control = request( "\n Do you want to change the number of production steps?" , "y/n" , '/' ) ;
      if( control == 'y' ) {
	printf( "  Insert the number of production steps: _" ) ;
	command_output = scanf( "%d%*c" , &production_steps ) ;
	if (command_output == 0) cout << "*** Warning: input steps: " << production_steps << endl ;
      }
    }

    printf( "\n  PRODUCTION STAGE . . .\n" ) ; fflush(0) ;
    lists_rebuilds = 0 ; // number of events of recalculation of the verlet lists
    initialize_storage_file( production_steps , seed ) ;     // open the file in which configurations are periodically saved
    if( SAVE.backup_interval == 0 ) SAVE.backup_interval = PRODUCTION_STEPS + 10 ;

    check_action( 0 , starting_time , total_step_div_10 ) ;
    if( strcmp( simulation_ensemble , "NVE" ) == 0 ) {
      for( current_step=1 ; current_step<=production_steps ; current_step++ ) {
	VelocityVerletIntegrationStep() ;
	check_action( current_step , starting_time , total_step_div_10 ) ;
      }
    }
  
    cout << "\tNumber of events of verlet lists recalculation : " << lists_rebuilds << endl ;

    cout << endl ;
    fflush( NULL ) ;
    close_storage_file() ;
  }
}
/*****************************************************************************************/

template <typename particle>
inline void configuration<particle>::VelocityVerletIntegrationStep(void) {
  double half_time_step = 0.5*time_step ;
  // positions update
  for( int i=0 ; i<numberOfPart ; i++ ) {
    particles[i]->position += ( ( velocities[i]->position + ( forces_t1[i]->position*half_time_step ) )*time_step ) ;
  }
  check_verlet_update() ;
  // velocities update
  calculate_forces_t2() ;
  for( int i=0 ; i<numberOfPart ; i++ ) {
    velocities[i]->position += ( ( forces_t2[i]->position + forces_t1[i]->position ) * half_time_step ) ;
  }
  support_pp2particle = forces_t1 ;
  forces_t1 = forces_t2 ;
  forces_t2 = support_pp2particle ;
  global_time += time_step ;
}
/****************************************************************************************/
