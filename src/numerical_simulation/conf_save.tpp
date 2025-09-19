/* this file contains methods to initialize output files and to save thermodynamic and per particle data along the simulations */
// MD VERSION
#include "../../lib/system_comm.h"
#include <limits>

template <typename particle>
void configuration<particle>::plot_equilibration_data(void) {
  cout << " ERROR: The function plot_equilibration_data() has still to be implemented!! " << endl ;
  exit( EXIT_FAILURE ) ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::save_options::copy_options( const save_options &other ) {
  PARTICLES.POS = other.PARTICLES.POS ;
  PARTICLES.POS_IMG = other.PARTICLES.POS_IMG ;
  PARTICLES.VEL = other.PARTICLES.VEL ;
  PARTICLES.FORCES = other.PARTICLES.FORCES ;
  PARTICLES.ENERGY = other.PARTICLES.ENERGY ;
  PARTICLES.VIRIAL = other.PARTICLES.VIRIAL ;
  PARTICLES.FORCES_SPLIT = other.PARTICLES.FORCES_SPLIT ;
  PARTICLES.MOL = other.PARTICLES.MOL ;
  PARTICLES.IDX = other.PARTICLES.IDX ;
  PARTICLES.TYPE = other.PARTICLES.TYPE ;
  POTENTIAL_ENERGY = other.POTENTIAL_ENERGY ;
  PARTIAL_POTENTIAL_ENERGY = other.PARTIAL_POTENTIAL_ENERGY ;
  KINETIC_ENERGY = other.KINETIC_ENERGY ;
  SECONDARY_POTENTIAL = other.SECONDARY_POTENTIAL ;
  TEMPERATURE = other.TEMPERATURE ;
  VOLUME = other.VOLUME ;
  MSD = other.MSD ;
  PRECISION10 = other.PRECISION10 ;
  VIRIAL = other.VIRIAL ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::save_options::clear( void ) {
  PARTICLES.POS = 0 ;
  PARTICLES.POS_IMG = 0 ;
  PARTICLES.VEL = 0 ;
  PARTICLES.FORCES = 0 ;
  PARTICLES.ENERGY = 0 ;
  PARTICLES.VIRIAL = 0 ;
  PARTICLES.FORCES_SPLIT = 0 ;
  PARTICLES.MOL = 0 ;
  PARTICLES.IDX = 0 ;
  PARTICLES.TYPE = 0 ;
  POTENTIAL_ENERGY = 0 ;
  PARTIAL_POTENTIAL_ENERGY = 0 ;
  KINETIC_ENERGY = 0 ;
  SECONDARY_POTENTIAL = 0 ;
  TEMPERATURE = 0 ;
  VOLUME = 0 ;
  MSD = 0 ;
  PRECISION10 = 12 ;
  VIRIAL = 0 ;
}
/****************************************************************************************/

template <typename particle>
configuration<particle>::save_options::save_options( void ) {
  PARTICLES.POS = 0 ;
  PARTICLES.POS_IMG = 0 ;
  PARTICLES.VEL = 0 ;
  PARTICLES.FORCES = 0 ;
  PARTICLES.ENERGY = 0 ;
  PARTICLES.VIRIAL = 0 ;
  PARTICLES.FORCES_SPLIT = 0 ;
  PARTICLES.MOL = 0 ;
  PARTICLES.IDX = 0 ;
  PARTICLES.TYPE = 0 ;
  POTENTIAL_ENERGY = 0 ;
  PARTIAL_POTENTIAL_ENERGY = 0 ;
  KINETIC_ENERGY = 0 ;
  SECONDARY_POTENTIAL = 0 ;
  TEMPERATURE = 0 ;
  VOLUME = 0 ;
  MSD = 0 ;
  PRECISION10 = 12 ;
  VIRIAL = 0 ;
}
/****************************************************************************************/

template <typename particle>
char *configuration<particle>::generate_new_directory_for_data(void) {
  int command_output;
  if( strcmp( _DIRECTORY_NAME, "" ) == 0 ) {
    time_t current;
    time(&current);
    struct tm *date_hour = localtime(&current);
    strftime( _DIRECTORY_NAME , 300 , "simulation_%I:%M_%d-%m-%Y" , date_hour );
  }
  char dir_creation[400];
  strcpy( dir_creation , "mkdir " );
  if( strcat( dir_creation , _DIRECTORY_NAME ) == NULL) cout << " Problems in generate_new_directory_for_data()\n";
  command_output = system( dir_creation );
  if (command_output < 0) cout << "*** ERROR: " << dir_creation << endl << "in method generate_new_directory_for_data()" << endl ;

  return _DIRECTORY_NAME;
}
/****************************************************************************************/

template <typename particle>
char *configuration<particle>::generate_equilibration_subdirectory(char *root) {
  int command_output;
  if( strcpy( root , _DIRECTORY_NAME ) == NULL) cout << " Problems in generate_equilibration_subdirectory\n";
  if( strcat( _DIRECTORY_NAME , "/equilibration_run" ) == NULL) cout << " Problems in generate_equilibration_subdirectory()\n";
  char dir_creation[300];
  strcpy( dir_creation , "mkdir " );
  if( strcat( dir_creation , _DIRECTORY_NAME ) == NULL) cout << " Problems in generate_equilibration_subdirectory()\n";
  command_output = system( dir_creation );
  if (command_output < 0) cout << "*** ERROR: " << dir_creation << endl << "in method generate_equilibration_subdirectory()" << endl ;

  return root;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle> :: initialize_storage_file ( int mcsteps , double seed ) {
  char file_path[300] ;
  // initialization of the supporting variables to compute the contribution of the different kind of interactions
  if( SAVE.PARTICLES.FORCES == 1 && SAVE.PARTICLES.FORCES_SPLIT == 1 ) {
    SAVE.force_contributions = new particle ** [ INTERACTIONS.number ] ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
      SAVE.force_contributions[k] = new particle * [ numberOfPart ] ;
      for( int i=0 ; i<numberOfPart ; i++ ) SAVE.force_contributions[k][i] = new particle ;
    }
  } else SAVE.force_contributions = NULL ;
  if( SAVE.PARTICLES.ENERGY == 1 ) {
    SAVE.per_part_energy = new double * [ INTERACTIONS.number+1 ] ;   // the last is the total
    for( int k=0 ; k<INTERACTIONS.number+1 ; k++ ) SAVE.per_part_energy[k] = new double [ numberOfPart ] ;
  } else SAVE.per_part_energy = NULL ;
  if( SAVE.PARTICLES.VIRIAL == 1 ) {
    SAVE.per_part_virial = new double * [ INTERACTIONS.number+1 ] ;   // the last is the total
    for( int k=0 ; k<INTERACTIONS.number+1 ; k++ ) SAVE.per_part_virial[k] = new double [ numberOfPart ] ;
  } else SAVE.per_part_virial = NULL ;

  strcpy( file_path , _DIRECTORY_NAME ) ;
  if( strcat( file_path , "/interactions_data.dat" ) == NULL) cout << " Problems in initialize_storage_file\n" ;
  ofstream _int_data ;
  _int_data.open( file_path ) ;
  _int_data << "total_steps " << mcsteps << " drand48_seed " << seed << endl ;
  for( int i=0 ; i<INTERACTIONS.number ; i++ ) {
    _int_data << endl << "INTERACTION_POTENTIAL " << i << " : cut-off " << INTERACTIONS.type[i]->cutoff << " ; shift " << INTERACTIONS.type[i]->SHIFT << " ; verlet_radii " << INTERACTIONS.type[i]->verlet_list->r_verlet << " , " << INTERACTIONS.type[i]->verlet_list->delta_verlet << " ; list_updated " << INTERACTIONS.type[i]->verlet_list->must_be_updated << " ; list_displacements " << INTERACTIONS.type[i]->verlet_list->disp_on << " ; DOUBLE_PAIRS " << INTERACTIONS.type[i]->verlet_list->DOUBLE_PAIRS << " ." << endl << INTERACTIONS.type[i]->ostr() << endl ;
  }
  _int_data.close() ;

  if( SAVE.PARTICLES.POS == 1 || SAVE.PARTICLES.VEL == 1 || SAVE.PARTICLES.FORCES == 1 || SAVE.PARTICLES.ENERGY == 1 || SAVE.PARTICLES.VIRIAL == 1 ) {
    if( strcmp( SAVE.configuration_format , "xyz" ) == 0 ) {
      strcpy( file_path , _DIRECTORY_NAME ) ;
      if( strcat( file_path , "/particles.dat" ) == NULL) cout << " Problems in initialize_storage_file\n" ;
      _particles_file.open( file_path ) ;
      _particles_file << setprecision(12) ;
      _particles_file << "# Columns_meaning:" ;
      if( SAVE.PARTICLES.POS == 1 ) {
	_particles_file << " X/SIGMA Y" ;
	if( DIMENSION == 3 ) _particles_file << " Z" ;
	if( SAVE.PARTICLES.POS_IMG == 1 ) {
	  _particles_file << " imgX/box_side imgY" ;
	  if( DIMENSION == 3 ) _particles_file << " imgZ" ;
	}
      }
      if( SAVE.PARTICLES.VEL == 1 ) {
	_particles_file << " Vx/sqrt(EPSILON/MASS) Vy" ;
	if( DIMENSION == 3 ) _particles_file << " Vz" ;
      }
      if( SAVE.PARTICLES.FORCES == 1 ) {
	if( SAVE.PARTICLES.FORCES_SPLIT == 1 ) {
	  for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
	    _particles_file << " Fx_" << k << "/(EPSILON/SIGMA) Fy_" << k ;
	    if( DIMENSION == 3 ) _particles_file << " Fz_" << k ;
	  }
	} else if( SAVE.PARTICLES.FORCES_SPLIT == 0 ) {
	  _particles_file << " Fx/(EPSILON/SIGMA) Fy" ;
	  if( DIMENSION == 3 ) _particles_file << " Fz" ;
	}
      }
      if( SAVE.PARTICLES.ENERGY == 1 ) {
	_particles_file << " energy/EPSILON" ;
      }
      if( SAVE.PARTICLES.VIRIAL == 1 ) {
	_particles_file << " virial/EPSILON" ;
      }
      _particles_file << endl << endl ;

    } else if( strcmp( SAVE.configuration_format , "sph" ) == 0 ) {
      strcpy( file_path , _DIRECTORY_NAME ) ;
      if( strcat( file_path , "/particles.sph" ) == NULL ) cout << " Problems in initialize_storage_file\n" ;
      _particles_file.open( file_path ) ;
      _particles_file << setprecision(12) ;
    }
  }


  if( SAVE.POTENTIAL_ENERGY == 1 ||
      SAVE.PARTIAL_POTENTIAL_ENERGY == 1 ||
      SAVE.KINETIC_ENERGY == 1 ||
      SAVE.SECONDARY_POTENTIAL == 1 ||
      SAVE.TEMPERATURE == 1 || SAVE.VOLUME == 1 || SAVE.VIRIAL == 1 ) {

    strcpy( file_path , _DIRECTORY_NAME ) ;
    if( strcat( file_path , "/thermoinfo.dat" ) == NULL) cout << " Problems in initialize_storage_file\n" ;
    _thermoinfo_file.open( file_path ) ;
    _thermoinfo_file << "# step time" ;

    if( SAVE.PARTIAL_POTENTIAL_ENERGY == 1 ) {
      for( int i=0 ; i<INTERACTIONS.number ; i++ ) _thermoinfo_file << " potential_" << i << "(EPSILON)" ;
    }
    if( SAVE.POTENTIAL_ENERGY == 1 ) _thermoinfo_file << " total_potential_en(EPSILON)" ;
    if( SAVE.KINETIC_ENERGY == 1 ) _thermoinfo_file << " kinetic_en(EPSILON)" ;
    if( SAVE.SECONDARY_POTENTIAL == 1 ) _thermoinfo_file << " secondary_potential_energy(EPSILON)" ;
    if( SAVE.TEMPERATURE == 1 ) _thermoinfo_file << " temperature(EPSILON)" ;
    if( SAVE.VOLUME == 1 ) {
      _thermoinfo_file << " boxside_x/(SIGMA) boxside_y" ;
      if( DIMENSION == 3 ) _thermoinfo_file << " boxside_z" ;
      _thermoinfo_file << " volume(SIGMA^" << DIMENSION << ")" ;
    }
    if( SAVE.VIRIAL == 1 ) _thermoinfo_file << " virial_pressure(EPSILON/SIGMA^" << DIMENSION << ")" ;
    _thermoinfo_file << endl ;
  }
  if( SAVE.MSD == 1 ) {
    for( int i=0 ; i<N_MSDs ; i++ ) MSD[i]->initialize_MSDfiles( _DIRECTORY_NAME ) ;
  }
}
/****************************************************************************************/

template <typename particle>
inline void configuration<particle> :: check_action ( int current_step , int starting_time , int screen_print_interval ) {
  if( current_step % screen_print_interval == 0 ) {
    int hour , min , sec , elapsed_time = time(0)-starting_time ;
    sec = elapsed_time % 60 ;
    min = ( (elapsed_time-sec) / 60 ) % 60 ;
    hour = (elapsed_time-sec-60*min) / 3600 ;
    printf( "\r   step %d ; elapsed time: %d:%d:%d", current_step , hour , min , sec ) ;
    fflush(0) ;
  }
  bool SAVE_DATA = ( current_step % SAVE.global_interval == 0 ) ;
  bool DUMP_CONF = ( current_step % SAVE.configuration_interval == 0 ) ;
  bool BACKUP = ( current_step % SAVE.backup_interval == 0 ) ;

  if( DEVICE == GPU ) {
    bool SYNC_positions = 0 ;
    bool SYNC_velocities = 0 ;
    bool SYNC_forces = 0 ;
    if( BACKUP ) SYNC_positions = 1 , SYNC_velocities = 1 ;
    if( DUMP_CONF ) {
      SYNC_positions = ( SYNC_positions || SAVE.PARTICLES.POS || SAVE.PARTICLES.POS_IMG ) ;
      SYNC_velocities = ( SYNC_velocities || SAVE.PARTICLES.VEL ) ;
      SYNC_forces = ( SAVE.PARTICLES.FORCES || SAVE.PARTICLES.FORCES_SPLIT ) ;
    }
    if( SAVE_DATA ) {
      SYNC_positions = ( SYNC_positions ) ;
      SYNC_velocities = ( SYNC_velocities || SAVE.KINETIC_ENERGY || SAVE.TEMPERATURE ) ;
    }
    copyParticlesToHost( SYNC_positions , SYNC_velocities , SYNC_forces ) ;
  }

  if( SAVE_DATA ) save_data( current_step ) ;
  if( DUMP_CONF ) save_configuration( current_step ) ;
  if( BACKUP ) write_backup_data( current_step ) ;
  if( SAVE.MSD == 1 ) {
    for( int i=0 ; i<N_MSDs ; i++ ) if( current_step == MSD[i]->saving_step ) MSD[i]->save_msd_conf( global_time , box_sides ) ;
  }
}
/****************************************************************************************/

template <typename particle>
inline void configuration<particle>::save_configuration( int step ) {
  // Preliminary computations :
  if( SAVE.PARTICLES.FORCES == 1 && SAVE.PARTICLES.FORCES_SPLIT == 1 ) {
    particle **true_forces_t2 = forces_t2 ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
      for( int i=0 ; i<numberOfPart ; i++ ) SAVE.force_contributions[k][i]->clear() ;
      forces_t2 = SAVE.force_contributions[k] ;
      INTERACTIONS.type[k]->compute_forces_contribution() ;
    }
    forces_t2 = true_forces_t2 ;
  }
  if( SAVE.PARTICLES.ENERGY == 1 ) {
    for( int i=0 ; i<numberOfPart ; i++ ) SAVE.per_part_energy[INTERACTIONS.number][i] = 0 ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) INTERACTIONS.type[k]->compute_per_part_energy_contribution( SAVE.per_part_energy[k] ) ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) for( int i=0 ; i<numberOfPart ; i++ ) SAVE.per_part_energy[INTERACTIONS.number][i] += SAVE.per_part_energy[k][i] ;
  }
  if( SAVE.PARTICLES.VIRIAL == 1 ) {
    for( int i=0 ; i<numberOfPart ; i++ ) SAVE.per_part_virial[INTERACTIONS.number][i] = 0 ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) INTERACTIONS.type[k]->compute_per_part_virial_contribution( SAVE.per_part_virial[k] ) ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) for( int i=0 ; i<numberOfPart ; i++ ) SAVE.per_part_virial[INTERACTIONS.number][i] += SAVE.per_part_virial[k][i] ;
  }

  // Saving actions :
  if( strcmp( SAVE.configuration_format , "xyz" ) == 0 ) {
    if( SAVE.PARTICLES.POS == 1 || SAVE.PARTICLES.VEL == 1 || SAVE.PARTICLES.FORCES == 1 || SAVE.PARTICLES.ENERGY == 1 || SAVE.PARTICLES.VIRIAL == 1 ) dump_xyz( _particles_file , 12 ) ;

  } else if ( strcmp( SAVE.configuration_format , "sph" ) == 0 ) {
    if( SAVE.PARTICLES.POS == 1 ) dump_sph( _particles_file , 12 ) ;

  } else if ( strcmp( SAVE.configuration_format , "lmp" ) == 0 ) {
    if( SAVE.PARTICLES.POS == 1 ) dump_lmp( _particles_file , 12 ) ;
  }
}
/****************************************************************************************/

template <typename particle>
inline void configuration<particle>::save_data( int step ) {
  if( SAVE.POTENTIAL_ENERGY == 1 ||
      SAVE.PARTIAL_POTENTIAL_ENERGY == 1 ||
      SAVE.KINETIC_ENERGY == 1 ||
      SAVE.SECONDARY_POTENTIAL == 1 ||
      SAVE.TEMPERATURE == 1 || SAVE.VOLUME == 1 || SAVE.VIRIAL == 1 ) {

    _thermoinfo_file << step << ' ' << global_time ;
    _thermoinfo_file << setprecision(SAVE.PRECISION10) ;
    if( SAVE.PARTIAL_POTENTIAL_ENERGY == 1 ) {
      for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
	if( INTERACTIONS.type[k]->partial_energy.time != global_time ) {
	  INTERACTIONS.type[k]->partial_energy.val = INTERACTIONS.type[k]->compute_energy_contribution() ;
	  INTERACTIONS.type[k]->partial_energy.time = global_time ;
	}
	_thermoinfo_file << ' ' << INTERACTIONS.type[k]->partial_energy.val ;
      }
      if( POTENTIAL_ENERGY.time != global_time ) {
	POTENTIAL_ENERGY.val = 0 ;
	for( int k=0 ; k<INTERACTIONS.number ; k++ ) POTENTIAL_ENERGY.val += INTERACTIONS.type[k]->partial_energy.val ;
	POTENTIAL_ENERGY.time = global_time ;
      }
    }
    if( SAVE.POTENTIAL_ENERGY == 1 ) _thermoinfo_file << ' ' << POTENTIAL_ENERGY.val ;
    if( SAVE.KINETIC_ENERGY == 1 ) _thermoinfo_file << ' ' << compute_kinetic_energy() ;
    if( SAVE.SECONDARY_POTENTIAL == 1 ) _thermoinfo_file << ' ' << calculate_secondary_potential() ;
    if( SAVE.TEMPERATURE == 1 ) _thermoinfo_file << ' ' << ( 2. * compute_kinetic_energy() / DIMENSION / numberOfPart ) ;
    if( SAVE.VOLUME == 1 ) _thermoinfo_file << ' ' << box_sides.position << ' ' << box_sides.volume() ;
    if( SAVE.VIRIAL == 1 ) _thermoinfo_file << ' ' << compute_virial_pressure() ;
    _thermoinfo_file << endl ;
  }
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::close_storage_file(void) {
  if( SAVE.PARTICLES.POS == 1 || SAVE.PARTICLES.VEL == 1 || SAVE.PARTICLES.FORCES == 1 || SAVE.PARTICLES.ENERGY == 1 || SAVE.PARTICLES.VIRIAL == 1 ) {
    _particles_file.flush() ;
    _particles_file.close() ;
  }
  if( SAVE.POTENTIAL_ENERGY == 1 ||
      SAVE.PARTIAL_POTENTIAL_ENERGY == 1 ||
      SAVE.KINETIC_ENERGY == 1 ||
      SAVE.SECONDARY_POTENTIAL == 1 ||
      SAVE.TEMPERATURE == 1 || SAVE.VOLUME == 1 || SAVE.VIRIAL == 1 ) {
    _thermoinfo_file.flush() ;
    _thermoinfo_file.close() ;
  }
  if( SAVE.MSD == 1 ) {
    for( int i=0 ; i<N_MSDs ; i++ ) {
      MSD[i]->_MSD_file.flush() ;
      MSD[i]->_MSD_file.close() ;
      MSD[i]->_saving_steps_file.close() ;
    }
  }
  if( SAVE.PARTICLES.FORCES == 1 && SAVE.PARTICLES.FORCES_SPLIT == 1 ) {
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
      for( int i=0 ; i<numberOfPart ; i++ ) delete SAVE.force_contributions[k][i] ;
      delete [] SAVE.force_contributions[k] ;
    }
    delete [] SAVE.force_contributions ;
  }
  if( SAVE.PARTICLES.ENERGY == 1 ) {
    for( int k=0 ; k<INTERACTIONS.number+1 ; k++ ) delete [] SAVE.per_part_energy[k] ;
    delete [] SAVE.per_part_energy ;
  }
  if( SAVE.PARTICLES.VIRIAL == 1 ) {
    for( int k=0 ; k<INTERACTIONS.number+1 ; k++ ) delete [] SAVE.per_part_virial[k] ;
    delete [] SAVE.per_part_virial ;
  }
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::set_saving_options( const char *MODE ) {
  static save_options OPTIONS ;

  if( strcmp( MODE, "OFF" ) == 0 ) {
    OPTIONS.copy_options( SAVE ) ;
    SAVE.clear() ;
  }else if( strcmp( MODE, "ON" ) == 0 ) {
    SAVE.copy_options( OPTIONS ) ;
  }else if( strcmp( MODE, "EQUILIBRATION" ) == 0 ) {
    SAVE.copy_options( EQ_SAVE_OPT ) ;
  }else if( strcmp( MODE, "PRODUCTION" ) == 0 ) {
    SAVE.copy_options( PROD_SAVE_OPT ) ;
  }else{
    cout << "ERROR: incorrect input for set_saving_options()" << endl;
  }
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::save_simulation_conditions( const char *file_name ) {
  int path_length = 0 ;
  while( _DIRECTORY_NAME[path_length] != '\0' ) path_length++ ;
  while( file_name[path_length] != '\0' ) path_length++ ;
  path_length += 200 ;
  char *command = (char *)calloc( path_length , sizeof(char) ) ;
  strcpy( command , _DIRECTORY_NAME ) ;
  strcat( command , "/" ) ;
  strcat( command , file_name ) ;
  ofstream conditions_file ;
  conditions_file.open( command ) ;
  free( command ) ;

  conditions_file << "DEVICE " << ( (DEVICE == GPU) ? "gpu" : "cpu" ) << endl ;
  conditions_file << "SIMULATION_MODE " << simulation_mode << endl ;
  conditions_file << "ENSEMBLE " << simulation_ensemble ;
  conditions_file << endl ;
  if( strcmp( simulation_ensemble, "NVE" ) == 0 ) {
    conditions_file << "ENERGY_PER_PARTICLE " << energyOfSystem/numberOfPart << endl ;
  } else {
    conditions_file << "KbT " << k_T << endl ;
  }
  conditions_file << "SEED " << seed << endl ;
  conditions_file << "SIGMA " << INTERACTIONS.SIGMA << endl ;
  conditions_file << "EPSILON " << INTERACTIONS.EPSILON << endl ;
  if( particle_mass > 0 ) conditions_file << "PARTICLE_MASS " << particle_mass << endl ;
  else conditions_file << "PARTICLE_MASS either 1 or in the initial configuration file" << endl ;
  conditions_file << "SIDES " << box_sides.position << endl ;
  conditions_file << "INTEGRATION_TIME_STEP " << time_step << endl ;
  conditions_file << "EQUILIBRATION_STEPS " << EQUILIBRATION_STEPS << endl ;
  conditions_file << "PRODUCTION_STEPS " << PRODUCTION_STEPS << endl ;
  conditions_file << "NUMBER_OF_PARTICLES " << numberOfPart << endl ;
  conditions_file << "SAVING_INTERVALS " << SAVE.global_interval << ' ' << SAVE.configuration_interval << ' ' << SAVE.backup_interval << endl ;
  if( INITIAL_CONF_BY_FILE == 1 ) conditions_file << "INITIAL_CONFIGURATION_BY_FILE " << configuration_file_name << endl ;
  if( INITIAL_VELOCITIES_BY_FILE == 1 ) conditions_file << "INITIAL_VELOCITIES_BY_FILE " << velocities_file_name << endl ;
  conditions_file << "____ POTENTIALS CONSTANTS " << endl ;
  conditions_file << "SIGMA : " << INTERACTIONS.SIGMA << endl ;
  conditions_file << "EPSILON : " << INTERACTIONS.EPSILON << endl ;
  conditions_file << INTERACTIONS << endl ;

  conditions_file.close() ;
}
/****************************************************************************************/

template <typename particle>
inline void configuration<particle>::write_backup_data( int current_step ) {
  char file_path[300] ;
  // saving of positions
  strcpy( file_path , _DIRECTORY_NAME ) ;
  strcat( file_path , "/positions_backup.dat" ) ;
  _positions_backup.open( file_path ) ;

  _positions_backup << setprecision(std::numeric_limits<double>::digits10 + 1) ;
  _positions_backup.seekp( 0 ) ;
  _positions_backup << "# step " << current_step << endl ;
  _positions_backup << "# N " << numberOfPart << endl ;
  _positions_backup << "# box " << box_sides.position << endl ;
  if ( strcmp( SAVE.configuration_format , "sph" ) == 0 ) {
    _positions_backup << "# format sph" << endl ;
    _positions_backup << numberOfPart << " " << global_time << endl ;
    _positions_backup << setprecision(std::numeric_limits<double>::digits10 + 1) << box_sides.position << endl ;
    if( DIMENSION ==3 ) for( int i=0 ; i<numberOfPart ; i++ ) _positions_backup << (char)(particles[i]->type-1+'a') << " " << particles[i]->dump_position_variables( std::numeric_limits<double>::digits10 + 1 , "sph" , NULL ) << " " << particles[i]->radius << endl ;
    else if( DIMENSION == 2 ) for( int i=0 ; i<numberOfPart ; i++ ) _positions_backup << (char)(particles[i]->type-1+'a') << " " << particles[i]->dump_position_variables( std::numeric_limits<double>::digits10 + 1 , "sph" , NULL ) << " 0.0 " << particles[i]->radius << endl ;

  } else {
    for(int i=0; i<numberOfPart; i++) {
      _positions_backup << particles[i]->dump_position_variables( std::numeric_limits<double>::digits10 + 1 , "all" ) << " " << particles[i]->periodic_box << endl ;
    }
  }

  // saving of velocities
  if( strcmp( simulation_mode , "MD" ) == 0 ) {
    if( velocities != NULL ) {
      strcpy( file_path , _DIRECTORY_NAME ) ;
      strcat( file_path , "/velocities_backup.dat" ) ;
      _velocities_backup.open( file_path ) ;

      _velocities_backup << setprecision(std::numeric_limits<double>::digits10 + 1) ;
      _velocities_backup.seekp( 0 ) ;
      _velocities_backup << "# N " << numberOfPart << endl ;
      _velocities_backup << "# step " << current_step << endl ;
      for(int i=0; i<numberOfPart; i++) {
	_velocities_backup << velocities[i]->position << endl ;
      }
      _velocities_backup << endl ;
      _velocities_backup.close() ;
    } else {
      cout << " *** ERROR: The velocities vector is not allocated, but simulation mode is MD !" << endl ;
      exit( EXIT_FAILURE ) ;
    }
  }

  // saving of information on bonds list and charges
  for( int i=0; i<numberOfPart; i++ ) {
    _positions_backup << i+1 << ' ' << particles[i]->valence << ' ' << particles[i]->charge << endl ;
    if( particles[i]->valence > 0 ) {
      for( int j=0 ; j<particles[i]->valence-1 ; j++ )  _positions_backup << particles[i]->bonded_monomer[j] + 1 << ' ' ;
      _positions_backup << particles[i]->bonded_monomer[ particles[i]->valence-1 ] + 1 << endl ;
    }
  }
  _positions_backup << endl ;

  _positions_backup.close() ;
}
/****************************************************************************************/
