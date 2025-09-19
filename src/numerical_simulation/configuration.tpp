#include "interactions.h"
#include "msd.h"
#include <climits>
#include <limits>

/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/
template <typename particle>
configuration<particle>::interactions_array::interactions_array() {
  number = 0 ;
  type = NULL ;
  SIGMA = 0.0 , EPSILON = 0.0 ;
}
/********************************************************************************************************************/

template <typename particle>
configuration<particle>::interactions_array::~interactions_array() {
  for( int i=0 ; i<number ; i++ ) delete type[i] ;
  free( type ) ;
  type = NULL ;
  number = 0 ;
}
/********************************************************************************************************************/

template <typename particle>
void configuration<particle>::interactions_array::operator=( const interactions_array &other ) {
  number = other.number ;
  type = (interaction<particle> **)calloc( number , sizeof(interaction<particle> *) ) ;
  for(int i=0; i<number; i++) {
    type[i] = other.type[i]->copy() ;
    type[i]->cutoff = other.type[i]->cutoff ;
    type[i]->SHIFT = other.type[i]->SHIFT ;
    type[i]->verlet_list = new verlet<particle>( *( other.type[i]->verlet_list ) ) ;     // copy of verlet list
  }
}
/********************************************************************************************************************/

/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/
template <typename particle>
configuration<particle>::configuration() {
  strcpy( simulation_ensemble, "0" ) ;
  strcpy( simulation_mode , "00" ) ;
  seed = -1 ;
  DIMENSION = -1 ;
  k_T = 1.0 ;
  energyOfSystem = 0.0 ;
  time_step = 1 ;   // The time-step is initialized to 1
  global_time = 0 ;
  EQUILIBRATION_STEPS = 0 ;
  PRODUCTION_STEPS = 0 ;
  box_side = 0 ;
  box_sides = particle() ;
  midside = particle() ;
  // default options saves only configurations and energy
  particle_mass = 0.0 ;
  SAVE.clear() ;
  SAVE.PARTICLES.POS = 1 ;
  SAVE.PARTICLES.POS_IMG = 1 ;
  SAVE.POTENTIAL_ENERGY = 1 ;
  SAVE.configuration_format = new char [4] ;
  sprintf( SAVE.configuration_format , "%s" , "xyz" ) ;
  numberOfPart = 0 ;
  particles_length = 0 ;
  particles = NULL ;
  velocities = NULL ;
  forces_t1 = NULL ;
  forces_t2 = NULL ;
  monomers_number = 0 ;
  strcpy( _DIRECTORY_NAME , "" ) ;
  inertia_matrix.rebuild( 3 , 3 ) ;
  N_MSDs = 0 ;
  MSD = NULL ;
  // GPU data structs
  GPU_arrays.d_idx = NULL ;
  GPU_arrays.d_x = NULL ;
  GPU_arrays.d_y = NULL ;
  GPU_arrays.d_z = NULL ;
  GPU_arrays.d_ix = NULL ;
  GPU_arrays.d_iy = NULL ;
  GPU_arrays.d_iz = NULL ;
  GPU_arrays.d_vx = NULL ;
  GPU_arrays.d_vy = NULL ;
  GPU_arrays.d_vz = NULL ;
  GPU_arrays.d_fx = NULL ;
  GPU_arrays.d_fy = NULL ;
  GPU_arrays.d_fz = NULL ;
  GPU_arrays.d_compute_vec = NULL ;
  GPU_cells.d_sorted_int = NULL ;
  GPU_cells.d_sorted_double = NULL ;
  GPU_cells.m_dCellStart = NULL ;
  GPU_cells.m_dCellEnd = NULL ;
  GPU_cells.m_dCellParticleIndex = NULL ;
  GPU_cells.m_dCellParticleHash = NULL ;
}
/********************************************************************************************************************/

template <typename particle>
configuration<particle>::configuration( const configuration &other_config ) {
  SAVE.copy_options( other_config.SAVE ) ;
  SAVE.global_interval = other_config.SAVE.global_interval ;
  SAVE.configuration_interval = other_config.SAVE.configuration_interval ;
  SAVE.backup_interval = other_config.SAVE.backup_interval ;
  sprintf( SAVE.configuration_format , "%s" , other_config.SAVE.configuration_format ) ;
  INTERACTIONS.SIGMA = other_config.INTERACTIONS.SIGMA ;
  INTERACTIONS.EPSILON = other_config.INTERACTIONS.EPSILON ;
  box_side = other_config.box_side ;
  box_sides = other_config.box_sides ;
  midside = other_config.midside ;
  k_T = other_config.k_T ;
  energyOfSystem = other_config.energyOfSystem ;
  POTENTIAL_ENERGY = other_config.POTENTIAL_ENERGY ;
  KINETIC_ENERGY = other_config.KINETIC_ENERGY ;
  time_step = other_config.time_step ;
  global_time = other_config.global_time ;
  particle_mass = other_config.particle_mass ;
  EQUILIBRATION_STEPS = other_config.EQUILIBRATION_STEPS ;
  PRODUCTION_STEPS = other_config.PRODUCTION_STEPS ;
  numberOfPart = other_config.numberOfPart ;
  monomers_number = other_config.monomers_number ;
  virial_press = other_config.virial_press ;
  INTERACTIONS = other_config.INTERACTIONS ; // assignment by copy
  N_MSDs = 0 ;
  MSD = NULL ;

  particles_length = numberOfPart ;
  particles = (particle **)calloc( particles_length , sizeof(particle *) ) ;
  if( particles == NULL ) {
    cout << "\n Allocation error 1 in configuration constructor" ;
    exit(EXIT_FAILURE) ;
  }
  forces_t1 = (particle **)calloc( particles_length , sizeof(particle *) ) ;
  if( forces_t1 == NULL ) {
    cout << "\n Allocation error 1.4 in configuration constructor" ;
    exit( EXIT_FAILURE ) ;
  }
  forces_t2 = (particle **)calloc( particles_length , sizeof(particle *) ) ;
  if( forces_t2 == NULL ) {
    cout << "\n Allocation error 1.5 in configuration constructor" ;
    exit( EXIT_FAILURE ) ;
  }
  // copy of positions and forces
  for( int i=0 ; i<numberOfPart ; i++ ) {
    particles[i] = new particle ;
    *( particles[i] ) = *( other_config.particles[i] ) ;
    if( other_config.forces_t1[i] != NULL ) {
      forces_t1[i] = new particle ;
      *( forces_t1[i] ) = *( other_config.forces_t1[i] ) ;
    }
    if( other_config.forces_t2[i] != NULL ) {
      forces_t2[i] = new particle ;
      *( forces_t2[i] ) = *( other_config.forces_t2[i] ) ;
    }
  }
  // copy of velocities
  if( other_config.velocities != NULL ) {
    velocities = (particle **)calloc( particles_length , sizeof(particle *) ) ;
    if( velocities == NULL ) {
      cout << "\n Allocation error 1.3 in configuration constructor" ;
      exit( EXIT_FAILURE ) ;
    }
    for( int i=0 ; i<numberOfPart ; i++ ) {
      velocities[i] = new particle ;
      *( velocities[i] ) = *( other_config.velocities[i] ) ;
    }
  }

}
/********************************************************************************************************************/

template <typename particle>
configuration< particle >::~configuration() {
  // deletion of the array of interactions and verlet lists is automatically done within the destructor of INTERACTIONS

          // deletion of the arrays of particles, velocities and forces
  for( int i=0; i<particles_length; i++ ) {
    delete particles[i] ;
    particles[i] = NULL ;
  }
  free( particles ) ;
  particles = NULL ;
  if( velocities != NULL ) {
    for( int i=0 ; i<numberOfPart ; i++ ) {
      delete velocities[i] ;
      velocities[i] = NULL ;
      delete forces_t1[i] ;
      forces_t1[i] = NULL ;
      delete forces_t2[i] ;
      forces_t2[i] = NULL ;
    }
    free( velocities ) ;
    velocities = NULL ;
    free( forces_t1 ) ;
    forces_t1 = NULL ;
    free( forces_t2 ) ;
    forces_t2 = NULL ;
  }

  particles_length = 0 ;
  numberOfPart = 0 ;
  monomers_number = 0 ;

          // deletion of the MSD structures
  if( N_MSDs != 0 ) {
    for( int i=0 ; i<N_MSDs ; i++ ) delete MSD[i] ;
    delete MSD ;
    MSD = NULL ;
    N_MSDs = 0 ;
  }

  // free mem on GPU
  if( DEVICE == GPU ) {
    freeDeviceArrays() ;
    freeDeviceCellList() ;
  }
}
/****************************************************************************************/
/****************************************************************************************/

/***********************    BLOCK 1   ***************************************************/
/****************************************************************************************/
template <typename particle>
void configuration<particle>::initialize_configuration( const char *input_file_name ) {
  ifstream config_file ;
  if( input_file_name == NULL ) {
    config_file.open( "input_file.dat" ) ;
  }else{
    config_file.open( input_file_name ) ;
  }
  if( config_file.fail() ) {
    cout << "\n  MD_config_file non trovato, apertura fallita" << endl ;
    exit( EXIT_FAILURE ) ;

  } else {
    interpret_input_script( config_file ) ;
    config_file.close() ;
  }

  if( strcmp( simulation_mode , "MD" ) == 0 ) {
    // generation of the velocities
    if( INITIAL_VELOCITIES_BY_FILE == 1 ) generate_velocities_by_file() ;
    else generate_MBvelocities( k_T ) ;

    // generation of the array containing forces
    forces_t1 = new particle* [ particles_length ] ;
    check_memalloc( forces_t1 , "\n Allocation error for forces_t1 in configuration constructor" ) ;
    forces_t2 = new particle* [ particles_length ] ;
    check_memalloc( forces_t2 , "\n Allocation error for forces_t2 in configuration constructor" ) ;
    for( int i=0 ; i<numberOfPart ; i++ ) {
      forces_t1[i] = new particle ;
      forces_t2[i] = new particle ;
    }
  }

  if( DEVICE == GPU )  initialize_GPU_data() ;

  compute_potential_energy() ;           // initial potential energy calculation

  // If MD NVE the kinetic energy is changed in order to fix the initial total energy
  if( strcmp( simulation_mode , "MD" ) == 0 ) {
    if( strcmp( simulation_ensemble, "NVE" ) == 0 ) {
      KINETIC_ENERGY.val = energyOfSystem - POTENTIAL_ENERGY.val ;
      if( KINETIC_ENERGY.val < 0 ) {
	cout << endl << "***ERROR: This fixed energy initial configuration has a negative kinetic energy" << endl;
	exit( EXIT_FAILURE ) ;
      }
      // mean square velocity will be set to fit the initial total energy, rescaling velocities
      set_kinetic_energy( KINETIC_ENERGY.val ) ;
    }
    CoM_velocity_reset() ;
    //    total_angularMomentum_reset() ;
    cout << "  Velocity of center of mass after reset (System units): " << compute_mean_velocity().position << endl << endl ;
    if( DEVICE == GPU )  copyParticlesToDevice( 0 , 1 , 0 ) ;
  }

  // initial forces are computed and stored in forces_t1 to start velocity verlet
  if( strcmp( simulation_mode , "MD" ) == 0 ) {
    if( DEVICE == GPU )  calculate_forces_gpu() ;
    else {
      calculate_forces_t2() ;
      support_pp2particle = forces_t1 ;
      forces_t1 = forces_t2 ;
      forces_t2 = support_pp2particle ;
    }
  }

  // creation of the folder containing all saved data
  if( generate_new_directory_for_data() == NULL ) cout << " Problems in configuration constructor\n" ;
}
/****************************************************************************************/

template <typename particle>
bool configuration<particle>::check_init_conf_set( void ) {
  if( numberOfPart == 0 || ( box_side == 0 && box_sides.volume() == 0 ) ) {
    cout << "ERROR: You must generate the initial configuration and the simulation box before specifying the interaction potential !" << endl ;
    cout << numberOfPart << " " << box_side << " " << box_sides.position << endl ;
    exit( EXIT_FAILURE ) ;
  }
  else return 1 ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::interpret_input_script( ifstream &config_file ) {
  double density = -1.0 ;
  record *file_rec = NULL ;
  file_rec = new record ;
  INITIAL_VELOCITIES_BY_FILE = 0 ;

  while( file_rec->getrecord( config_file ) ) {
    file_rec->split() ;
    if( file_rec->words_number > 0 ) {
      if( file_rec->word[0][0] != '#' ) {

	if( strcmp( file_rec->word[0] , "DIMENSION" ) == 0 ) sscanf( file_rec->word[1] , "%d" , &DIMENSION ) ;
	else if( strcmp( file_rec->word[0] , "ENSEMBLE" ) == 0 ) {
	  strcpy( simulation_ensemble , file_rec->word[1] ) ;
	  if( strcmp( simulation_ensemble , "NVE" ) != 0 ) {
	    cout << " *** ERROR : The simulation ensemble " << simulation_ensemble << " do not exist or was not yet coded !" << endl ;
	    exit( EXIT_FAILURE ) ;
	  }
	  strcpy( simulation_mode , "MD" ) ;
	}

	else if( strcmp( file_rec->word[0] , "DEVICE" ) == 0 ) {
	  if( strcmp( file_rec->word[1] , "GPU" ) == 0 ||
	      strcmp( file_rec->word[1] , "gpu" ) == 0 ||
	      strcmp( file_rec->word[1] , "1" ) == 0 ) DEVICE = GPU ;
	  else if( strcmp( file_rec->word[1] , "CPU" ) == 0 ||
		   strcmp( file_rec->word[1] , "cpu" ) == 0 ||
		   strcmp( file_rec->word[1] , "0" ) == 0 ) DEVICE = CPU ;
	  if( DEVICE == GPU ) GPU_query( &GPU_PROPS ) ;
	}

	else if( strcmp( file_rec->word[0] , "SIGMA" ) == 0 ) sscanf( file_rec->word[1] , "%lf" , &INTERACTIONS.SIGMA ) ;
	else if( strcmp( file_rec->word[0] , "EPSILON" ) == 0 ) sscanf( file_rec->word[1] , "%lf" , &INTERACTIONS.EPSILON ) ;
	else if( strcmp( file_rec->word[0] , "PARTICLE_MASS" ) == 0 ) sscanf( file_rec->word[1] , "%lf" , &particle_mass ) ;
	else if( strcmp( file_rec->word[0] , "SIDE") == 0 ) {
	  sscanf( file_rec->word[1] , "%lf" , &box_side ) ;
	  box_sides.set_equal_comp( box_side ) ;
	  midside.position = (box_sides.position * 0.5) ;
	}
	else if( strcmp( file_rec->word[0], "SIDES") == 0 ) {
	  box_sides.read_by_string( file_rec->word+1 ) ;
	  midside.position = (box_sides.position * 0.5) ;
	  if( box_sides.min_comp() == box_sides.max_comp() ) {
	    box_side = box_sides.min_comp() ;
	  }else{
	    box_side = nanl("") ;
	  }
	}
	else if( strcmp( file_rec->word[0], "KbT" ) == 0 ) sscanf( file_rec->word[1] , "%lf" , &k_T ) ;
	else if( strcmp( file_rec->word[0], "DENSITY" ) == 0 ) sscanf( file_rec->word[1] , "%lf" , &density ) ;

	else if( strcmp( file_rec->word[0], "INITIAL_CONFIGURATION" ) == 0 ) {
	  if( strcmp( file_rec->word[1], "file" ) == 0 ) {
	    sscanf( file_rec->word[3] , "%s" , configuration_file_name ) ;
	              // positions at time 0 are generated, giving the file name (configuration file) and the file format
	    generate_config_by_file( configuration_file_name , file_rec->word[2] ) ;
	    if( strcmp( file_rec->word[2] , "xyz" ) == 0 ) {
	      if( file_rec->words_number > 4 ) {
		sscanf( file_rec->word[4] , "%lf" , &energyOfSystem ) ;
		energyOfSystem *= (double)numberOfPart ;             // The input value is the total energy per particle
	      }else{
		cout << "ERROR: INITIAL_CONFIGURATION needs a value for the initial energy per particle !" << endl ;
		exit( EXIT_FAILURE ) ;
	      }
	    }

	  } else {
	    if( file_rec->words_number > 3 ) {
	      // generation of positions' array
	      sscanf( file_rec->word[2] , "%d" , &numberOfPart ) ;
	      particles_length = numberOfPart ;
	      particles = (particle **)calloc( particles_length , sizeof(particle *) ) ;
	      check_memalloc( particles , "\n Allocation error 1 during the generation of the initial configuration" ) ;
	      if( box_side == 0 && box_sides.volume() == 0 ) {
		if( density > 0 ) {
		  box_side = pow( (double)numberOfPart/density , 1./DIMENSION ) ;
		  box_sides.set_equal_comp( box_side ) ;
		  midside.position = (box_sides.position * 0.5) ;
		} else {
		  cout << "*** ERROR (initial configuration): Please, define the box size or the particles density!" << endl ;
		  exit( EXIT_FAILURE ) ;
		}
	      }
	      if( strcmp( file_rec->word[1] , "random" ) == 0 ) {
		sscanf( file_rec->word[3] , "%lf" , &energyOfSystem ) ;
		energyOfSystem *= (double)numberOfPart ;             // The input value is the total energy per particle
		generate_ControlledEnergy_randomConfig() ;
	      } else if( strcmp( file_rec->word[1] , "fcc" ) == 0 || strcmp( file_rec->word[1] , "bcc" ) == 0 ) {
		sscanf( file_rec->word[3] , "%lf" , &energyOfSystem ) ;
		energyOfSystem *= (double)numberOfPart ;             // The input value is the total energy per particle
		generate_ideal_lattice_3D( file_rec->word[1] ) ;
	      } else if( strcmp( file_rec->word[1] , "hardspheres" ) == 0 ) {
		if( file_rec->words_number < 5 ) { cout << "*** ERROR : invalid number of options for INITIAL_CONFIGURATION" << endl ; exit(EXIT_FAILURE) ; }
		sscanf( file_rec->word[3] , "%lf" , &energyOfSystem ) ;
		energyOfSystem *= (double)numberOfPart ;             // The input value is the total energy per particle
		double radius ;
		sscanf( file_rec->word[4] , "%lf" , &radius ) ;
		generate_HS_randomConfig( radius ) ;
	      } else if( strcmp( file_rec->word[1] , "lattice" ) == 0 ) {
		generate_ideal_lattice( file_rec->word[3] ) ;
	      } else {
		cout << endl << " You need to specify the configuration arrangement in the initialization in function main!" << endl;
	      }
	    } else {
	      cout << "ERROR: INITIAL_CONFIGURATION needs the values for the initial energy per particle and number of particles !" << endl ;
	      exit( EXIT_FAILURE ) ;
	    }
	  }
	}

	else if( strcmp( file_rec->word[0], "LENNARD_JONES" ) == 0 ) {
	  if( check_init_conf_set() ) {
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	    // pair potential definition

	    if( DEVICE == GPU ) {
	      INTERACTIONS.type[INTERACTIONS.number-1] = new lennard_jones_gpu<particle> ;
	      INTERACTIONS.type[INTERACTIONS.number-1]->ON_GPU = 1 ;
	    }else{
	      INTERACTIONS.type[INTERACTIONS.number-1] = new lennard_jones<particle> ;
	    }
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	    // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , numberOfPart ) ;
	  }
	}

	else if( strcmp( file_rec->word[0], "LENNARD_JONES_POLY" ) == 0 ) {
	  if( check_init_conf_set() ) {
	    INTERACTIONS.number ++ ;
	    INTERACTIONS.type = (interaction<particle> **)realloc( INTERACTIONS.type , INTERACTIONS.number*sizeof(interaction<particle> *) ) ;
	              // pair potential definition
	    if( DEVICE == GPU ) {
	      cout << "*** polydisperse Lennard-Jones was not yet coded for GPU!" << endl ;
	    }
	    INTERACTIONS.type[INTERACTIONS.number-1] = new lennard_jones_poly<particle> ;
	    INTERACTIONS.type[INTERACTIONS.number-1]->idx = INTERACTIONS.number-1 ;
	              // read-in of the constants, generation and initialization of the verlet list
	    INTERACTIONS.type[INTERACTIONS.number-1]->generate_by_record( file_rec , simulation_mode , this , numberOfPart ) ;
	  }
	}

	else if( strcmp( file_rec->word[0], "SECONDARY_POTENTIAL" ) == 0 ) {
	  cout << "ERROR: You still have to write the part of code about the computation of the secondary potential . . ." << endl ;
	  exit( EXIT_FAILURE ) ;
	    //	    sscanf( file_rec->word[1] , "%lf" , &SECONDARY_POTENTIAL_SHIFT ) ;
	}
	else if( strcmp( file_rec->word[0], "INTEGRATION_TIME_STEP" ) == 0 ) sscanf( file_rec->word[1] , "%lf" , &time_step ) ;
	else if( strcmp( file_rec->word[0], "EQUILIBRATION_STEPS" ) == 0 ) sscanf( file_rec->word[1] , "%d" , &EQUILIBRATION_STEPS ) ;
	else if( strcmp( file_rec->word[0], "PRODUCTION_STEPS" ) == 0 ) sscanf( file_rec->word[1] , "%d" , &PRODUCTION_STEPS ) ;
	else if( strcmp( file_rec->word[0], "SAVING_INTERVALS" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &SAVE.global_interval ) ;
	  sscanf( file_rec->word[2] , "%d" , &SAVE.configuration_interval ) ;
	  sscanf( file_rec->word[3] , "%d" , &SAVE.backup_interval ) ;
	}

	else if( strcmp( file_rec->word[0], "SAVE_FORMAT" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%s" , SAVE.configuration_format ) ;
	  if( strcmp( SAVE.configuration_format , "xyz" ) != 0 &&
	      strcmp( SAVE.configuration_format , "sph" ) != 0 ) {
	    cout << "*** ERROR : Unrecognized output format !" << endl ;
	    exit(EXIT_FAILURE) ;
	  }
	}
	else if( strcmp( file_rec->word[0], "SAVE_ATTRIBUTE" ) == 0 ) {
	  if( strcmp( file_rec->word[1], "particles" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.PARTICLES.POS ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.PARTICLES.POS ) ;
	  } else if( strcmp( file_rec->word[1], "part_img" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.PARTICLES.POS_IMG ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.PARTICLES.POS_IMG ) ;
	  } else if( strcmp( file_rec->word[1], "velocities" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.PARTICLES.VEL ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.PARTICLES.VEL ) ;
	  } else if( strcmp( file_rec->word[1], "forces" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.PARTICLES.FORCES ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.PARTICLES.FORCES ) ;
	  } else if( strcmp( file_rec->word[1], "forces_splitted" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.PARTICLES.FORCES ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.PARTICLES.FORCES ) ;
	    EQ_SAVE_OPT.PARTICLES.FORCES_SPLIT = EQ_SAVE_OPT.PARTICLES.FORCES ;
	    PROD_SAVE_OPT.PARTICLES.FORCES_SPLIT = PROD_SAVE_OPT.PARTICLES.FORCES ;
	  } else if( strcmp( file_rec->word[1], "energy" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.PARTICLES.ENERGY ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.PARTICLES.ENERGY ) ;
	  } else if( strcmp( file_rec->word[1], "virial" ) == 0 ) {
	    sscanf( file_rec->word[2] , "%d" , &EQ_SAVE_OPT.PARTICLES.VIRIAL ) ;
	    sscanf( file_rec->word[3] , "%d" , &PROD_SAVE_OPT.PARTICLES.VIRIAL ) ;
	  }
	}

	else if( strcmp( file_rec->word[0], "SAVE_PRECISION" ) == 0 ) {
	  int max_precision = std::numeric_limits<double>::digits10 ;
	  cout << " Size of double (bytes): " << sizeof(double) << endl ;
	  cout << " Max precision of double (decimal digits): " << max_precision << endl ;
	  if( strcmp( file_rec->word[1], "max" ) == 0 || strcmp( file_rec->word[1], "MAX" ) == 0 ) EQ_SAVE_OPT.PRECISION10 = max_precision + 1 ;
	  else sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.PRECISION10 ) ;
	  if( strcmp( file_rec->word[2], "max" ) == 0 || strcmp( file_rec->word[2], "MAX" ) == 0 ) PROD_SAVE_OPT.PRECISION10 = max_precision + 1 ;
	  else sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.PRECISION10 ) ;
	  cout << " Precision of Thermodynamic output (decimal digits): " << EQ_SAVE_OPT.PRECISION10 << " , " << PROD_SAVE_OPT.PRECISION10 << endl ;
	}
	else if( strcmp( file_rec->word[0], "SAVE_POTENTIAL_ENERGY" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.POTENTIAL_ENERGY ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.POTENTIAL_ENERGY ) ;
	}
	else if( strcmp( file_rec->word[0], "SAVE_PARTIAL_POTENTIAL_ENERGIES" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.PARTIAL_POTENTIAL_ENERGY ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.PARTIAL_POTENTIAL_ENERGY ) ;
	}
	else if( strcmp( file_rec->word[0], "SAVE_KINETIC_ENERGY" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.KINETIC_ENERGY ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.KINETIC_ENERGY ) ;
	}
	else if( strcmp( file_rec->word[0], "SAVE_SECONDARY_POTENTIAL" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.SECONDARY_POTENTIAL ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.SECONDARY_POTENTIAL ) ;
	}
	else if( strcmp( file_rec->word[0], "SAVE_TEMPERATURE" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.TEMPERATURE ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.TEMPERATURE ) ;
	}
	else if( strcmp( file_rec->word[0], "SAVE_VOLUME" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.VOLUME ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.VOLUME ) ;
	}
	else if( strcmp( file_rec->word[0], "SAVE_MSD" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.MSD ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.MSD ) ;
	  if( file_rec->words_number == 4 ) add_initialize_MSD( file_rec->word[3] ) ;
	  else add_initialize_MSD( file_rec->word[3] , file_rec->word[4] ) ;
	}
	else if( strcmp( file_rec->word[0], "SAVE_VIRIAL_PRESSURE" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &EQ_SAVE_OPT.VIRIAL ) ;
	  sscanf( file_rec->word[2] , "%d" , &PROD_SAVE_OPT.VIRIAL ) ;
	}
	else if( strcmp( file_rec->word[0], "INITIAL_VELOCITIES_BY_FILE" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%d" , &INITIAL_VELOCITIES_BY_FILE ) ;
	  sscanf( file_rec->word[2] , "%s" , velocities_file_name );
	}
	else if( strcmp( file_rec->word[0], "SAVING_DIRECTORY_NAME" ) == 0 ) {
	  sscanf( file_rec->word[1] , "%s" , _DIRECTORY_NAME );
	}
      }
    }
    delete file_rec ;
    file_rec = new record ;
  }
  delete file_rec ;
  file_rec = NULL ;

  // initialization of some variables taken in input script, printing some global information
  if( particle_mass > 0 ) for( int i=0; i<numberOfPart; i++ ) particles[i]->mass = particle_mass ;
  // Initialization of verlet lists pointers in particles
  for( int i=0; i<numberOfPart; i++ ) {
    particles[i]->list_ptr = (vlist_element<particle> **)calloc( INTERACTIONS.number , sizeof(vlist_element<particle> *) );
    if( check_memalloc( particles[i]->list_ptr , "\n list pointer in particles within initialization function" ) ) exit( EXIT_FAILURE ) ;
  }
  for( int k=0 ; k<INTERACTIONS.number ; k++ )  INTERACTIONS.type[k]->link_list2particles() ;
  if( DIMENSION != 2 && DIMENSION != 3 ) {
    cout << " *** Specify the geometrical dimension of the system!" << endl ;
    exit(EXIT_FAILURE) ;
  }
  cout << "kbT = " << k_T << endl ;
  midside.position = (box_sides.position * 0.5) ; // redundant, but midside is a pain... find the way to get rid of it
}
/****************************************************************************************/
/****************************************************************************************/

template <typename particle>
inline bool configuration<particle>::check_verlet_update(void) {
  bool flag = 0 ;
  verlet<particle> *verlet_list = NULL ;
  vlist_element<particle> *vl_element = NULL ;
  particle displacement ;
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
    verlet_list = INTERACTIONS.type[k]->verlet_list ;
    vl_element = verlet_list->element ;
            // In this version I trigger the update of all the verlet lists at the same time, if needed
    if( verlet_list->disp_on  &&  flag == 0 ) {
      double max_disp = verlet_list->delta_verlet * 0.5 ;
      for( int i=0 ; ( i<verlet_list->particles_in_list )&&( flag == 0 ) ; i++ ) {
	    // guadagno realmente cpu-time con questa cosa ?????????????????????????
            // If any particle has moved more than (rverlet-rcut)/2 I recalculate the verlet list
	displacement.position = vl_element[i].ptr->position - vl_element[i].cell_center->position ;
	if( displacement.position.norm() > max_disp ) flag = 1 ;
      }
    }
  }
  if( flag == 1 ) {
    lists_rebuilds ++ ;
    for( int i=0 ; i<numberOfPart ; i++ ) pbc( particles[i] ) ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
      verlet_list = INTERACTIONS.type[k]->verlet_list ;
      if( verlet_list->must_be_updated == 1 ) {
	verlet_list->clear() ;
	verlet_list->verlet_creation_bycell( simulation_mode ) ;
      }
    }
  }
  return flag ;
}
/****************************************************************************************/

template <typename particle>
inline bool configuration<particle>::check_verlet_update( particle *part ) {
  bool flag = 0 ;
  double max_disp ;
  verlet<particle> *verlet_list = NULL ;
  for( int k=0 ; k<INTERACTIONS.number && flag == 0  ; k++ ) {
    verlet_list = INTERACTIONS.type[k]->verlet_list ;
    if( verlet_list->disp_on  &&  flag == 0 ) {
      max_disp = verlet_list->delta_verlet * 0.5 ;
      if( ( part->position - part->list_ptr[k]->cell_center->position ).norm() > max_disp ) flag = 1 ;
    }
  }

  if( flag == 1 ) {
    lists_rebuilds ++ ;
    for( int i=0 ; i<numberOfPart ; i++ ) pbc( particles[i] ) ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
      verlet_list = INTERACTIONS.type[k]->verlet_list ;
      if( verlet_list->must_be_updated == 1 ) {
	verlet_list->clear() ;
	verlet_list->verlet_creation_bycell( simulation_mode ) ;
      }
    }
  }
  return flag ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::add_initialize_MSD( const char *mode , const char *steps_path ) {
  if( mode != NULL && strcmp( mode , "all" ) != 0 && strcmp( mode , "mols_cm" ) != 0 ) {
    cout << "ERROR: MSD initialization mode not recognized !!" << endl ;
    exit( EXIT_FAILURE ) ;
  }
  mean_squared_displacement<particle> **temp = NULL ;
  N_MSDs ++ ;
  temp = new mean_squared_displacement<particle> * [N_MSDs] ;
  for( int i=0 ; i<N_MSDs-1 ; i++ ) temp[i] = MSD[i] ;
  if( MSD != NULL ) delete [] MSD ;
  MSD = temp ;
  MSD[N_MSDs-1] = new mean_squared_displacement<particle> ;
  strcpy( MSD[N_MSDs-1]->_steps_path , steps_path ) ;

  if( mode == NULL || strcmp( mode , "all" ) == 0 ) {
    strcpy( MSD[N_MSDs-1]->_msd_filename , "msd_particles.dat" ) ;
    MSD[N_MSDs-1]->Nparts = numberOfPart ;
    MSD[N_MSDs-1]->parts = new particle * [ MSD[N_MSDs-1]->Nparts ] ;
    for( int i=0 ; i<MSD[N_MSDs-1]->Nparts ; i++ ) MSD[N_MSDs-1]->parts[i] = particles[i] ;

  }
};
/****************************************************************************************/
