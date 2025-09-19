#pragma once

#include <cuda_runtime.h>
#include "thrust/device_ptr.h"
#include "thrust/sort.h"
#include <thrust/reduce.h>
#include "configuration.cuh"
#include "Interactions/lennard_jones.cuh"
#include "../../lib/gpu_comm.cuh"

#include <ctime>
#include <fstream>
#include <iomanip>
#include <cstdlib>

#include "../../lib/system_comm.h"

using std::cout ;
using std::endl ;
using std::ifstream ;
using std::ofstream ;

#include "../../lib/matrix.h"
#include "../../lib/record.h"
#include "particles2D.h"
#include "particles3D.h"
#include "cells.h"
#include "verlet.h"

#define SQ(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define CPU 0
#define GPU 1

template <typename particle> class interaction ;
template <typename particle> class mean_squared_displacement ;

// DEFINITION OF THE CLASS CONFIGURATION, CONTAINING ALL THE INFORMATION ABOUT THE SYSTEM AND THE SIMULATION METHODS


template <typename particle>
class configuration {

 public :

  int seed ;                    // It can be used to initialize srand48()

  /*******************************************   SIMULATION BOX AND ENSEMBLES CONSTANTS   *******************************************/

  int DIMENSION ;
  double k_T ;
  double energyOfSystem ;               // total value, sum over all particles
  double box_side ;
  particle box_sides ; // If the box is not cubic    // AGGIORNARE LE FUNZIONI IN CUI RIENTRA
  particle midside ;
  inline void pbc( particle_3D *part ) ;   // imposes the periodic boundary conditions
  inline void pbc( particle_2D *part ) ;   // imposes the periodic boundary conditions
  char simulation_ensemble[15] ;
  char simulation_mode[3] ;  // It can be MC or MD and determines how the verlet lists are constructed
         // MD-ONLY VARIABLES
  double time_step , global_time ;           // global_time is the global time of the simulation, time_step is the integration step for MD
  double particle_mass ;

  /*******************************************   THERMODYNAMIC QUANTITIES (COMPUTATION NEEDED)   *******************************************/

  time_dep_variable<double> virial_press ;
  time_dep_variable<double> POTENTIAL_ENERGY , KINETIC_ENERGY ;               // total values, sum over all particles


  /*******************************************   MEAN SQUARED DISPLACEMENT (msd.h)   *******************************************/
  mean_squared_displacement<particle> **MSD = NULL ;
  int N_MSDs = 0 ;
  void add_initialize_MSD( const char *mode , const char *steps_path = "msd_saving_steps.dat" ) ;

                                               // READ CONFIGURATIONS BY FILES ( read_data.cpp / dump_data.cpp )
  void read_configuration( ifstream &data_file , char *format ) ;
  void read_xyz( ifstream &data_file ) ;
  void read_lammps_dump( ifstream &data_file ) ;
  void read_sph( ifstream &data_file ) ;
  void dump_configuration( ofstream &out_file , char *format , int Ndigits = 8 ) ;
  void dump_sph( ofstream &out_file , int Ndigits = 8 ) ;
  void dump_xyz( ofstream &out_file , int Ndigits = 8 ) ;
  void dump_lmp( ofstream &out_file , int Ndigits = 8 ) ;

  /*******************************************   INTERACTIONS ( folder Interactions )   *******************************************/

  class interactions_array {
  public:
    int number ;
    double SIGMA, EPSILON ;
    interaction<particle> **type ;

    void print_on_file( double delta , int steps ) ;
    //friend ostream & operator<< <particle>( ostream &out , const interactions_array &int_array ) ;
    friend ostream & operator<<( ostream &out , const interactions_array &int_array ) {
      for( int i=0 ; i<int_array.number ; i++ ) {
	out << int_array.type[i]->ostr() ;
      }
      return out ;
    };
    void operator=( const interactions_array &other ) ;
    interactions_array() ;
    ~interactions_array() ;
  } INTERACTIONS , SECONDARY_INTERACTIONS ;

  /*******************************************   CHECKS   *******************************************/

  struct check {
    int bonds_lists_number = 0 ;
    verlet<particle> **bonds_verlet_link = NULL ;
  } CHECKS ;
  int lists_rebuilds = 0 ;


  /*******************************************   SYSTEM VECTORS   *******************************************/

  int numberOfPart , monomers_number , particles_length ;   // number of all particles, monomers fene-bonded, charged monomers in the system
  particle **particles = NULL ;                       // It is a vector containing pointers to particle positions
  matrix inertia_matrix ;
  particle **velocities = NULL ;                      // It is a vector containing pointers to particle velocities
  particle **forces_t1 = NULL, **forces_t2 = NULL ;     // It is a vector containing total forces acting on each particle in the system
  particle **support_pp2particle = NULL ;
      // intended as calculated with unwrapped coordinates, used in ensemble MC_NVT_CMfixed, for the moment
      // CoM_displacement is the displacement of the centre of mass since the beginning of the simulation
  particle centre_of_mass , CoM_displacement ;
  double deltaEnergyOfCoM = 0.0 ;   // energy difference among the system with and without a constrained centre of mass in MC NVT ensemble
          // INPUT - OUTPUT
  int INITIAL_CONF_BY_FILE ;             // if I want to start from a configuration of stored positions
  int INITIAL_VELOCITIES_BY_FILE ;             // if I want to start from a configuration of stored velocities
  char configuration_file_name[300], velocities_file_name[300] ;

  int DEVICE = CPU ;
  cudaDeviceProp GPU_PROPS ;
  struct _GPU_arrays GPU_arrays ;
  struct _GPU_params GPU_params ;
  struct _GPU_cells GPU_cells ;

  struct save_options {      // I use this variables in the main() to set if storing this data in files produced
    struct PART_ATTR {
      int POS, POS_IMG, VEL, FORCES, FORCES_SPLIT, ENERGY, VIRIAL, MOL, IDX, TYPE;
    } PARTICLES;
    int POTENTIAL_ENERGY , PARTIAL_POTENTIAL_ENERGY , KINETIC_ENERGY , SECONDARY_POTENTIAL , MSD ;
    int VOLUME , VIRIAL , TEMPERATURE ;      // Data are saved if correspondent values are set to 1, 0 by default
    int PRECISION10 ;                          // number of significant digits that are used to save thermodynamic output
    particle ***force_contributions = NULL ;   // This is needed in case of forces_splitted are used, to split the several contribution of the forces
    double **per_part_energy = NULL , **per_part_virial = NULL ;  // This is needed if you want to dump the values of energy and/or virial for each particle
    int global_interval ;               // number of steps between successive savings of global quantities
    int configuration_interval ;        // number of steps between successive savings of positions and velocities
    int backup_interval ;               // number of steps between storage of positions and velocities for a restart of the simulation
    char *configuration_format ;
    void copy_options( const save_options &other ) ;
    void clear( void ) ;
    save_options( void ) ;
  } SAVE, EQ_SAVE_OPT, PROD_SAVE_OPT ;

  configuration( void ) ;
  configuration( const configuration &other_config ) ;   // copy constructor
  ~configuration() ;


  /*******************************************   BLOCK 1 - DEFINITIONS   *******************************************/

  void initialize_configuration( const char *input_file_name = NULL ) ;
  bool check_init_conf_set( void ) ;
  void interpret_input_script( ifstream &config_file ) ;


  /*******************************************   INITIAL POSITIONS SETTINGS ( initial_configuration.cpp )   *******************************************/

  void generate_ControlledEnergy_randomConfig( void ) ;     // this function generate a random configuration through controlled energy insertion for MD simulations
  void generate_HS_randomConfig( double radius = 0.5 ) ;  // this function generate a random configuration as if the particles were HS of determined radius
  void generate_ideal_lattice_3D( const char *ltype ) ;
  void generate_ideal_lattice( const char *ltype ) ;
  void generate_config_by_file( const char *file_name , const char *format = NULL ) ;


  /*******************************************   INITIAL VELOCITY SETTINGS ( velocities.h )   *******************************************/
  void generate_MBvelocities( double mean_square_velocity ) ; // this function generates a vector of velocity according to a
                                                            // Maxwell-Boltzmann distribution with mean value = compute_mean_velocity and mean square value = mean_square_velocity
  void generate_velocities_by_file( void ) ;
  inline void CoM_velocity_reset( void ) ;
  void total_angularMomentum_reset( const char protocol[] = "" ) ;  // It cancels the total angular momentum in the CM (protocol=RIGID|VTINVERSION|NULL)
                        // VTINVERSION = inverting the transverse component of velocities contributing to L ; RIGID = as if it were a rigid body : NULL = same as RIGID
  void set_kinetic_energy( double kin_en ) ;   // fix to kin_en the overall value of kinetic energy of system


  /*******************************************   BLOCK 2 - COMPUTATIONS   *******************************************/

  double compute_potential_energy( void ) ;            // calculates the total energy of the configuration
  double compute_potential_energy( particle *part ) ;      // calculates the energy of a single particle (all the amount of pair and molecular contributions are assigned)
  double compute_potential_energy( particle *part1 , particle *part2 ) ;   // calculates the energy of a single pair of particles (molecular contributions are assigned only if the pair belongs to the same molecule)
  double delta_energy_constrained_CoM( particle *displacement ) ;            // calculates the energy difference among the constrained and unconstrained CoM system after the displacement move of a particle in the MC NVT ensemble with fixed CoM integration scheme
  double compute_lambda_derivative_perpart( void ) ;            // calculates the lambda derivative of the potentials in the spring lattice fields with centre of mass constraint
  void calculate_forces_t2() ;              // calculates the vector of forces acting on particles and place it in forces_t2
  double compute_virial_pressure( void ) ;
  double compute_kinetic_energy() ;            // calculates the kinetic energy of the configuration
  double compute_packing_fraction() ;
  double calculate_secondary_potential( void ) ;
  particle compute_centre_of_mass( void ) ;
  particle compute_mean_position( void ) ;                    // center of mass in case of identical particles
  particle compute_mean_velocity( void ) ;                    // velocity of center of mass in case of identical particles
  posizione_3D<double> compute_CoM_angular_momentum( void ) ;                    // in case of identical particles of mass 1
  double compute_inertia_momentum( posizione_3D<double> axis ) ;                    // in case of identical particles of mass 1
  void compute_inertia_momentum( void ) ;
  double compute_max_bond_distance( void ) ;         // calculates the maximume stretching of bonds in the configuration


  /*******************************************   GPU IMPLEMENTATION   *******************************************/
  void initialize_GPU_data() ;
  void allocateDeviceBox() ;
  void allocateDeviceNeighHashes() ;
  void allocateDeviceArrays( bool VEL = 0 , bool FORCE = 0 ) ;
  void freeDeviceArrays() ;
  void copyParticlesToDevice( bool POS = 1 , bool VEL = 0 , bool FORCE = 0 ) ;
  void copyParticlesToHost( bool POS = 1 , bool VEL = 0 , bool FORCE = 0 , bool POS_IMG = 0 ) ;
  void allocateDeviceCellList( double cutoff ) ;
  void freeDeviceCellList() ;
  void buildDeviceCellList() ;


  /*******************************************   SIMULATION ENSEMBLES ( Integrators/ )   *******************************************/

  int EQUILIBRATION_STEPS ;
  int PRODUCTION_STEPS ;
  void MD( const char *questions ) ;    // Molecular Dynamics simulation of NVE, NVT_NHC, NVT_NHE, BROWNIAN ensembles. IF Questions == NO_STOP equilibration and production steps cannot be modified
  void MD_GPU( void ) ;    // Molecular Dynamics integration on GPU, NVE only for the moment.

 protected:

  /*******************************************   OUTPUT   *******************************************/

  char _DIRECTORY_NAME[300] ;
  ofstream _thermoinfo_file ;
  ofstream _particles_file ;
  ofstream _positions_backup ;
  ofstream _velocities_backup ;

  /*******************************************   CELL - VERLET METHOD   *******************************************/

public:
  inline bool check_verlet_update();  // This member checks if the verlet list has to be updated, triggering it if needed
  inline bool check_verlet_update( particle *part );  // same as previous method, only checks one particle
                                                  // VERLET METHOD
  inline void nearest_image( particle_3D *jptr , particle_3D *iptr , particle_3D *copy ) ;  // In a box with periodic boundary conditions this member returns the j-th part copy that is nearest to i-th part
  inline void nearest_image( particle_2D *jptr , particle_2D *iptr , particle_2D *copy ) ;  // In a box with periodic boundary conditions this member returns the j-th part copy that is nearest to i-th part


  /*******************************************   BLOCK 2 - DEFINITIONS   *******************************************/
protected:
  inline void check_action( int current_step , int starting_time , int screen_print_interval ) ;   // This function contains all the routines deciding if saving at each integration cycle


  /*******************************************   INTEGRATORS ( Integrators/ )   *******************************************/

  // MD NVE
  inline void VelocityVerletIntegrationStep( void ) ;
  inline void VelocityVerletIntegrationStep_GPU( void ) ;
  void calculate_forces_gpu( void ) ;


  /*******************************************   DATA SAVING ( data_saving.h )   *******************************************/

  void plot_equilibration_data( void ) ;       // produce plots of energy and number of particles versus MC-time
  char *generate_new_directory_for_data( void ) ;              // this function generate a new directory for saving data
  char *generate_equilibration_subdirectory( char *root ) ;        // this function generate a subdirectory in which are generated storage files for equilibration data. Returns Root-directory name
  void initialize_storage_file( int mcsteps, double seed ) ;  // opens the storage file and prints the headlines
  inline void save_configuration( int step ) ;                        // saves the configuration in positions' and velocities' files
  inline void save_data( int step ) ;                        // saves global quantities in file thermoinfo_file
  void set_saving_options( const char *MODE ) ;                     // This function allows to stop and resume the saving data process
  void close_storage_file( void ) ;
  void save_simulation_conditions( const char *file_name ) ;  // this function save in the simulation folder a copy of configuration file with a new name
  inline void write_backup_data( int current_step ) ;                    // This fuction update two files in which the program stores velocities and positions every (backup_interval) steps


  /*******************************************   PRINT DEBUG INFORMATION ( print_info.h )   *******************************************/

  void print_all_interactions( void ) ;
  void print_all_information( void ) ;
  void print_displacements( void ) ;
  void print_max_velocity( void ) ;
  double max_component_of_vector( particle **vec , int length ) ;
  double *min_part_distance( void ) ;
  double *max_part_distance( void ) ;

};


#include "configuration.tpp"
#include "conf_compute.tpp"

#include "configuration.cu"
#include "conf_compute.cu"
#include "Integrators/molecular_dynamics.cu"
#include "conf_kernels.cu"

#include "Integrators/molecular_dynamics.tpp"
#include "conf_box.tpp"
#include "conf_init.tpp"
#include "conf_velocities.tpp"
#include "conf_save.tpp"
#include "conf_print_info.tpp"
#include "conf_read.tpp"
#include "conf_dump.tpp"
