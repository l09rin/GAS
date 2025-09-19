/*****************************************************************************************/
#include "lennard_jones.cuh"

template <typename particle>
class lennard_jones_gpu : public interaction<particle> { // MD  // energy in units of EPSILON
  using interaction<particle>::verlet_list ;
  using interaction<particle>::cutoff ;
  using interaction<particle>::SHIFT ;
  using interaction<particle>::partial_energy ;
  using interaction<particle>::partial_virial ;
  using interaction<particle>::ON_GPU ;
  using interaction<particle>::idx ;

 public:
 private:
  particle force ;
  double multiplier = 0 , r6_inv = 0 , r2_inv = 0 , cutoff2 = 0 ;   // auxiliary variables
  inline double pair_potential( particle *i_part , particle *j_part , double distance ) ;
  inline particle pair_force( particle *i_part , particle *j_part , double distance ) ;
  dim3 block , grid ;
  struct _GPU_arrays *GPU_arrays = NULL ;
  struct _GPU_cells *GPU_cells = NULL ;

 public:
  double return_pair_potential( double distance ) ;
  double return_pair_potential( particle *part1 , particle *part2 ) ;
  particle return_pair_force( particle *i_part , particle *j_part ) ;
  void generate_by_record( record *file_rec , const char *simulation_kind , configuration<particle> *conf_ptr , int parts_number ) ;

 public:
  double compute_energy_contribution( void ) ;
  inline double compute_energy_contribution( particle *part ) ;
  double delta_energy_constrained_CoM( particle *disp ) ;
  double lambda_derivative_perpart( void ) ;
  void compute_forces_contribution( void ) ;
  double compute_virial_contribution( void ) ;
  void compute_per_part_energy_contribution( double *envec ) ;
  void compute_per_part_virial_contribution( double *virvec ) ;

  string ostr( void ) ;
  void build_bonds_network( double bonds_cutoff = 0.0 ) ;
  double *dump_parameters( void ) ;
  lennard_jones_gpu *copy( void ) ;
};




/************************   CONSTRUCTORS-DISTRUCTORS    *************************************/

/********************************************************************************************/

template <typename particle>
void lennard_jones_gpu<particle>::generate_by_record( record *file_rec , const char *simulation_kind , configuration<particle> *conf_ptr , int parts_number ) {
	              // CUTOFF
  sscanf( file_rec->word[1] , "%lf" , &cutoff ) ;
  cutoff2 = cutoff * cutoff ;
	              // SHIFT
  if( strcmp(file_rec->word[2], "auto") == 0 ) {
    SHIFT = 0 ;
    SHIFT -= pair_potential( NULL , NULL , cutoff ) ;
  }else{
    sscanf( file_rec->word[2] , "%lf" , &SHIFT ) ;
  }
	              // generation of the verlet list
  double verlet_ray ;
  sscanf( file_rec->word[3] , "%lf" , &verlet_ray ) ;
  // verlet_list = new verlet<particle>( conf_ptr , parts_number , cutoff + verlet_ray , verlet_ray ) ;
	              // non - generation of the verlet list
  verlet_list = new verlet<particle>( conf_ptr , 0 , 0.0 , 0.0 ) ;
  if( check_memalloc( verlet_list , "Allocation error of verlet_list in interaction generation!" ) ) exit( EXIT_FAILURE ) ;
  verlet_list->must_be_updated = 0 ; // Indeed it has been already done in the constr
  cout << "_____ Compute Forces" << endl ;
  if( conf_ptr->DIMENSION == 3 ) findOptimalGrid( grid , block ,
						  static_cast<void (*)(double*, double*, double*, double*, double*, double*, int*, int*, double)>(d_compute_forces_lj) ,
						  conf_ptr->numberOfPart , 0 ) ;
  else if( conf_ptr->DIMENSION == 2 ) findOptimalGrid( grid , block ,
						       static_cast<void (*)(double*, double*, double*, double*, int*, int*, double)>(d_compute_forces_lj) ,
						       conf_ptr->numberOfPart , 0 ) ;
  GPU_arrays = &(conf_ptr->GPU_arrays) ;
  GPU_cells = &(conf_ptr->GPU_cells) ;
};
/********************************************************************************************/

template <typename particle>
inline double lennard_jones_gpu<particle>::pair_potential( particle *i_part , particle *j_part , double distance ) {
    /***********   LENNARD-JONES POTENTIAL     offset energy Vlj(Rc=3)    **********/
  r2_inv = 1.0 / ( distance*distance ) ;
  r6_inv = r2_inv * r2_inv * r2_inv ;
  return 4.0 * ( r6_inv * r6_inv - r6_inv ) ;
};
/********************************************************************************************/

template <typename particle>
inline particle lennard_jones_gpu<particle>::pair_force( particle *i_part , particle *j_part , double distance ) {
    /***********   LENNARD-JONES POTENTIAL     **********/                                           // forces are in units of EPSILON/SIGMA
  r2_inv = 1.0 / ( distance*distance ) ;
  r6_inv = r2_inv * r2_inv * r2_inv ;
  multiplier = ( 48.0 * r6_inv * r6_inv - 24.0 * r6_inv ) * r2_inv ;
  force.position = i_part->position - j_part->position ;
  force.position *= multiplier ;

  return force ;
};
/********************************************************************************************/

template <typename particle>
double lennard_jones_gpu<particle>::return_pair_potential( double distance ) {
  return pair_potential( NULL , NULL , distance ) + SHIFT ;
};
/********************************************************************************************/

template <typename particle>
double lennard_jones_gpu<particle>::return_pair_potential( particle *part1 , particle *part2 ) {
  double energy = 0 ;
  double dist = 0 ;
  particle *image_part = NULL ;
  image_part = new particle ;

  verlet_list->god->nearest_image( part2 , part1 , image_part ) ;
  if( ( dist = part1->position.distance( image_part->position ) ) < cutoff ) {
    energy = ( pair_potential( part1 , image_part , dist ) + SHIFT ) ;
  }

  delete image_part ;
  return energy ;
};
/********************************************************************************************/

template <typename particle>
particle lennard_jones_gpu<particle>::return_pair_force( particle *i_part , particle *j_part ) {
  double dist = i_part->position.distance( j_part->position ) ;
  return pair_force( i_part , j_part , dist ) ;
};
/********************************************************************************************/

template <typename particle>
double lennard_jones_gpu<particle>::compute_energy_contribution( void ) {
  zero_vectors<<< grid , block >>>( GPU_arrays->d_compute_vec ) ;
  cudaDeviceSynchronize() ;

  d_compute_energy_lj<<< grid , block >>>( GPU_arrays->d_x , GPU_arrays->d_y , GPU_arrays->d_compute_vec , GPU_cells->m_dCellStart , GPU_cells->m_dCellEnd , cutoff2 , SHIFT ) ;
  cudaDeviceSynchronize() ;

  thrust::device_ptr<double> dev_ptr(GPU_arrays->d_compute_vec) ;
  double energy = thrust::reduce( thrust::device , dev_ptr , dev_ptr + verlet_list->god->numberOfPart , 0.0 , thrust::plus<double>()) ;
  return energy * 0.5 ;
};
/********************************************************************************************/

template <>
double lennard_jones_gpu<particle_3D>::compute_energy_contribution( void ) {
  zero_vectors<<< grid , block >>>( GPU_arrays->d_compute_vec ) ;
  cudaDeviceSynchronize() ;

  d_compute_energy_lj<<< grid , block >>>( GPU_arrays->d_x , GPU_arrays->d_y , GPU_arrays->d_z , GPU_arrays->d_compute_vec , GPU_cells->m_dCellStart , GPU_cells->m_dCellEnd , cutoff2 , SHIFT ) ;
  cudaDeviceSynchronize() ;

  thrust::device_ptr<double> dev_ptr(GPU_arrays->d_compute_vec) ;
  double energy = thrust::reduce( thrust::device , dev_ptr , dev_ptr + verlet_list->god->numberOfPart , 0.0 , thrust::plus<double>()) ;
  return energy * 0.5 ;
};
/********************************************************************************************/

template <typename particle>
inline double lennard_jones_gpu<particle>::compute_energy_contribution( particle *part ) {
  // This function computes the energy of ONE particle only
  double dist ;
  static double energy ;
  static neighbour_list<particle> *walker = NULL ;
  static particle image_part ;
  energy = 0 ;

  if( part->list_ptr[idx] ) {
    walker = part->list_ptr[idx]->neighbours ;
    while( walker != NULL ) {
      verlet_list->god->nearest_image( walker->ptr , part , &image_part ) ;
      if( ( dist = part->position.distance( image_part.position ) ) < cutoff ) {
	energy += ( pair_potential( part , &image_part , dist ) + SHIFT ) ;
      }
      walker = walker->next ;
    }
  }

  return energy ;
};
/********************************************************************************************/

template <typename particle>
double lennard_jones_gpu<particle>::delta_energy_constrained_CoM( particle *disp ) {
  // This function computes the energy correction due to a displacement of a particle with the constrain of keeping fixed the centre of mass
  // For a pair interaction potential it gives trivially 0 contribution
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
double lennard_jones_gpu<particle>::lambda_derivative_perpart( void ) {
  // This function computes the mean squared displacement from ideal lattice sites with the external harmonic spring potential
  // For a pair interaction potential it gives trivially 0 contribution
  return 0.0 ;
};
/********************************************************************************************/

template <typename particle>
void lennard_jones_gpu<particle>::compute_forces_contribution(void) {
  d_compute_forces_lj<<< grid , block >>>( GPU_arrays->d_x , GPU_arrays->d_y , GPU_arrays->d_fx , GPU_arrays->d_fy , GPU_cells->m_dCellStart , GPU_cells->m_dCellEnd , cutoff2 ) ;
  cudaDeviceSynchronize() ;
};
/********************************************************************************************/

template <>
void lennard_jones_gpu<particle_3D>::compute_forces_contribution(void) {
  d_compute_forces_lj<<< grid , block >>>( GPU_arrays->d_x , GPU_arrays->d_y , GPU_arrays->d_z , GPU_arrays->d_fx , GPU_arrays->d_fy , GPU_arrays->d_fz , GPU_cells->m_dCellStart , GPU_cells->m_dCellEnd , cutoff2 ) ;
  cudaDeviceSynchronize() ;
};
/********************************************************************************************/

template <typename particle>
double lennard_jones_gpu<particle>::compute_virial_contribution(void) {
  zero_vectors<<< grid , block >>>( GPU_arrays->d_compute_vec ) ;
  cudaDeviceSynchronize() ;

  d_compute_virial_lj<<< grid , block >>>( GPU_arrays->d_x , GPU_arrays->d_y , GPU_arrays->d_compute_vec , GPU_cells->m_dCellStart , GPU_cells->m_dCellEnd , cutoff2 ) ;
  cudaDeviceSynchronize() ;

  thrust::device_ptr<double> dev_ptr(GPU_arrays->d_compute_vec) ;
  double virial = thrust::reduce( thrust::device , dev_ptr , dev_ptr + verlet_list->god->numberOfPart , 0.0 , thrust::plus<double>()) ;
  return virial * 0.5 ;
};
/********************************************************************************************/

template <>
double lennard_jones_gpu<particle_3D>::compute_virial_contribution(void) {
  zero_vectors<<< grid , block >>>( GPU_arrays->d_compute_vec ) ;
  cudaDeviceSynchronize() ;

  d_compute_virial_lj<<< grid , block >>>( GPU_arrays->d_x , GPU_arrays->d_y , GPU_arrays->d_z , GPU_arrays->d_compute_vec , GPU_cells->m_dCellStart , GPU_cells->m_dCellEnd , cutoff2 ) ;
  cudaDeviceSynchronize() ;

  thrust::device_ptr<double> dev_ptr(GPU_arrays->d_compute_vec) ;
  double virial = thrust::reduce( thrust::device , dev_ptr , dev_ptr + verlet_list->god->numberOfPart , 0.0 , thrust::plus<double>()) ;
  return virial * 0.5 ;
};
/********************************************************************************************/

template <typename particle>
string lennard_jones_gpu<particle>::ostr( void ) {
  ostringstream st ;
  st << "LENNARD-JONES potential has no further parameters dependence." << endl ;
  return st.str() ;
};
/********************************************************************************************/

template <typename particle>
lennard_jones_gpu<particle> *lennard_jones_gpu<particle>::copy( void ) {
  lennard_jones_gpu *copy = new lennard_jones_gpu() ;
  return copy ;
};
/********************************************************************************************/

template <typename particle>
void lennard_jones_gpu<particle>::compute_per_part_energy_contribution( double *envec ) {
  cout << "*** ERROR : This method should still be implemented on GPU calculations!" << endl ;
  exit( EXIT_FAILURE ) ;
};
/********************************************************************************************/

template <typename particle>
void lennard_jones_gpu<particle>::compute_per_part_virial_contribution( double *virvec ) {
  cout << "*** ERROR : This method should still be implemented on GPU calculations!" << endl ;
  exit( EXIT_FAILURE ) ;
};
/********************************************************************************************/

template <typename particle>
void lennard_jones_gpu<particle>::build_bonds_network( double bonds_cutoff ) {
  cout << "*** WARNING : no method defined to build bonds_network" << endl ;
};
/********************************************************************************************/

template <typename particle>
double * lennard_jones_gpu<particle>::dump_parameters( void ) {
  double *params = NULL ;
  return params ;
};
/********************************************************************************************/
/********************************************************************************************/
