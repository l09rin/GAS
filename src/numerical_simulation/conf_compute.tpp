
template <typename particle>
double configuration<particle>::compute_potential_energy( void ) {
  POTENTIAL_ENERGY.val = 0 ;
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
    INTERACTIONS.type[k]->partial_energy.val = INTERACTIONS.type[k]->compute_energy_contribution() ;
    INTERACTIONS.type[k]->partial_energy.time = global_time ;
    POTENTIAL_ENERGY.val += INTERACTIONS.type[k]->partial_energy.val ;
  }
  POTENTIAL_ENERGY.time = global_time ;
  return POTENTIAL_ENERGY.val ;
}
/*****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_potential_energy( particle *part ) {   // MC mainly       // This requires verlet lists constructed in MC mode !
  double energy = 0 ;
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) energy += INTERACTIONS.type[k]->compute_energy_contribution( part ) ;
  return energy ;
}
/*****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_potential_energy( particle *part1 , particle *part2 ) {  // MC mainly
  double energy = 0 ;
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) energy += INTERACTIONS.type[k]->return_pair_potential( part1 , part2 ) ;
  if( strcmp( simulation_mode , "MC" ) == 0 ) return energy ;
  else {
    cout << " *** ERROR in function energy_compute( &part1 , &part2 ) : Verlet lists are constructed in MD mode, the single particle energy cannot be computed !" << endl ;
    exit( EXIT_FAILURE ) ;
  }
}
/*****************************************************************************************/

template <typename particle>
double configuration<particle>::delta_energy_constrained_CoM( particle *displacement ) {
  double energy = 0 ;
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) energy += INTERACTIONS.type[k]->delta_energy_constrained_CoM( displacement ) ;
  return energy ;
}
/*****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_lambda_derivative_perpart( void ) {
  double der = 0 ;
  for( int k=0 ; k<INTERACTIONS.number ; k++ )  der += INTERACTIONS.type[k]->lambda_derivative_perpart() ;
  return der ;
}
/*****************************************************************************************/

template <typename particle>
void configuration<particle>::calculate_forces_t2() {     // MD    // calculates the vector of forces acting on particles and place it in forces_t2
  if( strcmp( simulation_mode , "MC" ) == 0 ) {
    cout << " *** ERROR : function calculate_forces_t2() cannot be used in MC mode !" << endl ;
    exit( EXIT_FAILURE ) ;
  }
          // reset of forces_t2
  for( int i=0 ; i<numberOfPart ; i++ ) forces_t2[i]->clear() ;
          // computation of forces due to each type of interaction
  for( int k=0 ; k<INTERACTIONS.number ; k++ ) INTERACTIONS.type[k]->compute_forces_contribution() ;
}
/*****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_virial_pressure(void) {      // si può ottimizzare
  if( strcmp( simulation_mode , "MC" ) == 0 ) {
    cout << " *** ERROR : to compute the virial in Monte Carlo simulations you have to allocate the array of forces ! " << endl
	 << "             Remember that you should re-organize the information on particles, probably embedding forces in some way," << endl
	 << "             but you have also consider what this would imply for integration with velocity verlet algorithm . . ." << endl
	 << "             Maybe you have to store a variable in the functions that tells what force among two refers to the actual step . . ." << endl ;
    exit( EXIT_FAILURE ) ;
  }

  if( virial_press.time != global_time ) {
    virial_press.val = 0 ;
    for( int k=0 ; k<INTERACTIONS.number ; k++ ) {
      INTERACTIONS.type[k]->partial_virial.val = INTERACTIONS.type[k]->compute_virial_contribution() ;
      INTERACTIONS.type[k]->partial_virial.time = global_time ;
      virial_press.val += INTERACTIONS.type[k]->partial_virial.val ;
    }
    virial_press.val /= ( DIMENSION * box_sides.volume() ) ;
    virial_press.time = global_time ;
  }
  return virial_press.val ;
}
/****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_kinetic_energy() {           //MD
  if( strcmp( simulation_mode , "MC" ) == 0 ) {
    cout << " *** ERROR : function compute_kinetic_energy() cannot be used in MC mode !" << endl ;
    exit( EXIT_FAILURE ) ;
  }

  if( KINETIC_ENERGY.time != global_time ) {
    KINETIC_ENERGY.val = 0 ;

    if( DEVICE == GPU ) {
      zero_vectors<<< GPU_params.grid , GPU_params.block >>>( GPU_arrays.d_compute_vec ) ;
      cudaDeviceSynchronize() ;

      if( DIMENSION == 2 ) compute_v2<<< GPU_params.grid , GPU_params.block >>>( GPU_arrays.d_vx , GPU_arrays.d_vy , GPU_arrays.d_compute_vec ) ;
      else compute_v2<<< GPU_params.grid , GPU_params.block >>>( GPU_arrays.d_vx , GPU_arrays.d_vy , GPU_arrays.d_vz , GPU_arrays.d_compute_vec ) ;
      cudaDeviceSynchronize() ;
      thrust::device_ptr<double> dev_ptr(GPU_arrays.d_compute_vec) ;
      KINETIC_ENERGY.val = 0.5 * thrust::reduce( thrust::device , dev_ptr , dev_ptr + numberOfPart , 0.0 , thrust::plus<double>()) ;
    } else {

      for( int i=0 ; i<numberOfPart ; i++ ) {
	KINETIC_ENERGY.val += velocities[i]->position.square_norm() ;
      }
      KINETIC_ENERGY.val /= 2.0 ;
    }
    KINETIC_ENERGY.time = global_time ;
  }
  return KINETIC_ENERGY.val ;
}
/*****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_packing_fraction() {
  double packfrac = 0.0 ;
  for( int i=0 ; i<numberOfPart ; i++ ) packfrac += ( particles[i]->radius * particles[i]->radius ) ;
  packfrac *= ( M_PI / box_sides.volume() ) ;
  return packfrac ;
}
/*****************************************************************************************/

template <>
double configuration<particle_3D>::compute_packing_fraction() {
  double packfrac = 0.0 ;
  for( int i=0 ; i<numberOfPart ; i++ ) packfrac += ( particles[i]->radius * particles[i]->radius * particles[i]->radius ) ;
  packfrac *= ( 4. / 3 * M_PI / box_sides.volume() ) ;
  return packfrac ;
}
/*****************************************************************************************/

template <typename particle>
double configuration<particle>::calculate_secondary_potential(void) {                   //MC
  // double energy = 0 ;
  // double dist = 0 ;
  // struct neighbour_list *walker = NULL ;
  // particle image_part ;

  cout << " Rewrite this function (secondary potential calculation)" << endl ;
  exit( EXIT_FAILURE ) ;

  /* for(int i=0; i<numberOfPart-1; i++) {                   // attenzione al numberOfPart-1 */
  /*   walker = verlet_list->neighbours[i] ; */
  /*   while( walker != NULL ) {                 //for each i-particle I sum to the total the interaction energy with near j-particle if j>i (see verlet list construction) in order to avoid double counting */
  /*     if( (walker->near) > i ) {                //and if Ri-Rj<Rcutoff */
  /* 	image_part = nearest_image( walker->near , i ) ;      // I take the coordinates of the nearest copy of the neighbour */
  /* 	if( ( dist = particles[i]->distance( image_part ) ) < R_cutoff ) energy += secondary_potential( particles[i] , &image_part , dist ) ; */
  /*     } */
  /*     walker = walker->next ; */
  /*   } */
  /* } */
  /* return energy ; */
}
/****************************************************************************************/

template <typename particle>
particle configuration<particle>::compute_mean_position(void) {
  particle mean_pos ;
  for( int i=0; i<numberOfPart; i++ ) mean_pos.position += particles[i]->position ;
  mean_pos.position /= (double)numberOfPart ;

  return mean_pos ;
}
/****************************************************************************************/

template <typename particle>
particle configuration<particle>::compute_centre_of_mass(void) {
  centre_of_mass.position *= 0.0 ;
  double mass_tot = 0 ;
  for( int i=0; i<numberOfPart; i++ ) {
    centre_of_mass.position += ( particles[i]->unwrapped( box_sides.position ) * particles[i]->mass ) ;
    mass_tot += particles[i]->mass ;
  }
  centre_of_mass.position /= mass_tot ;

  return centre_of_mass ;
}
/****************************************************************************************/

template <typename particle>
particle configuration<particle>::compute_mean_velocity( void ) {
  particle mean_vel ;
  for( int i=0 ; i<numberOfPart ; i++ ) mean_vel.position += velocities[i]->position ;
  mean_vel.position /= (double)numberOfPart ;

  return mean_vel ;
}
/****************************************************************************************/

template <typename particle>
struct posizione_3D<double> configuration<particle>::compute_CoM_angular_momentum( void ) {    // in case of identical particles of mass = 1
  particle CoM ;
  posizione_3D<double> ang_mom { 0.0 , 0.0 , 0.0 } ;
  CoM = compute_mean_position() ;
  for( int i=0 ; i<numberOfPart ; i++ ) ang_mom -= ( CoM.position - particles[i]->position ).vectorial_product( velocities[i]->position ) ;

  return ang_mom ;
}
/****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_inertia_momentum( posizione_3D<double> axis ) {    // in case of identical particles of mass = 1
  double I = 0 ;
  particle CoM ;
  CoM = compute_mean_position() ;
  axis /= axis.norm() ;
  for( int i=0; i<numberOfPart; i++ ) {
    I += axis.vectorial_product( particles[i]->position - CoM.position ).vectorial_product( particles[i]->position - CoM.position ).norm() ;
  }

  return I ;
}
/****************************************************************************************/

template <typename particle>
void configuration<particle>::compute_inertia_momentum( void ) {    // in case of identical particles of mass = 1
  particle CoM ;
  posizione_3D<double> r ;
  CoM = compute_mean_position() ;
  for( int i=0; i<3; i++ ) {
    for( int j=0; j<3; j++ ) {
      inertia_matrix.element[i][j] = 0 ;
    }
  }
  for( int i=0; i<numberOfPart; i++ ) {
    r = ( particles[i]->position - CoM.position ) ;
    inertia_matrix.element[0][0] += particles[i]->mass * ( r.y * r.y + r.z * r.z ) ;
    inertia_matrix.element[0][1] -= particles[i]->mass * r.x * r.y ;
    inertia_matrix.element[0][2] -= particles[i]->mass * r.x * r.z ;
    inertia_matrix.element[1][0] -= particles[i]->mass * r.y * r.x ;
    inertia_matrix.element[1][1] += particles[i]->mass * ( r.x * r.x + r.z * r.z ) ;
    inertia_matrix.element[1][2] -= particles[i]->mass * r.y * r.z ;
    inertia_matrix.element[2][0] -= particles[i]->mass * r.z * r.x ;
    inertia_matrix.element[2][1] -= particles[i]->mass * r.z * r.y ;
    inertia_matrix.element[2][2] += particles[i]->mass * ( r.y * r.y + r.x * r.x ) ;
  }
}
/****************************************************************************************/

template <typename particle>
double configuration<particle>::compute_max_bond_distance( void ) {
  neighbour_list<particle> *walker=NULL ;
  particle *image_part = NULL ;
  image_part = new particle ;
  double max_bonding_distance = 0, dist = 0 ;

  for( int l=0 ; l<CHECKS.bonds_lists_number ; l++ ) {
    for( int i=0 ; i<CHECKS.bonds_verlet_link[l]->particles_in_list ; i++ ) {                 // bonding interactions
      walker = CHECKS.bonds_verlet_link[l]->element[i].neighbours ;
      while( walker != NULL ) {
	nearest_image( walker->ptr , CHECKS.bonds_verlet_link[l]->element[i].ptr , image_part ) ;      // I take the coordinates of the nearest copy of the neighbour. The list only one element for each couple of particles
	if( ( dist = particles[i]->position.distance( image_part->position ) ) > max_bonding_distance ) max_bonding_distance = dist ;
	walker = walker->next ;
      }
    }
  }

  delete image_part ;
  return max_bonding_distance ;
}
/*****************************************************************************************/
