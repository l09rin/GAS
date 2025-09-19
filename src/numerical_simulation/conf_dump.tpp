/* method to dump per particle information */

template <typename particle>
inline void dump_lammps_atom( ifstream &data_file , int pi , particle **particles , particle ***velocities ) ;

/********************************************************************************************************************/

template <typename particle>
void configuration<particle>::dump_configuration( ofstream &out_file , char *format , int Ndigits ) {
  if ( strcmp( format , "xyz" ) == 0 ) dump_xyz( out_file , Ndigits ) ;

  else if ( strcmp( format , "lmp" ) == 0 ) dump_lmp( out_file , Ndigits ) ;

  else if ( strcmp( format , "sph" ) == 0 ||
	    strcmp( format , "flf" ) == 0 ) dump_sph( out_file , Ndigits ) ;

  else {
    cout << "*** ERROR: unrecognized configuration file format!" << endl ;
    exit( EXIT_FAILURE ) ;
  }
}
/********************************************************************************************************************/

template <typename particle>
void configuration<particle>::dump_sph( ofstream &out_file , int Ndigits ) {
  out_file << numberOfPart << " " << global_time << endl ;
  out_file << setprecision(Ndigits) << box_sides.position << endl ;
  if( DIMENSION ==3 ) {
    if( strcmp( simulation_ensemble , "MC_NVT_CMfixed" ) == 0 ) {
      particle unconstrainedCM_pos ;
      for( int i=0 ; i<numberOfPart ; i++ ) {
	unconstrainedCM_pos.position = particles[i]->position ;
	particles[i]->position -= CoM_displacement.position ;
	out_file << (char)(particles[i]->type-1+'a') << " " << particles[i]->dump_position_variables( Ndigits , "sph" , NULL ) << " " << particles[i]->radius << endl ;
	particles[i]->position = unconstrainedCM_pos.position ;
      }
    } else for( int i=0 ; i<numberOfPart ; i++ ) out_file << (char)(particles[i]->type-1+'a') << " " << particles[i]->dump_position_variables( Ndigits , "sph" , NULL ) << " " << particles[i]->radius << endl ;
  } else if( DIMENSION == 2 ) {
    if( strcmp( simulation_ensemble , "MC_NVT_CMfixed" ) == 0 ) {
      particle unconstrainedCM_pos ;
      for( int i=0 ; i<numberOfPart ; i++ ) {
	unconstrainedCM_pos.position = particles[i]->position ;
	particles[i]->position -= CoM_displacement.position ;
	out_file << (char)(particles[i]->type-1+'a') << " " << particles[i]->dump_position_variables( Ndigits , "sph" , NULL ) << " 0.0 " << particles[i]->radius << endl ;
	particles[i]->position = unconstrainedCM_pos.position ;
      }
    } else for( int i=0 ; i<numberOfPart ; i++ ) out_file << (char)(particles[i]->type-1+'a') << " " << particles[i]->dump_position_variables( Ndigits , "sph" , NULL ) << " 0.0 " << particles[i]->radius << endl ;
  }
}
/********************************************************************************************************************/

template <typename particle>
void configuration<particle>::dump_xyz( ofstream &out_file , int Ndigits ) {
  out_file << setprecision(Ndigits) ;
  out_file << "# timestep " << global_time << endl ;
  out_file << "# N " << numberOfPart << endl ;
  out_file << "# box " << box_sides.position << endl ;
  for( int i=0 ; i<numberOfPart ; i++ ) {
    if( SAVE.PARTICLES.POS == 1 ) {
      if( strcmp( simulation_ensemble , "MC_NVT_CMfixed" ) == 0 ) {
	particle constrainedCM_pos( *(particles[i]) ) ;
	constrainedCM_pos.position -= CoM_displacement.position ;
	out_file << constrainedCM_pos.dump_position_variables( 12 , "all" ) << " " ;
      } else out_file << particles[i]->dump_position_variables( 12 , "all" ) << " " ;
    }
    if( SAVE.PARTICLES.POS_IMG == 1 ) out_file << particles[i]->periodic_box << " " ;
    if( SAVE.PARTICLES.VEL == 1 ) out_file << velocities[i]->position << " " ;
    if( SAVE.PARTICLES.FORCES == 1 ) {
      if( SAVE.PARTICLES.FORCES_SPLIT == 1 ) {
	for( int k=0 ; k<INTERACTIONS.number ; k++ ) out_file << SAVE.force_contributions[k][i]->position << " " ;
      } else out_file << forces_t1[i]->position << " " ;
    }
    if( SAVE.PARTICLES.ENERGY == 1 ) out_file << SAVE.per_part_energy[INTERACTIONS.number][i] << " " ;
    if( SAVE.PARTICLES.VIRIAL == 1 ) {
      for( int k=0 ; k<INTERACTIONS.number ; k++ ) out_file << SAVE.per_part_virial[k][i] << " " ;
      out_file << SAVE.per_part_virial[INTERACTIONS.number][i] << " " ;
    }
    out_file << endl ;
  }
  out_file << endl ;
}
/********************************************************************************************************************/

template <typename particle>
void configuration<particle>::dump_lmp( ofstream &out_file , int Ndigits ) {
  cout << "*** ERROR: dump_lmp() method has not yet been coded!" << endl ;
  exit(1);
}
