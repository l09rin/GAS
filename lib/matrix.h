/**
 * @file matrix.h
 * @brief Declaration of the Matrix class for linear algebra operations.
 *
 * Provides a column vector and Matrix classes with basic operations.
 *
 * @author Giovanni Del Monte
 */


#pragma once
#include <vector>
#include <iostream>

using std::cout ;
using std::endl ;
using std::ostream ;

/**
 * @class col_vector
 * @brief A simple dynamic column vector class.
 *
 * Manages a 1D array of doubles with basic constructors,
 * assignment, and utility functions for vector operations.
 */
class col_vector {
 public:
  /// Dimension (number of elements) in the vector.
  int D = 0 ;
  /// Pointer to dynamically allocated vector elements.
  double *element = NULL ;

  /**
   * @brief Default constructor.
   *
   * Initializes an empty (0-length) vector.
   */
  col_vector( void ) {
    element = NULL ;
    D = 0 ;
  };
  /**
   * @brief Construct a vector with N elements.
   * @param N Number of elements.
   */
  col_vector( int N ) {
    element = new double[N] ;
    D = N ;
  };
  /**
   * @brief Copy constructor.
   * @param other Vector to copy.
   */
  col_vector( const col_vector &other ) {
    D = other.D ;
    element = new double[D] ;
    for( int i=0 ; i<D ; i++ ) element[i] = other.element[i] ;
  };
  /**
   * @brief Destructor.
   *
   * Frees allocated memory.
   */
  ~col_vector() {
    delete [] element ;
    element = NULL ;
    D = 0 ;
  };

  void cartesian_base( int e ) ;
  void operator=( const col_vector &other ) ;
  void clear( void ) ;
  friend ostream & operator<< ( ostream &out, const col_vector &vec ) ;
};

/**
 * @class matrix
 * @brief A simple 2D dynamic matrix class.
 *
 * Provides storage and manipulation of matrices
 * with basic constructors, destructors, and operators.
 */
class matrix {
 public:
  /// Number of rows in the matrix.
  int rows ;
  /// Number of columns in the matrix.
  int columns ;
  /// Pointer to 2D array of elements.
  double **element = NULL ;

  matrix( void ) ;
  matrix( int rowsN, int columnsN ) ;
  matrix( const matrix &other_mat ) ;
  ~matrix() ;

  void rebuild( int rowsN, int columnsN ) ;
  double determinant( void ) ;
  double minor_determinant( int i , int j ) ;
  matrix *construct_minor( int i , int j ) ;
  matrix compute_inverse( void ) ;
  col_vector solve_cramer( const col_vector constant_vec ) ;
  friend ostream & operator<< ( ostream &out, const col_vector &vec ) ;
  matrix dot_prod( const matrix &B ) ;
  void clear( void ) ;
};
