/**
 * @file matrix.cpp
 * @brief Implementation of col_vector and matrix classes and related operations.
 *
 * This file defines basic linear algebra utilities, including:
 *  - Construction/destruction and rebuilding of dynamic matrices
 *  - Determinant and minor computation (Laplace expansion)
 *  - Cramer's rule solver and inverse computation
 *  - Matrix product and stream output helpers
 */

#include "matrix.h"

/*============================================**
                 CLASS COL_VECTOR
**============================================*/

/**
 * @brief Set vector to a Cartesian basis vector.
 * @param e Index of the basis element to set to 1 (others set to 0).
 */
void col_vector::cartesian_base( int e ) {
  for( int i=0 ; i<D ; i++ ) {
    if( i == e ) element[i] = 1 ;
    else element[i] = 0 ;
  }
}


/**
 * @brief Assignment operator.
 * @param other Vector to assign from.
 */
void col_vector::operator=( const col_vector &other ) {
  if( D != other.D ) {
    delete [] element ;
    D = other.D ;
    element = new double[D] ;
  } else {
    for( int i=0 ; i<D ; i++ ) element[i] = other.element[i] ;
  }
}


/**
 * @brief Clear vector values (set all elements to 0).
 */
void clear( void ) {
  for( int i=0; i<D; i++ ) element[i] = 0 ;
}

/**
 * @brief Stream output operator for printing.
 * @param out Output stream.
 * @param vec Vector to print.
 * @return Reference to the output stream.
 */
ostream & operator<< ( ostream &out, const col_vector &vec ) {
  out << "(\t" ;
  for( int i=0 ; i<vec.D ; i++ ) out << vec.element[i] << "\t" ;
  out << ")" ;
  return out ;
}



/*============================================**
                 CLASS MATRIX
**============================================*/

/**
 * @brief Default constructor for matrix.
 *
 * Creates an empty (0x0) matrix with no allocated memory.
 */
matrix::matrix( void ) {
  rows = 0 ;
  columns = 0 ;
  element = NULL ;
}

/**
 * @brief Construct a matrix with given dimensions.
 *
 * Allocates memory for an n_rows x n_columns matrix.
 * All elements are initialized to zero.
 *
 * @param n_rows Number of rows.
 * @param n_columns Number of columns.
 */
matrix::matrix( int rowsN , int columnsN ) {
  rows = rowsN ;
  columns = columnsN ;
  if( rows != 0 ) {
    element = new double* [rows] ;
    for( int i=0 ; i<rows ; i++ ) element[i] = new double [columns] ;
  }else{
    element = NULL ;
  }
}

/**
 * @brief Copy constructor for matrix.
 *
 * Performs a deep copy of the given matrix.
 *
 * @param other Matrix to copy.
 */
matrix::matrix( const matrix &other_mat ) {
  rows = other_mat.rows ;
  columns = other_mat.columns ;
  if( rows != 0 ) {
    element = new double* [rows] ;
    for( int i=0; i<rows; i++ ) {
      element[i] = new double [columns] ;
      for( int j=0; j<columns; j++ ) element[i][j] = other_mat.element[i][j] ;
    }
  }else{
    element = NULL ;
  }
}

/**
 * @brief Destructor for matrix.
 *
 * Frees allocated memory.
 */
matrix::~matrix() {
  if( rows != 0 ) {
    for( int i=0; i<rows; i++ ) delete[] element[i] ;
  }
  delete[] element ;
}


/**
 * @brief Rebuild matrix dimensions.
 *
 * Deallocates existing memory and allocates a new
 * matrix with specified size. New elements are initialized to zero.
 *
 * @param n_rows New number of rows.
 * @param n_columns New number of columns.
 */
void matrix::rebuild( int rowsN, int columnsN ) {
  if( rows != 0 ) {
    for( int i=0; i<rows; i++ ) delete[] element[i] ;
  }
  delete[] element ;
  rows = rowsN ;
  columns = columnsN ;
  if( rows != 0 ) {
    element = new double* [rows] ;
    for( int i=0 ; i<rows ; i++ ) element[i] = new double [columns] ;
  }else{
    element = NULL ;
  }
}


/**
 * @brief Determinant of the matrix.
 *
 * calculation of the determinant with the
 * usual recursive method of minors.
 *
 * @return determinant (double).
 */
double matrix::determinant( void ) {
  double determinante = 0 ;
  if( rows != columns ) return 0 ;
  if( rows == 1 ) {
    return element[0][0] ;
  }else{
    for( int i=0; i<rows; i++ ) {
      matrix *minor = NULL ;
      minor = construct_minor( i , 0 ) ;
      determinante += (1-2*(i%2)) * element[i][0] * minor->determinant() ;
      delete minor ;
    }
    return determinante ;
  }
}


/**
 * @brief Minor of the matrix.
 *
 * returns a given minor of the matrix
 * by eliminating the row i and column j.
 *
 * @param i Row index.
 * @param j Column index.
 * @return minor Matrix minor(i,j).
 */
matrix *matrix::construct_minor( int i, int j ) {
  matrix *minor = NULL ;
  minor = new matrix( rows-1 , columns-1 ) ;
  for( int l=0; l<i; l++ ) {
    for( int k=0; k<j; k++ ) {
      minor->element[l][k] = element[l][k] ;
    }
    for( int k=j+1; k<columns; k++ ) {
      minor->element[l][k-1] = element[l][k] ;
    }
  }
  for( int l=i+1; l<rows; l++ ) {
    for( int k=0; k<j; k++ ) {
      minor->element[l-1][k] = element[l][k] ;
    }
    for( int k=j+1; k<columns; k++ ) {
      minor->element[l-1][k-1] = element[l][k] ;
    }
  }

  return minor ;
}


/**
 * @brief Determinant of the Minor (i,j).
 *
 * returns the determinant of the minor constructed
 * by eliminating the row i and column j.
 *
 * @param i Row index.
 * @param j Column index.
 * @return det Minor determinant.
 */
double matrix::minor_determinant( int i, int j ) {
  matrix *minor = NULL ;
  minor = construct_minor( i, j ) ;
  double det = minor->determinant() ;
  delete minor ;
  return det ;
}


/**
 * @brief Cramer solution method.
 *
 * It solves the linear system M*x=b,
 * with the Cramer method.
 *
 * @param constant_vec b vector.
 * @return solution x vector, solution.
 */
col_vector matrix::solve_cramer( const col_vector constant_vec ) {
  double det = determinant() ;
  col_vector solution( constant_vec.D ) ;
  if( columns != rows || rows != constant_vec.D || det == 0 ) {
    cout << "The system is not solvable !" ;
    exit( EXIT_FAILURE ) ;
  }
  matrix appoggio( *this ) ;
  for( int j=0 ; j<constant_vec.D ; j++ ) {
    for( int i=0 ; i<constant_vec.D ; i++ ) {
      appoggio.element[i][j] = constant_vec.element[i] ;
    }
    solution.element[j] = appoggio.determinant() / det ;
    for( int i=0 ; i<constant_vec.D ; i++ ) {
      appoggio.element[i][j] = element[i][j] ;
    }
  }

  return solution ;
}


/**
 * @brief Inverse matrix.
 *
 * returns the inverse of a non-singular matrix.
 *
 * @return matrix inverse.
 */
matrix matrix::compute_inverse( void ) {
  double det = determinant() ;
  if( columns != rows || det == 0 ) {
    cout << "The inverse matrix cannot be computed !" ;
    exit( EXIT_FAILURE ) ;
  }
  matrix inverse( rows , columns ) ;
  col_vector solution( rows ) , cartesian_vec( rows ) ;

  for( int j=0 ; j<columns ; j++ ) {
    cartesian_vec.cartesian_base( j ) ;
    solution = solve_cramer( cartesian_vec ) ;
    for( int i=0 ; i<rows ; i++ ) inverse.element[i][j] = solution.element[i] ;
  }

  return inverse ;
}


/**
 * @brief Output stream.
 *
 * convert the matrix to string output.
 *
 * @return ostream object.
 */
ostream & operator<< ( ostream &out, const matrix &mat ) {
  for( int i=0 ; i<mat.rows ; i++ ) {
    out << "|\t" ;
    for( int j=0 ; j<mat.columns ; j++ ) out << mat.element[i][j] << "\t" ;
    out << "|" << endl ;
  }
  return out ;
};


/**
 * @brief Outer product (rows by columns).
 *
 * returns the matrix M*B result of the outer product rows by columns, if possible.
 *
 * @param B B matrix.
 * @return result M*B dot product.
 */
matrix matrix::dot_prod( const matrix &B ) {
  if( columns != B.rows ) {
    cout << "The outer product cannot be computed, the number of rows and columns of the two matrices does not match !" << endl ;
    exit( EXIT_FAILURE ) ;
  } else {
    matrix result( rows , B.columns ) ;
    for( int i=0 ; i<result.rows ; i++ ) {
      for( int j=0 ; j<result.columns ; j++ ) {
	result.element[i][j] = 0 ;
	for( int k=0 ; k<columns ; k++ ) result.element[i][j] += ( element[i][k] * B.element[k][j] ) ;
      }
    }

    return result ;
  }
}

/**
 * @brief Clear matrix values.
 *
 * Sets all elements to zero, preserving dimensions.
 */
void matrix::clear(void) {
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < columns; j++) element[i][j] = 0.0;
}
