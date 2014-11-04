//---------------------------------------------------------------------------//
//!
//! \file   Utility_TetrahedronHelpers_def.hpp
//! \author Alex Robinson, Eli Moll
//! \brief  Tetrahedron helper function template definitions
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_TETRAHEDRON_HELPERS_DEF_HPP
#define UTILITY_TETRAHEDRON_HELPERS_DEF_HPP

// Trilinos Includes
#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_BLAS_types.hpp>

// FRENSIE Includes
#include "Utility_ContractException.hpp"

namespace Utility{

// Calculate tetrahedron barycentric transform matrix                        
template<typename Matrix>                                                      
void calculateBarycentricTransformMatrix( const double vertex_a[3],    
					  const double vertex_b[3],    
					  const double vertex_c[3],    
					  const double reference_vertex[3],    
					  Matrix& matrix ) 
{
  // Make sure the matrix is valid
  testPrecondition( matrix.numRows() == 3 );
  testPrecondition( matrix.numCols() == 3 );
  
  matrix( 0, 0 ) = vertex_a[0] - reference_vertex[0];
  matrix( 0, 1 ) = vertex_b[0] - reference_vertex[0];
  matrix( 0, 2 ) = vertex_c[0] - reference_vertex[0];
  matrix( 1, 0 ) = vertex_a[1] - reference_vertex[1];
  matrix( 1, 1 ) = vertex_b[1] - reference_vertex[1];
  matrix( 1, 2 ) = vertex_c[1] - reference_vertex[1];
  matrix( 2, 0 ) = vertex_a[2] - reference_vertex[2];
  matrix( 2, 1 ) = vertex_b[2] - reference_vertex[2];
  matrix( 2, 2 ) = vertex_c[2] - reference_vertex[2];

  Teuchos::SerialDenseSolver<typename Matrix::ordinalType,
			     typename Matrix::scalarType> solver;
  solver.setMatrix( Teuchos::rcp<Matrix>( &matrix, false ) );
  
  int return_value = solver.invert();

  // Make sure the tet is valid
  testPrecondition( return_value == 0 );
}

// Determine if a point is in a given tet                        
template<typename Matrix>                                                      
bool isPointInTet( const double point[3],    
		   Matrix& matrix ) 
{
  // Make sure the matrix is valid
  testPrecondition( matrix.numRows() == 3 );
  testPrecondition( matrix.numCols() == 3 );
  
  Teuchos::SerialDenseMatrix<int,double> point_vector(3,1);
  Teuchos::SerialDenseMatrix<int,double> barycentric_location_vector(3,1);
  point_vector( 0, 0 ) = point[0];
  point_vector( 1, 0 ) = point[1];
  point_vector( 2, 0 ) = point[2];
  
  double tolerance = -1e-12;
  
  barycentric_location_vector.multiply(
             Teuchos::NO_TRANS,Teuchos::NO_TRANS, 1.0,matrix,point_vector,1.0);
  
  if ( barycentric_location_vector( 0, 0 ) <= tolerance ||
       barycentric_location_vector( 0, 1 ) <= tolerance ||
       barycentric_location_vector( 0, 2 ) <= tolerance )
  {
    bool point_in_tet = false;
    return point_in_tet;
  }
  else
  {
    bool point_in_tet = true;
    return point_in_tet;
  };
}

} // end Utility namespace

#endif // end UTILITY_TETRAHEDRON_HELPERS_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_TetrahedronHelpers_def.hpp
//---------------------------------------------------------------------------//
