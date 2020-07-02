//---------------------------------------------------------------------------//
//!
//! \file   Utility_PQLAQuadrature.hpp
//! \author Philip Britt
//! \brief  PQLA Direction Quadrature handler declaration
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_PQLA_QUADRATURE
#define UTILITY_PQLA_QUADRATURE

// FRENSIE includes
#include "Utility_Vector.hpp"
#include "Utility_Tuple.hpp"

namespace Utility{

class PQLAQuadrature
{

  public:

  //! Constructor
  PQLAQuadrature(unsigned quadrature_order);

  //! Destructor
  ~PQLAQuadrature()
  { /* ... */ }

  //! Find which triangle bin a direction vector is in
  unsigned findTriangleBin( const double[3]& direction) const;

  //! Find which triangle bin a direction vector is in
  unsigned findTriangleBin( const double x_direction, const double y_direction, const double z_direction) const;

  private:

  //! Take lower bounding plane indices of direction vector to form triangle index
  unsigned calculatePositiveTriangleBinIndex(const unsigned i_x, const unsigned i_y, const unsigned i_z) const;

  //! Take direction signs to calculate secondary index
  unsigned findSecondaryIndex(const bool x_sign, const bool y_sign, const bool z_sign) const;

  //! Quadrature order
  unsigned d_quadrature_order;

  //! Vector that holds the planes for the x,y,z direction
  std::vector<double> d_planes;

  //! Tuple that contains all relevant bin information (planes, nodes)
  std::vector<std::tuple<double[3], std::vector<double[3]>>> d_triangle_parameters

};

} // end Utility namespace

#endif // end UTILITY_PQLA_QUADRATURE

//---------------------------------------------------------------------------//
// end Utility_PQLADiscetization.hpp
//---------------------------------------------------------------------------//