//---------------------------------------------------------------------------//
//!
//! \file   tstSpatialCoordinateSystemType.cpp
//! \author Alex Robinson
//! \brief  Spatial coordinate system type unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>
#include <sstream>
#include <string>

// Trilinos Includes
#include <Teuchos_UnitTestHarness.hpp>

// FRENSIE Includes
#include "Utility_UnitTestHarnessExtensions.hpp"
#include "Utility_SpatialCoordinateSystemType.hpp"

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Check that the spatial coordinate system enums can be converted to strings
TEUCHOS_UNIT_TEST( SpatialCoordinateSystemType,
                   conertSpatialCoordinateystemTypeEnumToString )
{
  std::string coord_system_name =
    Utility::convertSpatialCoordinateSystemTypeEnumToString(
                                Utility::CARTESIAN_SPATIAL_COORDINATE_SYSTEM );

  TEST_EQUALITY_CONST( coord_system_name,
                       "Cartesian Spatial Coordinate System" );

  coord_system_name =
    Utility::convertSpatialCoordinateSystemTypeEnumToString(
                              Utility::CYLINDRICAL_SPATIAL_COORDINATE_SYSTEM );

  TEST_EQUALITY_CONST( coord_system_name,
                       "Cylindrical Spatial Coordinate System" );

  coord_system_name =
    Utility::convertSpatialCoordinateSystemTypeEnumToString(
                                Utility::SPHERICAL_SPATIAL_COORDINATE_SYSTEM );

  TEST_EQUALITY_CONST( coord_system_name,
                       "Spherical Spatial Coordinate System" );
}

//---------------------------------------------------------------------------//
// end tstSpatialCoordinateSystemType.cpp
//---------------------------------------------------------------------------//
