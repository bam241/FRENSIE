//---------------------------------------------------------------------------//
//!
//! \file   Ray.cpp
//! \author Alex Robinson
//! \brief  Ray class definition
//!
//---------------------------------------------------------------------------//

// FACEMC Includes
#include "Ray.hpp"
#include "DirectionHelpers.hpp"
#include "ContractException.hpp"

namespace FACEMC{

// Constructor
Ray::Ray( const double x_position,
	  const double y_position,
	  const double z_position,
	  const double x_direction,
	  const double y_direction,
	  const double z_direction )
  : PrintableObject( "Ray" ),
    d_position( {x_position, y_position, z_position} ),
    d_direction( {x_direction, y_direction, z_direction} )
{
  // Make sure the position is valid
  testPrecondition( !ST::isnaninf( d_position[0] ) );
  testPrecondition( !ST::isnaninf( d_position[1] ) );
  testPrecondition( !ST::isnaninf( d_position[2] ) );
  // Make sure the direction is a unit vector
  testPrecondition( validDirection( x_direction, y_direction, z_direction ) );
}

// Constructor
Ray::Ray( const double position[3],
	  const double direction[3] )
  : PrintableObject( "Ray" ),
    d_position( {position[0], position[1], position[2]} ), // deep copy
    d_direction( {direction[0], direction[1], direction[2]} ) // deep copy
{
  // Make sure the position and direction are valid
  testPrecondition( !ST::isnaninf( d_position[0] ) );
  testPrecondition( !ST::isnaninf( d_position[1] ) );
  testPrecondition( !ST::isnaninf( d_position[2] ) );
  // Make sure the direction is a unit vector
  testPrecondition( validDirection(direction[0], direction[1], direction[2]) );
}

// Return the x position of the ray
double Ray::getXPosition() const
{
  return d_position[0];
}

// Return the y position of the ray
double Ray::getYPosition() const
{
  return d_position[1];
}

// Return the z position of the ray
double Ray::getZPosition() const
{
  return d_position[2];
}

// Return the position of the ray
const double* Ray::getPosition() const
{
  return d_position;
}

// Return the x direction of the ray
double Ray::getXDirection() const
{
  return d_direction[0];
}

// Return the y direction of the ray
double Ray::getYDirection() const
{
  return d_direction[1];
}

// Return the z direction of the ray
double Ray::getZDirection() const
{
  return d_direction[2];
}

// Return the direction of the ray
const double* Ray::getDirection() const
{
  return d_direction;
}

// Advance the head along its direction by the requested distance
void Ray::advanceHead( const double distance )
{
  // Make sure the distance is valid
  testPrecondition( !ST::isnaninf( distance ) );

  d_position[0] += d_direction[0]*distance;
  d_position[1] += d_direction[1]*distance;
  d_position[2] += d_direction[2]*distance;
}

// Print method implementation
void Ray::print( std::ostream& os ) const
{
  os.precision( 16 );
  os << "Position: {" << d_position[0] << "," << d_position[1] << ","
     << d_position[2] << "}" << std::endl;
  os << "Direction: {" << d_direction[0] << "," 
     << d_direction[1] << ","
     << d_direction[2] << "}" << std::endl;
}

} // end FACEMC namespace

//---------------------------------------------------------------------------//
// end Ray.cpp
//---------------------------------------------------------------------------//
