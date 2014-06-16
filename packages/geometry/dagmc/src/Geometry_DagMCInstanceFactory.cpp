//---------------------------------------------------------------------------//
//!
//! \file   Geometry_DagMCInstanceFactory.hpp
//! \author Alex Robinson
//! \brief  Dagmc instance factory class definition.
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <map>

// FRENSIE Includes
#include "Geometry_DagMCInstanceFactory.hpp"
#include "Geometry_DagMCHelpers.hpp"
#include "Geometry_DagMCProperties.hpp"
#include "Utility_ContractException.hpp"
#include "Utility_ExceptionTestMacros.hpp"

namespace Geometry{

// Initialize DagMC
void DagMCInstanceFactory::initializeDagMC( 
				       const Teuchos::ParameterList& geom_rep )
{
  // Validate the geometry representation
  DagMCInstanceFactory::validateGeometryRep( geom_rep );

  // Get the CAD file name
  std::string cad_file_name = geom_rep.get<std::string>( "CAD File" );

  // Get the face tolerance
  double facet_tol = 1e-3; // default value

  if( geom_rep.isParameter( "Facet Tolerance" ) )
    facet_tol = geom_rep.get<double>( "Facet Tolerance" );

  // Get the property synonyms
  std::vector<std::string> properties( 9 );

  if( geom_rep.isParameter( "Termination Cell Synonym" ) )
    properties[0] = geom_rep.get<std::string>( "Termination Cell Synonym" );   
  else
    properties[0] = DagMCProperties::termination_cell_property;

  if( geom_rep.isParameter( "Material Synonym" ) )
    properties[1] = geom_rep.get<std::string>( "Material Synonym" );
  else
    properties[1] = DagMCProperties::material_property;

  if( geom_rep.isParameter( "Density Synonym" ) )
    properties[2] = geom_rep.get<std::string>( "Density Synonym" );
  else
    properties[2] = DagMCProperties::density_property;

  if( geom_rep.isParameter( "Estimator Synonym" ) )
    properties[3] = geom_rep.get<std::string>( "Estimator Synonym" );
  else
    properties[3] = DagMCProperties::estimator_property;

  if( geom_rep.isParameter( "Cell Track Length Flux Estimator Synonym" ) )
    properties[4] = geom_rep.get<std::string>( "Cell Track Length Flux Estimator Synonym" );
  else
    properties[4] = DagMCProperties::cell_track_length_flux_property;	     
  
  if( geom_rep.isParameter( "Cell Collision Flux Estimator Synonym" ) )
    properties[5] = geom_rep.get<std::string>( "Cell Collision Flux Estimator Synonym" );
  else
    properties[5] = DagMCProperties::cell_collision_flux_property;

  if( geom_rep.isParameter( "Cell Pulse Height Estimator Synonym" ) )
    properties[6] = geom_rep.get<std::string>( "Cell Pulse Height Estimator Synonym" );
  else
    properties[6] = DagMCProperties::cell_pulse_height_property;

  if( geom_rep.isParameter( "Surface Flux Estimator Synonym" ) )
    properties[7] = geom_rep.get<std::string>( "Surface Flux Estimator Synonym" );
  else
    properties[7] = DagMCProperties::surface_flux_property;

  if( geom_rep.isParameter( "Surface Current Estimator Synonym" ) )
    properties[8] = geom_rep.get<std::string>( "Surface Current Estimator Synonym" );
  else
    properties[8] = DagMCProperties::surface_current_property;

  // Initialize DagMC
  Geometry::initializeDagMC( cad_file_name, properties, facet_tol );

  // Print the unused parameters
  geom_rep.unused( std::cout );
}

// Validate a geometry representation
void DagMCInstanceFactory::validateGeometryRep( 
				       const Teuchos::ParameterList& geom_rep )
{
  testPrecondition( geom_rep.get<std::string>( "Handler" ) == "DagMC" );
  
  TEST_FOR_EXCEPTION( !geom_rep.isParameter( "CAD File" ),
		      InvalidGeometryRepresentation,
		      "Error: The CAD file needs to be specified!" );
}

} // end Geometry namespace

//---------------------------------------------------------------------------//
// end Geometry_DagMCInstanceFactory.cpp
//---------------------------------------------------------------------------//
