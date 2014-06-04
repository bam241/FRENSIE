//---------------------------------------------------------------------------//
//!
//! \file   Utility_DirectionalDistributionFactory.hpp
//! \author Alex Robinson
//! \brief  Directional distribution factory class declaration
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_DIRECTIONAL_DISTRIBUTION_FACTORY_HPP
#define UTILITY_DIRECTIONAL_DISTRIBUTION_FACTORY_HPP

// Std Lib Includes
#include <stdexcept>

// Trilinos Includes
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

// FRENSIE Includes
#include "Utility_DirectionalDistribution.hpp"
#include "Utility_OneDDistributionEntryConverterDB.hpp"
#include "Utility_Axis.hpp"

namespace Utility{

//! The directional distribution factory class
class DirectionalDistributionFactory
{

public:

  //! Create the directional distribution represented by the parameter list
  static Teuchos::RCP<DirectionalDistribution>
  createDistribution( const Teuchos::ParameterList& distribution_rep );

private:

  // Validate a distribution representation
  static void validateDistributionRep( 
			      const Teuchos::ParameterList& distribution_rep );

  // Validate the axis name
  static void validateAxisName( const std::string& axis_name );

  // Constructor
  DirectionalDistributionFactory();
};

//! The invalid directional distribution representation error
class InvalidDirectionalDistributionRepresentation : public std::logic_error
{

public:

  InvalidDirectionalDistributionRepresentation( const std::string& what_arg )
    : std::logic_error( what_arg )
  { /* ... */ }
};

#endif // end UTILITY_DIRECTIONAL_DISTRIBUTION_FACTORY_HPP

//---------------------------------------------------------------------------//
// end Utility_DirectionalDistributionFactory.hpp
//---------------------------------------------------------------------------//
