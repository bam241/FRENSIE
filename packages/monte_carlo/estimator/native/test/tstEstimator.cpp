//---------------------------------------------------------------------------//
//!
//! \file   tstEstimator.cpp
//! \author Alex Robinson
//! \brief  Estimator unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>

// Trilinos Includes
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

// FRENSIE Includes
#include "MonteCarlo_UnitTestHarnessExtensions.hpp"
#include "MonteCarlo_Estimator.hpp"
#include "MonteCarlo_EnergySpaceResponseFunction.hpp"
#include "MonteCarlo_PhaseSpaceDimension.hpp"
#include "MonteCarlo_PhotonState.hpp"
#include "Utility_UniformDistribution.hpp"

//---------------------------------------------------------------------------//
// Testing Structs.
//---------------------------------------------------------------------------//
class TestEstimator : public MonteCarlo::Estimator
{
public:
  TestEstimator( const unsigned long long id,
		 const double multiplier )
    : MonteCarlo::Estimator( id, multiplier )
  { /* ... */ }

  ~TestEstimator()
  { /* ... */ }

  void print( std::ostream& os ) const
  { 
    printEstimatorResponseFunctionNames( os );
    printEstimatorBins( os );
  }

  void enableThreadSupport( const unsigned num_threads )
  { /* ... */ }

  void commitHistoryContribution()
  { /* ... */ }

  // Allow public access to the estimator protected member functions
  using MonteCarlo::Estimator::DimensionValueMap;
  using MonteCarlo::Estimator::setHasUncommittedHistoryContribution;
  using MonteCarlo::Estimator::unsetHasUncommittedHistoryContribution;
  using MonteCarlo::Estimator::assignBinBoundaries;
  using MonteCarlo::Estimator::getMultiplier;
  using MonteCarlo::Estimator::getResponseFunctionName;
  using MonteCarlo::Estimator::evaluateResponseFunction;
  using MonteCarlo::Estimator::isPointInEstimatorPhaseSpace;
  using MonteCarlo::Estimator::calculateBinIndex;
  using MonteCarlo::Estimator::processMoments;
};

//---------------------------------------------------------------------------//
// Global Testing Variables.
//---------------------------------------------------------------------------//
TestEstimator estimator( 0ull, 1.0 );

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Check that the estimator id can be returned
TEUCHOS_UNIT_TEST( Estimator, getId )
{
  TEST_EQUALITY_CONST( estimator.getId(), 0ull );
}

//---------------------------------------------------------------------------//
// Check that the multiplier can be returned
TEUCHOS_UNIT_TEST( Estimator, getMultiplier )
{
  TEST_EQUALITY_CONST( estimator.getMultiplier(), 1.0 );
}

//---------------------------------------------------------------------------//
// Check that the energy bins can be set
TEUCHOS_UNIT_TEST( Estimator, setEnergyBinBoundaries )
{
  Teuchos::Array<double> energy_bin_boundaries( 7 );
  energy_bin_boundaries[0] = 0.0;
  energy_bin_boundaries[1] = 1e-1;
  energy_bin_boundaries[2] = 1e-1;
  energy_bin_boundaries[3] = 1.0;
  energy_bin_boundaries[4] = 10.0;
  energy_bin_boundaries[5] = 10.0;
  energy_bin_boundaries[6] = 20.0;

  estimator.setBinBoundaries<MonteCarlo::ENERGY_DIMENSION>( energy_bin_boundaries);

  TEST_EQUALITY_CONST( estimator.getNumberOfBins( MonteCarlo::ENERGY_DIMENSION ), 
		       6u );
}

//---------------------------------------------------------------------------//
// Check that cosine bins can be set
TEUCHOS_UNIT_TEST( Estimator, setCosineBinBoundaries )
{
  Teuchos::Array<double> cosine_bin_boundaries( 4 );
  cosine_bin_boundaries[0] = -1.0;
  cosine_bin_boundaries[1] = -1.0/3.0;
  cosine_bin_boundaries[2] = 1.0/3.0;
  cosine_bin_boundaries[3] = 1.0;
  
  estimator.setBinBoundaries<MonteCarlo::COSINE_DIMENSION>( cosine_bin_boundaries);
  
  TEST_EQUALITY_CONST( estimator.getNumberOfBins( MonteCarlo::COSINE_DIMENSION ), 
		       3u );
}

//---------------------------------------------------------------------------//
// Check that the time bins can be set
TEUCHOS_UNIT_TEST( Estimator, setTimeBinBoundaries )
{
  Teuchos::Array<double> time_bin_boundaries( 4 );
  time_bin_boundaries[0] = 0.0;
  time_bin_boundaries[1] = 1e3;
  time_bin_boundaries[2] = 1e5;
  time_bin_boundaries[3] = 1e7;

  estimator.setBinBoundaries<MonteCarlo::TIME_DIMENSION>( time_bin_boundaries );
  
  TEST_EQUALITY_CONST( estimator.getNumberOfBins( MonteCarlo::TIME_DIMENSION ), 
		       3u );
}

//---------------------------------------------------------------------------//
// Check that the collision number bins can be set
TEUCHOS_UNIT_TEST( Estimator, setCollisionNumberBins )
{
  Teuchos::Array<unsigned> collision_number_bins( 4 );
  collision_number_bins[0] = 0u;
  collision_number_bins[1] = 1u;
  collision_number_bins[2] = 2u;
  collision_number_bins[3] = std::numeric_limits<unsigned>::max();

  estimator.setBinBoundaries<MonteCarlo::COLLISION_NUMBER_DIMENSION>( 
						       collision_number_bins );

  TEST_EQUALITY_CONST( estimator.getNumberOfBins( MonteCarlo::COLLISION_NUMBER_DIMENSION ), 4u );
}

//---------------------------------------------------------------------------//
// Check that the total number of bins can be returned
TEUCHOS_UNIT_TEST( Estimator, getNumberOfBins )
{
  TEST_EQUALITY_CONST( estimator.getNumberOfBins(), 216u );
}

//---------------------------------------------------------------------------//
// Check that the response functions can be set
TEUCHOS_UNIT_TEST( Estimator, setResponseFunctions )
{
  Teuchos::Array<Teuchos::RCP<MonteCarlo::ResponseFunction> > 
    response_functions( 2 );
  
  Teuchos::RCP<Utility::OneDDistribution> energy_distribution(
			    new Utility::UniformDistribution( 0.0, 10, 1.0 ) );

  response_functions[0].reset( new MonteCarlo::EnergySpaceResponseFunction( 
						       0,
						       "uniform_energy",
						       energy_distribution ) );
  response_functions[1] = MonteCarlo::ResponseFunction::default_response_function;

  estimator.setResponseFunctions( response_functions );

  TEST_EQUALITY_CONST( estimator.getNumberOfResponseFunctions(), 2u );
}

//---------------------------------------------------------------------------//
// Check that the response function names can be returned
TEUCHOS_UNIT_TEST( Estimator, getResponseFunctionNames )
{
  TEST_EQUALITY_CONST( estimator.getResponseFunctionName( 0u ),
		       "uniform_energy" );
  TEST_EQUALITY_CONST( estimator.getResponseFunctionName( 1u ),
		       "default" );
}

//---------------------------------------------------------------------------//
// Check that the particle types that can contribute to the estimator can
// be set
TEUCHOS_UNIT_TEST( Estimator, setParticleTypes )
{
  Teuchos::Array<MonteCarlo::ParticleType> particle_types( 2 );
  particle_types[0] = MonteCarlo::PHOTON;
  particle_types[1] = MonteCarlo::NEUTRON;

  estimator.setParticleTypes( particle_types );

  TEST_ASSERT( estimator.isParticleTypeAssigned( MonteCarlo::PHOTON ) );
  TEST_ASSERT( estimator.isParticleTypeAssigned( MonteCarlo::NEUTRON ) );
  TEST_ASSERT( !estimator.isParticleTypeAssigned( MonteCarlo::ADJOINT_PHOTON ) );
  TEST_ASSERT( !estimator.isParticleTypeAssigned( MonteCarlo::ADJOINT_NEUTRON ) );
}

//---------------------------------------------------------------------------//
// Check that the response functions can be evaluated
TEUCHOS_UNIT_TEST( Estimator, evaluateResponseFunction )
{
  MonteCarlo::PhotonState particle( 0ull );
  particle.setEnergy( 1.0 );

  double response_function_value = 
    estimator.evaluateResponseFunction( particle, 0u );

  TEST_EQUALITY_CONST( response_function_value, 1.0 );

  response_function_value = 
    estimator.evaluateResponseFunction( particle, 1u );

  TEST_EQUALITY_CONST( response_function_value, 1.0 );
}

//---------------------------------------------------------------------------//
// Check if a point is in the estimator phase space
TEUCHOS_UNIT_TEST( Estimator, isPointInEstimatorPhaseSpace )
{
  TestEstimator::DimensionValueMap dimension_values;
  dimension_values[MonteCarlo::ENERGY_DIMENSION] = Teuchos::any( 0.0 );
  dimension_values[MonteCarlo::COSINE_DIMENSION] = Teuchos::any( -1.0 );
  dimension_values[MonteCarlo::TIME_DIMENSION] = Teuchos::any( 0.0 );
  dimension_values[MonteCarlo::COLLISION_NUMBER_DIMENSION] = Teuchos::any( 0u );

  TEST_ASSERT( estimator.isPointInEstimatorPhaseSpace( dimension_values ) );

  dimension_values[MonteCarlo::ENERGY_DIMENSION] = Teuchos::any( 20.0 );
  dimension_values[MonteCarlo::COSINE_DIMENSION] = Teuchos::any( 1.0 );
  dimension_values[MonteCarlo::TIME_DIMENSION] = Teuchos::any( 1e7 );
  dimension_values[MonteCarlo::COLLISION_NUMBER_DIMENSION] = 
    Teuchos::any( std::numeric_limits<unsigned>::max() );

  TEST_ASSERT( estimator.isPointInEstimatorPhaseSpace( dimension_values ) );

  dimension_values[MonteCarlo::ENERGY_DIMENSION] = Teuchos::any( 21.0 );

  TEST_ASSERT( !estimator.isPointInEstimatorPhaseSpace( dimension_values ) );

  dimension_values[MonteCarlo::ENERGY_DIMENSION] = Teuchos::any( 20.0 );
  dimension_values[MonteCarlo::TIME_DIMENSION] = Teuchos::any( 2e7 );

  TEST_ASSERT( !estimator.isPointInEstimatorPhaseSpace( dimension_values ) );
}

//---------------------------------------------------------------------------//
// Check that the bin index for the desired response function can be 
// calculated
TEUCHOS_UNIT_TEST( Estimator, calculateBinIndex )
{
  TestEstimator::DimensionValueMap dimension_values;
  dimension_values[MonteCarlo::ENERGY_DIMENSION] = Teuchos::any( 0.0 );
  dimension_values[MonteCarlo::COSINE_DIMENSION] = Teuchos::any( -1.0 );
  dimension_values[MonteCarlo::TIME_DIMENSION] = Teuchos::any( 0.0 );
  dimension_values[MonteCarlo::COLLISION_NUMBER_DIMENSION] = Teuchos::any( 0u );
  
  unsigned bin_index = estimator.calculateBinIndex( dimension_values, 0u );

  TEST_EQUALITY_CONST( bin_index, 0u );

  bin_index = estimator.calculateBinIndex( dimension_values, 1u );

  TEST_EQUALITY_CONST( bin_index, 216u );
						    
  dimension_values[MonteCarlo::ENERGY_DIMENSION] = Teuchos::any( 10.0 );
  dimension_values[MonteCarlo::COSINE_DIMENSION] = Teuchos::any( 0.0 );
  dimension_values[MonteCarlo::TIME_DIMENSION] = Teuchos::any( 1e6 );
  dimension_values[MonteCarlo::COLLISION_NUMBER_DIMENSION] = Teuchos::any( 2u );
  
  bin_index = estimator.calculateBinIndex( dimension_values, 0u );

  TEST_EQUALITY_CONST( bin_index, 154u );

  bin_index = estimator.calculateBinIndex( dimension_values, 1u );

  TEST_EQUALITY_CONST( bin_index, 370u );

  dimension_values[MonteCarlo::ENERGY_DIMENSION] = Teuchos::any( 20.0 );
  dimension_values[MonteCarlo::COSINE_DIMENSION] = Teuchos::any( 1.0 );
  dimension_values[MonteCarlo::TIME_DIMENSION] = Teuchos::any( 1e7 );
  dimension_values[MonteCarlo::COLLISION_NUMBER_DIMENSION] = 
    Teuchos::any( std::numeric_limits<unsigned>::max() );

  bin_index = estimator.calculateBinIndex( dimension_values, 0u );

  TEST_EQUALITY_CONST( bin_index, 215u );

  bin_index = estimator.calculateBinIndex( dimension_values, 1u );

  TEST_EQUALITY_CONST( bin_index, 431u );
}

//---------------------------------------------------------------------------//
// Check if the estimator has an uncommitted history contribution
TEUCHOS_UNIT_TEST( Estimator, hasUncommittedHisotryContribution )
{
  TEST_ASSERT( !estimator.hasUncommittedHistoryContribution() );
  
  estimator.setHasUncommittedHistoryContribution();

  TEST_ASSERT( estimator.hasUncommittedHistoryContribution() );

  estimator.unsetHasUncommittedHistoryContribution();

  TEST_ASSERT( !estimator.hasUncommittedHistoryContribution() );
}

//---------------------------------------------------------------------------//
// Process the first and second moments
TEUCHOS_UNIT_TEST( Estimator, processMoments_two )
{
  MonteCarlo::Estimator::setNumberOfHistories( 100ull );

  double mean;
  double relative_error;

  Utility::Pair<double,double> moments( 100.0, 150.0 );

  estimator.processMoments( moments,
			    1.0,
			    mean,
			    relative_error );

  TEST_EQUALITY_CONST( mean, 1.0 );
  TEST_FLOATING_EQUALITY( relative_error, 0.070710678118655, 1e-14 );
}

//---------------------------------------------------------------------------//
// Process the first, second, third and fourth moments
TEUCHOS_UNIT_TEST( Estimator, processMoments_four )
{
  MonteCarlo::Estimator::setNumberOfHistories( 100ull );
  MonteCarlo::Estimator::setStartTime( 0.0 );
  MonteCarlo::Estimator::setEndTime( 1.0 );

  double mean;
  double relative_error;
  double variance_of_variance;
  double figure_of_merit;

  Utility::Quad<double,double,double,double> moments( 100.0,
						      150.0,
						      300.0,
						      800.0 );

  estimator.processMoments( moments,
			    1.0,
			    mean,
			    relative_error,
			    variance_of_variance,
			    figure_of_merit );
  
  TEST_EQUALITY_CONST( mean, 1.0 );
  TEST_FLOATING_EQUALITY( relative_error, 0.070710678118655, 1e-14 );
  TEST_EQUALITY_CONST( variance_of_variance, 0.07 );
  TEST_EQUALITY_CONST( figure_of_merit, 200.0 );
}

//---------------------------------------------------------------------------//
// end tstEstimator.cpp
//---------------------------------------------------------------------------//
