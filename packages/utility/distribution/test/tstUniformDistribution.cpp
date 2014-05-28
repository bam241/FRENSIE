//---------------------------------------------------------------------------//
//!
//! \file   tstUniformDistribution.cpp
//! \author Alex Robinson
//! \brief  Uniform distribution unit tests.
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>

// Trilinos Includes
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

// FRENSIE Includes
#include "Utility_UnitTestHarnessExtensions.hpp"
#include "Utility_OneDDistribution.hpp"
#include "Utility_UniformDistribution.hpp"
#include "Utility_RandomNumberGenerator.hpp"

Teuchos::RCP<Utility::OneDDistribution> distribution( 
			   new Utility::UniformDistribution( -1.0, 1.0, 2.0 ) );

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Check that the distribution can be evaluated
TEUCHOS_UNIT_TEST( UniformDistribution, evaluate )
{  
  TEST_EQUALITY_CONST( distribution->evaluate( -2.0 ), 0.0 );
  TEST_EQUALITY_CONST( distribution->evaluate( -1.0 ), 2.0 );
  TEST_EQUALITY_CONST( distribution->evaluate( 0.0 ), 2.0 );
  TEST_EQUALITY_CONST( distribution->evaluate( 1.0 ), 2.0 );
  TEST_EQUALITY_CONST( distribution->evaluate( -2.0 ), 0.0 );
}

//---------------------------------------------------------------------------//
// Check that the PDF can be evaluated
TEUCHOS_UNIT_TEST( UniformDistribution, evaluatePDF )
{
  TEST_EQUALITY_CONST( distribution->evaluatePDF( -2.0 ), 0.0 );
  TEST_EQUALITY_CONST( distribution->evaluatePDF( -1.0 ), 0.5 );
  TEST_EQUALITY_CONST( distribution->evaluatePDF( 0.0 ), 0.5 );
  TEST_EQUALITY_CONST( distribution->evaluatePDF( 1.0 ), 0.5 );
  TEST_EQUALITY_CONST( distribution->evaluatePDF( 2.0 ), 0.0 );
}

//---------------------------------------------------------------------------//
// Check that the distribution can be sampled
TEUCHOS_UNIT_TEST( UniformDistribution, sample )
{
  std::vector<double> fake_stream( 3 );
  fake_stream[0] = 0.0;
  fake_stream[1] = 0.5;
  fake_stream[2] = 1.0 - 1e-15;
  
  Utility::RandomNumberGenerator::setFakeStream( fake_stream );
  
  double sample = distribution->sample();
  TEST_EQUALITY_CONST( sample, -1.0 );

  sample = distribution->sample();
  TEST_EQUALITY_CONST( sample, 0.0 ); 

  sample = distribution->sample();
  TEST_FLOATING_EQUALITY( sample, 1.0, 1e-14 );

  Utility::RandomNumberGenerator::unsetFakeStream();
}

//---------------------------------------------------------------------------//
// Check that the sampling efficiency can be returned
TEUCHOS_UNIT_TEST( UniformDistribution, getSamplingEfficiency )
{
  TEST_EQUALITY_CONST( distribution->getSamplingEfficiency(), 1.0 );
}

//---------------------------------------------------------------------------//
// Check that the upper bound of the distribution independent variable can be 
// returned
TEUCHOS_UNIT_TEST( UniformDistribution, getUpperBoundOfIndepVar )
{
  TEST_EQUALITY_CONST( distribution->getUpperBoundOfIndepVar(), 1.0 );
}

//---------------------------------------------------------------------------//
// Check that the lower bound of the distribution independent variable can be
// returned
TEUCHOS_UNIT_TEST( UniformDistribution, getLowerBoundOfIndepVar )
{
  TEST_EQUALITY_CONST( distribution->getLowerBoundOfIndepVar(), -1.0 );
}

//---------------------------------------------------------------------------//
// Check that the distribution type can be returned
TEUCHOS_UNIT_TEST( UniformDistribution, getDistributionType )
{
  TEST_EQUALITY_CONST( distribution->getDistributionType(),
		       Utility::UNIFORM_DISTRIBUTION );
}

//---------------------------------------------------------------------------//
// end tstUniformDistribution.cpp
//---------------------------------------------------------------------------//
