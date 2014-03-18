//---------------------------------------------------------------------------//
//!
//! \file   Estimator.hpp
//! \author Alex Robinson
//! \brief  Estimator base class declaration
//!
//---------------------------------------------------------------------------//

#ifndef ESTIMATOR_HPP
#define ESTIMATOR_HPP

// Std Lib Includes
#include <string>
#include <set>
#include <map>

// Boost Includes
#include <boost/unordered_map.hpp>

// Trilinos Includes
#include <Teuchos_Array.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_any.hpp>

// FACEMC Includes
#include "PrintableObject.hpp"
#include "ParticleType.hpp"
#include "Tuple.hpp"
#include "ResponseFunction.hpp"
#include "ContractException.hpp"
#include "PhaseSpaceDimension.hpp"
#include "EstimatorDimensionDiscretization.hpp"

namespace FACEMC{

//! The estimator base class
class Estimator : public PrintableObject
{

protected:

  // Typedef for Teuchos::ScalarTraits
  typedef Teuchos::ScalarTraits<double> ST;

  // Typedef for tuple of estimator moments (1st,2nd,3rd,4th)
  typedef Quad<double,double,double,double> EstimatorMoments;

  // Typedef for the array of estimator moments 
  typedef Teuchos::Array<EstimatorMoments> EstimatorMomentsArray;

public:

  //! Set the number of particle histories that will be simulated
  static void setNumberOfHistories( const unsigned long long num_histories );
  
  //! Set the start time for the figure of merit calculation
  static void setStartTime( const double start_time );
  
  //! Set the end time for the figure of merit calculation
  static void setEndTime( const double end_time );

  //! Constructor
  Estimator( const unsigned long long id,
	     const double multiplier );

  //! Destructor
  virtual ~Estimator()
  { /* ... */ }

  //! Return the estimator id
  unsigned long long getId() const;

  //! Set the bin boundaries for a dimension of the phase space (floating pt)
  template<PhaseSpaceDimension dimension, typename DimensionType>
  void setBinBoundaries( const Teuchos::Array<DimensionType>& bin_boundaries );
  
  //! Return the number of bins for a dimension of the phase space
  unsigned getNumberOfBins( const PhaseSpaceDimension dimension ) const;

  //! Return the total number of bins
  unsigned getNumberOfBins() const;

  //! Set the response functions
  virtual void setResponseFunctions( 
   const Teuchos::Array<Teuchos::RCP<ResponseFunction> >& response_functions );

  //! Return the number of response functions
  unsigned getNumberOfResponseFunctions() const;

  //! Set the particle types that can contribute to the estimator
  virtual void setParticleTypes( 
			  const Teuchos::Array<ParticleType>& particle_types );

  //! Check if the particle type is assigned to the estimator
  bool isParticleTypeAssigned( const ParticleType particle_type ) const;

  //! Commit the contribution from the current history to the estimator
  virtual void commitHistoryContribution() = 0;

protected:

  //! Assign bin boundaries to an estimator dimension
  virtual void assignBinBoundaries( 
	const Teuchos::RCP<EstimatorDimensionDiscretization>& bin_boundaries );

  //! Return the estimator constant multiplier
  double getMultiplier() const;

  //! Return the response function name
  const std::string& getResponseFunctionName( 
				const unsigned response_function_index ) const;

  //! Print the estimator response function names
  void printEstimatorResponseFunctionNames( std::ostream& os ) const;

  //! Print the estimator bins
  void printEstimatorBins( std::ostream& os ) const;

  //! Print the estimator data stored in an array
  void printEstimatorBinData( 
			    std::ostream& os,
			    const EstimatorMomentsArray& estimator_moment_data,
			    const double norm_constant ) const;

  //! Print the total estimator data stored in an array
  void printEstimatorTotalData( 
		     std::ostream& os,
		     const EstimatorMomentsArray& total_estimator_moments_data,
		     const double norm_constant ) const;

  //! Evaluate the desired response function
  double evaluateResponseFunction( 
				const BasicParticleState& particle,
				const unsigned response_function_index ) const;

  //! Check if the point is in the estimator phase space
  bool isPointInEstimatorPhaseSpace( 
		 const PhaseSpace::DimensionValueMap& dimension_values ) const;
			        

  //! Check if the point is in the estimator phase space
  bool isPointInEstimatorPhaseSpace(
	       const PhaseSpace::DimensionValueArray& dimension_values ) const;

  //! Calculate the bin index for the desired response function
  unsigned calculateBinIndex( 
			 const PhaseSpace::DimensionValueMap& dimension_values,
			 const unsigned response_function_index ) const;

  //! Calculate the bin index for the desired response function
  unsigned calculateBinIndex( 
		       const PhaseSpace::DimensionValueArray& dimension_values,
		       const unsigned response_function_index ) const;

  //! Calculate the mean of a set of contributions
  double calculateMean( const double first_moment_contributions ) const;

  //! Calculate the relative error of a set of contributions
  double calculateRelativeError( 
			      const double first_moment_contributions,
			      const double second_moment_contributions ) const;

  //! Calculate the variance of the variance (VOV) of a set of contributions
  double calculateVOV( const double first_moment_contributions,
		       const double second_moment_contributions,
		       const double third_moment_contributions,
		       const double fourth_moment_contributions ) const;

  //! Calculate the figure of merit (FOM) of an estimator bin
  double calculateFOM( const double relative_error ) const;

private:

  // The number of particle histories that will be run
  static unsigned long long num_histories;

  // The start time used for the figure of merit calculation
  static double start_time;

  // The end time used for the figure of merit calculation
  static double end_time;

  // The estimator id
  unsigned long long d_id;

  // The constant multiplier for the estimator
  double d_multiplier;

  // The response functions
  Teuchos::Array<Teuchos::RCP<ResponseFunction> > d_response_functions;
  
  // The estimator phase space dimension bin boundaries map
  boost::unordered_map<PhaseSpaceDimension,
  		       Teuchos::RCP<EstimatorDimensionDiscretization> >
  d_dimension_bin_boundaries_map;

  // The estimator phase space dimension index step size map
  boost::unordered_map<PhaseSpaceDimension,unsigned>
  d_dimension_index_step_size_map;

  // The estimator phase space dimension ordering
  Teuchos::Array<PhaseSpaceDimension> d_dimension_ordering;

  // The particle types that this estimator will take contributions from
  std::set<ParticleType> d_particle_types;
};

// Return the estimator id
inline unsigned long long Estimator::getId() const
{
  return d_id;
}

// Return the estimator constant multiplier
inline double Estimator::getMultiplier() const
{
  return d_multiplier;
}

// Return the number of bins for a dimension of the phase space
inline unsigned Estimator::getNumberOfBins( 
				    const PhaseSpaceDimension dimension ) const
{
  if( d_dimension_bin_boundaries_map.count( dimension ) != 0 )
    return d_dimension_bin_boundaries_map.find(dimension)->second->getNumberOfBins();
  else
    return 1u;
}

// Return the total number of bins
inline unsigned Estimator::getNumberOfBins() const
{
  unsigned number_of_bins = 1u;
  
  for( unsigned i = 0u; i < d_dimension_ordering.size(); ++i )
    number_of_bins *= getNumberOfBins( d_dimension_ordering[i] );
  
  return number_of_bins;
}

// Return the number of response functions
inline unsigned Estimator::getNumberOfResponseFunctions() const
{
  return d_response_functions.size();
}

// Check if the particle type is assigned to the estimator
inline bool Estimator::isParticleTypeAssigned( 
					const ParticleType particle_type) const
{
  return d_particle_types.count( particle_type );
}

// Return the response function name
inline const std::string& Estimator::getResponseFunctionName( 
				 const unsigned response_function_index ) const
{
  // Make sure the response function index is valid
  testPrecondition( response_function_index < getNumberOfResponseFunctions() );

  return d_response_functions[response_function_index]->getName();
}

// Evaluate the desired response function
inline double Estimator::evaluateResponseFunction( 
				 const BasicParticleState& particle,
				 const unsigned response_function_index ) const
{
  // Make sure the response function index is valid
  testPrecondition( response_function_index < getNumberOfResponseFunctions() );
  
  return d_response_functions[response_function_index]->evaluate( particle );
}

} // end FACEMC namespace

//---------------------------------------------------------------------------//
// Template Includes
//---------------------------------------------------------------------------//

#include "Estimator_def.hpp"

//---------------------------------------------------------------------------//

#endif // end ESTIMATOR_HPP

//---------------------------------------------------------------------------//
// end Estimator.hpp
//---------------------------------------------------------------------------//