//---------------------------------------------------------------------------//
//!
//! \file   Utility_InterpolatedTabularBasicBivariateDistributionImplBase_def.hpp
//! \author Alex Robinson
//! \brief  The interpolated tabular two-dimensional dist. helper class defs.
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_INTERPOLATED_TABULAR_BASIC_BIVARIATE_DISTRIBUTION_IMPL_BASE_DEF_HPP
#define UTILITY_INTERPOLATED_TABULAR_BASIC_BIVARIATE_DISTRIBUTION_IMPL_BASE_DEF_HPP

// FRENSIE Includes
#include "Utility_PartiallyTabularBasicBivariateDistribution.hpp"
#include "Utility_FullyTabularBasicBivariateDistribution.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_ExceptionTestMacros.hpp"
#include "Utility_ContractException.hpp"

namespace Utility{

// Constructor
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase(
     const std::vector<PrimaryIndepQuantity>& primary_indep_grid,
     const std::vector<std::shared_ptr<const BaseUnivariateDistributionType> >&
     secondary_distributions,
     const double fuzzy_boundary_tol,
     const double relative_error_tol,
     const double error_tol )
  : BaseType( primary_indep_grid, secondary_distributions ),
    d_fuzzy_boundary_tol( fuzzy_boundary_tol ),
    d_relative_error_tol( relative_error_tol ),
    d_error_tol( error_tol )
{
  // Verify that the distribution data is valid
  this->verifyValidData( primary_indep_grid, secondary_distributions );
  this->verifyValidTolerances( fuzzy_boundary_tol,
                               relative_error_tol,
                               error_tol );
}

// Set the distribution
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
void UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::setDistribution(
     const std::vector<PrimaryIndepQuantity>& primary_indep_grid,
     const std::vector<std::shared_ptr<const BaseUnivariateDistributionType> >&
     secondary_distributions )
{
  // Set the distribution
  Distribution::setDistribution( primary_indep_grid, secondary_distributions );
  
  // Verify that the distribution data is valid
  this->verifyValidData( primary_indep_grid, secondary_distributions );
}

// Set the evaluation tolerances
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
void UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::setEvaluationTolerances(
                                               const double fuzzy_boundary_tol,
                                               const double relative_error_tol,
                                               const double error_tol )
{  
  d_fuzzy_boundary_tol = fuzzy_boundary_tol;
  d_relative_error_tol = relative_error_tol;
  d_error_tol = error_tol;

  this->verifyValidTolerances( fuzzy_boundary_tol,
                               relative_error_tol,
                               error_tol );
}

// Return the evaluation fuzzy bound tolerance
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
double UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::getFuzzyBoundTolerance() const
{
  return d_fuzzy_boundary_tol;
}

// Return the evaluation relative error tolerance
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
double UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::getRelativeErrorTolerance() const
{
  return d_relative_error_tol;
}

// Return the evaluation error tolerance
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
double UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::getErrorTolerance() const
{
  return d_error_tol;
}

// Verify that the distribution data is valid
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
void UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::verifyValidData(
     const std::vector<PrimaryIndepQuantity>& primary_indep_grid,
     const std::vector<std::shared_ptr<const BaseUnivariateDistributionType> >&
     secondary_distributions )
{
  // Check that there is at least one bin specified
  TEST_FOR_EXCEPTION( primary_indep_grid.size() <= 1,
                      Utility::BadBivariateDistributionParameter,
                      "The interpolated tabular basic bivariate distribution "
                      "cannot be created because at least one primary bin "
                      "needs to be specified!" );

  // Check that every univariate distribution is continuous
  typename std::vector<std::shared_ptr<const BaseUnivariateDistributionType> >::const_iterator bad_secondary_dist =
    std::find_if( secondary_distributions.begin(),
                  secondary_distributions.end(),
                  [](const std::shared_ptr<const BaseUnivariateDistributionType>& dist){ return !dist->isContinuous(); } );

  TEST_FOR_EXCEPTION( bad_secondary_dist != secondary_distributions.end(),
                      Utility::BadBivariateDistributionParameter,
                      "The interpolated tabular basic bivariate distribution "
                      "cannot be created because the secondary distribution "
                      "at index "
                      << std::distance( secondary_distributions.begin(), bad_secondary_dist ) <<
                      " is not continuous!" );

  // Check that every univariate distribution is compatible with the
  // interpolation type
  bad_secondary_dist =
    std::find_if( secondary_distributions.begin(),
                  secondary_distributions.end(),
                  [](const std::shared_ptr<const BaseUnivariateDistributionType>& dist){ return !dist->template isCompatibleWithInterpType<typename TwoDInterpPolicy::SecondaryBasePolicy>(); } );

  TEST_FOR_EXCEPTION( bad_secondary_dist != secondary_distributions.end(),
                      Utility::BadBivariateDistributionParameter,
                      "The interpolated tabular basic bivariate distribution "
                      "cannot be created because the secondary distribution "
                      "at index "
                      << std::distance( secondary_distributions.begin(), bad_secondary_dist ) <<
                      " is not compatible with " << TwoDInterpPolicy::name() <<
                      " interpolation!" );
}

// Verify that the tolerances are valid
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
void UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::verifyValidTolerances(
                                               const double fuzzy_boundary_tol,
                                               const double relative_error_tol,
                                               const double error_tol )
{
  TEST_FOR_EXCEPTION( fuzzy_boundary_tol < 0.0 || fuzzy_boundary_tol >= 1,
                      Utility::BadBivariateDistributionParameter,
                      "The fuzzy boundary tolerance must in the interval "
                      "[0, 1)!" );

  TEST_FOR_EXCEPTION( relative_error_tol < 0.0 || relative_error_tol >= 1,
                      Utility::BadBivariateDistributionParameter,
                      "The relative error tolerance must in the interval "
                      "[0, 1)!" );

  TEST_FOR_EXCEPTION( error_tol < 0.0 || error_tol >= 1,
                      Utility::BadBivariateDistributionParameter,
                      "The error tolerance must in the interval "
                      "[0, 1)!" );
}
 
// Evaluate the distribution
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
auto UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::evaluate(
                 const PrimaryIndepQuantity primary_indep_var_value,
                 const SecondaryIndepQuantity secondary_indep_var_value ) const
  -> DepQuantity
{
  return this->evaluateImpl<TwoDInterpPolicy,DepQuantity>(
                                    primary_indep_var_value,
                                    secondary_indep_var_value,
                                    &BaseUnivariateDistributionType::evaluate);
}

// Evaluate the secondary conditional PDF
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
auto UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::evaluateSecondaryConditionalPDF(
                 const PrimaryIndepQuantity primary_indep_var_value,
                 const SecondaryIndepQuantity secondary_indep_var_value ) const
  -> InverseSecondaryIndepQuantity
{
  return this->evaluateImpl<TwoDInterpPolicy,InverseSecondaryIndepQuantity>(
                                primary_indep_var_value,
                                secondary_indep_var_value,
                                &BaseUnivariateDistributionType::evaluatePDF );
}

// Evaluate the distribution using the desired evaluation method
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
template<typename LocalTwoDInterpPolicy,
         typename ReturnType,
         typename EvaluationMethod>
inline ReturnType UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::evaluateImpl(
                        const PrimaryIndepQuantity primary_indep_var_value,
                        const SecondaryIndepQuantity secondary_indep_var_value,
                        EvaluationMethod evaluate ) const
{
  // Create the lower bound functor
  std::function<SecondaryIndepQuantity(const PrimaryIndepQuantity)>
    min_secondary_indep_var_functor = std::bind<SecondaryIndepQuantity>(
                                &ThisType::getLowerBoundOfSecondaryConditionalIndepVar,
                                std::cref( *this ),
                                std::placeholders::_1 );

  // Create the upper bound functor
  std::function<SecondaryIndepQuantity(const PrimaryIndepQuantity)>
    max_secondary_indep_var_functor = std::bind<SecondaryIndepQuantity>(
                                &ThisType::getUpperBoundOfSecondaryConditionalIndepVar,
                                std::cref( *this ),
                                std::placeholders::_1 );

  return this->evaluateImpl<LocalTwoDInterpPolicy,ReturnType,EvaluationMethod>(
              primary_indep_var_value,
              secondary_indep_var_value,
              min_secondary_indep_var_functor,
              max_secondary_indep_var_functor,
              evaluate );
}

// Evaluate the distribution using the desired evaluation method
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
template<typename LocalTwoDInterpPolicy,
         typename ReturnType,
         typename EvaluationMethod>
inline auto UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::evaluateImpl(
             const PrimaryIndepQuantity primary_indep_var_value,
             const SecondaryIndepQuantity secondary_indep_var_value,
             const std::function<SecondaryIndepQuantity(PrimaryIndepQuantity)>&
             min_secondary_indep_var_functor,
             const std::function<SecondaryIndepQuantity(PrimaryIndepQuantity)>&
             max_secondary_indep_var_functor,
             EvaluationMethod evaluate,
             unsigned max_number_of_iterations ) const
  -> ReturnType
{
  // Find the bin boundaries
  DistributionDataConstIterator lower_bin_boundary, upper_bin_boundary;
  
  this->findBinBoundaries( primary_indep_var_value,
                           lower_bin_boundary,
                           upper_bin_boundary );
  
  if( lower_bin_boundary != upper_bin_boundary )
  {
    return TwoDSamplePolicy::template evaluatePDF<LocalTwoDInterpPolicy,BaseUnivariateDistributionType,PrimaryIndepQuantity,SecondaryIndepQuantity,ReturnType>(
        primary_indep_var_value,
        secondary_indep_var_value,
        min_secondary_indep_var_functor,
        max_secondary_indep_var_functor,
        evaluate,
        lower_bin_boundary,
        upper_bin_boundary,
        d_fuzzy_boundary_tol,
        d_relative_error_tol,
        d_error_tol,
        max_number_of_iterations );
  }
  // Primary value is outside of the primary grid limits
  else
  {
    if( this->arePrimaryLimitsExtended() )
    {
      return ((*lower_bin_boundary->second).*evaluate)(secondary_indep_var_value);
    }
    else 
      return QuantityTraits<ReturnType>::zero();
  }
}

// Return a random sample from the secondary conditional PDF
/*! \details If the primary value provided is outside of the primary grid 
 * limits the appropriate limiting secondary distribution will be used to 
 * create the sample. The alternative to this behavior is to throw an exception
 * unless the distribution has been extended by calling the 
 * extendBeyondPrimaryIndepLimits method. Since this is a performance critical 
 * method we decided against this behavior.
 */
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
auto UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::sampleSecondaryConditional(
                     const PrimaryIndepQuantity primary_indep_var_value ) const
  -> SecondaryIndepQuantity
{
  // Create the sampling functor
  std::function<SecondaryIndepQuantity(const BaseUnivariateDistributionType&)>
    sampling_functor = std::bind<SecondaryIndepQuantity>(
                                       &BaseUnivariateDistributionType::sample,
                                       std::placeholders::_1 );

  return this->sampleImpl( primary_indep_var_value, sampling_functor );
}

// Return a random sample and record the number of trials
/*! \details If the primary value
 * provided is outside of the primary grid limits the appropriate limiting
 * secondary distribution will be used to create the sample. The alternative
 * to this behavior is to throw an exception unless the distribution has 
 * been extended by calling the extendBeyondPrimaryIndepLimits method. Since
 * this is a performance critical method we decided against this behavior.
 */
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
auto UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::sampleSecondaryConditionalAndRecordTrials(
                            const PrimaryIndepQuantity primary_indep_var_value,
                            DistributionTraits::Counter& trials ) const
  -> SecondaryIndepQuantity
{
  // Create the sampling functor
  std::function<SecondaryIndepQuantity(const BaseUnivariateDistributionType&)>
    sampling_functor = std::bind<SecondaryIndepQuantity>(
                        &BaseUnivariateDistributionType::sampleAndRecordTrials,
                        std::placeholders::_1,
                        std::ref( trials ) );

  return this->sampleImpl( primary_indep_var_value, sampling_functor );
}

// Sample from the distribution using the desired sampling functor
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
template<typename SampleFunctor>
inline auto UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::sampleDetailedImpl(
                            const PrimaryIndepQuantity primary_indep_var_value,
                            SampleFunctor sample_functor,
                            SecondaryIndepQuantity& raw_sample,
                            size_t& primary_bin_index ) const
  -> SecondaryIndepQuantity
{
  // Create the lower bound functor
  std::function<SecondaryIndepQuantity(const PrimaryIndepQuantity)>
    lower_bound_functor = std::bind<SecondaryIndepQuantity>(
                              &ThisType::getLowerBoundOfSecondaryConditionalIndepVar,
                              std::cref( *this ),
                              std::placeholders::_1 );

  // Create the upper bound functor
  std::function<SecondaryIndepQuantity(const PrimaryIndepQuantity)>
    upper_bound_functor = std::bind<SecondaryIndepQuantity>(
                              &ThisType::getUpperBoundOfSecondaryConditionalIndepVar,
                              std::cref( *this ),
                              std::placeholders::_1 );

  return this->sampleDetailedImpl( primary_indep_var_value,
                                   sample_functor,
                                   raw_sample,
                                   primary_bin_index,
                                   lower_bound_functor,
                                   upper_bound_functor );
}

// Sample from the distribution using the desired sampling functor
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
template<typename SampleFunctor>
inline auto UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::sampleDetailedImpl(
            const PrimaryIndepQuantity primary_indep_var_value,
            SampleFunctor sample_functor,
            SecondaryIndepQuantity& raw_sample,
            size_t& primary_bin_index,
            const std::function<SecondaryIndepQuantity(PrimaryIndepQuantity)>&
            min_secondary_indep_var_functor,
            const std::function<SecondaryIndepQuantity(PrimaryIndepQuantity)>&
            max_secondary_indep_var_functor ) const
  -> SecondaryIndepQuantity
{
  // Find the bin boundaries
  DistributionDataConstIterator lower_bin_boundary, upper_bin_boundary;
  
  this->findBinBoundaries( primary_indep_var_value,
                           lower_bin_boundary,
                           upper_bin_boundary );

  SecondaryIndepQuantity sample;
  
  if( lower_bin_boundary != upper_bin_boundary )
  {
    DistributionDataConstIterator sampled_bin_boundary;
    
    sample = TwoDSamplePolicy::template sampleDetailed<TwoDInterpPolicy,PrimaryIndepQuantity,SecondaryIndepQuantity>(
          sample_functor,
          min_secondary_indep_var_functor,
          max_secondary_indep_var_functor,
          primary_indep_var_value,
          lower_bin_boundary,
          upper_bin_boundary,
          sampled_bin_boundary,
          raw_sample );

    // Calculate the index of the primary bin boundary that will be used to
    // create the secondary conditional sample
    primary_bin_index = this->calculateBinIndex( sampled_bin_boundary );
  }
  else
  {
    if( this->arePrimaryLimitsExtended() )
    {
      raw_sample = sample_functor( *lower_bin_boundary->second );
      sample = raw_sample;

      // Calculate the index of the primary bin boundary that will be used to
      // create the secondary conditional sample
      primary_bin_index = this->calculateBinIndex( lower_bin_boundary );
    }
    else
    {
      THROW_EXCEPTION( std::logic_error,
                       "Sampling beyond the primary grid boundaries "
                       "cannot be done unless the grid has been extended ("
                       << primary_indep_var_value << " not in ["
                       << this->getLowerBoundOfPrimaryIndepVar() << ","
                       << this->getUpperBoundOfPrimaryIndepVar() << "])!" );
    }
  }
  return sample;
}

// Sample from the distribution using the desired sampling functor
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
template<typename SampleFunctor>
inline auto UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::sampleImpl(
                            const PrimaryIndepQuantity primary_indep_var_value,
                            SampleFunctor sample_functor ) const
  -> SecondaryIndepQuantity
{
  // Dummy variables
  SecondaryIndepQuantity dummy_raw_sample;
  size_t dummy_primary_bin_index;

  return this->sampleDetailedImpl( primary_indep_var_value,
                                   sample_functor,
                                   dummy_raw_sample,
                                   dummy_primary_bin_index );
}

// Sample from the distribution using the desired sampling functor
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
template<typename SampleFunctor>
inline auto UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::sampleImpl(
            const PrimaryIndepQuantity primary_indep_var_value,
            SampleFunctor sample_functor,
            const std::function<SecondaryIndepQuantity(PrimaryIndepQuantity)>&
              min_secondary_indep_var_functor,
            const std::function<SecondaryIndepQuantity(PrimaryIndepQuantity)>&
              max_secondary_indep_var_functor ) const
  -> SecondaryIndepQuantity
{
  SecondaryIndepQuantity dummy_raw_sample;
  size_t dummy_primary_bin_index;

  return this->sampleDetailedImpl( primary_indep_var_value,
                                   sample_functor,
                                   dummy_raw_sample,
                                   dummy_primary_bin_index,
                                   min_secondary_indep_var_functor,
                                   max_secondary_indep_var_functor );
}

// Return the upper bound of the conditional distribution
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
auto UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::getUpperBoundOfSecondaryConditionalIndepVar(
                     const PrimaryIndepQuantity primary_indep_var_value ) const
  -> SecondaryIndepQuantity
{
  DistributionDataConstIterator lower_bin_boundary, upper_bin_boundary;
  
  this->findBinBoundaries( primary_indep_var_value,
                           lower_bin_boundary,
                           upper_bin_boundary );

  // Check for a primary value outside of the primary grid limits
  if( lower_bin_boundary == upper_bin_boundary )
  {
    if( this->arePrimaryLimitsExtended() )
      return Utility::get<1>( *lower_bin_boundary )->getUpperBoundOfIndepVar();
    else
      return QuantityTraits<SecondaryIndepQuantity>::zero();
  }
  else
  {
    return TwoDInterpPolicy::calculateIntermediateGridLimit(
           Utility::get<0>( *lower_bin_boundary ),
           Utility::get<0>( *upper_bin_boundary ),
           primary_indep_var_value,
           Utility::get<1>( *lower_bin_boundary )->getUpperBoundOfIndepVar(),
           Utility::get<1>( *upper_bin_boundary )->getUpperBoundOfIndepVar() );
  }
}

// Return the lower bound of the conditional distribution
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
auto UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::getLowerBoundOfSecondaryConditionalIndepVar(
                     const PrimaryIndepQuantity primary_indep_var_value ) const
  -> SecondaryIndepQuantity
{
  DistributionDataConstIterator lower_bin_boundary, upper_bin_boundary;
  
  this->findBinBoundaries( primary_indep_var_value,
                           lower_bin_boundary,
                           upper_bin_boundary );

  // Check for a primary value outside of the primary grid limits
  if( lower_bin_boundary == upper_bin_boundary )
  {
    if( this->arePrimaryLimitsExtended() )
      return Utility::get<1>( *lower_bin_boundary )->getLowerBoundOfIndepVar();
    else
      return QuantityTraits<SecondaryIndepQuantity>::zero();
  }
  else
  {
    return TwoDInterpPolicy::calculateIntermediateGridLimit(
           Utility::get<0>( *lower_bin_boundary ),
           Utility::get<0>( *upper_bin_boundary ),
           primary_indep_var_value,
           Utility::get<1>( *lower_bin_boundary )->getLowerBoundOfIndepVar(),
           Utility::get<1>( *upper_bin_boundary )->getLowerBoundOfIndepVar() );
  }
}

// Test if the distribution is continuous in the primary dimension
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
bool UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::isPrimaryDimensionContinuous() const
{
  return true;
}

// Save the distribution to an archive
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
template<typename Archive>
void UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::save( Archive& ar, const unsigned version ) const
{
  // Save the base class first
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( BaseType );

  // Save the local member data
  ar & BOOST_SERIALIZATION_NVP( d_fuzzy_boundary_tol );
  ar & BOOST_SERIALIZATION_NVP( d_relative_error_tol );
  ar & BOOST_SERIALIZATION_NVP( d_error_tol );
}

// Load the distribution from an archive
template<typename TwoDInterpPolicy, typename TwoDSamplePolicy, typename Distribution>
template<typename Archive>
void UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase<TwoDInterpPolicy,TwoDSamplePolicy,Distribution>::load( Archive& ar, const unsigned version )
{
  // Load the base class first
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( BaseType );

  // Load the local member data
  ar & BOOST_SERIALIZATION_NVP( d_fuzzy_boundary_tol );
  ar & BOOST_SERIALIZATION_NVP( d_relative_error_tol );
  ar & BOOST_SERIALIZATION_NVP( d_error_tol );
}

} // end Utility namespace

#define __INTERPOLATED_TABULAR_BASIC_BIVARIATE_DIST_WITH_SAMPLE_POLICY__( DECL_TYPE, BivariateDistribution, TwoDSamplePolicy, ... ) \
  DECL_TYPE##_DISTRIBUTION_INST( BivariateDistribution<Utility::LinLinLin,TwoDSamplePolicy,__VA_ARGS__ > ); \
  DECL_TYPE##_DISTRIBUTION_INST( BivariateDistribution<Utility::LinLogLin,TwoDSamplePolicy,__VA_ARGS__ > ); \
  DECL_TYPE##_DISTRIBUTION_INST( BivariateDistribution<Utility::LinLogCosLin,TwoDSamplePolicy,__VA_ARGS__ > ); \
  DECL_TYPE##_DISTRIBUTION_INST( BivariateDistribution<Utility::LinLinLog,TwoDSamplePolicy,__VA_ARGS__ > ); \
  DECL_TYPE##_DISTRIBUTION_INST( BivariateDistribution<Utility::LinLogLog,TwoDSamplePolicy,__VA_ARGS__ > ); \
  DECL_TYPE##_DISTRIBUTION_INST( BivariateDistribution<Utility::LinLogCosLog,TwoDSamplePolicy,__VA_ARGS__ > ); \
  DECL_TYPE##_DISTRIBUTION_INST( BivariateDistribution<Utility::LogLinLin,TwoDSamplePolicy,__VA_ARGS__ > ); \
  DECL_TYPE##_DISTRIBUTION_INST( BivariateDistribution<Utility::LogLogLin,TwoDSamplePolicy,__VA_ARGS__ > ); \
  DECL_TYPE##_DISTRIBUTION_INST( BivariateDistribution<Utility::LogLogCosLin,TwoDSamplePolicy,__VA_ARGS__ > ); \
  DECL_TYPE##_DISTRIBUTION_INST( BivariateDistribution<Utility::LogLinLog,TwoDSamplePolicy,__VA_ARGS__ > ); \
  DECL_TYPE##_DISTRIBUTION_INST( BivariateDistribution<Utility::LogLogLog,TwoDSamplePolicy,__VA_ARGS__ > ); \
  DECL_TYPE##_DISTRIBUTION_INST( BivariateDistribution<Utility::LogLogCosLog,TwoDSamplePolicy,__VA_ARGS__ > )

#define ___INTERPOLATED_TABULAR_BASIC_BIVARIATE_DIST__( DECL_TYPE, BivariateDistribution, ... ) \
  __INTERPOLATED_TABULAR_BASIC_BIVARIATE_DIST_WITH_SAMPLE_POLICY__( DECL_TYPE, BivariateDistribution, Utility::Direct, __VA_ARGS__ ); \
  __INTERPOLATED_TABULAR_BASIC_BIVARIATE_DIST_WITH_SAMPLE_POLICY__( DECL_TYPE, BivariateDistribution, Utility::UnitBase, __VA_ARGS__ ); \
  __INTERPOLATED_TABULAR_BASIC_BIVARIATE_DIST_WITH_SAMPLE_POLICY__( DECL_TYPE, BivariateDistribution, Utility::CumulativePoints, __VA_ARGS__ ); \
  __INTERPOLATED_TABULAR_BASIC_BIVARIATE_DIST_WITH_SAMPLE_POLICY__( DECL_TYPE, BivariateDistribution, Utility::Correlated, __VA_ARGS__ ); \
  __INTERPOLATED_TABULAR_BASIC_BIVARIATE_DIST_WITH_SAMPLE_POLICY__( DECL_TYPE, BivariateDistribution, Utility::UnitBaseCorrelated, __VA_ARGS__ )

#define EXTERN_EXPLICIT_INTERPOLATED_TABULAR_BASIC_BIVARIATE_DIST( BivariateDistribution, ... ) \
  ___INTERPOLATED_TABULAR_BASIC_BIVARIATE_DIST__( EXTERN_EXPLICIT, BivariateDistribution, __VA_ARGS__ )

#define EXPLICIT_INTERPOLATED_TABULAR_BASIC_BIVARIATE_DIST( BivariateDistribution, ... ) \
  ___INTERPOLATED_TABULAR_BASIC_BIVARIATE_DIST__( EXPLICIT, BivariateDistribution, __VA_ARGS__ )

EXTERN_EXPLICIT_INTERPOLATED_TABULAR_BASIC_BIVARIATE_DIST( Utility::UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase, Utility::UnitAwareFullyTabularBasicBivariateDistribution<void,void,void> );
EXTERN_EXPLICIT_INTERPOLATED_TABULAR_BASIC_BIVARIATE_DIST( UnitAwareInterpolatedTabularBasicBivariateDistributionImplBase, Utility::UnitAwarePartiallyTabularBasicBivariateDistribution<void,void,void> );

#endif // end UTILITY_INTERPOLATED_TABULAR_BASIC_BIVARIATE_DISTRIBUTION_IMPL_BASE_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_InterpolatedTabularBasicBivariateDistributionImplBase_def.hpp
//---------------------------------------------------------------------------//
