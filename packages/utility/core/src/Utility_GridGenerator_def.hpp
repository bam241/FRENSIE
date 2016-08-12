//---------------------------------------------------------------------------//
//!
//! \file   Utility_LinearGridGenerator_def.hpp
//! \author Alex Robinson
//! \brief  Linear grid generator class template definitions
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_LINEAR_GRID_GENERATOR_DEF_HPP
#define UTILITY_LINEAR_GRID_GENERATOR_DEF_HPP

// Std Lib Includes
#include <deque>
#include <iterator>
#include <sstream>

// FRENSIE Includes
#include "Utility_ContractException.hpp"
#include "Utility_InterpolationPolicy.hpp"
#include "Utility_SortAlgorithms.hpp"
#include "Utility_ComparePolicy.hpp"
#include "Utility_ExceptionTestMacros.hpp"

namespace Utility{

// Constructor
template<typename InterpPolicy>
GridGenerator<InterpPolicy>::GridGenerator( const double convergence_tol,
					    const double absolute_diff_tol,
					    const double distance_tol )
  : d_convergence_tol( convergence_tol ),
    d_absolute_diff_tol( absolute_diff_tol ),
    d_distance_tol( distance_tol ),
    d_throw_exceptions( false ),
    d_os_warn( &std::cerr )
{
  // Make sure the convergence tolerance is valid
  testPrecondition( convergence_tol <= 1.0 );
  testPrecondition( convergence_tol > 0.0 );
  // Make sure the absolute difference tolerance is valid
  testPrecondition( absolute_diff_tol <= 1.0 );
  testPrecondition( absolute_diff_tol >= 0.0 );
  // Make sure the distance tolerance is valid
  testPrecondition( distance_tol <= 1.0 );
  testPrecondition( distance_tol >= 0.0 );
}

// Throw exception on dirty convergence
/*! \details "Dirty Convergence" has occured when the distance tolerance or
 * the absolute difference tolerance is reached before the convergence
 * tolerance. This type of convergence should be avoided because the grid
 * has not truely converged.
 */
template<typename InterpPolicy>
void GridGenerator<InterpPolicy>::throwExceptionOnDirtyConvergence()
{
  d_throw_exceptions = true;
}

// Warn on dirty convergence (default)
/*! \details "Dirty Convergence" has occured when the distance tolerance or
 * the absolute difference tolerance is reached before the convergence
 * tolerance. This type of convergence should be avoided because the grid
 * has not truely converged.
 */
template<typename InterpPolicy>
void GridGenerator<InterpPolicy>::warnOnDirtyConvergence(
                                                        std::ostream* os_warn )
{
  d_throw_exceptions = false;

  d_os_warn = os_warn;
}

// Check if an exception will be thrown on dirty convergence
template<typename InterpPolicy>
bool GridGenerator<InterpPolicy>::isExceptionThrownOnDirtyConvergence() const
{
  return d_throw_exceptions;
}

// Set the convergence tolerance
template<typename InterpPolicy>
void GridGenerator<InterpPolicy>::setConvergenceTolerance(
						 const double convergence_tol )
{
  // Make sure the convergence tolerance is valid
  testPrecondition( convergence_tol <= 1.0 );
  testPrecondition( convergence_tol > 0.0 );

  d_convergence_tol = convergence_tol;
}

// Get the convergence tolerance
template<typename InterpPolicy>
double GridGenerator<InterpPolicy>::getConvergenceTolerance() const
{
  return d_convergence_tol;
}

// Set the absolute difference tolerance
template<typename InterpPolicy>
void GridGenerator<InterpPolicy>::setAbsoluteDifferenceTolerance(
					       const double absolute_diff_tol )
{
  // Make sure the absolute difference tolerance is valid
  testPrecondition( absolute_diff_tol <= 1.0 );
  testPrecondition( absolute_diff_tol >= 0.0 );

  d_absolute_diff_tol = absolute_diff_tol;
}

// Get the absolute difference tolerance
template<typename InterpPolicy>
double GridGenerator<InterpPolicy>::getAbsoluteDifferenceTolerance() const
{
  return d_absolute_diff_tol;
}

// Set the distance tolerance
template<typename InterpPolicy>
void GridGenerator<InterpPolicy>::setDistanceTolerance(
						    const double distance_tol )
{
  // Make sure the distance tolerance is valid
  testPrecondition( distance_tol <= 1.0 );
  testPrecondition( distance_tol >= 0.0 );

  d_distance_tol = distance_tol;
}

// Get the distance tolerance
template<typename InterpPolicy>
double GridGenerator<InterpPolicy>::getDistanceTolerance() const
{
  return d_distance_tol;
}

// Generate the grid in place
/*! \details There must be at least two initial grid points given (the lower
 * grid boundary and the upper grid boundary). If there are discontinuities in
 * the function, the grid points just below and just above the discontinuity
 * should also be given to speed up the algorithm. The convergence tolerance
 * is used to determine if two consecutive grid points are acceptable - if
 * the relative error between the estimated value of the function (from
 * desired interpolation) at the midpoint between two grid points end the
 * actual value of the function at the midpoint is less that or equal to the
 * convergence tolerance, the two grid points are kept. Otherwise the midpoint
 * is inserted into the grid and the process is repeated. Do not process
 * the grid points before passing them into this function.
 */
template<typename InterpPolicy>
template<typename STLCompliantContainer, typename Functor>
void GridGenerator<InterpPolicy>::generateInPlace(
						STLCompliantContainer& grid,
						const Functor& function ) const
{
  // Make sure the container value type is a floating point type
  testStaticPrecondition( (boost::is_float<typename STLCompliantContainer::value_type>::value) );
  // Make sure at least 2 initial grid points have been given
  testPrecondition( grid.size() >= 2 );
  // Make sure the initial grid points are sorted
  testPrecondition( Sort::isSortedAscending( grid.begin(), grid.end(), true ));

  STLCompliantContainer evaluated_function;

  this->generateAndEvaluateInPlace( grid, evaluated_function, function );
}

// Generate the grid in place (return evaluated function on grid)
/*! \details There must be at least two initial grid points given (the lower
 * grid boundary and the upper grid boundary). If there are discontinuities in
 * the function, the grid points just below and just above the discontinuity
 * should also be given to speed up the algorithm. The convergence tolerance
 * is used to determine if two consecutive grid points are acceptable - if
 * the relative error between the estimated value of the function (from
 * desired interpolation) at the midpoint between two grid points end the
 * actual value of the function at the midpoint is less that or equal to the
 * convergence tolerance, the two grid points are kept. Otherwise the midpoint
 * is inserted into the grid and the process is repeated. Do not process
 * the grid points before passing them into this function.
 */
template<typename InterpPolicy>
template<typename STLCompliantContainerA,
	 typename STLCompliantContainerB,
	 typename Functor>
void GridGenerator<InterpPolicy>::generateAndEvaluateInPlace(
				    STLCompliantContainerA& grid,
				    STLCompliantContainerB& evaluated_function,
				    const Functor& function ) const
{
  // Make sure the container value type is a floating point type
  testStaticPrecondition( (boost::is_float<typename STLCompliantContainerA::value_type>::value) );
  testStaticPrecondition( (boost::is_float<typename STLCompliantContainerB::value_type>::value) );
  // Make sure at least 2 initial grid points have been given
  testPrecondition( grid.size() >= 2 );
  // Make sure the intial grid points are sorted
  testPrecondition( Sort::isSortedAscending( grid.begin(), grid.end(), true ));

  // Use a queue data structure to calculate the grid points
  std::deque<double> grid_queue( grid.begin(), grid.end() );

  // Clear the initial grid
  grid.clear();
  evaluated_function.clear();

  // Variables used to calculate the linearized grid
  double x0, x1, x_mid, y0, y1, y_mid_exact, y_mid_estimated;

  // Evaluate the first grid point
  x0 = grid_queue.front();
  grid_queue.pop_front();

  y0 = function( x0 );

  // Calculate the grid points
  while( !grid_queue.empty() )
  {
    x1 = grid_queue.front();

    x_mid = InterpPolicy::recoverProcessedIndepVar(
				     0.5*(InterpPolicy::processIndepVar(x0) +
					  InterpPolicy::processIndepVar(x1)) );

    y1 = function( x1 );
    y_mid_exact = function( x_mid );

    y_mid_estimated = InterpPolicy::interpolate( x0, x1, x_mid, y0, y1 );

    bool converged =
      this->hasGridConverged( x0, x_mid, x1, y_mid_estimated, y_mid_exact );

    // Keep the grid points
    if( converged )
    {
      grid.push_back( x0 );
      evaluated_function.push_back( y0 );

      x0 = x1;
      grid_queue.pop_front();

      y0 = y1;
    }
    // Refine the grid
    else
      grid_queue.push_front( x_mid );
  }

  // Add the last point to the linearized grid
  grid.push_back( x0 );
  evaluated_function.push_back( y0 );

  // Make sure the linearized grid has at least 2 points
  testPostcondition( grid.size() >= 2 );
  testPostcondition( grid.size() == evaluated_function.size() );
  // Make sure the linearized grid is sorted
  testPostcondition( Sort::isSortedAscending( grid.begin(), grid.end() ) );
}

// Generate the linearized grid
/*! \details There must be at least two initial grid points given (the lower
 * grid boundary and the upper grid boundary). If there are discontinuities in
 * the function, the grid points just below and just above the discontinuity
 * should also be given to speed up the algorithm. The convergence tolerance
 * is used to determine if two consecutive grid points are acceptable - if
 * the relative error between the estimated value of the function (from
 * lin-lin interpolation) at the midpoint between two grid points end the
 * actual value of the function at the midpoint is less that or equal to the
 * convergence tolerance, the two grid points are kept. Otherwise the midpoint
 * is inserted into the grid and the process is repeated. Do not process
 * the grid points before passing them into this function.
 */
template<typename InterpPolicy>
template<typename STLCompliantContainerA,
	 typename STLCompliantContainerB,
	 typename Functor>
void GridGenerator<InterpPolicy>::generate(
			     STLCompliantContainerA& grid,
			     const STLCompliantContainerB& initial_grid_points,
			     const Functor& function ) const
{
  // Make sure the container value type is a floating point type
  testStaticPrecondition( (boost::is_float<typename STLCompliantContainerA::value_type>::value) );
  testStaticPrecondition( (boost::is_float<typename STLCompliantContainerB::value_type>::value) );
  // Make sure at least 2 initial grid points have been given
  testPrecondition( initial_grid_points.size() >= 2 );
  // Make sure the intial grid points are sorted
  testPrecondition( Sort::isSortedAscending( initial_grid_points.begin(),
					     initial_grid_points.end(),
					     true ) );

  grid.assign( initial_grid_points.begin(), initial_grid_points.end() );

  STLCompliantContainerA evaluated_function;

  this->generateAndEvaluateInPlace( grid, evaluated_function, function );
}

// Generate the linearized grid
/*! \details There must be at least two initial grid points given (the lower
 * grid boundary and the upper grid boundary). If there are discontinuities in
 * the function, the grid points just below and just above the discontinuity
 * should also be given to speed up the algorithm. The convergence tolerance
 * is used to determine if two consecutive grid points are acceptable - if
 * the relative error between the estimated value of the function (from
 * lin-lin interpolation) at the midpoint between two grid points end the
 * actual value of the function at the midpoint is less that or equal to the
 * convergence tolerance, the two grid points are kept. Otherwise the midpoint
 * is inserted into the grid and the process is repeated. Do not process
 * the grid points before passing them into this function.
 */
template<typename InterpPolicy>
template<typename STLCompliantContainerA,
	 typename STLCompliantContainerB,
	 typename STLCompliantContainerC,
	 typename Functor>
void GridGenerator<InterpPolicy>::generateAndEvaluate(
			     STLCompliantContainerA& grid,
			     STLCompliantContainerB& evaluated_function,
			     const STLCompliantContainerC& initial_grid_points,
			     const Functor& function ) const
{
  // Make sure the container value type is a floating point type
  testStaticPrecondition( (boost::is_float<typename STLCompliantContainerA::value_type>::value) );
  testStaticPrecondition( (boost::is_float<typename STLCompliantContainerB::value_type>::value) );
  testStaticPrecondition( (boost::is_float<typename STLCompliantContainerC::value_type>::value) );
  // Make sure at least 2 initial grid points have been given
  testPrecondition( initial_grid_points.size() >= 2 );
  // Make sure the intial grid points are sorted
  testPrecondition( Sort::isSortedAscending( initial_grid_points.begin(),
					     initial_grid_points.end() ) );

  grid.assign( initial_grid_points.begin(), initial_grid_points.end() );

  this->generateAndEvaluateInPlace( grid, evaluated_function, function );
}

// Check for convergence
template<typename InterpPolicy>
bool GridGenerator<InterpPolicy>::hasGridConverged(
                                               const double lower_grid_point,
                                               const double mid_grid_point,
                                               const double upper_grid_point,
                                               const double y_mid_estimated,
                                               const double y_mid_exact ) const
                                            
{
  bool converged = false;

  // Calculate the convergence parameters
  double relative_error = Policy::relError( y_mid_exact, y_mid_estimated );

  double absolute_difference =
      Teuchos::ScalarTraits<double>::magnitude( y_mid_exact - y_mid_estimated);

  double relative_distance =
    Policy::relError( lower_grid_point, upper_grid_point );

  // Check if the distance tolerance was hit - dirty convergence
  if( relative_distance <= d_distance_tol &&
      relative_error > d_convergence_tol )
  {
    std::ostringstream oss;
    oss.precision( 18 );
    oss << "distance tolerance hit before convergence - "
        << "relError(x0,x1) = relError(" << lower_grid_point << ","
        << upper_grid_point << ") = " << relative_distance
        << ", relError(ym,ym_exact) = relError(" << y_mid_estimated
        << "," << y_mid_exact << ") = " << relative_error;

    if( d_throw_exceptions )
    {
      THROW_EXCEPTION( std::runtime_error, "Error: " << oss.str() );
    }
    else
    {
      converged = true;
      
      d_os_warn->precision( 18 );
      (*d_os_warn) << "Warning: " << oss.str() << std::endl;
    }
  }

  // Check if the absolute difference tolerance was hit - dirty convergence
  if( absolute_difference <= d_absolute_diff_tol &&
      relative_error > d_convergence_tol )
  {
    std::ostringstream oss;
    oss.precision( 18 );
    oss << "absolute difference tolerance hit before "
        << "convergence - x_mid=" << mid_grid_point << ", y_mid_exact="
        << y_mid_exact << ", y_mid_estimated="
        << y_mid_estimated << ", abs_diff=" << absolute_difference;

    if( d_throw_exceptions )
    {
      THROW_EXCEPTION( std::runtime_error,
                       "Error: " << oss.str() );
    }
    else
    {
      converged = true;
      
      d_os_warn->precision( 18 );
      (*d_os_warn) << "Warning: " << oss.str() << std::endl;
    }
  }

  // Check if the convergence tolerance was hit - clean convergence
  if( relative_error <= d_convergence_tol )
    converged = true;

  return converged;
}

} // end Utility namespace

#endif // end UTILITY_LINEAR_GRID_GENERATOR_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_LinearGridGenerator_def.hpp
//---------------------------------------------------------------------------//
