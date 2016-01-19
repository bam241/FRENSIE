//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_CompleteDopplerBroadenedPhotonEnergyDistribution_def.hpp
//! \author Alex Robinson
//! \brief  The complete Doppler broadened photon energy distribution def.
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_COMPLETE_DOPPLER_BROADENED_PHOTON_ENERGY_DISTRIBUTION_DEF_HPP
#define MONTE_CARLO_COMPLETE_DOPPLER_BROADENED_PHOTON_ENERGY_DISTRIBUTION_DEF_HPP

// FRENSIE Includes
#include "MonteCarlo_CompleteDopplerBroadenedPhotonEnergyDistribution.hpp"
#include "Utility_DiscreteDistribution.hpp"
#include "Utility_ContractException.hpp"

namespace MonteCarlo{

// Constructor
template<typename ComptonProfilePolicy>
CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::CompleteDopplerBroadenedPhotonEnergyDistribution(
		const Teuchos::Array<double>& endf_subshell_occupancies,
                const Teuchos::Array<SubshellType>& endf_subshell_order,
                const std::shared_ptr<const ComptonProfileSubshellConverter>&
                subshell_converter,
                const ElectronMomentumDistArray& electron_momentum_dist_array )
  : d_endf_subshell_occupancy_distribution(),
    d_endf_subshell_order(),
    d_endf_subshell_occupancies( endf_subshell_occupancies ),
    d_subshell_converter( subshell_converter ),
    d_electron_momentum_dist_array( electron_momentum_dist_array )
{
  // Make sure the shell interaction data is valid
  testPrecondition( endf_subshell_occupancies.size() > 0 );
  testPrecondition( endf_subshell_order.size() == 
		    endf_subshell_occupancies.size() );
  // Make sure that the subshell converter is valid
  testPrecondition( subshell_converter.get() );
  // Make sure the comptron profile array is valid
  testPrecondition( electron_momentum_dist_array.size() > 0 );
  testPrecondition( ComptonProfilePolicy::isValidProfile( electron_momentum_dist_array.front() ) );
  testPrecondition( ComptonProfilePolicy::isValidProfile( electron_momentum_dist_array.back() ) );
  
  // Create the ENDF subshell interaction distribution
  Teuchos::Array<double> dummy_indep_vals( endf_subshell_occupancies.size() );

  d_endf_subshell_occupancy_distribution.reset(
	      new Utility::DiscreteDistribution( dummy_indep_vals,
						 endf_subshell_occupancies ) );

  // Create the endf subshell order bimap
  for( unsigned i = 0; i < endf_subshell_order.size(); ++i )
  {
    d_endf_subshell_order.insert( SubshellOrderMapType::value_type( 
                                                 i, endf_subshell_order[i] ) );
  }

  // Check if a half (standard) or full profile is being used.
  if( electron_momentum_dist_array.front()->getLowerBoundOfMomentum() < 
      0.0*Utility::Units::mec_momentum )
    d_half_profiles = false;
  else
    d_half_profiles = true;
}

// Evaluate the distribution
template<typename ComptonProfilePolicy>
double CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::evaluate( 
				   const double incoming_energy,
				   const double outgoing_energy,
				   const double scattering_angle_cosine ) const
{
  // Make sure the incoming energy is valid
  testPrecondition( incoming_energy > 0.0 );
  // Make sure the outgoing energy is valid
  testPrecondition( outgoing_energy < incoming_energy );
  // Make sure the scattering angle is valid
  testPrecondition( scattering_angle_cosine >= -1.0 );
  testPrecondition( scattering_angle_cosine <= 1.0 );

  // The total double differential cross section
  double cross_section = 0.0;

  // Evaluate each subshell
  SubshellOrderMapType::const_iterator subshell_it = 
    d_endf_subshell_order.begin();

  while( subshell_it != d_endf_subshell_order.end() )
  {
    cross_section += this->evaluateSubshell( incoming_energy,
                                             outgoing_energy,
                                             scattering_angle_cosine,
                                             subshell_it->right );

    ++subshell_it;
  }

  // Make sure the cross section is valid
  testPostcondition( cross_section >= 0.0 );

  return cross_section;
}
    
// Evaluate the subshell distribution
template<typename ComptonProfilePolicy>
double CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::evaluateSubshell( 
				          const double incoming_energy,
					  const double outgoing_energy,
				          const double scattering_angle_cosine,
					  const SubshellType subshell ) const
{
  // Make sure the incoming energy is valid
  testPrecondition( incoming_energy > 0.0 );
  // Make sure the outgoing energy is valid
  testPrecondition( outgoing_energy < incoming_energy );
  // Make sure the scattering angle is valid
  testPrecondition( scattering_angle_cosine >= -1.0 );
  testPrecondition( scattering_angle_cosine <= 1.0 );

  // Get the subshell binding energy
  const double subshell_binding_energy = 
    this->getSubshellBindingEnergy( subshell );

  // The evaluated double differential cross section
  double cross_section;

  if( outgoing_energy < incoming_energy - subshell_binding_energy )
  {
    // Get the Compton profile for the subshell
    const ComptonProfile& compton_profile = 
      this->getComptonProfile( subshell );
      
    // Get the subshell occupancy
    const double subshell_occupancy = this->getSubshellOccupancy( subshell );
  
    // Calculate the electron momentum projection
    const ComptonProfile::MomentumQuantity electron_momentum_projection = 
      Utility::Units::mec_momentum*
      calculateElectronMomentumProjection( incoming_energy,
                                           outgoing_energy,
                                           scattering_angle_cosine);

    // Evaluate the Compton profile
    ComptonProfile::ProfileQuantity compton_profile_value = 
      ComptonProfilePolicy::evaluate( compton_profile, 
                                      electron_momentum_projection );

    // Evaluate the cross section
    const double multiplier = this->evaluateMultiplier(
                                                     incoming_energy,
                                                     scattering_angle_cosine );

    cross_section = multiplier*subshell_occupancy*compton_profile_value;
  }
  else
    cross_section = 0.0;

  // Make sure the cross section is valid
  testPostcondition( cross_section >= 0.0 );

  return cross_section;
}

// Evaluate the PDF
template<typename ComptonProfilePolicy>
double CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::evaluatePDF( 
				   const double incoming_energy,
				   const double outgoing_energy,
				   const double scattering_angle_cosine ) const
{
  return this->evaluate( incoming_energy, 
                         outgoing_energy,
                         scattering_angle_cosine )/
    this->evaluateIntegratedCrossSection( incoming_energy,
                                          scattering_angle_cosine,
                                          1e-3 );
}

// Evaluate the PDF
template<typename ComptonProfilePolicy>
double CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::evaluateSubshellPDF( 
					  const double incoming_energy,
					  const double outgoing_energy,
				          const double scattering_angle_cosine,
					  const SubshellType subshell ) const
{
  return this->evaluateSubshell( incoming_energy,
                                 outgoing_energy,
                                 scattering_angle_cosine,
                                 subshell )/
    this->evaluateSubshellIntegratedCrossSection( incoming_energy,
                                                  scattering_angle_cosine,
                                                  subshell,
                                                  1e-3 );
}

// Evaluate the integrated cross section (b/mu)
template<typename ComptonProfilePolicy>
double CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::evaluateIntegratedCrossSection( 
					  const double incoming_energy,
					  const double scattering_angle_cosine,
					  const double precision ) const
{
  // Make sure the incoming energy is valid
  testPrecondition( incoming_energy > 0.0 );
  // Make sure the scattering angle is valid
  testPrecondition( scattering_angle_cosine >= -1.0 );
  testPrecondition( scattering_angle_cosine <= 1.0 );

  double cross_section = 0.0;

  // Evaluate the integrated cross section for each subshell
  for( unsigned i = 0; i < d_endf_subshell_order.size(); ++i )
  {
    cross_section += this->evaluateSubshellIntegratedCrossSection(
                                                       incoming_energy,
                                                       scattering_angle_cosine,
                                                       precision );
  }
  
  // Make sure the integrated cross section is valid
  testPrecondition( cross_section >= 0.0 );
  
  return cross_section;
}

// Evaluate the integrated cross section (b/mu)
template<typename ComptonProfilePolicy>
double CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::evaluateSubshellIntegratedCrossSection( 
				          const double incoming_energy,
					  const double scattering_angle_cosine,
					  const SubshellType subshell,
					  const double precision ) const
{
  // Make sure the incoming energy is valid
  testPrecondition( incoming_energy > 0.0 );
  // Make sure the scattering angle is valid
  testPrecondition( scattering_angle_cosine >= -1.0 );
  testPrecondition( scattering_angle_cosine <= 1.0 );
  
  boost::function<double (double x)> double_diff_cs_wrapper = 
    boost::bind<double>( &CompleteDopplerBroadenedPhotonEnergyDistribution::evaluateSubshell,
                         boost::cref( *this ),
                         incoming_energy,
                         _1,
                         scattering_angle_cosine,
                         subshell );

  double abs_error, diff_cs;

  const double binding_energy = this->getSubshellBindingEnergy( subshell );

  Utility::GaussKronrodQuadratureSet quadrature_set( presicion );
  
  quadrature_set.integrateAdaptively<15>( double_diff_cs_wrapper,
                                          0.0,
                                          incoming_energy - binding_energy,
                                          diff_cs,
                                          abs_error );

  // Make sure that the differential cross section is valid
  testPostcondition( diff_cs > 0.0 );

  return diff_cs;
}

// Sample an outgoing energy from the distribution
template<typename ComptonProfilePolicy>
void CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::sample( 
				     const double incoming_energy,
				     const double scattering_angle_cosine,
				     double& outgoing_energy,
				     SubshellType& shell_of_interaction ) const
{
  unsigned trial_dummy;

  this->sampleAndRecordTrials( incoming_energy,
			       scattering_angle_cosine,
			       outgoing_energy,
			       shell_of_interaction,
			       trial_dummy );
}

// Sample an outgoing energy and record the number of trials
/*! \details The sampling of the Compton profile and the interaction subshell
 * are decoupled in this procedure.
 */
template<typename ComptonProfilePolicy>
void CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::sampleAndRecordTrials( 
				     const double incoming_energy,
				     const double scattering_angle_cosine,
				     double& outgoing_energy,
				     SubshellType& shell_of_interaction,
				     unsigned& trials ) const
{
  // Make sure the incoming energy is valid
  testPrecondition( incoming_energy > 0.0 );
  // Make sure the scattering angle cosine is valid
  testPrecondition( scattering_angle_cosine >= -1.0 );
  testPrecondition( scattering_angle_cosine <= 1.0 );

  // Record if a valid Doppler broadening energy is calculated
  bool valid_doppler_broadening = false;

  while( true )
  {
    // Sample the shell that is interacted with
    unsigned compton_subshell_index;
    double subshell_binding_energy;
  
    this->sampleInteractionSubshell( compton_subshell_index,
                                     subshell_binding_energy,
                                     shell_of_interaction );

    // Get the Compton profile for the sampled subshell
    const ComptonProfile& compton_profile = 
      *d_electron_momentum_distribution[compton_subshell_index];

    // Calculate the maximum outgoing photon energy
    double energy_max = incoming_energy - subshell_binding_energy;

    // Compton scattering can only occur if there is enough energy to release
    // the electron from its shell
    if( energy_max <= 0.0 )
      break;

    // Calculate the maximum electron momentum projection
    ComptonProfile::MomentumQuantity pz_max = Utility::Units::mec_momentum*
      calculateMaxElectronMomentumProjection( incoming_energy,
					      subshell_binding_energy,
					      scattering_angle_cosine );

    // Sample an electron momentum projection
    ComptonProfile::MomentumQuantity pz = 
      ComptonProfilePolicy::sample( compton_profile, pz_max );
    
    bool energetically_possible;

    outgoing_energy = calculateDopplerBroadenedEnergy(pz.value(),
						      incoming_energy,
						      scattering_angle_cosine,
						      energetically_possible );
    
    if( !energetically_possible || outgoing_energy < 0.0 )
      break;
    else
    {
      if( outgoing_energy == 0.0 )
	outgoing_energy = std::numeric_limits<double>::min();
      
      valid_doppler_broadening = true;

      break;
    }
  }

  // Increment the number of trials
  ++trials;

  // Reset the outgoing energy to the Compton line energy if the
  // Doppler broadening was unsuccessful
  if( !valid_doppler_broadening ) 
  {   
    outgoing_energy = calculateComptonLineEnergy( incoming_energy,
						  scattering_angle_cosine );
  }

  // Make sure the outgoing energy is valid
  testPostcondition( outgoing_energy <= incoming_energy );
  testPostcondition( outgoing_energy > 0.0 );
  // Make sure that the sampled subshell is valid
  testPostcondition( shell_of_interaction != UNKNOWN_SUBSHELL );
  testPostcondition( shell_of_interaction != INVALID_SUBSHELL );
}

// Check if the subshell is valid
template<typename ComptonProfilePolicy>
bool CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::isValidSubshell( 
                                            const SubshellType subshell ) const
{
  return d_endf_subshell_order.right.find( subshell ) != 
    d_endf_subshell_order.right.end();
}

// Return the occupancy of a subshell (default is the ENDF occupacy)
template<typename ComptonProfilePolicy>
double CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::getSubshellOccupancy( const SubshellType subshell ) const
{
  // Make sure the subshell is valid
  testPrecondition( this->isValidSubshell( subshell ) );

  unsigned endf_subshell_index = this->getENDFSubshellIndex( subshell );

  return d_endf_subshell_occupancies[old_subshell_index];
}

// Return the old subshell index corresponding to the subshell
template<typename ComptonProfilePolicy>
unsigned CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::getOldSubshellIndex(
                                            const SubshellType subshell ) const
{
  return d_subshell_converter->convertSubshellToIndex( subshell );
}

// Return the endf subshell index corresponding to the subshell
template<typename ComptonProfilePolicy>
unsigned CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::getENDFSubshellIndex(
                                            const SubshellType subshell ) const
{
  // Make sure the subshell is valid
  testPrecondition( this->isValidSubshell( subshell ) );
  
  return d_endf_subshell_order.right.find( subshell )->second;
}

// Return the subshell corresponding to the endf subshell index
template<typename ComptonProfilePolicy>
Subshell CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::getSubshell( 
                                     const unsigned endf_subshell_index ) const
{
  SubshellOrderMapType::left_map::const_iterator endf_subshell_index_it =
    d_endf_subshell_order.left.find( endf_subshell_index );

  // Make sure the index was found
  testPostcondition( endf_subshell_index_it != 
                     d_endf_subshell_order.left.end() );
}

// Return the Compton profile for a subshell
template<typename ComptonProfilePolicy>
const ComptonProfile& CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::getComptonProfile( 
                                           const SubshellType& subshell ) const
{
  // Make sure the subshell is valid
  testPrecondition( this->isValidSubshell( subshell ) );

  // Get the old subshell corresponding to the subshell type
  unsigned old_subshell_index = this->getOldSubshellIndex( subshell );

  return d_electron_momentum_distribution[old_subshell_index];
}
  
// Return the Compton profile for an old subshell index 
template<typename ComptonProfilePolicy>
const ComptonProfile& CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::getComptonProfile( 
                               const subshell_index& old_subshell_index ) const
{
  // Make sure the old subshell index is valid
  testPrecondition( old_subshell_index < 
                    d_electron_momentum_distribution.size() );

  return d_electron_momentum_distribution[old_subshell_index];
}

// Sample an ENDF subshell
template<typename ComptonProfilePolicy>
void CompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::sampleENDFInteractionSubshell( 
					    SubshellType& shell_of_interaction,
					    unsigned& shell_index ) const
{
  d_endf_subshell_occupancy_distribution->sampleAndRecordBinIndex(shell_index);
								  
  shell_of_interaction = d_endf_subshell_order[shell_index];
}

} // end MonteCarlo namespace

#endif // end MONTE_CARLO_COMPLETE_DOPPLER_BROADENED_PHOTON_ENERGY_DISTRIBUTION_DEF_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_CompleteDopplerBroadenedPhotonEnergyDistribution_def.hpp
//---------------------------------------------------------------------------//
