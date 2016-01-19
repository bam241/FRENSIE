//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_CoupledCompleteDopplerBroadenedPhotonEnergyDistribution_def.hpp
//! \author Alex Robinson
//! \brief  The coupled complete Doppler broadened photon energy dist. def.
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_COUPLED_COMPLETE_DOPPLER_BROADENED_PHOTON_ENERGY_DISTRIBUTION_DEF_HPP
#define MONTE_CARLO_COUPLED_COMPLETE_DOPPLER_BROADENED_PHOTON_ENERGY_DISTRIBUTION_DEF_HPP

// FRENSIE Includes
#include "MonteCarlo_CoupledCompleteDopplerBroadenedPhotonEnergyDistribution.hpp"
#include "MonteCarlo_PhotonKinematicsHelpers.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_ContractException.hpp"

namespace MonteCarlo{

// Constructor
/*! \details The Compton profile grids must be in me*c units (not atomic 
 * units). The Compton profiles must be in inverse me*c units (not inverse 
 * atomic units). Only half profiles should be provided (grid ranges from 0.0 
 * to 1.0).
 */
CoupledCompleteDopplerBroadenedPhotonEnergyDistribution::CoupledCompleteDopplerBroadenedPhotonEnergyDistribution(
                const Teuchos::Array<double>& subshell_binding_energies,
                const Teuchos::Array<double>& subshell_occupancies,
                const Teuchos::Array<SubshellType>& subshell_order,
                const std::shared_ptr<const ComptonProfileSubshellConverter>&
                subshell_converter,
                const ElectronMomentumDistArray& electron_momentum_dist_array )
  : CompleteDopplerBroadenedPhotonEnergyDistribution( subshell_occupancies,
						      subshell_order ),
    d_subshell_converter( subshell_converter ),
    d_subshell_binding_energies( subshell_binding_energies ),
{
  // Make sure the shell interaction data is valid
  testPrecondition( subshell_occupancies.size() > 0 );
  testPrecondition( subshell_order.size() == 
		    subshell_occupancies.size() );
  testPrecondition( subshell_binding_energies.size() ==
		    subshell_occupancies.size() );
}

// Return the binding energy of a subshell
template<typename ComptonProfilePolicy>
double CoupledCompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::getSubshellBindingEnergy( 
                                            const SubshellType subshell ) const
{
  // Make sure the subshell is valid
  testPrecondition( this->isValidSubshell( subshell ) );

  unsigned endf_subshell_index = this->getENDFSubshellIndex( subshell );

  return d_subshell_binding_energy[endf_subshell_index];
}

// Sample an interaction subshell
/*! \details The old subshell index used to select the Compton profile and
 * and the binding energy is the same as the subshell (i.e. they are coupled).
 */
template<typename ComptonProfilePolicy>
void CoupledCompleteDopplerBroadenedPhotonEnergyDistribution<ComptonProfilePolicy>::sampleInteractionSubshell( 
                                               unsigned& old_subshell_index,
                                               double& subshell_binding_energy,
                                               Subshell& subshell ) const
{
  subshell = this->sampleENDFInteractionSubshell();

  subshell_binding_energy = this->getSubshellBindingEnergy( subshell );
  
  old_subshell_index = this->getOldSubshellIndex( subshell );
}

} // end MonteCarlo namespace

#endif // end MONTE_CARLO_COUPLED_COMPLETE_DOPPLER_BROADENED_PHOTON_ENERGY_DISTRIBUTION_DEF_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_CoupledCompleteDopplerBroadenedPhotonEnergyDistribution_def.hpp
//---------------------------------------------------------------------------//
