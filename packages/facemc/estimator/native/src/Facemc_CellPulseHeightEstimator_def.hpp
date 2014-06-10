//---------------------------------------------------------------------------//
//!
//! \file   Facemc_CellPulseHeightEstimator_def.hpp
//! \author Alex Robinson
//! \brief  Cell pulse height estimator class definition
//!
//---------------------------------------------------------------------------//

#ifndef FACEMC_CELL_PULSE_HEIGHT_ESTIMATOR_DEF_HPP
#define FACEMC_CELL_PULSE_HEIGHT_ESTIMATOR_DEF_HPP

// Std Lib Includes
#include <iostream>

// FRENSIE Includes
#include "Utility_ContractException.hpp"

namespace Facemc{

// Constructor
template<typename ContributionMultiplierPolicy>
CellPulseHeightEstimator<
                       ContributionMultiplierPolicy>::CellPulseHeightEstimator(
       const Estimator::idType id,
       const double multiplier,
       const Teuchos::Array<CellPulseHeightEstimator::cellIdType>& entity_ids,
       const bool auto_register_with_dispatchers )			   
  : EntityEstimator<cellIdType>( id, multiplier, entity_ids ),
    ParticleGenerationEventObserver( id, 
				     entity_ids, 
				     auto_register_with_dispatchers ),
    d_total_energy_deposition_moments( 1 )
{
  // Set up the entity map
  for( unsigned i = 0; i < entity_ids.size(); ++i )
  {
    // Ignore duplicate entity ids
    if( d_cell_energy_deposition_map.count( entity_ids[i] ) == 0 )
      d_cell_energy_deposition_map[ entity_ids[i] ] = 0.0;
  }
}

// Set the response functions
template<typename ContributionMultiplierPolicy>
void CellPulseHeightEstimator<
                           ContributionMultiplierPolicy>::setResponseFunctions(
    const Teuchos::Array<Teuchos::RCP<ResponseFunction> >& response_functions )
{
  std::cerr << "Warning: Response functions cannot be set for pulse height "
	    << "estimators. The response functions requested for pulse height "
	    << "estimator " << this->getId() << " will be ignored."
	    << std::endl;
}

// Set the particle types that can contribute to the estimator
/*! \details Only photons and electrons can contribute to this estimator
 */
template<typename ContributionMultiplierPolicy>
void CellPulseHeightEstimator<ContributionMultiplierPolicy>::setParticleTypes( 
			   const Teuchos::Array<ParticleType>& particle_types )
{
  Teuchos::Array<ParticleType> valid_particle_types;
  
  bool warning_issued = false;
  
  for( unsigned i = 0; i < particle_types.size(); ++i )
  {
    if( particle_types[i] != PHOTON )
    {
      if( !warning_issued )
      {
	std::cerr << "Warning: Only photons and charged particles can "
		  << "contribute to pulse height estimators. The other "
		  << "particle types requested for pulse height estimator "
		  << this->getId() << " will be ignored."
		  << std::endl;
      
	warning_issued = true;
      }
    }
    else
      valid_particle_types.push_back( particle_types[i] );
  }
  
  Estimator::setParticleTypes( valid_particle_types );

  testPostcondition( !this->isParticleTypeAssigned( NEUTRON ) );
  testPostcondition( !this->isParticleTypeAssigned( ADJOINT_NEUTRON ) );
  testPostcondition( !this->isParticleTypeAssigned( ADJOINT_PHOTON ) );
}

// Add estimator contribution from a portion of the current history
/*! \detials This function is to be called at a surface crossing
 */
template<typename ContributionMultiplierPolicy>
void CellPulseHeightEstimator<
		  ContributionMultiplierPolicy>::addPartialHistoryContribution(
		    const ParticleState& particle,
		    const CellPulseHeightEstimator::cellIdType& cell_leaving,
		    const CellPulseHeightEstimator::cellIdType& cell_entering )
{
  
  if( this->isParticleTypeAssigned( particle.getParticleType() ) )
  {
    double contribution = particle.getWeight()*particle.getEnergy();
    
    // Subtract the contribution from the cell being exited
    if( d_cell_energy_deposition_map.count( cell_leaving ) != 0 )
      d_cell_energy_deposition_map[cell_leaving] -= contribution;
    
    // Add the contribution to the cell being entered
    if( d_cell_energy_deposition_map.count( cell_entering ) != 0 )
      d_cell_energy_deposition_map[cell_entering] += contribution;
  }
}

// Add estimator contribution from a portion of the current history
/*! \details This function is to be called at a source point
 */
template<typename ContributionMultiplierPolicy>
void CellPulseHeightEstimator<
		  ContributionMultiplierPolicy>::addPartialHistoryContribution(
						const ParticleState& particle )
{
  // Make sure the cell is in the energy deposition map
  testPrecondition( d_cell_energy_deposition_map.count( particle.getCell() )
		    != 0 );
  
  if( this->isParticleTypeAssigned( particle.getParticleType() ) )
  {
    d_cell_energy_deposition_map[particle.getCell()] += 
      particle.getWeight()*particle.getEnergy();
  }
}

// Add estimator contribution from a portion of the current history
template<typename ContributionMultiplierPolicy>
void CellPulseHeightEstimator<
		     ContributionMultiplierPolicy>::commitHistoryContribution()
{
  typename boost::unordered_map<cellIdType,double>::iterator 
    entity, end_entity;

  entity = d_cell_energy_deposition_map.begin();
  end_entity = d_cell_energy_deposition_map.end();

  double energy_deposition_in_all_cells = 0.0;
  
  unsigned bin_index;
  double bin_contribution;

  while( entity != end_entity )
  {
    d_dimension_values[ENERGY_DIMENSION] = Teuchos::any( entity->second );
    
    if( this->isPointInEstimatorPhaseSpace( d_dimension_values ) )
    {
      bin_index = this->calculateBinIndex( d_dimension_values, 0u );
      
      bin_contribution = calculateHistoryContribution( 
					      entity->second,
					      ContributionMultiplierPolicy() );
      
      this->commitHistoryContributionToBinOfEntity( entity->first,
						    bin_index,
						    bin_contribution );
      
      // Add the energy deposition in this cell to the total energy deposition
      energy_deposition_in_all_cells += entity->second;

      // Reset the energy deposition in this cell
      entity->second = 0.0;
    }

    ++entity;
  }

  // Store the total energy deposition in the dimension values map
  d_dimension_values[ENERGY_DIMENSION] = 
    Teuchos::any( energy_deposition_in_all_cells );
  
  // Determine the pulse bin for the combination of all cells
  if( this->isPointInEstimatorPhaseSpace( d_dimension_values ) )
  {
    bin_index = this->calculateBinIndex( d_dimension_values, 0u );

    bin_contribution = calculateHistoryContribution( 
					      energy_deposition_in_all_cells,
					      ContributionMultiplierPolicy() );

    // Compute the moment contributions
    d_total_energy_deposition_moments[bin_index].first += bin_contribution;
    
    bin_contribution *= bin_contribution;
    d_total_energy_deposition_moments[bin_index].second += bin_contribution;

    bin_contribution *= bin_contribution;
    d_total_energy_deposition_moments[bin_index].third += bin_contribution;

    bin_contribution *= bin_contribution;
    d_total_energy_deposition_moments[bin_index].fourth += bin_contribution;
  }
}

// Update the observer
template<typename ContributionMultiplierPolicy>
inline void CellPulseHeightEstimator<
              ContributionMultiplierPolicy>::updateFromParticleGenerationEvent(
						const ParticleState& particle )
{
  this->addPartialHistoryContribution( particle );
}

// Print the estimator data
template<typename ContributionMultiplierPolicy>
void CellPulseHeightEstimator<ContributionMultiplierPolicy>::print( 
						       std::ostream& os ) const
{
  os << "Cell Pulse Height Estimator: " << this->getId() << std::endl;

  this->printImplementation( os, "Cell" );

  os << "All Cells" << std::endl;
  os << "--------" << std::endl;

  this->printEstimatorBinData( os,
			       d_total_energy_deposition_moments,
			       this->getTotalNormConstant() );
}

// Assign bin boundaries to an estimator dimension
template<typename ContributionMultiplierPolicy>
void CellPulseHeightEstimator<
			    ContributionMultiplierPolicy>::assignBinBoundaries(
	 const Teuchos::RCP<EstimatorDimensionDiscretization>& bin_boundaries )
{
  if( bin_boundaries->getDimension() == ENERGY_DIMENSION )
  {
    EntityEstimator<cellIdType>::assignBinBoundaries( bin_boundaries );

    d_total_energy_deposition_moments.resize( this->getNumberOfBins() );
  }
  else
  {
    std::cerr << "Warning: " << bin_boundaries->getDimensionName()
	      << " bins cannot be set for pulse height estimators. The bins "
	      << "requested for pulse height estimator " << this->getId()
	      << " will be ignored." 
	      << std::endl;
  }
}	       

// Calculate the estimator contribution from the entire history
/*! \details The multiplier policy cannot be used directly since it only
 * operates on a particle state. The policy classes will instead be used as
 * tags for tag dispatching. This function simply returns 1.0.
 */
template<typename ContributionMultiplierPolicy>
inline double CellPulseHeightEstimator<
		   ContributionMultiplierPolicy>::calculateHistoryContribution(
					        const double energy_deposition,
						WeightMultiplier )
{
  return 1.0;
}

// Calculate the estimator contribution from the entire history
/*! \details The multiplier policy cannot be used directly since it only
 * operates on a particle state. The policy classes will instead be used as
 * tags for tag dispatching. This function returns the energy that has been
 * deposited in the cell(s) of interest for the entire history.
 */
template<typename ContributionMultiplierPolicy>
inline double CellPulseHeightEstimator<
		   ContributionMultiplierPolicy>::calculateHistoryContribution(
					        const double energy_deposition,
						WeightAndEnergyMultiplier )
{
  return energy_deposition;
}

} // end Facemc namespace

#endif // end FACEMC_CELL_PULSE_HEIGHT_ESTIMATOR_DEF_HPP

//---------------------------------------------------------------------------//
// end Facemc_CellPulseHeightEstimator_def.hpp
//---------------------------------------------------------------------------//
