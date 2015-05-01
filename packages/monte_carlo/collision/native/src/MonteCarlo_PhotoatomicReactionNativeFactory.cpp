//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_PhotoatomicReactionNativeFactory.cpp
//! \author Alex Robinson
//! \brief  The photoatomic reaction native data factory declaration
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <algorithm>
#include <limits>

// FRENSIE Includes
#include "MonteCarlo_PhotoatomicReactionNativeFactory.hpp"
#include "MonteCarlo_IncoherentPhotonScatteringDistributionNativeFactory.hpp"
#include "MonteCarlo_CoherentScatteringDistributionNativeFactory.hpp"
#include "MonteCarlo_IncoherentPhotoatomicReaction.hpp"
#include "MonteCarlo_SubshellIncoherentPhotoatomicReaction.hpp"
#include "MonteCarlo_CoherentPhotoatomicReaction.hpp"
#include "MonteCarlo_PairProductionPhotoatomicReaction.hpp"
#include "MonteCarlo_PhotoelectricPhotoatomicReaction.hpp"
#include "MonteCarlo_SubshellPhotoelectricPhotoatomicReaction.hpp"
#include "MonteCarlo_AbsorptionPhotoatomicReaction.hpp"
#include "MonteCarlo_SubshellType.hpp"
#include "Utility_TabularDistribution.hpp"
#include "Utility_ContractException.hpp"

namespace MonteCarlo{

// Create the total incoherent photoatomic reaction
/*! \details This will use the Waller-Hartree incoherent cross section.
 */
void PhotoatomicReactionNativeFactory::createTotalIncoherentReaction(
       const Data::ElectronPhotonRelaxationDataContainer& raw_photoatom_data,
       const Teuchos::ArrayRCP<const double>& energy_grid,
       const Teuchos::RCP<const Utility::HashBasedGridSearcher>& grid_searcher,
       Teuchos::RCP<PhotoatomicReaction>& incoherent_reaction,
       const bool use_doppler_broadening_data )
{
  // Make sure the energy grid is valid
  testPrecondition( raw_photoatom_data.getPhotonEnergyGrid().size() ==
		    energy_grid.size() );
  testPrecondition( Utility::Sort::isSortedAscending( energy_grid.begin(),
						      energy_grid.end() ) );

  // Extract the cross section
  Teuchos::ArrayRCP<double> incoherent_cross_section;
  incoherent_cross_section.assign( 
	   raw_photoatom_data.getWallerHartreeIncoherentCrossSection().begin(),
	   raw_photoatom_data.getWallerHartreeIncoherentCrossSection().end() );

  unsigned threshold_index = 
    raw_photoatom_data.getWallerHartreeIncoherentCrossSectionThresholdEnergyIndex();

  // Create the scattering distribution
  Teuchos::RCP<const IncoherentPhotonScatteringDistribution> distribution;

  if( use_doppler_broadening_data )
  {
    IncoherentPhotonScatteringDistributionNativeFactory::createAdvancedDopplerBroadenedIncoherentDistribution( 
							    raw_photoatom_data,
							    distribution,
							    3.0 );
  }
  else
  {
    IncoherentPhotonScatteringDistributionNativeFactory::createIncoherentDistribution( 
							    raw_photoatom_data,
							    distribution,
							    3.0 );
  }
  
  // Create the incoherent reaction
  incoherent_reaction.reset(
		      new IncoherentPhotoatomicReaction<Utility::LinLin,false>(
						     energy_grid,
						     incoherent_cross_section,
						     threshold_index,
						     grid_searcher,
						     distribution ) );
}

// Create the subshell incoherent photoatomic reactions
/*! \details This will use the impulse approximation incoherent cross sections
 */ 
void PhotoatomicReactionNativeFactory::createSubshellIncoherentReactions(
       const Data::ElectronPhotonRelaxationDataContainer& raw_photoatom_data,
       const Teuchos::ArrayRCP<const double>& energy_grid,
       const Teuchos::RCP<const Utility::HashBasedGridSearcher>& grid_searcher,
       Teuchos::Array<Teuchos::RCP<PhotoatomicReaction> >&
       subshell_incoherent_reactions,
       const bool use_doppler_broadening_data )
{
  // Make sure the energy grid is valid
  testPrecondition( raw_photoatom_data.getPhotonEnergyGrid().size() ==
		    energy_grid.size() );
  testPrecondition( Utility::Sort::isSortedAscending( energy_grid.begin(),
						      energy_grid.end() ) );

  subshell_incoherent_reactions.clear();

  Teuchos::RCP<PhotoatomicReaction> subshell_incoherent_reaction;

  std::set<unsigned>::const_iterator subshell_it = 
    raw_photoatom_data.getSubshells().begin();

  while( subshell_it != raw_photoatom_data.getSubshells().end() )
  {
    // Extract the cross section
    Teuchos::ArrayRCP<double> subshell_incoherent_cross_section;
    subshell_incoherent_cross_section.assign(
	   raw_photoatom_data.getImpulseApproxSubshellIncoherentCrossSection(*subshell_it).begin(),
	   raw_photoatom_data.getImpulseApproxSubshellIncoherentCrossSection(*subshell_it).end() );

    unsigned subshell_threshold_index = 
      raw_photoatom_data.getImpulseApproxSubshellIncoherentCrossSectionThresholdEnergyIndex(*subshell_it);

    // Create the occupation number distribution
    Teuchos::RCP<Utility::OneDDistribution> occupation_number(
       new Utility::TabularDistribution<Utility::LinLin>( 
	    raw_photoatom_data.getOccupationNumberMomentumGrid( *subshell_it ),
	    raw_photoatom_data.getOccupationNumber( *subshell_it ) ) );

    if( use_doppler_broadening_data )
    {
      // Create the Compton profile
      Teuchos::RCP<Utility::TabularOneDDistribution> compton_profile(
	new Utility::TabularDistribution<Utility::LinLin>(
	      raw_photoatom_data.getComptonProfileMomentumGrid( *subshell_it ),
	      raw_photoatom_data.getComptonProfile( *subshell_it ) ) );

      // Create the subshell incoherent reaction
      subshell_incoherent_reaction.reset(
	      new SubshellIncoherentPhotoatomicReaction<Utility::LinLin,false>(
		   energy_grid,
		   subshell_incoherent_cross_section,
		   subshell_threshold_index,
		   grid_searcher,
		   convertENDFDesignatorToSubshellEnum( *subshell_it ),
		   raw_photoatom_data.getSubshellBindingEnergy( *subshell_it ),
		   occupation_number,
		   compton_profile ) );
    }
    // Ignore Doppler broadening
    else
    {
      subshell_incoherent_reaction.reset(
	      new SubshellIncoherentPhotoatomicReaction<Utility::LinLin,false>(
		   energy_grid,
		   subshell_incoherent_cross_section,
		   subshell_threshold_index,
		   grid_searcher,
		   convertENDFDesignatorToSubshellEnum( *subshell_it ),
		   raw_photoatom_data.getSubshellBindingEnergy( *subshell_it ),
		   occupation_number ) );
    }

    subshell_incoherent_reactions.push_back( subshell_incoherent_reaction );

    ++subshell_it;
  }
}

// Create the coherent scattering photoatomic reaction
void PhotoatomicReactionNativeFactory::createCoherentReaction(
       const Data::ElectronPhotonRelaxationDataContainer& raw_photoatom_data,
       const Teuchos::ArrayRCP<const double>& energy_grid,
       const Teuchos::RCP<const Utility::HashBasedGridSearcher>& grid_searcher,
       Teuchos::RCP<PhotoatomicReaction>& coherent_reaction )
{
  // Make sure the energy grid is valid
  testPrecondition( raw_photoatom_data.getPhotonEnergyGrid().size() ==
		    energy_grid.size() );
  testPrecondition( Utility::Sort::isSortedAscending( energy_grid.begin(),
						      energy_grid.end() ) );

  Teuchos::ArrayRCP<double> coherent_cross_section;
  coherent_cross_section.assign(
	     raw_photoatom_data.getWallerHartreeCoherentCrossSection().begin(),
	     raw_photoatom_data.getWallerHartreeCoherentCrossSection().end() );

  unsigned threshold_index = 
    raw_photoatom_data.getWallerHartreeCoherentCrossSectionThresholdEnergyIndex();

  // Create the coherent scattering distribution
  Teuchos::RCP<const CoherentScatteringDistribution> distribution;

  CoherentScatteringDistributionNativeFactory::createEfficientCoherentDistribution(
					                    raw_photoatom_data,
							    distribution ); 
  
  // Create the coherent reaction
  coherent_reaction.reset(
		        new CoherentPhotoatomicReaction<Utility::LinLin,false>(
							energy_grid,
							coherent_cross_section,
							threshold_index,
							grid_searcher,
							distribution ) );
}

// Create the pair production photoatomic reaction
void PhotoatomicReactionNativeFactory::createPairProductionReaction(
       const Data::ElectronPhotonRelaxationDataContainer& raw_photoatom_data,
       const Teuchos::ArrayRCP<const double>& energy_grid,
       const Teuchos::RCP<const Utility::HashBasedGridSearcher>& grid_searcher,
       Teuchos::RCP<PhotoatomicReaction>& pair_production_reaction,
       const bool use_detailed_pair_production_data )
{
  // Make sure the energy grid is valid
  testPrecondition( raw_photoatom_data.getPhotonEnergyGrid().size() ==
		    energy_grid.size() );
  testPrecondition( Utility::Sort::isSortedAscending( energy_grid.begin(),
						      energy_grid.end() ) );

  Teuchos::ArrayRCP<double> pair_production_cross_section;
  pair_production_cross_section.assign(
	     raw_photoatom_data.getPairProductionCrossSection().begin(),
	     raw_photoatom_data.getPairProductionCrossSection().end() );

  unsigned threshold_index = 
    raw_photoatom_data.getPairProductionCrossSectionThresholdEnergyIndex();

  pair_production_reaction.reset(
		  new PairProductionPhotoatomicReaction<Utility::LinLin,false>(
					 energy_grid,
					 pair_production_cross_section,
					 threshold_index,
					 grid_searcher,
					 use_detailed_pair_production_data ) );
}

// Create the total photoelectric photoatomic reaction
void PhotoatomicReactionNativeFactory::createTotalPhotoelectricReaction(
       const Data::ElectronPhotonRelaxationDataContainer& raw_photoatom_data,
       const Teuchos::ArrayRCP<const double>& energy_grid,
       const Teuchos::RCP<const Utility::HashBasedGridSearcher>& grid_searcher,
       Teuchos::RCP<PhotoatomicReaction>& photoelectric_reaction )
{
  // Make sure the energy grid is valid
  testPrecondition( raw_photoatom_data.getPhotonEnergyGrid().size() ==
		    energy_grid.size() );
  testPrecondition( Utility::Sort::isSortedAscending( energy_grid.begin(),
						      energy_grid.end() ) );

  Teuchos::ArrayRCP<double> photoelectric_cross_section;
  photoelectric_cross_section.assign(
	     raw_photoatom_data.getPhotoelectricCrossSection().begin(),
	     raw_photoatom_data.getPhotoelectricCrossSection().end() );

  unsigned threshold_index = 
    raw_photoatom_data.getPhotoelectricCrossSectionThresholdEnergyIndex();

  // Create the total photoelectric reaction
  photoelectric_reaction.reset( 
		   new PhotoelectricPhotoatomicReaction<Utility::LinLin,false>(
						   energy_grid,
						   photoelectric_cross_section,
						   threshold_index,
						   grid_searcher ) );
}
					       
// Create the subshell photoelectric photoatomic reactions
void PhotoatomicReactionNativeFactory::createSubshellPhotoelectricReactions(
       const Data::ElectronPhotonRelaxationDataContainer& raw_photoatom_data,
       const Teuchos::ArrayRCP<const double>& energy_grid,
       const Teuchos::RCP<const Utility::HashBasedGridSearcher>& grid_searcher,
       Teuchos::Array<Teuchos::RCP<PhotoatomicReaction> >&
       subshell_photoelectric_reactions )
{
  // Make sure the energy grid is valid
  testPrecondition( raw_photoatom_data.getPhotonEnergyGrid().size() ==
		    energy_grid.size() );
  testPrecondition( Utility::Sort::isSortedAscending( energy_grid.begin(),
						      energy_grid.end() ) );

  subshell_photoelectric_reactions.clear();

  Teuchos::RCP<PhotoatomicReaction> subshell_photoelectric_reaction;

  std::set<unsigned>::const_iterator subshell_it = 
    raw_photoatom_data.getSubshells().begin();

  while( subshell_it != raw_photoatom_data.getSubshells().end() )
  {
    // Extract the cross section
    Teuchos::ArrayRCP<double> subshell_photoelectric_cross_section;
    subshell_photoelectric_cross_section.assign(
      raw_photoatom_data.getSubshellPhotoelectricCrossSection(*subshell_it).begin(),
      raw_photoatom_data.getSubshellPhotoelectricCrossSection(*subshell_it).end() );

    unsigned subshell_threshold_index = 
      raw_photoatom_data.getSubshellPhotoelectricCrossSectionThresholdEnergyIndex( *subshell_it );

    subshell_photoelectric_reaction.reset(
	   new SubshellPhotoelectricPhotoatomicReaction<Utility::LinLin,false>(
	       energy_grid,
	       subshell_photoelectric_cross_section,
	       subshell_threshold_index,
	       grid_searcher,
	       convertENDFDesignatorToSubshellEnum( *subshell_it ),
	       raw_photoatom_data.getSubshellBindingEnergy( *subshell_it ) ) );

    subshell_photoelectric_reactions.push_back( 
					     subshell_photoelectric_reaction );

    ++subshell_it;
  }
}

// Create the heating photoatomic reaction
void PhotoatomicReactionNativeFactory::createHeatingReaction(
       const Data::ElectronPhotonRelaxationDataContainer& raw_photoatom_data,
       const Teuchos::ArrayRCP<const double>& energy_grid,
       const Teuchos::RCP<const Utility::HashBasedGridSearcher>& grid_searcher,
       Teuchos::RCP<PhotoatomicReaction>& heating_reaction )
{
  // Make sure the energy grid is valid
  testPrecondition( raw_photoatom_data.getPhotonEnergyGrid().size() ==
		    energy_grid.size() );
  testPrecondition( Utility::Sort::isSortedAscending( energy_grid.begin(),
						      energy_grid.end() ) );

  Teuchos::ArrayRCP<double> heating_cross_section;
  heating_cross_section.assign(
	     raw_photoatom_data.getAveragePhotonHeatingNumbers().begin(),
	     raw_photoatom_data.getAveragePhotonHeatingNumbers().end() );

  unsigned threshold_index = 0u;

  // Create the heating reaction
  heating_reaction.reset(
		      new AbsorptionPhotoatomicReaction<Utility::LinLin,false>(
					      energy_grid,
					      heating_cross_section,
					      threshold_index,
					      grid_searcher,
					      HEATING_PHOTOATOMIC_REACTION ) );
}

} // end MonteCarlo namespace

//---------------------------------------------------------------------------//
// end MonteCarlo_PhotoatomicReactionNativeFactory.cpp
//---------------------------------------------------------------------------//
