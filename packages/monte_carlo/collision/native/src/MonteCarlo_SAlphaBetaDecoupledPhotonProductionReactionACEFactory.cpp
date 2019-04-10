//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_SAlphaBetaDecoupledPhotonProductionNuclearReactionACEFactory.cpp
//! \author Eli Moll
//! \brief  Nuclear reaction factory class declaration
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <limits>

// Trilinos Includes
#include <Teuchos_ScalarTraits.hpp>

// FRENSIE Includes
#include "MonteCarlo_DecoupledPhotonProductionReaction.hpp"
#include "MonteCarlo_SAlphaBetaDecoupledPhotonProductionReactionACEFactory.hpp"
#include "MonteCarlo_PhotonProductionNuclearScatteringDistributionACEFactory.hpp"
#include "MonteCarlo_DecoupledYieldBasedPhotonProductionReaction.hpp"
#include "MonteCarlo_DecoupledCrossSectionBasedPhotonProductionReaction.hpp"
#include "MonteCarlo_NuclearScatteringDistribution.hpp"
#include "MonteCarlo_NeutronAbsorptionReaction.hpp"
#include "Utility_ExceptionTestMacros.hpp"
#include "Utility_ContractException.hpp"
#include "Utility_InterpolationPolicy.hpp"

namespace MonteCarlo{

// Constructor
/*! \details All blocks from the ACE file will be stored. 
 * \param[in] raw_nuclide_data The necessary data blocks will be extracted 
 * using the data extractor.
 */
SAlphaBetaDecoupledPhotonProductionReactionACEFactory::SAlphaBetaDecoupledPhotonProductionReactionACEFactory( 
		 const std::string& table_name,
		 const double atomic_weight_ratio,
		 const double temperature,
		 const Teuchos::ArrayRCP<const double>& energy_grid,
		 const Data::XSSNeutronDataExtractor& raw_nuclide_data,
		 const Data::XSSSabDataExtractor& sab_nuclide_data )
		 : SAlphaBetaNuclearReactionACEFactory( table_name,
		                              atomic_weight_ratio,
		                              temperature,
		                              energy_grid,
		                              raw_nuclide_data,
		                              sab_nuclide_data )
{ 
  // Create the scattering distribution factory
  PhotonProductionNuclearScatteringDistributionACEFactory 
    photon_production_dist_factory( table_name,
			                              atomic_weight_ratio,
			                              raw_nuclide_data );

  // Extract the required blocks
  Teuchos::ArrayView<const double> mtrp_block = 
    raw_nuclide_data.extractMTRPBlock();
  Teuchos::ArrayView<const double> lsigp_block = 
    raw_nuclide_data.extractLSIGPBlock();
  Teuchos::ArrayView<const double> sigp_block = 
    raw_nuclide_data.extractSIGPBlock();
  
  // Create a map of the reaction types and their table ordering
  boost::unordered_map<unsigned,unsigned> reaction_ordering;
  SAlphaBetaDecoupledPhotonProductionReactionACEFactory::createReactionOrderingMap( mtrp_block,
							reaction_ordering );

  // Parse the SIGP data and create the necessary data maps
  boost::unordered_map<unsigned,Teuchos::ArrayView<const double> > yield_energy_map;
  boost::unordered_map<unsigned,Teuchos::ArrayView<const double> > yield_values_map;
  boost::unordered_map<unsigned,Teuchos::ArrayRCP<double> > xs_based_map;
  boost::unordered_map<unsigned,unsigned> threshold_energy_map;
  boost::unordered_map<unsigned,NuclearReactionType> base_reaction_type_map;
  
  SAlphaBetaDecoupledPhotonProductionReactionACEFactory::parseSIGP( 
                                                     lsigp_block,
                                                     sigp_block,
                                                     reaction_ordering,
                                                     yield_energy_map,
                                                     yield_values_map,
                                                     xs_based_map,
                                                     threshold_energy_map,
                                                     base_reaction_type_map );

  // Construct a map of required base reaction classes
  boost::unordered_map<NuclearReactionType,Teuchos::RCP<NuclearReaction> > base_reaction_map;
  
  SAlphaBetaDecoupledPhotonProductionReactionACEFactory::constructBaseReactionMap(
                                                    base_reaction_type_map,
                                                    base_reaction_map,
                                                    yield_energy_map );
                                                    
  // Construct a map of photon MT numbers to yield distributions
  SAlphaBetaDecoupledPhotonProductionReactionACEFactory::constructMTPYieldDistributions(
                                                        yield_energy_map,
                                                        yield_values_map );

  // Construct a map of photon MT numbers to yield distributions
  SAlphaBetaDecoupledPhotonProductionReactionACEFactory::constructMTYieldArrays(
                                                      base_reaction_type_map,
                                                      yield_energy_map );

  // Construct the total neutron cross section
  SAlphaBetaDecoupledPhotonProductionReactionACEFactory::createTotalReaction(
                                   raw_nuclide_data.extractTotalCrossSection(),
                                   energy_grid,
                                   temperature );

  // Create the yield based photon production reactions
  initializeYieldBasedPhotonProductionReactions( base_reaction_type_map,
	                                               temperature,
	                                               yield_energy_map,
	                                               base_reaction_map,
	                                               photon_production_dist_factory );
	     
	// Create the cross section based photon production reactions
  initializeCrossSectionBasedPhotonProductionReactions( base_reaction_type_map,
	                                                      temperature,
	                                                      threshold_energy_map,
	                                                      xs_based_map,
                                                        energy_grid,
                                                        photon_production_dist_factory );
}

// Create the photon production reactions 
void SAlphaBetaDecoupledPhotonProductionReactionACEFactory::createPhotonProductionReactions( 
      boost::unordered_map<unsigned,Teuchos::RCP<DecoupledPhotonProductionReaction> >&
      photon_production_reactions ) const
{
  photon_production_reactions.insert( d_photon_production_reactions.begin(),
			       d_photon_production_reactions.end() );
}

// Create the total reaction for weight normalization
void SAlphaBetaDecoupledPhotonProductionReactionACEFactory::createTotalReaction(
                       const Teuchos::ArrayView<const double>& total_xs_block,
                       const Teuchos::ArrayRCP<const double>& energy_grid,
                       const double temperature )
{
  Teuchos::ArrayRCP<double> total_cross_section;
  total_cross_section.deepCopy(total_xs_block);

  d_total_reaction.reset( new NeutronAbsorptionReaction( N__TOTAL_REACTION,
							 temperature,
							 0.0,
							 0u,
							 energy_grid,
							 total_cross_section ) );
}

// Create the reaction type ordering map
void SAlphaBetaDecoupledPhotonProductionReactionACEFactory::createReactionOrderingMap( 
        const Teuchos::ArrayView<const double>& mtrp_block,
        boost::unordered_map<unsigned,unsigned>& reaction_ordering )

{  
  unsigned reaction;
  
  for( unsigned i = 0; i < mtrp_block.size(); ++i )
  {
    reaction = static_cast<unsigned>( mtrp_block[i] );
    
    reaction_ordering[reaction] = i;
  }
}

// Parse the SIGP to create the yield energy map, yield values, xs map, and 
//   threshold energy map
void SAlphaBetaDecoupledPhotonProductionReactionACEFactory::parseSIGP(
  const Teuchos::ArrayView<const double>& lsigp_block,
  const Teuchos::ArrayView<const double>& sigp_block,
  const boost::unordered_map<unsigned,unsigned>& reaction_ordering,
  boost::unordered_map<unsigned,Teuchos::ArrayView<const double> >& yield_energy_map,
  boost::unordered_map<unsigned,Teuchos::ArrayView<const double> >& yield_values_map,
  boost::unordered_map<unsigned,Teuchos::ArrayRCP<double> >& xs_based_map,
  boost::unordered_map<unsigned,unsigned>& threshold_energy_map,
  boost::unordered_map<unsigned,NuclearReactionType>& base_reaction_type_map )
{
  boost::unordered_map<unsigned,unsigned>::const_iterator
    reaction, end_reaction;
  reaction = reaction_ordering.begin();
  end_reaction = reaction_ordering.end();

  unsigned cs_index;
  unsigned cs_array_size;
  unsigned energy_array_size;
  
  while( reaction != end_reaction )
  {
    cs_index = static_cast<unsigned>( lsigp_block[reaction->second] ) - 1u;
       
	  if ( static_cast<unsigned>( sigp_block[cs_index] ) == 13u )
	  { 
      cs_array_size = static_cast<unsigned>( sigp_block[cs_index + 2u] );

      Teuchos::ArrayRCP<double>& cross_section = 
	                                    xs_based_map[reaction->first];

      cross_section.deepCopy( sigp_block( cs_index + 3u, cs_array_size ) );
      
      threshold_energy_map[reaction->first] = 
                      static_cast<unsigned>( sigp_block[cs_index + 1u] ) ;
                      
      base_reaction_type_map[reaction->first] = 
	          convertUnsignedToNuclearReactionType( reaction->first/1000u );
	  }
	  else if ( static_cast<unsigned>( sigp_block[cs_index] ) == 12u  ||
	            static_cast<unsigned>( sigp_block[cs_index] ) == 16u )
	  {
	    std::cout << " " << std::endl;
	    std::cout << "SIGP Block Interpolation Regions" << std::endl;
	    std::cout << sigp_block(cs_index + 1u, 3) << std::endl;
	  
	    TEST_FOR_EXCEPTION( static_cast<unsigned>( sigp_block[cs_index + 2u] ) != 0,
	                        std::runtime_error,
	                        "Error: multiple interpolation regions were defined "
	                        "in the ACE table for MT " << reaction->first <<
	                        "." );
	                        
	    energy_array_size =  static_cast<unsigned>( sigp_block[cs_index + 3u] );
	                     
	    yield_energy_map[reaction->first] = 
	                              sigp_block( cs_index + 4u, energy_array_size );
	                             
	    yield_values_map[reaction->first] =
	          sigp_block( cs_index + 4u + energy_array_size, energy_array_size );
	          
	    base_reaction_type_map[reaction->first] = 
	          convertUnsignedToNuclearReactionType( reaction->first/1000u );
	  }
	  else
	  {
	    THROW_EXCEPTION( std::runtime_error,
	                     "Error: MFTYPE was found to be " << 
	                     static_cast<unsigned>( sigp_block[cs_index + 2u] ) <<
	                     " which is not equal to 12, 13, or 16 (only allowable"
	                     " entries).");
	  }
    ++reaction;
  }
}

// Create the base reaction map
void SAlphaBetaDecoupledPhotonProductionReactionACEFactory::constructBaseReactionMap(
  boost::unordered_map<unsigned,NuclearReactionType>& base_reaction_type_map,
  boost::unordered_map<NuclearReactionType,Teuchos::RCP<NuclearReaction> >& base_reaction_map,
  boost::unordered_map<unsigned,Teuchos::ArrayView<const double> >& yield_energy_map )
{
  boost::unordered_map<unsigned,Teuchos::ArrayView<const double> >::const_iterator
    reaction, end_reaction;
  reaction = yield_energy_map.begin();
  end_reaction = yield_energy_map.end();
  
  Teuchos::RCP<NuclearReaction> base_reaction;
  
  while( reaction != end_reaction )
  {
    this->getReactionFromReactionType( 
                          base_reaction_type_map.find(reaction->first)->second,
                          base_reaction );
    base_reaction_map[base_reaction_type_map.find(reaction->first)->second] = 
                                                                 base_reaction;
    
    ++reaction;
  }
}

// Construct a map of photon MT numbers to yield distributions
void SAlphaBetaDecoupledPhotonProductionReactionACEFactory::constructMTPYieldDistributions(
	     const boost::unordered_map<unsigned,Teuchos::ArrayView<const double> >& yield_energy_map,
	     const boost::unordered_map<unsigned,Teuchos::ArrayView<const double> >& yield_values_map )
{
  boost::unordered_map<unsigned,Teuchos::ArrayView<const double> >::const_iterator
    iter_reaction, end_reaction;
  iter_reaction = yield_energy_map.begin();
  end_reaction = yield_energy_map.end();
  
  while( iter_reaction != end_reaction )
  {
    unsigned reaction_type = iter_reaction->first;
    
    std::shared_ptr<Utility::OneDDistribution> tabular_yield_pointer( 
        new Utility::TabularDistribution<Utility::LinLin>( 
                              yield_energy_map.find(reaction_type)->second,
                              yield_values_map.find(reaction_type)->second ) );
                                          
    d_mtp_yield_distributions_map[reaction_type] = tabular_yield_pointer;
		  
	  ++iter_reaction;
  }
}

// Construct a map of base reaction types to yield distribution arrays
void SAlphaBetaDecoupledPhotonProductionReactionACEFactory::constructMTYieldArrays(
       const boost::unordered_map<unsigned,NuclearReactionType>& base_reaction_type_map,
       const boost::unordered_map<unsigned,Teuchos::ArrayView<const double> >& yield_energy_map )
{
  boost::unordered_map<unsigned,Teuchos::ArrayView<const double> >::const_iterator
    iter_reaction, end_reaction;
  iter_reaction = yield_energy_map.begin();
  end_reaction = yield_energy_map.end();
  
  while( iter_reaction != end_reaction )
  {
    unsigned reaction_type = iter_reaction->first;
    
    NuclearReactionType base_reaction_type = 
                            base_reaction_type_map.find(reaction_type)->second;
                                          
    d_mt_yield_distributions[base_reaction_type].push_back(
                  d_mtp_yield_distributions_map.find(reaction_type)->second );
		  
	  ++iter_reaction;
  }
}

// Initialize the yield based photon production reactions
void SAlphaBetaDecoupledPhotonProductionReactionACEFactory::initializeYieldBasedPhotonProductionReactions( 
       const boost::unordered_map<unsigned,NuclearReactionType>& base_reaction_type_map,
	     const double temperature,
	     const boost::unordered_map<unsigned,Teuchos::ArrayView<const double> >& yield_energy_map,
	     const boost::unordered_map<NuclearReactionType,Teuchos::RCP<NuclearReaction> >& base_reaction_map,
	     PhotonProductionNuclearScatteringDistributionACEFactory photon_production_dist_factory )	
{
  boost::unordered_map<unsigned,Teuchos::ArrayView<const double> >::const_iterator
    iter_reaction, end_reaction;
  iter_reaction = yield_energy_map.begin();
  end_reaction = yield_energy_map.end();
  
  Teuchos::RCP<NuclearScatteringDistribution<NeutronState,PhotonState> > 
    photon_production_distribution;
  
  while( iter_reaction != end_reaction )
  {
    unsigned reaction_type = iter_reaction->first;
   
    Teuchos::RCP<DecoupledPhotonProductionReaction>& reaction = 
                    d_photon_production_reactions[reaction_type];
    
    photon_production_dist_factory.createScatteringDistribution(
						     reaction_type,
						     photon_production_distribution );
    
    reaction.reset( new DecoupledYieldBasedPhotonProductionReaction(
		  base_reaction_type_map.find(reaction_type)->second,
		  reaction_type,
		  temperature,
		  d_mt_yield_distributions[base_reaction_type_map.find(reaction_type)->second],
		  d_mtp_yield_distributions_map[reaction_type],
		  base_reaction_map.find(base_reaction_type_map.find(reaction_type)->second)->second,
		  photon_production_distribution,
		  d_total_reaction ) );  
		  
	  ++iter_reaction;
  }

}

// Initialize the yield based photon production reactions
void SAlphaBetaDecoupledPhotonProductionReactionACEFactory::initializeCrossSectionBasedPhotonProductionReactions( 
  const boost::unordered_map<unsigned,NuclearReactionType>& base_reaction_type_map,
  const double temperature,
  const boost::unordered_map<unsigned,unsigned>& threshold_energy_map,
  const boost::unordered_map<unsigned,Teuchos::ArrayRCP<double> >& xs_based_map,
  const Teuchos::ArrayRCP<const double>& energy_grid,
  PhotonProductionNuclearScatteringDistributionACEFactory photon_production_dist_factory )			
{
  boost::unordered_map<unsigned,unsigned>::const_iterator
    iter_reaction, end_reaction;
  iter_reaction = threshold_energy_map.begin();
  end_reaction = threshold_energy_map.end();
  
  Teuchos::RCP<NuclearScatteringDistribution<NeutronState,PhotonState> >
    photon_production_distribution;
  
  while( iter_reaction != end_reaction )
  {
    unsigned reaction_type = iter_reaction->first;
    
    Teuchos::RCP<DecoupledPhotonProductionReaction>& reaction = 
                    d_photon_production_reactions[reaction_type];
    
    photon_production_dist_factory.createScatteringDistribution(
						     reaction_type,
						     photon_production_distribution );		     
						     
    reaction.reset( new DecoupledCrossSectionBasedPhotonProductionReaction(
		  base_reaction_type_map.find(reaction_type)->second,
		  reaction_type,
		  temperature,
		  threshold_energy_map.find(reaction_type)->second,
		  energy_grid,
		  xs_based_map.find(reaction_type)->second,
		  photon_production_distribution,
		  d_total_reaction,
		  Teuchos::Array<std::shared_ptr<Utility::OneDDistribution> >() ) ); 
		  
	  ++iter_reaction;  
  }
}

} // end MonteCarlo namespace

//---------------------------------------------------------------------------//
// end MonteCarlo_NuclearReactionACEFactory.cpp
//---------------------------------------------------------------------------//