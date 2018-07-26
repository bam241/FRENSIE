//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_ParticleLeavingCellEventLocalDispatcher.hpp
//! \author Alex Robinson
//! \brief  Particle leaving cell event dispatcher class declaration.
//!
//---------------------------------------------------------------------------//

// FRENSIE Includes
#include "FRENSIE_Archives.hpp"
#include "MonteCarlo_ParticleLeavingCellEventLocalDispatcher.hpp"
#include "Utility_DesignByContract.hpp"

namespace MonteCarlo{

// Constructor
ParticleLeavingCellEventLocalDispatcher::ParticleLeavingCellEventLocalDispatcher(
                                      const Geometry::Model::EntityId cell_id )
  : BaseType( cell_id )
{ /* ... */ }

// Dispatch the new event to the observers
void ParticleLeavingCellEventLocalDispatcher::dispatchParticleLeavingCellEvent(
                                 const ParticleState& particle,
                                 const Geometry::Model::EntityId cell_leaving )
{
  // Make sure the cell being entered is valid
  testPrecondition( cell_leaving == this->getEntityId() );

  ObserverSet& observer_set =
    this->getObserverSet( particle.getParticleType() );
  
  ObserverSet::iterator it = observer_set.begin();

  while( it != observer_set.end() )
  {
    it->second->updateFromParticleLeavingCellEvent( particle, cell_leaving );

    ++it;
  }
}

} // end MonteCarlo namespace

EXPLICIT_CLASS_SERIALIZE_INST( MonteCarlo, ParticleLeavingCellEventLocalDispatcher );

//---------------------------------------------------------------------------//
// end MonteCarlo_ParticleLeavingCellEventLocalDispatcher.cpp
//---------------------------------------------------------------------------//
