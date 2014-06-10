//---------------------------------------------------------------------------//
//!
//! \file   Facemc_StandardSurfaceEstimator.hpp
//! \author Alex Robinson
//! \brief  Standard surface estimator class declaration.
//!
//---------------------------------------------------------------------------//

#ifndef FACEMC_STANDARD_SURFACE_ESTIMATOR_HPP
#define FACEMC_STANDARD_SURFACE_ESTIMATOR_HPP

// FRENSIE Includes
#include "Facemc_StandardEntityEstimator.hpp"
#include "Geometry_ModuleTraits.hpp"

namespace Facemc{

//! The standard surface estimator base class
class StandardSurfaceEstimator : public StandardEntityEstimator<Geometry::ModuleTraits::InternalSurfaceHandle>
{

public:

  //! Typedef for the surface id type
  typedef Geometry::ModuleTraits::InternalSurfaceHandle surfaceIdType;

  //! Set the angle cosine cutoff value
  static void setAngleCosineCutoff( const double angle_cosine_cutoff );

  //! Constructor
  StandardSurfaceEstimator( const Estimator::idType id,
			    const double multiplier,
			    const Teuchos::Array<surfaceIdType>& surface_ids,
			    const Teuchos::Array<double>& surface_areas );

  //! Destructor
  virtual ~StandardSurfaceEstimator()
  { /* ... */ }

  //! Set the particle types that can contribute to the estimator
  void setParticleTypes( const Teuchos::Array<ParticleType>& particle_types );

  //! Add estimator contribution from a portion of the current history
  virtual void addPartialHistoryContribution( 
					   const ParticleState& particle,
					   const surfaceIdType surface_crossed,
					   const double angle_cosine ) = 0;

protected:

  //! Get the angle cosine cutoff value
  static double getAngleCosineCutoff();

private:

  // Angle cosine cutoff value (default = 0.01)
  static double angle_cosine_cutoff;
};

// Get the angle cosine cutoff value
inline double StandardSurfaceEstimator::getAngleCosineCutoff()
{
  return StandardSurfaceEstimator::angle_cosine_cutoff;
}

} // end Facemc namespace

#endif // end FACEMC_STANDARD_SURFACE_ESTIMATOR_HPP

//---------------------------------------------------------------------------//
// end Facemc_StandardSurfaceEstimator.hpp
//---------------------------------------------------------------------------//