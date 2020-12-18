//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_WeightWindow.hpp
//! \author Philip Britt
//! \brief  Weight window class declaration
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_WEIGHT_WINDOW_HPP
#define MONTE_CARLO_WEIGHT_WINDOW_HPP

// Std Lib Includes
#include <memory>

// FRENSIE Includes
#include "MonteCarlo_PopulationControl.hpp"


namespace MonteCarlo{

//! An actual weight window object.
struct WeightWindow{
  double upper_weight;
  double survival_weight;
  double lower_weight;

  // Serialize the data
  template<typename Archive>
  void serialize( Archive& ar, const unsigned version )
  { 
    ar & BOOST_SERIALIZATION_NVP( upper_weight );
    ar & BOOST_SERIALIZATION_NVP( survival_weight );
    ar & BOOST_SERIALIZATION_NVP( lower_weight );
  }
};

//! The weight window base class
class WeightWindowBase: public PopulationControl
{

public:

  //! Constructor
  WeightWindowBase();

  //! Destructor
  ~WeightWindowBase()
  { /* ... */ }

  void checkParticleWithPopulationController( ParticleState& particle, 
                                              ParticleBank& bank) const;

  void setMaxSplit( const unsigned max_split_integer );

  virtual const WeightWindow& getWeightWindow( const ParticleState& particle ) const = 0;

  virtual bool isParticleInWeightWindowDiscretization( const ParticleState& particle ) const = 0;

private:

  // Declare the boost serialization access object as a friend
  friend class boost::serialization::access;

  // Serialize the data
  template<typename Archive>
  void serialize( Archive& ar, const unsigned version )
  {     
    // Serialize the base class data
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( PopulationControl ); 
  }

};

} // end MonteCarlo namespace

BOOST_SERIALIZATION_CLASS_VERSION( WeightWindowBase, MonteCarlo, 0 );
BOOST_SERIALIZATION_ASSUME_ABSTRACT_CLASS( WeightWindowBase, MonteCarlo );
EXTERN_EXPLICIT_CLASS_SERIALIZE_INST( MonteCarlo, WeightWindowBase );

#endif // end MONTE_CARLO_WEIGHT_WINDOW_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_WeightWindow.hpp
//---------------------------------------------------------------------------//
