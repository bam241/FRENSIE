//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_ScreenedRutherfordElasticElectronScatteringDistribution.hpp
//! \author Luke Kersting
//! \brief  The screened Rutherford elastic electron scattering distribution base class
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_SCREENED_RUTHERFORD_ELASTIC_ELECTRON_SCATTERING_DISTRIBUTION_HPP
#define MONTE_CARLO_SCREENED_RUTHERFORD_ELASTIC_ELECTRON_SCATTERING_DISTRIBUTION_HPP

// Std Lib Includes
#include <limits>

// Boost Includes
#include <boost/function.hpp>
#include <boost/bind.hpp>

// Trilinos Includes
#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

// FRENSIE Includes
#include "MonteCarlo_ElectronState.hpp"
#include "MonteCarlo_ParticleBank.hpp"
#include "Utility_TabularDistribution.hpp"
#include "Utility_TabularOneDDistribution.hpp"
#include "MonteCarlo_ElectronScatteringDistribution.hpp"
#include "MonteCarlo_AdjointElectronScatteringDistribution.hpp"
#include "MonteCarlo_CutoffElasticElectronScatteringDistribution.hpp"

namespace MonteCarlo{

//! The scattering distribution base class
class ScreenedRutherfordElasticElectronScatteringDistribution : public ElectronScatteringDistribution,
    public AdjointElectronScatteringDistribution
{

public:

  //! Typedef for the array of energy dependent screened rutherford paramters
  //! (first = energy, second = Moliere screening constant,
  //!  third = normalization constant
  typedef Teuchos::Array<Utility::Trip<double,double,double> > ParameterArray;

  typedef Teuchos::RCP<const CutoffElasticElectronScatteringDistribution>
            ElasticDistribution;

  //! Constructor
  ScreenedRutherfordElasticElectronScatteringDistribution(
    const ElasticDistribution& elastic_cutoff_distribution,
    const int atomic_number );

  //! Destructor
  virtual ~ScreenedRutherfordElasticElectronScatteringDistribution()
  { /* ... */ }

  //! Evaluate the distribution
  double evaluate( const double incoming_energy,
                   const double scattering_angle_cosine ) const;

  //! Evaluate the distribution
  double evaluate( const double incoming_energy,
                   const double scattering_angle_cosine,
                   const double eta ) const;

  //! Evaluate the PDF
  double evaluatePDF( const double incoming_energy,
                      const double scattering_angle_cosine ) const;

  //! Evaluate the PDF
  double evaluatePDF( const double incoming_energy,
                      const double scattering_angle_cosine,
                      const double eta ) const;

  //! Evaluate the CDF
  double evaluateCDF( const double incoming_energy,
                      const double scattering_angle_cosine ) const;
  //! Evaluate the CDF
  double evaluateCDF( const double incoming_energy,
                      const double scattering_angle_cosine,
                      const double eta ) const;

  //! Sample an outgoing energy and direction from the distribution
  void sample( const double incoming_energy,
               double& outgoing_energy,
               double& scattering_angle_cosine ) const;

  //! Sample an outgoing energy and direction and record the number of trials
  void sampleAndRecordTrials( const double incoming_energy,
                              double& outgoing_energy,
                              double& scattering_angle_cosine,
                              unsigned& trials ) const;

  //! Randomly scatter the electron
  void scatterElectron( ElectronState& electron,
                        ParticleBank& bank,
                        Data::SubshellType& shell_of_interaction ) const;

  //! Randomly scatter the adjoint electron
  void scatterAdjointElectron( AdjointElectronState& adjoint_electron,
                               ParticleBank& bank,
                               Data::SubshellType& shell_of_interaction ) const;

  //! Evaluate Moliere's atomic screening constant at the given electron energy
  double evaluateMoliereScreeningConstant( const double energy ) const;

  //! Evaluate the integrated PDF
  double evaluateIntegratedPDF( const double incoming_energy ) const;

  //! Evaluate the integrated PDF
  double evaluateIntegratedPDF( const double incoming_energy,
                                const double eta ) const;

protected:

   //! Sample an outgoing direction from the distribution
  void sampleAndRecordTrialsImpl( const double incoming_energy,
                                  double& scattering_angle_cosine,
                                  unsigned& trials ) const;

private:

  // The change scattering angle cosine below which the screened Rutherford distribution is used
  static double s_cutoff_delta_mu;

  // The scattering angle cosine above which the screened Rutherford distribution is used
  static double s_cutoff_mu;

  // The fine structure constant (fsc) squared
  static double s_fine_structure_const_squared;

  // A parameter for moliere's screening factor  (1/2*(fsc/0.885)**2)
  static double s_screening_param1;

  // Atomic number (Z) of the target atom
  int d_atomic_number;

  // Atomic number (Z) of the target atom to the 2/3 power (Z^2/3)
  double d_Z_two_thirds_power;

  // A parameter for moliere's screening factor (3.76*fsc**2*Z**2)
  double d_screening_param2;

  // Cutoff elastic scattering distribution
  ElasticDistribution d_elastic_cutoff_distribution;

  // Screened Rutherford energy depended paramters: Moliere's screening constant and normalization constant
  ParameterArray d_screened_rutherford_parameters;
};

} // end MonteCarlo namespace

#endif // end MONTE_CARLO_SCREENED_RUTHERFORD_ELASTIC_ELECTRON_SCATTERING_DISTRIBUTION_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_ScreenedRutherfordElasticElectronScatteringDistribution.hpp
//---------------------------------------------------------------------------//
