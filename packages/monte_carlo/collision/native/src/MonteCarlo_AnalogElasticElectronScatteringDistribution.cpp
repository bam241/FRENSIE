//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_AnalogElasticElectronScatteringDistribution.cpp
//! \author Luke Kersting
//! \brief  The analog elastic electron scattering distribution definition
//!
//---------------------------------------------------------------------------//

// FRENSIE Includes
#include "MonteCarlo_AnalogElasticElectronScatteringDistribution.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_SearchAlgorithms.hpp"
#include "Utility_DirectionHelpers.hpp"
#include "Utility_KinematicHelpers.hpp"
#include "Utility_PhysicalConstants.hpp"

namespace MonteCarlo{


// Initialize static member data

// The change in scattering angle cosine below which the screened Rutherford distribution is used
double AnalogElasticElectronScatteringDistribution::s_cutoff_delta_mu = 1.0e-6;

// The scattering angle cosine above which the screened Rutherford distribution is used
double AnalogElasticElectronScatteringDistribution::s_cutoff_mu = 0.999999;

// The fine structure constant squared
double AnalogElasticElectronScatteringDistribution::s_fine_structure_const_squared=
        Utility::PhysicalConstants::fine_structure_constant *
        Utility::PhysicalConstants::fine_structure_constant;

// A parameter for moliere's screening factor  (1/2*(fsc/0.885)**2)
double AnalogElasticElectronScatteringDistribution::s_screening_param1 =
      AnalogElasticElectronScatteringDistribution::s_fine_structure_const_squared/
      ( 2.0*0.885*0.885 );

// Constructor
AnalogElasticElectronScatteringDistribution::AnalogElasticElectronScatteringDistribution(
    const std::shared_ptr<TwoDDist>& elastic_cutoff_distribution,
    const int atomic_number,
    const bool linlinlog_interpolation_mode_on )
  : d_elastic_cutoff_distribution( elastic_cutoff_distribution ),
    d_atomic_number( atomic_number ),
    d_linlinlog_interpolation_mode_on( linlinlog_interpolation_mode_on ),
    d_Z_two_thirds_power( pow( atomic_number, 2.0/3.0 ) ),
    d_screening_param2( 3.76*s_fine_structure_const_squared*
                              d_atomic_number*d_atomic_number )
{
  // Make sure the array is valid
  testPrecondition( d_elastic_cutoff_distribution.use_count() > 0 );
  // Make sure the atomic number is valid
  testPrecondition( d_atomic_number > 0 );
  testPrecondition( d_atomic_number <= 100u );
}

// Evaluate the distribution at the given energy and scattering angle cosine
//! \details Because the scattering angle cosine is very close to one, precision will be lost.
double AnalogElasticElectronScatteringDistribution::evaluate(
        const double incoming_energy,
        const double scattering_angle_cosine ) const
{
  // Make sure the energy and angle are valid
  testPrecondition( incoming_energy > 0.0 );
  testPrecondition( scattering_angle_cosine >= -1.0 );
  testPrecondition( scattering_angle_cosine <= 1.0 );

  if ( scattering_angle_cosine > s_cutoff_mu )
  {
    double eta = evaluateMoliereScreeningConstant( incoming_energy );

    // evaluate on the screened Rutherford distribution
    return this->evaluateScreenedRutherford( 
                    incoming_energy,
                    scattering_angle_cosine,
                    eta );
  }
  else
  {
    // evaluate on the cutoff distribution
    return d_elastic_cutoff_distribution->evaluateExact( incoming_energy,
                                                         scattering_angle_cosine );
  }
}

// Evaluate the screened Rutherford distribution given energy, eta, and scattering angle cosine
//! \details Because the scattering angle cosine is very close to one, precision will be lost.
double AnalogElasticElectronScatteringDistribution::evaluateScreenedRutherford(
        const double incoming_energy,
        const double scattering_angle_cosine,
        const double eta ) const
{
  // Make sure the energy, eta and angle are valid
  testPrecondition( incoming_energy > 0.0 );
  testPrecondition( scattering_angle_cosine >= s_cutoff_mu );
  testPrecondition( scattering_angle_cosine <= 1.0 );
  testPrecondition( eta > 0.0 );

  double cutoff_pdf =
    d_elastic_cutoff_distribution->evaluateExact( incoming_energy, s_cutoff_mu );

  double pdf = evaluateScreenedRutherfordPDF( incoming_energy,
                                              scattering_angle_cosine,
                                              eta,
                                              cutoff_pdf );

  return pdf;
}

// Evaluate the PDF at the given energy and scattering angle cosine
//! \details This PDF is normalize to equal 1 when integrated from mu = -1.0 to mu = s_cutoff_mu (0.999999)
double AnalogElasticElectronScatteringDistribution::evaluatePDF(
        const double incoming_energy,
        const double scattering_angle_cosine ) const
{
  // Make sure the energy, eta and angle are valid
  testPrecondition( incoming_energy > 0.0 );
  testPrecondition( scattering_angle_cosine >= -1.0 );
  testPrecondition( scattering_angle_cosine <= 1.0 );

  if ( scattering_angle_cosine > s_cutoff_mu )
  {
    double eta = evaluateMoliereScreeningConstant( incoming_energy );

    // evaluate PDF on the screened Rutherford distribution
    return this->evaluateScreenedRutherfordPDF(
                    incoming_energy,
                    scattering_angle_cosine,
                    eta );
  }
  else
  {
    // evaluate PDF on the cutoff distribution
    return d_elastic_cutoff_distribution->evaluateSecondaryConditionalPDFExact(
            incoming_energy,
            scattering_angle_cosine );
  }
}

// Evaluate the screened Rutherford PDF at the given energy and scattering angle cosine
//! \details This screened Rutherford pdf uses the normalized cutoff pdf value at s_cutoff_mu (0.999999)
double AnalogElasticElectronScatteringDistribution::evaluateScreenedRutherfordPDF(
        const double incoming_energy,
        const double scattering_angle_cosine,
        const double eta ) const
{
  // Make sure the energy, eta and angle are valid
  testPrecondition( incoming_energy > 0.0 );
  testPrecondition( scattering_angle_cosine >= -1.0 );
  testPrecondition( scattering_angle_cosine <= 1.0 );
  testPrecondition( eta > 0.0 );

  double cutoff_pdf =
    d_elastic_cutoff_distribution->evaluateSecondaryConditionalPDFExact(
            incoming_energy,
            s_cutoff_mu );

  return evaluateScreenedRutherfordPDF( incoming_energy,
                                        scattering_angle_cosine,
                                        eta,
                                        cutoff_pdf );
}

// Evaluate the screened Rutherford PDF at the given energy and scattering angle cosine
double AnalogElasticElectronScatteringDistribution::evaluateScreenedRutherfordPDF(
        const double incoming_energy,
        const double scattering_angle_cosine,
        const double eta,
        const double norm_factor ) const
{
  // Make sure the energy, eta and angle are valid
  testPrecondition( incoming_energy > 0.0 );
  testPrecondition( scattering_angle_cosine >= -1.0 );
  testPrecondition( scattering_angle_cosine <= 1.0 );
  testPrecondition( eta > 0.0 );
  testPrecondition( norm_factor > 0.0 );

  double delta_mu = 1.0 - scattering_angle_cosine;

  return norm_factor*
            ( s_cutoff_delta_mu + eta )*( s_cutoff_delta_mu + eta )/(
            ( delta_mu + eta )*( delta_mu + eta ) );
}

// Evaluate the CDF
//! \details This CDF is normalize to 1 at s_cutoff_mu (0.999999)
double AnalogElasticElectronScatteringDistribution::evaluateCDF(
        const double incoming_energy,
        const double scattering_angle_cosine ) const
{
  // Make sure the energy and angle are valid
  testPrecondition( incoming_energy > 0.0 );
  testPrecondition( scattering_angle_cosine >= -1.0 );
  testPrecondition( scattering_angle_cosine <= 1.0 );

  if ( scattering_angle_cosine > s_cutoff_mu )
  {
    double eta = this->evaluateMoliereScreeningConstant( incoming_energy );

    // evaluate CDF on the screened Rutherford distribution
    return this->evaluateScreenedRutherfordCDF(
                    incoming_energy,
                    scattering_angle_cosine,
                    eta );
  }
  else
  {
    // evaluate CDF on the cutoff distribution
    return d_elastic_cutoff_distribution->evaluateSecondaryConditionalCDFExact(
            incoming_energy,
            scattering_angle_cosine );
  }
}

// Evaluate the screened Rutherford CDF
//! \details This CDF is normalize to the elastic cutoff cdf at s_cutoff_mu (0.999999)
double AnalogElasticElectronScatteringDistribution::evaluateScreenedRutherfordCDF(
        const double incoming_energy,
        const double scattering_angle_cosine,
        const double eta ) const
{
  // Make sure the energy, eta and angle cosine are valid
  testPrecondition( incoming_energy > 0.0 );
  testPrecondition( scattering_angle_cosine >= s_cutoff_mu );
  testPrecondition( scattering_angle_cosine <= 1.0 );
  testPrecondition( eta > 0.0 );

  double delta_mu = 1.0 - scattering_angle_cosine;

  double cutoff_pdf =
    d_elastic_cutoff_distribution->evaluateSecondaryConditionalPDFExact(
            incoming_energy,
            s_cutoff_mu );

  return 1.0 + evaluateScreenedRutherfordCDF( incoming_energy,
                                              scattering_angle_cosine,
                                              eta,
                                              cutoff_pdf );
}

// Evaluate the screened Rutherford CDF
double AnalogElasticElectronScatteringDistribution::evaluateScreenedRutherfordCDF(
        const double incoming_energy,
        const double scattering_angle_cosine,
        const double eta,
        const double norm_factor ) const
{
  // Make sure the energy, eta and angle cosine are valid
  testPrecondition( incoming_energy > 0.0 );
  testPrecondition( scattering_angle_cosine >= s_cutoff_mu );
  testPrecondition( scattering_angle_cosine <= 1.0 );
  testPrecondition( eta > 0.0 );
  testPrecondition( norm_factor > 0.0 );

  double delta_mu = 1.0 - scattering_angle_cosine;

  return norm_factor*( scattering_angle_cosine - s_cutoff_mu )*
                ( s_cutoff_delta_mu + eta )/( delta_mu + eta );
}

// Sample an outgoing energy and direction from the distribution
void AnalogElasticElectronScatteringDistribution::sample(
         const double incoming_energy,
         double& outgoing_energy,
         double& scattering_angle_cosine ) const
{
  // The outgoing energy is always equal to the incoming energy
  outgoing_energy = incoming_energy;

  unsigned trial_dummy;

  // Sample an outgoing direction
  this->sampleAndRecordTrialsImpl( incoming_energy,
                                   scattering_angle_cosine,
                                   trial_dummy );
}

// Sample an outgoing energy and direction and record the number of trials
void AnalogElasticElectronScatteringDistribution::sampleAndRecordTrials(
        const double incoming_energy,
        double& outgoing_energy,
        double& scattering_angle_cosine,
        unsigned& trials ) const
{
  // The outgoing energy is always equal to the incoming energy
  outgoing_energy = incoming_energy;

  // Sample an outgoing direction
  this->sampleAndRecordTrialsImpl( incoming_energy,
                                   scattering_angle_cosine,
                                   trials );
}

// Randomly scatter the electron
void AnalogElasticElectronScatteringDistribution::scatterElectron(
         ElectronState& electron,
         ParticleBank& bank,
         Data::SubshellType& shell_of_interaction ) const
{
  double scattering_angle_cosine;

  unsigned trial_dummy;

  // Sample an outgoing direction
  this->sampleAndRecordTrialsImpl( electron.getEnergy(),
                                   scattering_angle_cosine,
                                   trial_dummy );

  shell_of_interaction =Data::UNKNOWN_SUBSHELL;

  // Set the new direction
  electron.rotateDirection( scattering_angle_cosine,
                            this->sampleAzimuthalAngle() );
}

// Randomly scatter the adjoint electron
void AnalogElasticElectronScatteringDistribution::scatterAdjointElectron(
         AdjointElectronState& adjoint_electron,
         ParticleBank& bank,
         Data::SubshellType& shell_of_interaction ) const
{
  double scattering_angle_cosine;

  unsigned trial_dummy;

  // Sample an outgoing direction
  this->sampleAndRecordTrialsImpl( adjoint_electron.getEnergy(),
                                   scattering_angle_cosine,
                                   trial_dummy );

  shell_of_interaction = Data::UNKNOWN_SUBSHELL;

  // Set the new direction
  adjoint_electron.rotateDirection( scattering_angle_cosine,
                                    this->sampleAzimuthalAngle() );
}

// Evaluate Moliere's atomic screening constant (modified by Seltzer) at the given electron energy
double AnalogElasticElectronScatteringDistribution::evaluateMoliereScreeningConstant(
                                              const double energy ) const
{
  // get the energy-momentum**2 of the electron in units of electron_rest_mass_energy
  double electron_energy_momentum_squared =
           Utility::calculateDimensionlessRelativisticMomentumSquared(
                          Utility::PhysicalConstants::electron_rest_mass_energy,
                          energy );

  // get the velocity of the electron divided by the speed of light beta = v/c
  double beta_squared = Utility::calculateDimensionlessRelativisticSpeedSquared(
           Utility::PhysicalConstants::electron_rest_mass_energy,
           energy );

  double screening_param3 = 1.0/beta_squared*sqrt( energy/
            ( energy + Utility::PhysicalConstants::electron_rest_mass_energy) );

 // Calculate Moliere's atomic screening constant
 return s_screening_param1 * 1.0/electron_energy_momentum_squared *
        d_Z_two_thirds_power * ( 1.13 + d_screening_param2*screening_param3 );
}

// Sample an outgoing direction from the distribution
void AnalogElasticElectronScatteringDistribution::sampleAndRecordTrialsImpl(
                                                const double incoming_energy,
                                                double& scattering_angle_cosine,
                                                unsigned& trials ) const
{
  // Make sure the incoming energy is valid
  testPrecondition( incoming_energy > 0.0 );

  // Increment the number of trials
  ++trials;

  // Find the lower and upper distribution bins
  double upper_energy, lower_energy;

  // Find the bin boundaries
  TwoDDist::DistributionType::const_iterator lower_bin, upper_bin;

  d_elastic_cutoff_distribution->findBinBoundaries( incoming_energy,
                                                    lower_bin,
                                                    upper_bin );

  // Get a random number
  double random_number =
            Utility::RandomNumberGenerator::getRandomNumber<double>();

  if ( lower_bin->first == incoming_energy )
  {
    sampleBin( lower_bin, random_number, scattering_angle_cosine );
  }  
  else if ( upper_bin->first == incoming_energy )
  {
    sampleBin( upper_bin, random_number, scattering_angle_cosine );
  }
  else if ( lower_bin != upper_bin )
  {
    // Sample lower bin
    double lower_angle;
    sampleBin( lower_bin, random_number, lower_angle );

    // Sample upper bin
    double upper_angle;
    sampleBin( upper_bin, random_number, upper_angle );

    if ( d_linlinlog_interpolation_mode_on )
    {
      // LinLinLog interpolation between energy bins
      scattering_angle_cosine = Utility::LinLog::interpolate(
                                  lower_bin->first,
                                  upper_bin->first,
                                  incoming_energy,
                                  lower_angle,
                                  upper_angle );
    }
    else
    {
      // LinLinLin interpolation between energy bins
      scattering_angle_cosine = Utility::LinLin::interpolate(
                                  lower_bin->first,
                                  upper_bin->first,
                                  incoming_energy,
                                  lower_angle,
                                  upper_angle );
    }
  }
  else
  {
    scattering_angle_cosine = 1.0;
  }

  // Make sure the scattering angle cosine is valid
  testPostcondition( scattering_angle_cosine >= -1.0 );
  testPostcondition( scattering_angle_cosine <= 1.0 );
}

// Sample an outgoing direction from the given distribution
/*! \details Due to roundoff error, that algorithm used to calculate the
 *  scattering angle cosine can sometimes return a number slightly greater that
 *  1.0. If this is the case, the scattering angle cosine is set to 1.0.
 */
void AnalogElasticElectronScatteringDistribution::sampleBin(
        const TwoDDist::DistributionType::const_iterator& distribution_bin,
        const double random_number,
        double& scattering_angle_cosine ) const
{
  // Get the bin energy
  double energy = distribution_bin->first;
  // Get the scattering constant at the bin energy
  double eta = evaluateMoliereScreeningConstant( energy );
  // Get maximum CDF at the bin energy
  double max_cdf = evaluateScreenedRutherfordCDF( energy, 1.0, eta );

  if ( max_cdf*random_number > 1.0 ) // Sample screened Rutherford
  {
    // calculated a reapeated variable
    double var = s_cutoff_delta_mu*random_number;

    // calculate the screened Rutherford scattering angle
    scattering_angle_cosine = std::min( 1.0,
        ( var*( 1.0 + eta ) + eta*s_cutoff_mu )/
        ( var + eta ) );

    // Make sure the scattering angle cosine is valid
    testPostcondition( scattering_angle_cosine >= s_cutoff_mu );
  }
  else // Sample Cutoff
  {
    // calculate the cutoff distribution scattering angle
    scattering_angle_cosine =
      distribution_bin->second->sampleWithRandomNumber( max_cdf*random_number );

    // Make sure the scattering angle cosine is valid
    testPostcondition( scattering_angle_cosine <= s_cutoff_mu );
  }
  // Make sure the scattering angle cosine is valid
  testPostcondition( scattering_angle_cosine <= 1.0 );
}

} // end MonteCarlo namespace

//---------------------------------------------------------------------------//
// end MonteCarlo_AnalogElasticElectronScatteringDistribution.cpp
//---------------------------------------------------------------------------//
