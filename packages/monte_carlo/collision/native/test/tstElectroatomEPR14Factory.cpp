//---------------------------------------------------------------------------//
//!
//! \file   tstElectroatomACEFactory.cpp
//! \author Luke Kersting
//! \brief  Electroatom factory using EPR14 data unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>

// Trilinos Includes
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_RCP.hpp>

// FRENSIE Includes
#include "MonteCarlo_ElectroatomACEFactory.hpp"
#include "MonteCarlo_AtomicRelaxationModelFactory.hpp"
#include "MonteCarlo_BremsstrahlungAngularDistributionType.hpp"
#include "MonteCarlo_SimulationProperties.hpp"
#include "Data_ACEFileHandler.hpp"
#include "Data_XSSEPRDataExtractor.hpp"
#include "Utility_InterpolationPolicy.hpp"
#include "Utility_PhysicalConstants.hpp"
#include "Utility_UnitTestHarnessExtensions.hpp"

//---------------------------------------------------------------------------//
// Testing Variables
//---------------------------------------------------------------------------//

Teuchos::RCP<Data::XSSEPRDataExtractor> xss_data_extractor;
Teuchos::RCP<MonteCarlo::AtomicRelaxationModel> relaxation_model;
std::string electroatom_name;
double atomic_weight;

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Check that a basic electroatom can be created
TEUCHOS_UNIT_TEST( ElectroatomACEFactory, createElectroatom_basic )
{
  MonteCarlo::SimulationProperties properties;
  properties.setBremsstrahlungAngularDistributionFunction( MonteCarlo::DIPOLE_DISTRIBUTION );
  properties.setElasticCutoffAngleCosine( 1.0 );
  properties.setAtomicRelaxationModeOff( MonteCarlo::ELECTRON );
  properties.setNumberOfElectronHashGridBins( 100 );

  Teuchos::RCP<MonteCarlo::Electroatom> atom;

  atom.reset();
  MonteCarlo::ElectroatomACEFactory::createElectroatom( *xss_data_extractor,
                                                        electroatom_name,
                                                        atomic_weight,
                                                        relaxation_model,
                                                        properties, 
                                                        atom );

  // Test the electroatom properties
  TEST_EQUALITY_CONST( atom->getAtomName(), "82000.14p" );
  TEST_EQUALITY_CONST( atom->getAtomicNumber(), 82 );
  TEST_FLOATING_EQUALITY( atom->getAtomicWeight(), 207.1999470456033, 1e-12 );

  double energy, cross_section;

  // Test that the total cross section can be returned
  energy = 2.0e-3;
  cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 4.806193787852e8, 1e-12 );

  energy = 4.0e-4;
  cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 1.278128947846e9, 1e-12 );

  energy = 9.0e-5;
  cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 2.411603154884e9, 1e-12 );


  // Test that the absorption cross section can be returned
  cross_section = atom->getAbsorptionCrossSection( 1.e-2 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 2.0e-3 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 4.0e-4 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 9.0e-5 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  MonteCarlo::ElectroatomicReactionType reaction;

  // Test that the atomic excitation cross section can be returned
  reaction = MonteCarlo::ATOMIC_EXCITATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.96517e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 6.22682e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 9.0e-5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.16042e9, 1e-12 );

  // Test that the bremsstrahlung cross section can be returned
  reaction = MonteCarlo::BREMSSTRAHLUNG_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 9.258661418255e3, 1e-12 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 8.914234996439e3, 1e-12 );

  cross_section = atom->getReactionCrossSection( 9.0e-5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 7.249970966838e3, 1e-12 );

  // Test that the decoupled elastic cross section can be returned
  reaction = MonteCarlo::DECOUPLED_ELASTIC_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction);
  TEST_FLOATING_EQUALITY( cross_section, 2.100574153670e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction);
  TEST_FLOATING_EQUALITY( cross_section, 4.436635458458e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 9.0e-5, reaction);
  TEST_FLOATING_EQUALITY( cross_section, 8.887469904554e8, 1e-12 );

  // Test that the cutoff elastic cross section can be returned
  reaction = MonteCarlo::CUTOFF_ELASTIC_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  // Test that there is no total electroionization
  reaction = MonteCarlo::TOTAL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  // Test that the K subshell electroionization cross section can be returned
  reaction = MonteCarlo::K_SUBSHELL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 8.9754e-2, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.250673571307e-1, 1e-12 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 3.64919e1, 1e-12 );

  // Test that the P3 subshell electroionization cross section can be returned
  reaction = MonteCarlo::P3_SUBSHELL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.06530e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.774672656160e7, 1e-12 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.82234e5, 1e-12 );

  atom.reset();
}

//---------------------------------------------------------------------------//
/* Check that a electroatom with detailed 2BS photon angular distribution
 * data can be created
 */
TEUCHOS_UNIT_TEST( ElectroatomACEFactory, createElectroatom_detailed_brem )
{
  MonteCarlo::SimulationProperties properties;
  properties.setBremsstrahlungAngularDistributionFunction( MonteCarlo::TWOBS_DISTRIBUTION );
  properties.setElasticCutoffAngleCosine( 1.0 );
  properties.setAtomicRelaxationModeOn( MonteCarlo::ELECTRON );
  properties.setNumberOfElectronHashGridBins( 100 );

  Teuchos::RCP<MonteCarlo::Electroatom> atom;

  MonteCarlo::ElectroatomACEFactory::createElectroatom( *xss_data_extractor,
                                                        electroatom_name,
                                                        atomic_weight,
                                                        relaxation_model,
                                                        properties, 
                                                        atom );

  // Test the electroatom properties
  TEST_EQUALITY_CONST( atom->getAtomName(), "82000.14p" );
  TEST_EQUALITY_CONST( atom->getAtomicNumber(), 82 );
  TEST_FLOATING_EQUALITY( atom->getAtomicWeight(), 207.1999470456033, 1e-12 );

  double energy, cross_section;

  // Test that the total cross section can be returned
  energy = 2.0e-3;
  cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 4.806193787852e8, 1e-12 );

  energy = 4.0e-4;
  cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 1.278128947846e9, 1e-12 );

  energy = 9.0e-5;
  cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 2.411603154884e9, 1e-12 );

  // Test that the absorption cross section can be returned
  cross_section = atom->getAbsorptionCrossSection( 1.0e-2 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 2.0e-3 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 4.0e-4 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 9.0e-5 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  MonteCarlo::ElectroatomicReactionType reaction;

  // Test that the atomic excitation cross section can be returned
  reaction = MonteCarlo::ATOMIC_EXCITATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.96517e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 6.22682e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 9.0e-5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.16042e9, 1e-12 );

  // Test that the bremsstrahlung cross section can be returned
  reaction = MonteCarlo::BREMSSTRAHLUNG_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 9.258661418255e3, 1e-12 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 8.914234996439e3, 1e-12 );

  cross_section = atom->getReactionCrossSection( 9.0e-5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 7.249970966838e3, 1e-12 );

  // Test that the decoupled elastic cross section can be returned
  reaction = MonteCarlo::DECOUPLED_ELASTIC_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction);
  TEST_FLOATING_EQUALITY( cross_section, 2.100574153670e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction);
  TEST_FLOATING_EQUALITY( cross_section, 4.436635458458e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 9.0e-5, reaction);
  TEST_FLOATING_EQUALITY( cross_section, 8.887469904554e8, 1e-12 );

  // Test that the cutoff elastic cross section can be returned
  reaction = MonteCarlo::CUTOFF_ELASTIC_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  // Test that there is no total electroionization
  reaction = MonteCarlo::TOTAL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  // Test that the K subshell electroionization cross section can be returned
  reaction = MonteCarlo::K_SUBSHELL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 8.9754e-2, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.250673571307e-1, 1e-12 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 3.64919e1, 1e-12 );

  // Test that the P3 subshell electroionization cross section can be returned
  reaction = MonteCarlo::P3_SUBSHELL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.06530e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.774672656160e7, 1e-12 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.82234e5, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that a electroatom without elastic subshell data can be created
TEUCHOS_UNIT_TEST( ElectroatomACEFactory, createElectroatom_no_elastic )
{
  MonteCarlo::SimulationProperties properties;
  properties.setBremsstrahlungAngularDistributionFunction( MonteCarlo::DIPOLE_DISTRIBUTION );
  properties.setElasticCutoffAngleCosine( 1.0 );
  properties.setAtomicRelaxationModeOn( MonteCarlo::ELECTRON );
  properties.setNumberOfElectronHashGridBins( 100 );
  properties.setElasticModeOff();

  Teuchos::RCP<MonteCarlo::Electroatom> atom;

  MonteCarlo::ElectroatomACEFactory::createElectroatom( *xss_data_extractor,
                                                        electroatom_name,
                                                        atomic_weight,
                                                        relaxation_model,
                                                        properties,
                                                        atom );

  // Test the electroatom properties
  TEST_EQUALITY_CONST( atom->getAtomName(), "82000.14p" );
  TEST_EQUALITY_CONST( atom->getAtomicNumber(), 82 );
  TEST_FLOATING_EQUALITY( atom->getAtomicWeight(), 207.1999470456033, 1e-12 );

  // Test that the total cross section can be returned
  double energy = 2.0e-3;
  double cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 4.806193787852e8-2.100574153670e8, 1e-12 );

  energy = 4.0e-4;
  cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 1.278128947846e9-4.436635458458e8, 1e-12 );

  energy = 9.0e-5;
  cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 2.411603154884e9-8.887469904554e8, 1e-12 );


  // Test that the absorption cross section can be returned
  cross_section = atom->getAbsorptionCrossSection( 1.e-2 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 2.0e-3 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 4.0e-4 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 9.0e-5 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  MonteCarlo::ElectroatomicReactionType reaction;

  // Test that the atomic excitation cross section can be returned
  reaction = MonteCarlo::ATOMIC_EXCITATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.96517e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 6.22682e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 9.0e-5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.16042e9, 1e-12 );

  // Test that the bremsstrahlung cross section can be returned
  reaction = MonteCarlo::BREMSSTRAHLUNG_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 9.258661418255e3, 1e-12 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 8.914234996439e3, 1e-12 );

  cross_section = atom->getReactionCrossSection( 9.0e-5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 7.249970966838e3, 1e-12 );

  // Test that the decoupled elastic cross section can be returned
  reaction = MonteCarlo::DECOUPLED_ELASTIC_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  // Test that the cutoff elastic cross section can be returned
  reaction = MonteCarlo::CUTOFF_ELASTIC_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  // Test that there is no total electroionization
  reaction = MonteCarlo::TOTAL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  // Test that the K subshell electroionization cross section can be returned
  reaction = MonteCarlo::K_SUBSHELL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 8.9754e-2, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.250673571307e-1, 1e-12 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 3.64919e1, 1e-12 );

  // Test that the P3 subshell electroionization cross section can be returned
  reaction = MonteCarlo::P3_SUBSHELL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.06530e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.774672656160e7, 1e-12 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.82234e5, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that a electroatom without electroionization subshell data can be created
TEUCHOS_UNIT_TEST( ElectroatomACEFactory, createElectroatom_no_electroionization )
{
  MonteCarlo::SimulationProperties properties;
  properties.setBremsstrahlungAngularDistributionFunction( MonteCarlo::DIPOLE_DISTRIBUTION );
  properties.setElasticCutoffAngleCosine( 1.0 );
  properties.setAtomicRelaxationModeOn( MonteCarlo::ELECTRON );
  properties.setNumberOfElectronHashGridBins( 100 );
  properties.setElectroionizationModeOff();

  Teuchos::RCP<MonteCarlo::Electroatom> atom;

  MonteCarlo::ElectroatomACEFactory::createElectroatom( *xss_data_extractor,
                                                        electroatom_name,
                                                        atomic_weight,
                                                        relaxation_model,
                                                        properties,
                                                        atom );

  // Test the electroatom properties
  TEST_EQUALITY_CONST( atom->getAtomName(), "82000.14p" );
  TEST_EQUALITY_CONST( atom->getAtomicNumber(), 82 );
  TEST_FLOATING_EQUALITY( atom->getAtomicWeight(), 207.1999470456033, 1e-12 );

  // Test that the total cross section can be returned
  double energy = 2.0e-3;
  double cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 4.806193787852e8-7.403570475672e7, 1e-12 );

  energy = 4.0e-4;
  cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 1.278128947846e9-2.117744877649e8, 1e-12 );

  energy = 9.0e-5;
  cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 2.411603154884e9-3.624289144575e8, 1e-12 );


  // Test that the absorption cross section can be returned
  cross_section = atom->getAbsorptionCrossSection( 1.e-2 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 2.0e-3 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 4.0e-4 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 9.0e-5 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  MonteCarlo::ElectroatomicReactionType reaction;

  // Test that the atomic excitation cross section can be returned
  reaction = MonteCarlo::ATOMIC_EXCITATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.96517e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 6.22682e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 9.0e-5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.16042e9, 1e-12 );

  // Test that the bremsstrahlung cross section can be returned
  reaction = MonteCarlo::BREMSSTRAHLUNG_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 9.258661418255e3, 1e-12 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 8.914234996439e3, 1e-12 );

  cross_section = atom->getReactionCrossSection( 9.0e-5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 7.249970966838e3, 1e-12 );

  // Test that the decoupled elastic cross section can be returned
  reaction = MonteCarlo::DECOUPLED_ELASTIC_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction);
  TEST_FLOATING_EQUALITY( cross_section, 2.100574153670e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction);
  TEST_FLOATING_EQUALITY( cross_section, 4.436635458458e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 9.0e-5, reaction);
  TEST_FLOATING_EQUALITY( cross_section, 8.887469904554e8, 1e-12 );

  // Test that the cutoff elastic cross section can be returned
  reaction = MonteCarlo::CUTOFF_ELASTIC_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  // Test that there is no total electroionization
  reaction = MonteCarlo::TOTAL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  // Test that the K subshell electroionization cross section can be returned
  reaction = MonteCarlo::K_SUBSHELL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  // Test that the P3 subshell electroionization cross section can be returned
  reaction = MonteCarlo::P3_SUBSHELL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );
}

//---------------------------------------------------------------------------//
// Check that a electroatom without bremsstrahlung data can be created
TEUCHOS_UNIT_TEST( ElectroatomACEFactory, createElectroatom_no_bremsstrahlung )
{
  MonteCarlo::SimulationProperties properties;
  properties.setBremsstrahlungAngularDistributionFunction( MonteCarlo::DIPOLE_DISTRIBUTION );
  properties.setElasticCutoffAngleCosine( 1.0 );
  properties.setAtomicRelaxationModeOn( MonteCarlo::ELECTRON );
  properties.setNumberOfElectronHashGridBins( 100 );
  properties.setBremsstrahlungModeOff();

  Teuchos::RCP<MonteCarlo::Electroatom> atom;
 
  MonteCarlo::ElectroatomACEFactory::createElectroatom( *xss_data_extractor,
                                                        electroatom_name,
                                                        atomic_weight,
                                                        relaxation_model,
                                                        properties,
                                                        atom );


  // Test the electroatom properties
  TEST_EQUALITY_CONST( atom->getAtomName(), "82000.14p" );
  TEST_EQUALITY_CONST( atom->getAtomicNumber(), 82 );
  TEST_FLOATING_EQUALITY( atom->getAtomicWeight(), 207.1999470456033, 1e-12 );

  // Test that the total cross section can be returned
  double energy = 2.0e-3;
  double cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 4.806193787852e8-9.258661418255e3, 1e-12 );

  energy = 4.0e-4;
  cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 1.278128947846e9-8.914234996439e3, 1e-12 );

  energy = 9.0e-5;
  cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 2.411603154884e9-7.249970966838e3, 1e-12 );


  // Test that the absorption cross section can be returned
  cross_section = atom->getAbsorptionCrossSection( 1.e-2 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 2.0e-3 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 4.0e-4 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 9.0e-5 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  MonteCarlo::ElectroatomicReactionType reaction;

  // Test that the atomic excitation cross section can be returned
  reaction = MonteCarlo::ATOMIC_EXCITATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.96517e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 6.22682e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 9.0e-5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.16042e9, 1e-12 );

  // Test that the bremsstrahlung cross section can be returned
  reaction = MonteCarlo::BREMSSTRAHLUNG_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  // Test that the decoupled elastic cross section can be returned
  reaction = MonteCarlo::DECOUPLED_ELASTIC_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction);
  TEST_FLOATING_EQUALITY( cross_section, 2.100574153670e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction);
  TEST_FLOATING_EQUALITY( cross_section, 4.436635458458e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 9.0e-5, reaction);
  TEST_FLOATING_EQUALITY( cross_section, 8.887469904554e8, 1e-12 );

  // Test that the cutoff elastic cross section can be returned
  reaction = MonteCarlo::CUTOFF_ELASTIC_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  // Test that there is no total electroionization
  reaction = MonteCarlo::TOTAL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  // Test that the K subshell electroionization cross section can be returned
  reaction = MonteCarlo::K_SUBSHELL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 8.9754e-2, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.250673571307e-1, 1e-12 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 3.64919e1, 1e-12 );

  // Test that the P3 subshell electroionization cross section can be returned
  reaction = MonteCarlo::P3_SUBSHELL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.06530e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.774672656160e7, 1e-12 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.82234e5, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that a electroatom without atomic excitation data can be created
TEUCHOS_UNIT_TEST( ElectroatomACEFactory, createElectroatom_no_atomic_excitation )
{
  MonteCarlo::SimulationProperties properties;
  properties.setBremsstrahlungAngularDistributionFunction( MonteCarlo::DIPOLE_DISTRIBUTION );
  properties.setElasticCutoffAngleCosine( 1.0 );
  properties.setAtomicRelaxationModeOn( MonteCarlo::ELECTRON );
  properties.setNumberOfElectronHashGridBins( 100 );
  properties.setAtomicExcitationModeOff();

  Teuchos::RCP<MonteCarlo::Electroatom> atom;

  MonteCarlo::ElectroatomACEFactory::createElectroatom( *xss_data_extractor,
                                                        electroatom_name,
                                                        atomic_weight,
                                                        relaxation_model,
                                                        properties,
                                                        atom );

  // Test the electroatom properties
  TEST_EQUALITY_CONST( atom->getAtomName(), "82000.14p" );
  TEST_EQUALITY_CONST( atom->getAtomicNumber(), 82 );
  TEST_FLOATING_EQUALITY( atom->getAtomicWeight(), 207.1999470456033, 1e-12 );

  // Test that the total cross section can be returned
  double energy = 2.0e-3;
  double cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 4.806193787852e8-1.965170e8, 1e-12 );

  energy = 4.0e-4;
  cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 1.278128947846e9-6.22682e08, 1e-12 );

  energy = 9.0e-5;
  cross_section = atom->getTotalCrossSection( energy );
  TEST_FLOATING_EQUALITY( cross_section, 2.411603154884e9-1.16042e09, 1e-12 );

  // Test that the absorption cross section can be returned
  cross_section = atom->getAbsorptionCrossSection( 1.e-2 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 2.0e-3 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 4.0e-4 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getAbsorptionCrossSection( 9.0e-5 );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  MonteCarlo::ElectroatomicReactionType reaction;

  // Test that the atomic excitation cross section can be returned
  reaction = MonteCarlo::ATOMIC_EXCITATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  // Test that the bremsstrahlung cross section can be returned
  reaction = MonteCarlo::BREMSSTRAHLUNG_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 9.258661418255e3, 1e-12 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 8.914234996439e3, 1e-12 );

  cross_section = atom->getReactionCrossSection( 9.0e-5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 7.249970966838e3, 1e-12 );

  // Test that the decoupled elastic cross section can be returned
  reaction = MonteCarlo::DECOUPLED_ELASTIC_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction);
  TEST_FLOATING_EQUALITY( cross_section, 2.100574153670e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction);
  TEST_FLOATING_EQUALITY( cross_section, 4.436635458458e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 9.0e-5, reaction);
  TEST_FLOATING_EQUALITY( cross_section, 8.887469904554e8, 1e-12 );

  // Test that the cutoff elastic cross section can be returned
  reaction = MonteCarlo::CUTOFF_ELASTIC_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 4.0e-4, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction);
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  // Test that there is no total electroionization
  reaction = MonteCarlo::TOTAL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  // Test that the K subshell electroionization cross section can be returned
  reaction = MonteCarlo::K_SUBSHELL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_EQUALITY_CONST( cross_section, 0.0 );

  cross_section = atom->getReactionCrossSection( 8.9754e-2, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.250673571307e-1, 1e-12 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 3.64919e1, 1e-12 );

  // Test that the P3 subshell electroionization cross section can be returned
  reaction = MonteCarlo::P3_SUBSHELL_ELECTROIONIZATION_ELECTROATOMIC_REACTION;
  cross_section = atom->getReactionCrossSection( 1.0e-5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.06530e8, 1e-12 );

  cross_section = atom->getReactionCrossSection( 2.0e-3, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.774672656160e7, 1e-12 );

  cross_section = atom->getReactionCrossSection( 1.0e5, reaction );
  TEST_FLOATING_EQUALITY( cross_section, 1.82234e5, 1e-12 );
}

//---------------------------------------------------------------------------//
// Custom setup
//---------------------------------------------------------------------------//
UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_SETUP_BEGIN();

std::string test_ace_file_name, test_ace_table_name;

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_COMMAND_LINE_OPTIONS()
{
  clp().setOption( "test_ace_file",
                   &test_ace_file_name,
                   "Test ACE file name" );
  clp().setOption( "test_ace_table",
                   &test_ace_table_name,
                   "Test ACE table name" );
}

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_DATA_INITIALIZATION()
{
  {
    // Create a file handler and data extractor
    Teuchos::RCP<Data::ACEFileHandler> ace_file_handler(
        new Data::ACEFileHandler( test_ace_file_name,
                                  test_ace_table_name,
                                  1u ) );
    xss_data_extractor.reset(
        new Data::XSSEPRDataExtractor( ace_file_handler->getTableNXSArray(),
                                       ace_file_handler->getTableJXSArray(),
                                       ace_file_handler->getTableXSSArray() ) );

    MonteCarlo::AtomicRelaxationModelFactory::createAtomicRelaxationModel(
           *xss_data_extractor,
           relaxation_model,
           1e-3,
           1e-5,
           true );

    electroatom_name = test_ace_table_name;
    atomic_weight = ace_file_handler->getTableAtomicWeightRatio()*
      Utility::PhysicalConstants::neutron_rest_mass_amu;
  }
}

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_SETUP_END();

//---------------------------------------------------------------------------//
// end tstElectroatomACEFactory.cpp
//---------------------------------------------------------------------------//