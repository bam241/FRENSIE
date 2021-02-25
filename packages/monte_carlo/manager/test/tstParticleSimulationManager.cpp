//---------------------------------------------------------------------------//
//!
//! \file   tstParticleSimulationManager.cpp
//! \author Alex Robinson
//! \brief  The particle simulation manager unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>
#include <memory>
#include <csignal>
#include <functional>
#include <map>

// Boost Includes
#include <boost/filesystem.hpp>

// FRENSIE Includes
#include "MonteCarlo_ParticleSimulationManagerFactory.hpp"
#include "MonteCarlo_StandardParticleSource.hpp"
#include "MonteCarlo_StandardParticleSourceComponent.hpp"
#include "MonteCarlo_StandardAdjointParticleSourceComponent.hpp"
#include "MonteCarlo_StandardParticleDistribution.hpp"
#include "Data_ScatteringCenterPropertiesDatabase.hpp"
#include "Geometry_InfiniteMediumModel.hpp"
#include "Utility_UnitTestHarnessWithMain.hpp"
#include "ArchiveTestHelpers.hpp"
#include "FRENSIE_config.hpp"




//import PyFrensie.Geometry as Geometry
//import PyFrensie.Geometry.DagMC as DagMC

//#include "PyFrensie_PythonTypeTraits.hpp"
#include "Geometry_InfiniteMediumNavigator.hpp"
#include "Geometry_InfiniteMediumModel.hpp"
#include "Geometry_DagMCModelProperties.hpp"
#include "Geometry_DagMCModel.hpp"
#include "Geometry_DagMCNavigator.hpp"
#include "Geometry_Exceptions.hpp"
#include "Utility_SerializationHelpers.hpp"
#include "Utility_DesignByContract.hpp"


//import PyFrensie.Utility.Mesh as Mesh
#include "Utility_StructuredHexMesh.hpp"
//import PyFrensie.Utility.Distribution as Distribution
#include "Utility_HistogramDistribution.hpp"
//import PyFrensie.Utility.DirectionDiscretization as DirectionDiscretization
#include "Utility_PQLAQuadrature.hpp"

//import PyFrensie.MonteCarlo.Collision as Collision
//import PyFrensie.MonteCarlo.ActiveRegion as ActiveRegion
#include "MonteCarlo_StandardParticleDistribution.hpp"
#include "MonteCarlo_StandardParticleSourceComponent.hpp"
#include "MonteCarlo_StandardParticleSource.hpp"
//import PyFrensie.MonteCarlo.Event as Event
#include "MonteCarlo_CellCollisionFluxEstimator.hpp"
#include "MonteCarlo_MeshTrackLengthFluxEstimator.hpp"
//import PyFrensie.MonteCarlo.Manager as Manager
#include "MonteCarlo_ParticleSimulationManager.hpp"
#include "MonteCarlo_StandardParticleSimulationManager.hpp"
#include "MonteCarlo_ParticleSimulationManagerFactory.hpp"

//import PyFrensie.Data as Data
//import PyFrensie.Data.Native as Native
#include "Data_ScatteringCenterPropertiesDatabase.hpp"

#include "MonteCarlo_WeightImportanceMesh.hpp"
#include "MonteCarlo_ParticleType.hpp"
#include "MonteCarlo_ScatteringCenterDefinitionDatabase.hpp"



//---------------------------------------------------------------------------//
// Testing Types
//---------------------------------------------------------------------------//

using boost::units::si::kelvin;
using boost::units::cgs::cubic_centimeter;
using Utility::Units::MeV;

//---------------------------------------------------------------------------//
// Testing Variables
//---------------------------------------------------------------------------//

std::string test_scattering_center_database_name;

std::shared_ptr<MonteCarlo::ScatteringCenterDefinitionDatabase>
scattering_center_definition_database;

std::shared_ptr<MonteCarlo::MaterialDefinitionDatabase>
material_definition_database;

std::shared_ptr<const Geometry::Model> unfilled_model;

std::shared_ptr<const MonteCarlo::ParticleDistribution> particle_distribution;

int threads;

std::shared_ptr<MonteCarlo::ParticleSimulationManager> global_manager;

//---------------------------------------------------------------------------//
// Testing functions
//---------------------------------------------------------------------------//
// void (*default_signal_handler)( int );

// extern "C" void custom_signal_handler( int signal )
// {
//   if( global_manager )
//     global_manager->signalHandler( signal );
// }

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Check that history details can be returned
FRENSIE_UNIT_TEST( ParticleSimulationManager, get_history_details )
{
  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::NEUTRON_PHOTON_MODE );
    properties->setNumberOfHistories( 5 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardNeutronSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
  }

  FRENSIE_CHECK_EQUAL( manager->getNextHistory(), 0 );
  FRENSIE_CHECK_EQUAL( manager->getNumberOfRendezvous(), 0 );
  FRENSIE_CHECK_EQUAL( manager->getRendezvousBatchSize(), 5 );
  FRENSIE_CHECK_EQUAL( manager->getBatchSize(), 5 );

  manager.reset();

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::NEUTRON_PHOTON_MODE );
    properties->setNumberOfHistories( 1005 );
    properties->setMinNumberOfRendezvous( 10 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardNeutronSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
  }

  FRENSIE_CHECK_EQUAL( manager->getNextHistory(), 0 );
  FRENSIE_CHECK_EQUAL( manager->getNumberOfRendezvous(), 0 );
  FRENSIE_CHECK_EQUAL( manager->getRendezvousBatchSize(), 100 );
  FRENSIE_CHECK_EQUAL( manager->getBatchSize(), 100 );

  manager.reset();

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::NEUTRON_PHOTON_MODE );
    properties->setNumberOfHistories( 1005 );
    properties->setMinNumberOfRendezvous( 10 );
    properties->setMaxRendezvousBatchSize( 50 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardNeutronSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
  }

  FRENSIE_CHECK_EQUAL( manager->getNextHistory(), 0 );
  FRENSIE_CHECK_EQUAL( manager->getNumberOfRendezvous(), 0 );
  FRENSIE_CHECK_EQUAL( manager->getRendezvousBatchSize(), 50 );
  FRENSIE_CHECK_EQUAL( manager->getBatchSize(), 50 );

  manager.reset();

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::NEUTRON_PHOTON_MODE );
    properties->setNumberOfHistories( 1000000 );
    properties->setMinNumberOfRendezvous( 5 );
    properties->setMaxRendezvousBatchSize( 100000 );
    properties->setMinNumberOfBatchesPerRendezvous( 10 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardNeutronSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
  }

  FRENSIE_CHECK_EQUAL( manager->getNextHistory(), 0 );
  FRENSIE_CHECK_EQUAL( manager->getNumberOfRendezvous(), 0 );
  FRENSIE_CHECK_EQUAL( manager->getRendezvousBatchSize(), 100000 );
  FRENSIE_CHECK_EQUAL( manager->getBatchSize(), 10000 );

  manager.reset();

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::NEUTRON_PHOTON_MODE );
    properties->setNumberOfHistories( 1000000 );
    properties->setMinNumberOfRendezvous( 5 );
    properties->setMaxRendezvousBatchSize( 100000 );
    properties->setMinNumberOfBatchesPerRendezvous( 10 );
    properties->setMaxBatchSize( 5000 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardNeutronSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
  }

  FRENSIE_CHECK_EQUAL( manager->getNextHistory(), 0 );
  FRENSIE_CHECK_EQUAL( manager->getNumberOfRendezvous(), 0 );
  FRENSIE_CHECK_EQUAL( manager->getRendezvousBatchSize(), 100000 );
  FRENSIE_CHECK_EQUAL( manager->getBatchSize(), 5000 );

  manager.reset();

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::NEUTRON_PHOTON_MODE );
    properties->setNumberOfHistories( 0 );
    properties->setMinNumberOfRendezvous( 5 );
    properties->setMaxRendezvousBatchSize( 1000000 );
    properties->setMinNumberOfBatchesPerRendezvous( 10 );
    properties->setMaxBatchSize( 50000 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardNeutronSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
  }

  FRENSIE_CHECK_EQUAL( manager->getNextHistory(), 0 );
  FRENSIE_CHECK_EQUAL( manager->getNumberOfRendezvous(), 0 );
  FRENSIE_CHECK_EQUAL( manager->getRendezvousBatchSize(), 1000000 );
  FRENSIE_CHECK_EQUAL( manager->getBatchSize(), 50000 );

  manager.reset();
}

//---------------------------------------------------------------------------//
// Check that a particle simulation manager can rename the simulation
FRENSIE_UNIT_TEST( ParticleSimulationManager, setSimulationName_default )
{
  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::NEUTRON_PHOTON_MODE );
    properties->setNumberOfHistories( 5 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardNeutronSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
  }

  FRENSIE_CHECK_EQUAL( manager->getSimulationName(), "test_sim" );

  manager->setSimulationName( "test_sim_2" );

  FRENSIE_CHECK_EQUAL( manager->getSimulationName(), "test_sim_2" );
  FRENSIE_CHECK( boost::filesystem::exists( "test_sim_2_rendezvous.xml" ) );
}

//---------------------------------------------------------------------------//
// Check that a particle simulation manager can rename the simulation
FRENSIE_UNIT_TEST( ParticleSimulationManager, setSimulationName_single_file )
{
  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::NEUTRON_PHOTON_MODE );
    properties->setNumberOfHistories( 5 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardNeutronSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
    manager->useSingleRendezvousFile();
  }

  FRENSIE_CHECK_EQUAL( manager->getSimulationName(), "test_sim" );

  manager->setSimulationName( "test_sim_2" );

  FRENSIE_CHECK_EQUAL( manager->getSimulationName(), "test_sim_2" );
  FRENSIE_CHECK( boost::filesystem::exists( "test_sim_2_rendezvous.xml" ) );
}

//---------------------------------------------------------------------------//
// Check that a particle simulation manager can rename the simulation
FRENSIE_UNIT_TEST( ParticleSimulationManager, setSimulationName_multiple_files )
{
  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::NEUTRON_PHOTON_MODE );
    properties->setNumberOfHistories( 5 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardNeutronSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
    manager->useMultipleRendezvousFiles();
  }

  FRENSIE_CHECK_EQUAL( manager->getSimulationName(), "test_sim" );

  manager->setSimulationName( "test_sim_2" );

  FRENSIE_CHECK_EQUAL( manager->getSimulationName(), "test_sim_2" );
  FRENSIE_CHECK( boost::filesystem::exists( "test_sim_2_rendezvous_0.xml" ) );
}

//---------------------------------------------------------------------------//
// Check that a particle simulation manager can change the archive
// type
FRENSIE_UNIT_TEST( ParticleSimulationManager, setSimulationArchiveType )
{
  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::NEUTRON_PHOTON_MODE );
    properties->setNumberOfHistories( 5 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardNeutronSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
    manager->useMultipleRendezvousFiles();
  }

  FRENSIE_CHECK_EQUAL( manager->getSimulationArchiveType(), "xml" );

  manager->setSimulationArchiveType( "txt" );

  FRENSIE_CHECK_EQUAL( manager->getSimulationArchiveType(), "txt" );
  FRENSIE_CHECK( boost::filesystem::exists( "test_sim_rendezvous_0.txt" ) );
}

//---------------------------------------------------------------------------//
// Check that the simulation name and archive can be changed simultaneously
FRENSIE_UNIT_TEST( ParticleSimulationManager, setSimulationNameAndArchiveType )
{
  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::NEUTRON_PHOTON_MODE );
    properties->setNumberOfHistories( 5 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardNeutronSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
    manager->useMultipleRendezvousFiles();
  }

  FRENSIE_CHECK_EQUAL( manager->getSimulationName(), "test_sim" );
  FRENSIE_CHECK_EQUAL( manager->getSimulationArchiveType(), "xml" );

  manager->setSimulationNameAndArchiveType( "test_sim_2", "txt" );

  FRENSIE_CHECK_EQUAL( manager->getSimulationName(), "test_sim_2" );
  FRENSIE_CHECK_EQUAL( manager->getSimulationArchiveType(), "txt" );
  FRENSIE_CHECK( boost::filesystem::exists( "test_sim_2_rendezvous_0.txt" ) );
}

//---------------------------------------------------------------------------//
// Check that the geometry model can be returned
FRENSIE_UNIT_TEST( ParticleSimulationManager, getModel )
{
  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager;

  std::shared_ptr<const MonteCarlo::FilledGeometryModel> model;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::NEUTRON_PHOTON_MODE );
    properties->setNumberOfHistories( 5 );

    model.reset( new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardNeutronSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
  }

  FRENSIE_CHECK( &manager->getModel() == model.get() );
}

//---------------------------------------------------------------------------//
// Check that the source can be returned
FRENSIE_UNIT_TEST( ParticleSimulationManager, getSource )
{
  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager;

  std::shared_ptr<MonteCarlo::ParticleSource> source;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::NEUTRON_PHOTON_MODE );
    properties->setNumberOfHistories( 5 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardNeutronSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
  }

  FRENSIE_CHECK( &manager->getSource() == source.get() );
}

//---------------------------------------------------------------------------//
// Check that the event handler can be returned
FRENSIE_UNIT_TEST( ParticleSimulationManager, getEventHandler )
{
  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager;

  std::shared_ptr<MonteCarlo::EventHandler> event_handler;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::NEUTRON_PHOTON_MODE );
    properties->setNumberOfHistories( 5 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardNeutronSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    event_handler.reset( new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
  }

  FRENSIE_CHECK( &manager->getEventHandler() == event_handler.get() );
}

//---------------------------------------------------------------------------//
// Check that a simulation can be run
FRENSIE_UNIT_TEST( ParticleSimulationManager, runSimulation_history_wall )
{
  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::PHOTON_MODE );
    properties->setNumberOfHistories( 5 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardPhotonSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
  }

  FRENSIE_REQUIRE_NO_THROW( manager->runSimulation() );

  FRENSIE_CHECK_EQUAL( manager->getNextHistory(), 5 );
  FRENSIE_CHECK_EQUAL( manager->getNumberOfRendezvous(), 2 );
}

//---------------------------------------------------------------------------//
// Check that a simulation can be run
FRENSIE_UNIT_TEST( ParticleSimulationManager, runSimulation_wall_time )
{
  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::PHOTON_MODE );
    properties->setSimulationWallTime( 0.5 );
    properties->setMaxRendezvousBatchSize( 10 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardElectronSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
  }

  FRENSIE_REQUIRE_NO_THROW( manager->runSimulation() );

  FRENSIE_CHECK( manager->getNextHistory() > 0 );
  FRENSIE_CHECK( manager->getNumberOfRendezvous() > 0 );
}

//---------------------------------------------------------------------------//
// Check that a particle simulation manager can handle a signal
#ifdef HAVE_FRENSIE_OPENMP
FRENSIE_UNIT_TEST( ParticleSimulationManager, runInterruptibleSimulation )
{
  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::PHOTON_MODE );
    properties->setMaxRendezvousBatchSize( 100 );
    properties->setMaxBatchSize( 10 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardNeutronSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
  }

  #pragma omp parallel num_threads( 2 )
  {
    if( Utility::OpenMPProperties::getThreadId() == 0 )
      manager->runInterruptibleSimulation();
    else
    {
      std::shared_ptr<Utility::Timer> timer =
        Utility::OpenMPProperties::createTimer();

      timer->start();

      while( timer->elapsed().count() < 0.2 );

      timer->stop();
      timer.reset();

      // Terminate the simulation (it is set up to run indefinitely unless it
      // receives an interput signal)
      std::raise( SIGINT );
    }
  }

  FRENSIE_CHECK( manager->getNextHistory() > 0 );
  FRENSIE_CHECK( manager->getNumberOfRendezvous() > 0 );
}
#endif // end HAVE_FRENSIE_OPENMP

//---------------------------------------------------------------------------//
// Check that a particle simulation summary can be printed
FRENSIE_UNIT_TEST( ParticleSimulationManager, printSimulationSummary )
{
  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::PHOTON_MODE );
    properties->setNumberOfHistories( 5 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardPhotonSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
  }

  manager->runSimulation();

  FRENSIE_REQUIRE_NO_THROW( manager->printSimulationSummary( std::cout ) );
}

//---------------------------------------------------------------------------//
// Check that a particle simulation summary can be logged
FRENSIE_UNIT_TEST( ParticleSimulationManager, logSimulationSummary )
{
  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::PHOTON_MODE );
    properties->setNumberOfHistories( 5 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardPhotonSourceComponent(
                                                     0,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

    factory.reset(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              "xml",
                                                              threads ) );

    manager = factory->getManager();
  }

  manager->runSimulation();

  FRENSIE_REQUIRE_NO_THROW( manager->logSimulationSummary() );
}

//---------------------------------------------------------------------------//
// Check that a particle simulation can be restarted
FRENSIE_DATA_UNIT_TEST_DECL( ParticleSimulationManager, restart_basic )
{
  FETCH_FROM_TABLE( std::string, archive_type );
  FETCH_FROM_TABLE( uint32_t, source_id );

  uint64_t next_history;
  uint64_t rendezvous_number;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::PHOTON_MODE );
    properties->setSimulationWallTime( 0.25 );
    properties->setMaxRendezvousBatchSize( 10 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardPhotonSourceComponent(
                                                     source_id,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              archive_type,
                                                              threads ) );

    std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager =
      factory->getManager();
    manager->useMultipleRendezvousFiles();

    FRENSIE_REQUIRE_NO_THROW( manager->runSimulation() );

    next_history = manager->getNextHistory();
    rendezvous_number = manager->getNumberOfRendezvous();
  }

  std::string archive_name( "test_sim_rendezvous_" );
  archive_name += Utility::toString( rendezvous_number - 1 );
  archive_name += ".";
  archive_name += archive_type;

  std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

  FRENSIE_REQUIRE_NO_THROW( factory.reset( new MonteCarlo::ParticleSimulationManagerFactory( archive_name, (unsigned)threads ) ) );

  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager =
    factory->getManager();

  FRENSIE_REQUIRE_NO_THROW( manager->runSimulation() );

  FRENSIE_CHECK( manager->getNextHistory() > next_history );
  FRENSIE_CHECK( manager->getNumberOfRendezvous() > rendezvous_number );
}

FRENSIE_DATA_UNIT_TEST_INST( ParticleSimulationManager, restart_basic )
{
  COLUMNS()         << "archive_type" << "source_id" ;
  NEW_ROW( "xml" )  <<    "xml"       <<    0;
  NEW_ROW( "txt" )  <<    "txt"       <<    1;
  NEW_ROW( "bin" )  <<    "bin"       <<    2;
#ifdef HAVE_FRENSIE_HDF5
  NEW_ROW( "h5fa" ) <<    "h5fa"      <<    3;
#endif
}

//---------------------------------------------------------------------------//
// Check that a particle simulation manager can be restarted
FRENSIE_DATA_UNIT_TEST_DECL( ParticleSimulationManager, restart_add_histories )
{
  FETCH_FROM_TABLE( std::string, archive_type );
  FETCH_FROM_TABLE( uint32_t, source_id );

  uint64_t next_history;
  uint64_t rendezvous_number;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::PHOTON_MODE );
    properties->setSimulationWallTime( 0.25 );
    properties->setMaxRendezvousBatchSize( 10 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardPhotonSourceComponent(
                                                     source_id,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              archive_type,
                                                              threads ) );

    std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager =
      factory->getManager();
    manager->useSingleRendezvousFile();

    FRENSIE_REQUIRE_NO_THROW( manager->runSimulation() );

    next_history = manager->getNextHistory();
    rendezvous_number = manager->getNumberOfRendezvous();
  }

  std::string archive_name( "test_sim_rendezvous." );
  archive_name += archive_type;

  std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

  FRENSIE_REQUIRE_NO_THROW( factory.reset( new MonteCarlo::ParticleSimulationManagerFactory( archive_name, (uint64_t)5, (unsigned)threads ) ) );

  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager =
    factory->getManager();

  FRENSIE_REQUIRE_NO_THROW( manager->runSimulation() );

  FRENSIE_CHECK_EQUAL( manager->getNextHistory(), next_history+5 );
  FRENSIE_CHECK( manager->getNumberOfRendezvous() > rendezvous_number );
}

FRENSIE_DATA_UNIT_TEST_INST( ParticleSimulationManager, restart_add_histories )
{
  COLUMNS()         << "archive_type" << "source_id" ;
  NEW_ROW( "xml" )  <<    "xml"       <<    0;
  NEW_ROW( "txt" )  <<    "txt"       <<    1;
  NEW_ROW( "bin" )  <<    "bin"       <<    2;
#ifdef HAVE_FRENSIE_HDF5
  NEW_ROW( "h5fa" ) <<    "h5fa"      <<    3;
#endif
}

//---------------------------------------------------------------------------//
// Check that a particle simulation manager can be restarted
FRENSIE_DATA_UNIT_TEST_DECL( ParticleSimulationManager, restart_new_wall_time )
{
  FETCH_FROM_TABLE( std::string, archive_type );
  FETCH_FROM_TABLE( uint32_t, source_id );

  uint64_t next_history;
  uint64_t rendezvous_number;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::PHOTON_MODE );
    properties->setSimulationWallTime( 0.25 );
    properties->setMaxRendezvousBatchSize( 10 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardPhotonSourceComponent(
                                                     source_id,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              archive_type,
                                                              threads ) );

    std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager =
      factory->getManager();
    manager->useMultipleRendezvousFiles();

    FRENSIE_REQUIRE_NO_THROW( manager->runSimulation() );

    next_history = manager->getNextHistory();
    rendezvous_number = manager->getNumberOfRendezvous();
  }

  std::string archive_name( "test_sim_rendezvous_" );
  archive_name += Utility::toString( rendezvous_number - 1 );
  archive_name += ".";
  archive_name += archive_type;

  std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

  FRENSIE_REQUIRE_NO_THROW( factory.reset( new MonteCarlo::ParticleSimulationManagerFactory( archive_name, 0.1, (unsigned)threads ) ) );

  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager =
    factory->getManager();

  FRENSIE_REQUIRE_NO_THROW( manager->runSimulation() );

  FRENSIE_CHECK( manager->getNextHistory() > next_history );
  FRENSIE_CHECK( manager->getNumberOfRendezvous() > rendezvous_number );
}

FRENSIE_DATA_UNIT_TEST_INST( ParticleSimulationManager, restart_new_wall_time )
{
  COLUMNS()         << "archive_type" << "source_id" ;
  NEW_ROW( "xml" )  <<    "xml"       <<    0;
  NEW_ROW( "txt" )  <<    "txt"       <<    1;
  NEW_ROW( "bin" )  <<    "bin"       <<    2;
#ifdef HAVE_FRENSIE_HDF5
  NEW_ROW( "h5fa" ) <<    "h5fa"      <<    3;
#endif
}

//---------------------------------------------------------------------------//
// Check that a particle simulation manager can be restarted
FRENSIE_DATA_UNIT_TEST_DECL( ParticleSimulationManager,
                             restart_add_histories_new_wall_time )
{
  FETCH_FROM_TABLE( std::string, archive_type );
  FETCH_FROM_TABLE( uint32_t, source_id );

  uint64_t next_history;
  uint64_t rendezvous_number;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::PHOTON_MODE );
    properties->setSimulationWallTime( 0.25 );
    properties->setMaxRendezvousBatchSize( 10 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardPhotonSourceComponent(
                                                     source_id,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              archive_type,
                                                              threads ) );

    std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager =
      factory->getManager();
    manager->useMultipleRendezvousFiles();

    FRENSIE_REQUIRE_NO_THROW( manager->runSimulation() );

    next_history = manager->getNextHistory();
    rendezvous_number = manager->getNumberOfRendezvous();
  }

  std::string archive_name( "test_sim_rendezvous_" );
  archive_name += Utility::toString( rendezvous_number - 1 );
  archive_name += ".";
  archive_name += archive_type;

  std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

  FRENSIE_REQUIRE_NO_THROW( factory.reset( new MonteCarlo::ParticleSimulationManagerFactory( archive_name, (uint64_t)5, 0.1, (unsigned)threads ) ) );

  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager =
    factory->getManager();

  FRENSIE_REQUIRE_NO_THROW( manager->runSimulation() );

  FRENSIE_CHECK_EQUAL( manager->getNextHistory(), next_history+5 );
  FRENSIE_CHECK( manager->getNumberOfRendezvous() > rendezvous_number );
}

FRENSIE_DATA_UNIT_TEST_INST( ParticleSimulationManager,
                             restart_add_histories_new_wall_time )
{
  COLUMNS()         << "archive_type" << "source_id" ;
  NEW_ROW( "xml" )  <<    "xml"       <<    0;
  NEW_ROW( "txt" )  <<    "txt"       <<    1;
  NEW_ROW( "bin" )  <<    "bin"       <<    2;
#ifdef HAVE_FRENSIE_HDF5
  NEW_ROW( "h5fa" ) <<    "h5fa"      <<    3;
#endif
}

//---------------------------------------------------------------------------//
// Check that a particle simulation manager can be restarted
FRENSIE_DATA_UNIT_TEST_DECL( ParticleSimulationManager,
                             restart_updated_props )
{
  FETCH_FROM_TABLE( std::string, archive_type );
  FETCH_FROM_TABLE( uint32_t, source_id );

  uint64_t next_history;
  uint64_t rendezvous_number;

  {
    std::shared_ptr<MonteCarlo::SimulationProperties> properties(
                                        new MonteCarlo::SimulationProperties );
    properties->setParticleMode( MonteCarlo::PHOTON_MODE );
    properties->setSimulationWallTime( 0.25 );
    properties->setMaxRendezvousBatchSize( 10 );

    std::shared_ptr<const MonteCarlo::FilledGeometryModel> model(
                               new MonteCarlo::FilledGeometryModel(
                                        test_scattering_center_database_name,
                                        scattering_center_definition_database,
                                        material_definition_database,
                                        properties,
                                        unfilled_model,
                                        false ) );

    std::shared_ptr<MonteCarlo::ParticleSource> source;

    {
      std::shared_ptr<MonteCarlo::ParticleSourceComponent>
        source_component( new MonteCarlo::StandardPhotonSourceComponent(
                                                     source_id,
                                                     1.0,
                                                     unfilled_model,
                                                     particle_distribution ) );

      source.reset( new MonteCarlo::StandardParticleSource( {source_component} ) );
    }

    std::shared_ptr<MonteCarlo::EventHandler> event_handler(
                                 new MonteCarlo::EventHandler( *properties ) );

    std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory(
            new MonteCarlo::ParticleSimulationManagerFactory( model,
                                                              source,
                                                              event_handler,
                                                              properties,
                                                              "test_sim",
                                                              archive_type,
                                                              threads ) );

    std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager =
      factory->getManager();
    manager->useSingleRendezvousFile();

    FRENSIE_REQUIRE_NO_THROW( manager->runSimulation() );

    next_history = manager->getNextHistory();
    rendezvous_number = manager->getNumberOfRendezvous();
  }

  std::string archive_name( "test_sim_rendezvous." );
  archive_name += archive_type;

  MonteCarlo::SimulationGeneralProperties updated_properties;
  updated_properties.setNumberOfHistories( 16 );
  updated_properties.setMinNumberOfRendezvous( 2 );
  updated_properties.setMaxRendezvousBatchSize( 100 );
  updated_properties.setMinNumberOfBatchesPerRendezvous( 2 );
  updated_properties.setMaxBatchSize( 10 );
  updated_properties.setNumberOfSnapshotsPerBatch( 3 );
  updated_properties.setSimulationWallTime( 1.0 );

  std::unique_ptr<MonteCarlo::ParticleSimulationManagerFactory> factory;

  FRENSIE_REQUIRE_NO_THROW( factory.reset( new MonteCarlo::ParticleSimulationManagerFactory( archive_name, updated_properties, (unsigned)threads ) ) );

  std::shared_ptr<MonteCarlo::ParticleSimulationManager> manager =
    factory->getManager();

  FRENSIE_REQUIRE_NO_THROW( manager->runSimulation() );

  FRENSIE_CHECK_EQUAL( manager->getNextHistory(), next_history+16 );
  FRENSIE_CHECK( manager->getNumberOfRendezvous() > rendezvous_number );
  FRENSIE_CHECK_EQUAL( manager->getRendezvousBatchSize(), 8 );
  FRENSIE_CHECK_EQUAL( manager->getBatchSize(), 4 );
}

FRENSIE_DATA_UNIT_TEST_INST( ParticleSimulationManager, restart_updated_props )
{
  COLUMNS()         << "archive_type" << "source_id" ;
  NEW_ROW( "xml" )  <<    "xml"       <<    0;
  NEW_ROW( "txt" )  <<    "txt"       <<    1;
  NEW_ROW( "bin" )  <<    "bin"       <<    2;
#ifdef HAVE_FRENSIE_HDF5
  NEW_ROW( "h5fa" ) <<    "h5fa"      <<    3;
#endif
}

//---------------------------------------------------------------------------//
// Check that a 
FRENSIE_UNIT_TEST( ParticleSimulationManager, i_choose_thename )
{

  boost::filesystem::path database_path = "";
  // Load the database
  const Data::ScatteringCenterPropertiesDatabase database( database_path );

  double num_particles = 1e2;

  
  std::shared_ptr<Geometry::DagMCModel>   forward_model;
// Forward model (must be seperate due to estimator-geometry simultaneous declaration)
  Geometry::DagMCModelProperties* forward_model_properties = new Geometry::DagMCModelProperties( "test_dagmc_geom_file_name" );
  forward_model_properties->setMaterialPropertyName("material");
  forward_model_properties->setDensityPropertyName("density");
  forward_model_properties->setTerminationCellPropertyName("termination->cell");
  forward_model_properties->setCellCollisionFluxName("cell->c->flux");
  forward_model_properties->useFastIdLookup();

  forward_model.reset( new Geometry::DagMCModel( *forward_model_properties));


// Needs to be the same for all 3 meshes to avoid split/terminate on birth
  double mesh_increment = 50.0;
  double x_distance = 5000;
  double y_distance = 1000;
  double z_distance = 1000;
//Form the mesh for the entire geometry
  double x0 = -2500.0;
  std::vector<double> x_planes;
  for (int i=0; i <  (int)((x_distance/mesh_increment) + 1); i++ ) {
    x_planes.push_back(i * mesh_increment + x0);
  }
  double y0 = -500.0;
  double z0 = -500.0;
  std::vector<double> y_planes;
  std::vector<double> z_planes;
  for (int i=0; i <  (int)((y_distance/mesh_increment) + 1); i++ ) {
    y_planes.push_back(i * mesh_increment + y0);
    z_planes.push_back(i * mesh_increment + z0);
  }

// entire geometry mesh
  using StructuredHexMesh::StructuredHexMesh;
  std::shared_ptr<StructuredHexMesh> geometry_mesh =
      std::make_shared<StructuredHexMesh>(x_planes, y_planes,
                                                             z_planes);

  // Form the mesh for the source and response function
  double x0_src = -2350;
  double x0_resp = 2300;

  double x_planes_src[2] = [ x0_src, x0_src + 50 ];
  double x_planes_resp[2] = [ x0_resp, x0_resp + 50 ];

  y0 = -300;
  z0 = -300;
  y_planes.clear();
  z_planes.clear();
  for (i=0; i<13; i++){
    y_planes.push_back( i*mesh_increment + y0 );
    z_planes.push_back( i*mesh_increment + z0 );
  }

// source/detector mesh objects
  std::shared_ptr<StructuredHexMesh> source_mesh =
      std::make_shared<StructuredHexMesh>(x_planes_src, y_planes, z_planes);
  std::shared_ptr<StructuredHexMesh> response_mesh =
      std::make_shared<StructuredHexMesh>(x_planes_resp, y_planes, z_planes);
  size_t number_of_mesh_elements = source_mesh.getNumberOfElements();

// Direction Discretization
  int PQLA_quadrature_order = 2;
  using Utility::PQLAQuadrature;
  std::shared_pr<PQLAQuadrature> direction_discretization =
      std::make_shared<PQLAQuadrature>(PQLA_quadrature_order);
  size_t = number_of_direction_elements =
      direction_discretization.getNumberOfTriangles();

// Non-importance sampled forward source (for initial forward run)
        // Raw energy histogram in MeV
  std::vector<double> raw_forward_energy_distribution_bounds = {
      7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0};
  std::vector<double> raw_forward_energy_distribution_values = {1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.};
//Use later for source
  using Utility::HistogramDistribution;
  std::shared_ptr<HistogramDistribution>
      raw_forward_energy_distribution = std::make_shared<HistogramDistribution>(
          raw_forward_energy_distribution_bounds,
          raw_forward_energy_distribution_values);

  // Raw mesh distribution
  std::vector<double> raw_forward_mesh_distribution_bounds;
  for (int i= 0; i < number_of_mesh_elements; i++ )
    raw_forward_mesh_distribution_bounds.push_back(i);
  std::vector<double> raw_forward_mesh_distribution_values;
  for (int i= 0; i < number_of_mesh_elements-1; i++ )
    raw_forward_mesh_distribution_values.push_back(1.0);

//Use later for source
  std::shared_ptr<HistogramDistribution>
      raw_forward_mesh_distribution = std::make_shared<HistogramDistribution>(
          raw_forward_mesh_distribution_bounds,
          raw_forward_mesh_distribution_values);

// Form INITIAL particle distributions for source
  using MonteCarlo::IndependentPhaseSpaceDimensionDistribution;
  std::shared_ptr<IndependentPhaseSpaceDimensionDistribution>
      forward_energy_distribution =
          std::make_shared<IndependentEnergyDimensionDistribution>(
              raw_forward_energy_distribution);
  std::shared_ptr<IndependentPhaseSpaceDimensionDistribution>
    forward_mesh_distribution =
        std::make_shared<IndependentEnergyDimensionDistribution>(
            raw_forward_mesh_distribution);

  using MonteCarlo::StandardParticleDistribution;
  std::shared_ptr<StandardParticleDistribution>
      forward_particle_distribution =
          std::make_shared<StandardParticleDistribution>(
              "Initial forward source distribution");
  forward_particle_distribution->setDimensionDistribution(
      forward_energy_distribution);
  forward_particle_distribution->setMeshIndexDimensionDistribution(
      forward_mesh_distribution, source_mesh);
  forward_particle_distribution
      ->constructDimensionDistributionDependencyTree();

// Set up source
  using MonteCarlo::StandardParticleSource;
  using MonteCarlo::StandardPhotonSourceComponent;

  std::shared_ptr<StandardPhotonSourceComponent> forward_particle_source_component;

  forward_particle_source_component = std::make_shared<StandardPhotonSourceComponent>(
      1, 1, forward_model, forward_particle_distribution);
  std::shared_ptr<StandardParticleSource> forward_source = std::make_shared<
      StandardParticleSource>(forward_particle_source_component);

  // initialize forward weight-importance hybrids (all are equal to 1)
  using MonteCarlo::WeightImportanceMesh;

  std::shared_ptr<WeightImportanceMesh> forward_weight_importance_mesh =
      std::make_shared<WeightImportanceMesh>();
  forward_weight_importance_mesh->setMesh(geometry_mesh);
  forward_weight_importance_mesh->setDirectionDiscretization(
      Event.ObserverDirectionDimensionDiscretization.PQLA, 2, True);
  std::vector<double> geometry_mesh_observer_energy_discretization = {
      0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
      5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0};

  std::vector<double> forward_detector_energy_discretization = {0.0, 0.5, 1.0,
                                                                1.5};
  forward_weight_importance_mesh->setEnergyDiscretization(
      geometry_mesh_observer_energy_discretization);

  using Utility::Mesh::ElementHandle;
  std::unordered_map<ElementHandle, std::vector<double>> weight_importance_dictionary;
  for (int i =0; i < geometry_mesh.getNumberOfElements(); i++) {
    std::vector<double> local_weight_importance_vector;
    for (int j = 0; j < 32; j++) { //direction_index
          for (int u = 0; u < geometry_mesh_observer_energy_discretization.size()-1; u++) { //energy_discretization_index
            local_weight_importance_vector.push_back(1);
          }
    }
    weight_importance_dictionary[i] = local_weight_importance_vector;
  }
  forward_weight_importance_mesh.setWeightImportanceMap( weight_importance_dictionary );

  
        
// Set the simulation properties
  using MonteCarlo::SimulationProperties;
  share_ptr<SimulationProperties> simulation_properties =
      std::make_shared<SimulationProperties>();

  // Simulate photons only
  simulation_properties->setParticleMode(MonteCarlo::ParticleType::PHOTON);

// Set the number of histories to run and the number of rendezvous
  simulation_properties->setNumberOfHistories(num_particles);
  simulation_properties->setMinNumberOfRendezvous(10);
  simulation_properties->setNumberOfSnapshotsPerBatch(10);
  simulation_properties->setNumberOfPhotonHashGridBins(100);

//Set up the materials
  std::shared_ptr<Data::ScatteringCenterPropertiesDatabase> database =
      std::make_shared<Data::ScatteringCenterPropertiesDatabase>(db_path);

  // Extract the properties for H from the database
  Data::AtomProperties& H_properties = database->getAtomProperties(Data.ZAID(1000));

  // Extract the properties for Pb from the database
  Data::AtomProperties& Pb_properties = database->getAtomProperties(Data.ZAID(82000));

// Extract the properties for K from the database
  Data::AtomProperties& K_properties = database->getAtomProperties(Data.ZAID(19000));

// Extract the properties for Ge from the database
  Data::AtomProperties& Ge_properties = database->getAtomProperties(Data.ZAID(32000));

// Set the definition for H, Pb, K, Ge for this simulation
  std::shared_ptr<MonteCarlo::ScatteringCenterDefinitionDatabase>
      scattering_center_definitions =
          std::make_shared<MonteCarlo::ScatteringCenterDefinitionDatabase>();
  MonteCarlo::ScatteringCenterDefinition& H_definition =
      scattering_center_definitions->createDefinition("H", Data.ZAID(1000));
  MonteCarlo::ScatteringCenterDefinition& Pb_definition =
      scattering_center_definitions->createDefinition("Pb", Data.ZAID(82000));
   MonteCarlo::ScatteringCenterDefinition& K_definition =
      scattering_center_definitions->createDefinition("K", Data.ZAID(19000));
   MonteCarlo::ScatteringCenterDefinition& Ge_definition =
      scattering_center_definitions->createDefinition("Ge", Data.ZAID(32000));

   auto data_file_type = Data::PhotoatomicDataProperties::Native_EPR_FILE;
   unsigned file_version = 0;

   H_definition.setPhotoatomicDataProperties(
       H_properties.getSharedPhotoatomicDataProperties(data_file_type,
                                                       file_version));

   Pb_definition.setPhotoatomicDataProperties(
       Pb_properties.getSharedPhotoatomicDataProperties(data_file_type,
                                                        file_version));

   K_definition.setPhotoatomicDataProperties(
       K_properties.getSharedPhotoatomicDataProperties(data_file_type,
                                                       file_version));

   Ge_definition.setPhotoatomicDataProperties(
       Ge_properties.getSharedPhotoatomicDataProperties(data_file_type,
                                                        file_version));

   // Set the definition for materials
   std::shared_ptr<MonteCarlo::MaterialDefinitionDatabase> material_definitions = std::make_shared<MaterialDefinitionDatabase>();
   material_definitions->addDefinition("H", 1, {"H"}, {1.0});
   material_definitions->addDefinition("Pb", 2, {"Pb"}, {1.0});
   material_definitions->addDefinition("K", 3, {"K"}, {1.0});
   material_definitions->addDefinition("Ge", 4, {"Ge"}, {1.0});

  std::shared_ptr<const MonteCarlo::FilledGeometryModel> filled_model;

  filled_model.reset(new MonteCarlo::FilledGeometryModel(
      db_path, scattering_center_definitions, material_definitions,
      simulation_properties, model, True));

  // Set up the event handler
  MonteCarlo::EventHandler event_handler = Event.EventHandler(model, simulation_properties);


  // Detector Collision Estimator (main estimator, not VR producing estimator)
  event_handler.getEstimator(1).setDiscretization<MonteCarlo::ENERGY_DIMENSION>(
      forward_detector_energy_discretization);

  // Mesh estimator (weight importance mesh producing estimator)
  std::shared_ptr<MonteCarlo::WeightMultipliedMeshTrackLengthFluxEstimator>
      forward_mesh_estimator = std::make_shared<
          MonteCarlo::WeightMultipliedMeshTrackLengthFluxEstimator>(
          2, 1.0, geometry_mesh);

  forward_mesh_estimator->setDirectionDiscretization(
      Event.ObserverDirectionDimensionDiscretization.PQLA, 2, False);
  forward_mesh_estimator->setEnergyDiscretization(
      geometry_mesh_observer_energy_discretization);
  forward_mesh_estimator->setParticleTypes([MonteCarlo.PHOTON]);
  event_handler.addEstimator(forward_mesh_estimator);

  // Detector mesh estimator (for adjoint source biasing)
  std::shared_ptr<MonteCarlo::WeightMultipliedMeshTrackLengthFluxEstimator>
    forward_mesh_detector_estimator = std::make_shared<
        MonteCarlo::WeightMultipliedMeshTrackLengthFluxEstimator>(
        3, 1.0, detector_mesh);

  forward_mesh_detector_estimator->setDirectionDiscretization(
      Event.ObserverDirectionDimensionDiscretization.PQLA, 2, False);
  forward_mesh_detector_estimator->setEnergyDiscretization(
      forward_detector_energy_discretization);
  forward_mesh_detector_estimator->setParticleTypes([MonteCarlo.PHOTON]);
  event_handler.addEstimator(forward_mesh_detector_estimator);


  // Set up the simulation manager
  std::shared_ptr<MonteCarlo::ParticleSimulationManagerFactory>factory = std::make_shared<MonteCarlo::ParticleSimulationManagerFactory>(
      filled_model, forward_source, event_handler, simulation_properties,
      sim_name + "_forward", "xml", threads);

// Create the simulation manager
//factory.setPopulationControl(forward_weight_importance_mesh);
  manager = factory.getManager();
  manager.useSingleRendezvousFile();

//Run the simulation manager.runInterruptibleSimulation();
}

//---------------------------------------------------------------------------//
// Custom setup
//---------------------------------------------------------------------------//
FRENSIE_CUSTOM_UNIT_TEST_SETUP_BEGIN();

FRENSIE_CUSTOM_UNIT_TEST_COMMAND_LINE_OPTIONS()
{
  ADD_STANDARD_OPTION_AND_ASSIGN_VALUE( "test_database",
                                        test_scattering_center_database_name, "",
                                        "Test scattering center database name "
                                        "with path" );
  ADD_STANDARD_OPTION_AND_ASSIGN_VALUE( "threads",
                                        threads, 1,
                                        "Number of threads to use" );
}

FRENSIE_CUSTOM_UNIT_TEST_INIT()
{
  {
    // Determine the database directory
    boost::filesystem::path database_path =
      test_scattering_center_database_name;

    // Load the database
    const Data::ScatteringCenterPropertiesDatabase database( database_path );

    const Data::AtomProperties& h_properties =
      database.getAtomProperties( 1001 );

    const Data::NuclideProperties& h1_properties =
      database.getNuclideProperties( 1001 );

    // Set the sattering center definitions
    scattering_center_definition_database.reset(
                          new MonteCarlo::ScatteringCenterDefinitionDatabase );

    MonteCarlo::ScatteringCenterDefinition& h_definition =
      scattering_center_definition_database->createDefinition( "H1 @ 293.6K", 1001 );

    h_definition.setPhotoatomicDataProperties(
          h_properties.getSharedPhotoatomicDataProperties(
                       Data::PhotoatomicDataProperties::Native_EPR_FILE, 0 ) );

    h_definition.setAdjointPhotoatomicDataProperties(
          h_properties.getSharedAdjointPhotoatomicDataProperties(
                Data::AdjointPhotoatomicDataProperties::Native_EPR_FILE, 0 ) );

    h_definition.setElectroatomicDataProperties(
          h_properties.getSharedElectroatomicDataProperties(
                     Data::ElectroatomicDataProperties::Native_EPR_FILE, 0 ) );

    h_definition.setAdjointElectroatomicDataProperties(
          h_properties.getSharedAdjointElectroatomicDataProperties(
              Data::AdjointElectroatomicDataProperties::Native_EPR_FILE, 0 ) );

    h_definition.setNuclearDataProperties(
          h1_properties.getSharedNuclearDataProperties(
                                         Data::NuclearDataProperties::ACE_FILE,
                                         7,
                                         2.53010E-08*MeV,
                                         true ) );

    material_definition_database.reset(
                                  new MonteCarlo::MaterialDefinitionDatabase );

    material_definition_database->addDefinition( "H1 @ 293.6K", 1,
                                                 {"H1 @ 293.6K"}, {1.0} );
  }

  unfilled_model.reset(
            new Geometry::InfiniteMediumModel( 1, 1, -1.0/cubic_centimeter ) );

  {
    std::shared_ptr<MonteCarlo::StandardParticleDistribution>
      tmp_particle_distribution( new MonteCarlo::StandardParticleDistribution( "test dist" ) );

    particle_distribution = tmp_particle_distribution;
  }
}

FRENSIE_CUSTOM_UNIT_TEST_SETUP_END();

//---------------------------------------------------------------------------//
// end tstParticleSimulationManager.cpp
//---------------------------------------------------------------------------//
