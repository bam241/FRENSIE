//---------------------------------------------------------------------------//
//!
//! \file   tstACEElectroatomicDataProperties.cpp
//! \author Alex Robinson
//! \brief  ACEElectroatomicDataProperties class unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <string>
#include <memory>
#include <iostream>

// FRENSIE Includes
#include "Data_ACEElectroatomicDataProperties.hpp"
#include "Utility_UnitTestHarnessWithMain.hpp"
#include "ArchiveTestHelpers.hpp"

//---------------------------------------------------------------------------//
// Testing Types
//---------------------------------------------------------------------------//

using Utility::Units::amu;

typedef TestArchiveHelper::TestArchives TestArchives;

//---------------------------------------------------------------------------//
// Testing Variables
//---------------------------------------------------------------------------//
std::unique_ptr<const Data::ElectroatomicDataProperties> properties;

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Check that the atom can be returned
FRENSIE_UNIT_TEST( ACEElectroatomicDataProperties, atom )
{
  FRENSIE_CHECK_EQUAL( properties->atom(), Data::H_ATOM );
}

//---------------------------------------------------------------------------//
// Check that the file type can be returned
FRENSIE_UNIT_TEST( ACEElectroatomicDataProperties, fileType )
{
  FRENSIE_CHECK_EQUAL( properties->fileType(),
                       Data::ElectroatomicDataProperties::ACE_EPR_FILE );

  Data::ACEElectroatomicDataProperties local_properties( 1.0*amu,
                                                         "electroatomic_data/h_data.txt",
                                                         10,
                                                         "1000.03e" );

  FRENSIE_CHECK_EQUAL( local_properties.fileType(),
                       Data::ElectroatomicDataProperties::ACE_FILE );
}

//---------------------------------------------------------------------------//
// Check that the atomic weight can be returned
FRENSIE_UNIT_TEST( ACEElectroatomicDataProperties, atomicWeight )
{
  FRENSIE_CHECK_EQUAL( properties->atomicWeight(), 1.0*amu );
}

//---------------------------------------------------------------------------//
// Check that the file path can be returned
FRENSIE_UNIT_TEST( ACEElectroatomicDataProperties, filePath )
{
  FRENSIE_CHECK_EQUAL( properties->filePath().string(),
                       "electroatomic_data/h_data.txt" );
}

//---------------------------------------------------------------------------//
// Check that the file start line can be returned
FRENSIE_UNIT_TEST( ACEElectroatomicDataProperties, fileStartLine )
{
  FRENSIE_CHECK_EQUAL( properties->fileStartLine(), 10 );
}

//---------------------------------------------------------------------------//
// Check that the file version can be returned
FRENSIE_UNIT_TEST( ACEElectroatomicDataProperties, fileVersion )
{
  FRENSIE_CHECK_EQUAL( properties->fileVersion(), 12 );
}

//---------------------------------------------------------------------------//
// Check that the table name can be returned
FRENSIE_UNIT_TEST( ACEElectroatomicDataProperties, tableName )
{
  FRENSIE_CHECK_EQUAL( properties->tableName(), "1000.12p" );
}

//---------------------------------------------------------------------------//
// Check that the properties can be cloned
FRENSIE_UNIT_TEST( ACEElectroatomicDataProperties, clone )
{
  std::unique_ptr<const Data::ElectroatomicDataProperties>
    properties_clone( properties->clone() );

  FRENSIE_REQUIRE( properties_clone.get() != NULL );
  FRENSIE_CHECK( properties_clone.get() != properties.get() );
}

//---------------------------------------------------------------------------//
// Check that the properties can be archived
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( ACEElectroatomicDataProperties,
                                   archive,
                                   TestArchives )
{
  FETCH_TEMPLATE_PARAM( 0, RawOArchive );
  FETCH_TEMPLATE_PARAM( 1, RawIArchive );

  typedef typename std::remove_pointer<RawOArchive>::type OArchive;
  typedef typename std::remove_pointer<RawIArchive>::type IArchive;

  std::string archive_base_name( "test_ace_electroatomic_data_properties" );
  std::ostringstream archive_ostream;

  // Create and archive some properties
  {
    std::unique_ptr<OArchive> oarchive;

    createOArchive( archive_base_name, archive_ostream, oarchive );

    Data::ACEElectroatomicDataProperties local_properties( 2.0*amu,
                                                           "electroatomic_data/he_data.txt",
                                                           2,
                                                           "2000.03e" );

    std::shared_ptr<const Data::ElectroatomicDataProperties>
      shared_properties( properties->clone() );

    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP( local_properties ) );
    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP( shared_properties ) );
  }

  // Copy the archive ostream to an istream
  std::istringstream archive_istream( archive_ostream.str() );

  // Load the archived distributions
  std::unique_ptr<IArchive> iarchive;

  createIArchive( archive_istream, iarchive );

  Data::ACEElectroatomicDataProperties
    local_properties( 10.0*amu, "dummy", 100000, "1000.00e" );

  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP( local_properties ) );
  FRENSIE_CHECK_EQUAL( local_properties.atom(), Data::He_ATOM );
  FRENSIE_CHECK_EQUAL( local_properties.atomicNumber(), 2 );
  FRENSIE_CHECK_EQUAL( local_properties.atomicWeight(), 2.0*amu );
  FRENSIE_CHECK_EQUAL( local_properties.filePath().string(),
                       "electroatomic_data/he_data.txt" );
  FRENSIE_CHECK_EQUAL( local_properties.fileStartLine(), 2 );
  FRENSIE_CHECK_EQUAL( local_properties.fileVersion(), 3 );
  FRENSIE_CHECK_EQUAL( local_properties.tableName(), "2000.03e" );

  std::shared_ptr<const Data::ElectroatomicDataProperties>
    shared_properties;

  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP( shared_properties ) );
  FRENSIE_CHECK_EQUAL( shared_properties->atom(), Data::H_ATOM );
  FRENSIE_CHECK_EQUAL( shared_properties->atomicNumber(), 1 );
  FRENSIE_CHECK_EQUAL( shared_properties->atomicWeight(), 1.0*amu );
  FRENSIE_CHECK_EQUAL( shared_properties->filePath().string(),
                       "electroatomic_data/h_data.txt" );
  FRENSIE_CHECK_EQUAL( shared_properties->fileStartLine(), 10 );
  FRENSIE_CHECK_EQUAL( shared_properties->fileVersion(), 12 );
  FRENSIE_CHECK_EQUAL( shared_properties->tableName(), "1000.12p" );
}

//---------------------------------------------------------------------------//
// Custom setup
//---------------------------------------------------------------------------//
FRENSIE_CUSTOM_UNIT_TEST_SETUP_BEGIN();

FRENSIE_CUSTOM_UNIT_TEST_INIT()
{
  properties.reset( new Data::ACEElectroatomicDataProperties( 1.0*amu,
                                                              "electroatomic_data/h_data.txt",
                                                              10,
                                                              "1000.12p" ) );
}

FRENSIE_CUSTOM_UNIT_TEST_SETUP_END();

//---------------------------------------------------------------------------//
// end tstACEElectroatomicDataProperties.cpp
//---------------------------------------------------------------------------//
