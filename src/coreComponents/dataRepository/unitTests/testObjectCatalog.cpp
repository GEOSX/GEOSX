/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#define OBJECTCATALOGVERBOSE 2

// Source includes
#include "dataRepository/ObjectCatalog.hpp"
#include "common/Logger.hpp"
#include "managers/initialization.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace LvArray;
using namespace geosx;

//START_SPHINX_BASE
class Base
{
public:
  Base( int & junk, double const & junk2 )
  {
    GEOSX_LOG( "calling Base constructor with arguments (" << junk << " "
                                                           << junk2 << ")" );
  }

  virtual ~Base()
  {
    GEOSX_LOG( "calling Base destructor" );
  }

  using CatalogInterface =
    dataRepository::CatalogInterface< Base, int &, double const & >;
  static CatalogInterface::CatalogType &
  GetCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  virtual std::string
  getCatalogName() = 0;
};
//STOP_SPHINX

//START_SPHINX_DERIVED1
class Derived1 : public Base
{
public:
  Derived1( int & junk, double const & junk2 ) :
    Base( junk, junk2 )
  {
    GEOSX_LOG( "calling Derived1 constructor with arguments (" << junk << " "
                                                               << junk2 << ")" );
  }

  ~Derived1()
  {
    GEOSX_LOG( "calling Derived1 destructor" );
  }
  static std::string
  CatalogName()
  {
    return "derived1";
  }
  std::string
  getCatalogName()
  {
    return CatalogName();
  }
};
REGISTER_CATALOG_ENTRY( Base, Derived1, int &, double const & )
//STOP_SPHINX

//START_SPHINX_DERIVED2
class Derived2 : public Base
{
public:
  Derived2( int & junk, double const & junk2 ) :
    Base( junk, junk2 )
  {
    GEOSX_LOG( "calling Derived2 constructor with arguments (" << junk << " "
                                                               << junk2 << ")" );
  }

  ~Derived2()
  {
    GEOSX_LOG( "calling Derived2 destructor" );
  }
  static std::string
  CatalogName()
  {
    return "derived2";
  }
  std::string
  getCatalogName()
  {
    return CatalogName();
  }
};
REGISTER_CATALOG_ENTRY( Base, Derived2, int &, double const & )
//STOP_SPHINX

//START_SPHINX_TEST
TEST( testObjectCatalog, testRegistration )
{
  GEOSX_LOG( "EXECUTING MAIN" );
  int junk = 1;
  double junk2 = 3.14;

  // allocate a new Derived1 object
  std::unique_ptr< Base > derived1 =
    Base::CatalogInterface::Factory( "derived1", junk, junk2 );

  // allocate a new Derived2 object
  std::unique_ptr< Base > derived2 =
    Base::CatalogInterface::Factory( "derived2", junk, junk2 );

  EXPECT_STREQ( derived1->getCatalogName().c_str(),
                Derived1::CatalogName().c_str() );

  EXPECT_STREQ( derived2->getCatalogName().c_str(),
                Derived2::CatalogName().c_str() );
  GEOSX_LOG( "EXITING MAIN" );
}
//STOP_SPHINX

int
main( int argc, char ** argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
