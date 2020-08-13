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

#include <gtest/gtest.h>
#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"

using namespace geosx;
using namespace dataRepository;

template< typename T >
class WrapperSetGet : public ::testing::Test
{
public:
  WrapperSetGet() :
    m_group( "root", nullptr ),
    m_wrapper( "wrapper", &m_group ),
    m_wrapperBase( m_wrapper )
  {}

  void
  testSizedFromParent( int const value )
  {
    {
      Wrapper< T > * rval = m_wrapper.setSizedFromParent( value );
      EXPECT_EQ( value, m_wrapper.sizedFromParent() );
      EXPECT_EQ( rval, &m_wrapper );
    }

    {
      WrapperBase * rval = m_wrapperBase.setSizedFromParent( value );
      EXPECT_EQ( value, m_wrapperBase.sizedFromParent() );
      EXPECT_EQ( rval, &m_wrapperBase );
    }
  }

  void
  testRestartFlags( RestartFlags const value )
  {
    {
      Wrapper< T > * rval = m_wrapper.setRestartFlags( value );
      EXPECT_EQ( value, m_wrapper.getRestartFlags() );
      EXPECT_EQ( rval, &m_wrapper );
    }

    {
      WrapperBase * rval = m_wrapperBase.setRestartFlags( value );
      EXPECT_EQ( value, m_wrapperBase.getRestartFlags() );
      EXPECT_EQ( rval, &m_wrapperBase );
    }
  }

  void
  testPlotLevel( PlotLevel const value )
  {
    {
      Wrapper< T > * rval = m_wrapper.setPlotLevel( value );
      EXPECT_EQ( value, m_wrapper.getPlotLevel() );
      EXPECT_EQ( rval, &m_wrapper );
    }

    {
      WrapperBase * rval = m_wrapperBase.setPlotLevel( value );
      EXPECT_EQ( value, m_wrapperBase.getPlotLevel() );
      EXPECT_EQ( rval, &m_wrapperBase );
    }
  }

  void
  testInputFlags( InputFlags const value )
  {
    {
      Wrapper< T > * rval = m_wrapper.setInputFlag( value );
      EXPECT_EQ( value, m_wrapper.getInputFlag() );
      EXPECT_EQ( rval, &m_wrapper );
    }

    {
      WrapperBase * rval = m_wrapperBase.setInputFlag( value );
      EXPECT_EQ( value, m_wrapperBase.getInputFlag() );
      EXPECT_EQ( rval, &m_wrapperBase );
    }
  }

  void
  testDescription( std::string const & value )
  {
    {
      Wrapper< T > * rval = m_wrapper.setDescription( value );
      EXPECT_EQ( value, m_wrapper.getDescription() );
      EXPECT_EQ( rval, &m_wrapper );
    }

    {
      WrapperBase * rval = m_wrapperBase.setDescription( value );
      EXPECT_EQ( value, m_wrapperBase.getDescription() );
      EXPECT_EQ( rval, &m_wrapperBase );
    }
  }

private:
  Group m_group;
  Wrapper< T > m_wrapper;
  WrapperBase & m_wrapperBase;
};

using WrapperSetGetTypes = ::testing::Types< int, array1d< real64 > >;

TYPED_TEST_SUITE( WrapperSetGet, WrapperSetGetTypes, );

TYPED_TEST( WrapperSetGet, SizedFromParent )
{
  this->testSizedFromParent( true );
  this->testSizedFromParent( false );
}

TYPED_TEST( WrapperSetGet, RestartFlags )
{
  this->testRestartFlags( RestartFlags::NO_WRITE );
  this->testRestartFlags( RestartFlags::WRITE_AND_READ );
}

TYPED_TEST( WrapperSetGet, PlotLevel )
{
  this->testPlotLevel( PlotLevel::LEVEL_0 );
  this->testPlotLevel( PlotLevel::LEVEL_1 );
}

TYPED_TEST( WrapperSetGet, InputFlag )
{
  this->testInputFlags( InputFlags::OPTIONAL );
  this->testInputFlags( InputFlags::REQUIRED );
}

TYPED_TEST( WrapperSetGet, Description )
{
  this->testDescription( "First description." );
  this->testDescription( "Second description." );
}
