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

// Source includes
#include "common/DataTypes.hpp"
#include "managers/initialization.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"
#include "dataRepository/ConduitRestart.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "utils.hpp"

// TPL includes
#include <gtest/gtest.h>

// System includes
#include <random>

namespace geosx
{
namespace dataRepository
{
namespace testing
{
template< typename T >
class SingleWrapperTest : public ::testing::Test
{
public:
  virtual void
  SetUp() override
  {
    m_group = new Group( m_groupName, nullptr );
    m_groupSize = rand( 0, 100 );
    m_group->resize( m_groupSize );

    m_wrapper = m_group->registerWrapper< T >( m_wrapperName );
    m_wrapperSizedFromParent = rand( 0, 100 );
    m_wrapper->setSizedFromParent( m_wrapperSizedFromParent );
  }

  virtual void
  TearDown() override
  {
    delete m_group;
  }

  void
  test()
  {
    T value;
    fill( value, 100 );

    // Set the value
    m_wrapper->reference() = value;

    // Write out the tree
    m_group->prepareToWrite();
    writeTree( m_fileName );
    m_group->finishWriting();

    // Delete geosx tree and reset the conduit tree.
    delete m_group;
    rootConduitNode.reset();

    // Load in the tree
    loadTree( m_fileName );
    m_group = new Group( m_groupName, nullptr );
    m_wrapper = m_group->registerWrapper< T >( m_wrapperName );
    m_group->loadFromConduit();

    // Compare metadata
    EXPECT_EQ( m_group->size(), m_groupSize );
    EXPECT_EQ( m_wrapper->sizedFromParent(), m_wrapperSizedFromParent );

    // Compare the values
    compare( value, m_wrapper->reference() );
  }

private:
  std::string const m_groupName = "root";
  std::string const m_wrapperName = "wrapper";
  std::string const m_fileName = "testRestartBasic_SingleWrapperTest";

  Group * m_group;
  int m_groupSize;
  Wrapper< T > * m_wrapper;
  int m_wrapperSizedFromParent;
};

using TestTypes =
  ::testing::Types< int,
                    double,
                    R1Tensor,
                    std::pair< int, R1Tensor >,  // This should be passed to conduit via an external
                    // pointer but currently we're packing it.
                    std::pair< std::string, double >,
                    std::string,
                    std::vector< int >,
                    // std::vector< std::string > bufferOps currently can't pack this
                    array1d< double >,
                    array1d< std::string >,
                    array1d< array1d< double > >,
                    array1d< array1d< std::string > >,
                    array2d< double >,
                    array2d< std::string >,
                    array2d< double, RAJA::PERM_JI >,
                    array2d< std::string, RAJA::PERM_JI >,
                    array3d< double >,
                    array3d< std::string >,
                    array3d< double, RAJA::PERM_KJI >,
                    array3d< std::string, RAJA::PERM_IKJ >,
                    SortedArray< int >,
                    SortedArray< std::string >,
                    map< std::string, int >,
                    unordered_map< std::string, int >,
                    map< long, int >,
                    unordered_map< long, int > >;
TYPED_TEST_SUITE( SingleWrapperTest, TestTypes, );

TYPED_TEST( SingleWrapperTest, WriteAndRead )
{
  this->test();
}

}  // namespace testing
}  // namespace dataRepository
}  // namespace geosx

int
main( int argc, char * argv[] )
{
  testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
