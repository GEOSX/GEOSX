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
#include "dataRepository/wrapperHelpers.hpp"
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
void
checkNoSizeMethod( T const & var )
{
  EXPECT_EQ( wrapperHelpers::size( var ), 1 );
}

template< typename T >
void
checkSizeMethod( T const & var )
{
  EXPECT_EQ( wrapperHelpers::size( var ), var.size() );
}

template< typename T, int USD >
void
checkAverageOverSecondDimHelper( ArrayView< T const, 2, USD > const & input,
                                 ArrayView< T const, 1 > const & output )
{
  for( localIndex i = 0; i < input.size( 0 ); ++i )
  {
    T value {};
    for( localIndex j = 0; j < input.size( 1 ); ++j )
    {
      value += input( i, j );
    }

    value /= input.size( 1 );
    EXPECT_EQ( value, output( i ) );
  }
}

template< typename T, int USD >
void
checkAverageOverSecondDimHelper( ArrayView< T const, 3, USD > const & input,
                                 ArrayView< T const, 2 > const & output )
{
  for( localIndex i = 0; i < input.size( 0 ); ++i )
  {
    for( localIndex k = 0; k < input.size( 2 ); ++k )
    {
      T value {};
      for( localIndex j = 0; j < input.size( 1 ); ++j )
      {
        value += input( i, j, k );
      }

      value /= input.size( 1 );
      EXPECT_EQ( value, output( i, k ) );
    }
  }
}

template< typename T, int USD >
void
checkAverageOverSecondDimHelper( ArrayView< T const, 4, USD > const & input,
                                 ArrayView< T const, 3 > const & output )
{
  for( localIndex i = 0; i < input.size( 0 ); ++i )
  {
    for( localIndex k = 0; k < input.size( 2 ); ++k )
    {
      for( localIndex l = 0; l < input.size( 3 ); ++l )
      {
        T value {};
        for( localIndex j = 0; j < input.size( 1 ); ++j )
        {
          value += input( i, j, k, l );
        }

        value /= input.size( 1 );
        EXPECT_EQ( value, output( i, k, l ) );
      }
    }
  }
}

TEST( wrapperHelpers, size )
{
  checkNoSizeMethod( int() );
  checkNoSizeMethod( long() );
  checkNoSizeMethod( double() );
  checkNoSizeMethod( R1Tensor() );

  checkSizeMethod( std::string( "hello" ) );
  checkSizeMethod( std::set< int > { 1, 2, 3 } );
  checkSizeMethod( std::vector< int > { 1, 2, 3 } );
  checkSizeMethod( std::array< int, 5 > {} );
  checkSizeMethod( array1d< int >( 5 ) );
  checkSizeMethod( array2d< int >( 5, 5 ) );
}

/// TODO: dataPtr
/// TODO: resize
/// TODO: resizeDefault
/// TODO: resizeDimensions
/// TODO: byteSizeOfElement
/// TODO: byteSize
/// TODO: numElementsFromByteSize
/// TODO: reserve
/// TODO: capacity

using ArrayTypes = ::testing::Types< array2d< int, RAJA::PERM_IJ >,
                                     array2d< double, RAJA::PERM_JI >,
                                     array3d< double, RAJA::PERM_IJK >,
                                     array3d< int, RAJA::PERM_KJI >,
                                     array3d< R1Tensor, RAJA::PERM_JIK >,
                                     array4d< int, RAJA::PERM_IJKL >,
                                     array4d< int, RAJA::PERM_LKJI >,
                                     array4d< int, RAJA::PERM_KJLI > >;

template< typename ARRAY >
class PopulateMCArray : public ::testing::Test
{
public:
  static_assert( traits::is_array< ARRAY >, "T must be an LvArray::Array!" );
  using T = typename ARRAY::value_type;

  void
  test()
  {
    using ConduitType = typename conduitTypeInfo< T >::type;
    constexpr int numComponentsPerValue = conduitTypeInfo< T >::numConduitValues;

    fill( m_array, 20 );

    conduit::Node node;
    wrapperHelpers::populateMCArray( m_array.toViewConst(), node, {} );

    localIndex const numComponents = node.number_of_children();
    EXPECT_EQ( numComponents,
               numComponentsPerValue * m_array.size() / m_array.size( 0 ) );

    std::vector< conduit::DataArray< ConduitType > > valuesFromNode;
    for( localIndex i = 0; i < numComponents; ++i )
    {
      EXPECT_TRUE( node.child( i ).is_data_external() );
      conduit::DataArray< ConduitType > const conduitArray = node.child( i ).value();
      ASSERT_EQ( conduitArray.number_of_elements(), m_array.size( 0 ) );
      valuesFromNode.push_back( conduitArray );
    }

    for( localIndex i = 0; i < m_array.size( 0 ); ++i )
    {
      int curComponent = 0;
      LvArray::forValuesInSlice(
        m_array[i],
        [&valuesFromNode, &curComponent, i]( T const & value ) {
          for( int j = 0; j < numComponentsPerValue; ++j )
          {
            auto const valueOfComponent =
              *wrapperHelpers::internal::getPointerToComponent( value, j );
            EXPECT_EQ( valueOfComponent, valuesFromNode[curComponent][i] );
            ++curComponent;
          }
        } );
    }
  }

private:
  ARRAY m_array;
};

TYPED_TEST_SUITE( PopulateMCArray, ArrayTypes, );
TYPED_TEST( PopulateMCArray, test )
{
  this->test();
}

template< typename ARRAY >
class AverageOverSecondDim : public ::testing::Test
{
public:
  static_assert( traits::is_array< ARRAY >, "T must be an LvArray::Array!" );
  using T = typename ARRAY::value_type;

  static constexpr int NDIM = ARRAY::NDIM;
  static_assert( NDIM >= 2, "Cannot average over a 1D array!" );

  void
  test()
  {
    fill( m_array, 20 );

    std::unique_ptr< Array< T, NDIM - 1 > > const output =
      wrapperHelpers::averageOverSecondDim( m_array.toViewConst() );

    ASSERT_EQ( output->size( 0 ), m_array.size( 0 ) );

    for( int i = 2; i < NDIM; ++i )
    {
      ASSERT_EQ( output->size( i - 1 ), m_array.size( i ) );
    }

    checkAverageOverSecondDimHelper( m_array.toViewConst(), output->toViewConst() );
  }

private:
  ARRAY m_array;
};

TYPED_TEST_SUITE( AverageOverSecondDim, ArrayTypes, );
TYPED_TEST( AverageOverSecondDim, test )
{
  this->test();
}

}  // namespace testing
}  // namespace dataRepository
}  // end namespace geosx

int
main( int argc, char * argv[] )
{
  testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
