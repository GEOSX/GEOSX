/// Source includes
#include "common/DataTypes.hpp"

/// TPL includes
#include <gtest/gtest.h>

/// System includes
#include <random>

namespace geosx
{
namespace dataRepository
{
namespace testing
{
int
rand( int const min, int const max )
{
  static std::mt19937_64 gen;
  return std::uniform_int_distribution< int >( min, max )( gen );
}

template< typename T >
void
fill( T & val, localIndex )
{
  val = rand( -100, 100 );
}

void
fill( R1Tensor & val, localIndex )
{
  for( int i = 0; i < 3; ++i )
  {
    val[i] = rand( -100, 100 );
  }
}

void
fill( std::string & val, localIndex )
{
  int const num = rand( -100, 100 );
  val = std::to_string( num ) +
    std::string(
          " The rest of this is to avoid any small string optimizations. " ) +
    std::to_string( 2 * num );
}

template< typename T, typename U >
void
fill( std::pair< T, U > val, localIndex maxSize )
{
  fill( val.first, maxSize );
  fill( val.second, maxSize );
}

template< typename T >
void
fill( std::vector< T > & val, localIndex const maxSize )
{
  val.resize( rand( 0, maxSize ) );
  for( T & v : val )
  {
    fill( v, maxSize );
  }
}

template< typename T, int NDIM, typename PERMUTATION >
void
fill( Array< T, NDIM, PERMUTATION > & val, localIndex const maxSize )
{
  localIndex dims[NDIM];
  for( int i = 0; i < NDIM; ++i )
  {
    dims[i] = rand( 1, maxSize );
  }

  val.resize( NDIM, dims );

  for( localIndex i = 0; i < val.size(); ++i )
  {
    fill( val.data()[i], 1 );
  }
}

template< typename T >
void
fill( SortedArray< T > & val, localIndex const maxSize )
{
  int const nVals = rand( 0, maxSize );
  for( int i = 0; i < nVals; ++i )
  {
    T v;
    fill( v, maxSize );
    val.insert( v );
  }
}

template< typename K, typename V, typename SORTED >
void
fill( mapBase< K, V, SORTED > & val, localIndex const maxSize )
{
  int const nVals = rand( 0, maxSize );
  for( int i = 0; i < nVals; ++i )
  {
    K k;
    V v;
    fill( k, maxSize );
    fill( v, maxSize );
    val[k] = v;
  }
}

template< typename T >
void
compare( T const & val, T const & valFromFile )
{
  EXPECT_EQ( val, valFromFile );
}

template< typename T, int NDIM, typename PERMUTATION >
void
compare( Array< T, NDIM, PERMUTATION > const & val,
         Array< T, NDIM, PERMUTATION > const & valFromFile )
{
  ASSERT_EQ( val.size(), valFromFile.size() );
  for( localIndex i = 0; i < val.size(); ++i )
  {
    compare( val.data()[i], valFromFile.data()[i] );
  }
}

template< typename T >
void
compare( SortedArray< T > const & val, SortedArray< T > const & valFromFile )
{
  ASSERT_EQ( val.size(), valFromFile.size() );
  for( localIndex i = 0; i < val.size(); ++i )
  {
    compare( val[i], valFromFile[i] );
  }
}

}  // namespace testing
}  // namespace dataRepository
}  // namespace geosx
