#include "fileIO/timeHistory/TimeHistHDF.hpp"
#include "mainInterface/initialization.hpp"

#include <gtest/gtest.h>

using namespace geosx;

TEST( testHDFIO_parallel, SingleValueHistory )
{
  GEOSX_MARK_FUNCTION;
  string filename( "single_value_parallel" );
  HistoryMetadata spec( "Time History", 1, std::type_index( typeid(real64)));

  int rank = MpiWrapper::commRank( );

  HDFHistIO io( filename, spec );
  io.init( true );
  real64 val = 0.0;
  for( localIndex tidx = 0; tidx < 100; ++tidx )
  {
    val += 0.5 * (rank+1);
    buffer_unit_type * buffer = io.getBufferHead( );
    memcpy( buffer, &val, sizeof(real64));
  }

  io.write( );
  io.compressInFile();
}


TEST( testHDFIO_parallel, DynamicHistory )
{
  GEOSX_MARK_FUNCTION;
  string filename( "dynamic_parallel" );
  HistoryMetadata singleItemSpec( "some_real", 1, std::type_index( typeid(real64)));

  int rnk = MpiWrapper::commRank( );

  HDFHistIO io( filename, singleItemSpec );
  io.init( true );
  real64 val = 0.0;
  for( localIndex tidx = 0; tidx < 4; ++tidx )
  {
    val += 0.5 * ( rnk + 1 );
    if( ( rnk + tidx ) % 2 == 0 )
    {
      io.updateCollectingCount( 1 );
      buffer_unit_type * buffer = io.getBufferHead( );
      memcpy( buffer, &val, sizeof(real64));
    }
    else
    {
      io.updateCollectingCount( 0 );
      io.getBufferHead( );
    }
  }
  io.write( );
  val = 0.0;
  for( localIndex tidx = 0; tidx < 8; ++tidx )
  {
    val += 0.5 * ( rnk + 1 );
    if( ( rnk + tidx ) % 2 == 1 )
    {
      io.updateCollectingCount( 1 );
      buffer_unit_type * buffer = io.getBufferHead( );
      memcpy( buffer, &val, sizeof(real64) );
    }
    else
    {
      io.updateCollectingCount( 0 );
      io.getBufferHead( );
    }
  }
  io.write( );
  io.compressInFile();
}

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  geosx::basicSetup( ac, av );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
