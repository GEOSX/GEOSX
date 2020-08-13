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
#include "dataRepository/ConduitRestart.hpp"
#include "dataRepository/Wrapper.hpp"

// TPL includes
#include <gtest/gtest.h>

namespace geosx
{
namespace dataRepository
{
template< typename T >
Wrapper< array1d< T > > *
createArrayView( Group * parent,
                 const string & name,
                 int sfp,
                 const array1d< T > & data )
{
  Wrapper< array1d< T > > * view = parent->registerWrapper< array1d< T > >( name );
  view->setSizedFromParent( sfp );

  /* Resize the array */
  view->resize( data.size() );

  /* Check that the Wrapper size and byteSize return the proper values */
  EXPECT_EQ( view->size(), data.size() );

  /* Set the data */
  array1d< T > & view_data = view->reference();
  for( int i = 0; i < view->size(); i++ )
  {
    view_data[i] = data[i];
  }

  return view;
}

template< typename T >
void
checkArrayView( const Wrapper< array1d< T > > * view, int sfp, const array1d< T > & data )
{
  EXPECT_EQ( view->sizedFromParent(), sfp );
  EXPECT_EQ( view->size(), data.size() );
  arrayView1d< T const > const & view_data = view->reference();
  for( int i = 0; i < view->size(); i++ )
  {
    EXPECT_EQ( view_data[i], data[i] );
  }
}

template< typename T >
Wrapper< array2d< T > > *
createArray2dView( Group * parent,
                   const string & name,
                   int sfp,
                   const array2d< T > & data )
{
  Wrapper< array2d< T > > * view = parent->registerWrapper< array2d< T > >( name );
  view->setSizedFromParent( sfp );

  /* Resize the array */
  localIndex dims[2];
  dims[0] = data.size( 0 );
  dims[1] = data.size( 1 );
  view->resize( 2, dims );

  /* Check that the Wrapper size and byteSize return the proper values */
  EXPECT_EQ( view->size(), data.size() );

  /* Set the data */
  array2d< T > & view_data = view->reference();
  for( int i = 0; i < dims[0]; i++ )
  {
    for( int j = 0; j < dims[1]; j++ )
    {
      view_data[i][j] = data[i][j];
    }
  }

  return view;
}

template< typename T >
void
checkArray2dView( const Wrapper< array2d< T > > * view,
                  int sfp,
                  const array2d< T > & data )
{
  EXPECT_EQ( view->sizedFromParent(), sfp );
  EXPECT_EQ( view->size(), data.size() );

  arrayView2d< T const > const & view_data = view->reference();
  for( int i = 0; i < data.size( 0 ); i++ )
  {
    for( int j = 0; j < data.size( 1 ); j++ )
    {
      EXPECT_EQ( view_data[i][j], data[i][j] );
    }
  }
}

template< typename T >
Wrapper< SortedArray< T > > *
createSetView( Group * parent,
               const string & name,
               localIndex sfp,
               const SortedArray< T > & data )
{
  Wrapper< SortedArray< T > > * view = parent->registerWrapper< SortedArray< T > >( name );
  view->setSizedFromParent( int( sfp ) );

  /* Insert the data */
  view->reference().insert( data.begin(), data.end() );

  /* Check that the Wrapper size and byteSize return the proper values */
  EXPECT_EQ( view->size(), data.size() );

  return view;
}

template< typename T >
void
checkSetView( const Wrapper< SortedArray< T > > * view,
              localIndex sfp,
              const SortedArray< T > & data )
{
  EXPECT_EQ( view->sizedFromParent(), sfp );
  EXPECT_EQ( view->size(), data.size() );
  SortedArrayView< T const > const & view_data = view->reference();
  for( int i = 0; i < view->size(); i++ )
  {
    EXPECT_EQ( view_data[i], data[i] );
  }
}

Wrapper< string > *
createStringView( Group * parent,
                  const string & name,
                  int sfp,
                  const string & str )
{
  Wrapper< string > * view = parent->registerWrapper< string >( name );
  view->setSizedFromParent( sfp );

  /* Set the data */
  view->reference() = str;

  /* Check that the Wrapper size and byteSize return the proper values */
  EXPECT_EQ( static_cast< uint >( view->size() ), str.size() );

  return view;
}

void
checkStringView( const Wrapper< string > * view, const int sfp, const string & str )
{
  EXPECT_EQ( view->sizedFromParent(), sfp );
  EXPECT_EQ( view->reference(), str );
}

Wrapper< string_array > *
createStringArrayView( Group * parent,
                       const string & name,
                       int sfp,
                       const string_array & arr )
{
  Wrapper< string_array > * view = parent->registerWrapper< string_array >( name );
  view->setSizedFromParent( sfp );

  view->resize( static_cast< localIndex >( arr.size() ) );

  EXPECT_EQ( static_cast< uint >( view->size() ), arr.size() );

  string_array & view_data = view->reference();
  for( localIndex i = 0; i < arr.size(); ++i )
  {
    view_data[i] = arr[i];
  }

  return view;
}

void
checkStringArrayView( const Wrapper< string_array > * view,
                      const int sfp,
                      const string_array & arr )
{
  EXPECT_EQ( view->sizedFromParent(), sfp );
  EXPECT_EQ( view->size(), arr.size() );
  arrayView1d< string const > const & view_data = view->reference();
  for( int i = 0; i < view->size(); i++ )
  {
    EXPECT_EQ( view_data[i], arr[i] );
  }
}

template< typename T >
Wrapper< T > *
createScalarView( Group * parent,
                  const string & name,
                  int sfp,
                  const T & value )
{
  Wrapper< T > * view = parent->registerWrapper< T >( name );
  view->setSizedFromParent( sfp );

  /* Set the data */
  view->reference() = value;

  /* Check that the Wrapper size and byteSize return the proper values */
  EXPECT_EQ( view->size(), 1 );

  return view;
}

template< typename T >
void
checkScalarView( const Wrapper< T > * view, int sfp, const T value )
{
  EXPECT_EQ( view->sizedFromParent(), sfp );
  EXPECT_EQ( view->reference(), value );
}

TEST( testRestartExtended, testRestartExtended )
{
  const string path = "testRestartExtended";
  const int group_size = 44;
  int sfp = 55;

  /* Create a new Group directly below the conduit root. */
  Group * root = new Group( std::string( "data" ), nullptr );
  root->resize( group_size );

  /* Create a new globalIndex_array Wrapper. */
  string view_globalIndex_name = "globalIndex";
  int view_globalIndex_sfp = sfp++;
  int view_globalIndex_size = 100;
  globalIndex_array view_globalIndex_data( view_globalIndex_size );
  for( int i = 0; i < view_globalIndex_size; i++ )
  {
    view_globalIndex_data[i] = i * i * i;
  }
  createArrayView( root,
                   view_globalIndex_name,
                   view_globalIndex_sfp,
                   view_globalIndex_data );

  /* Create a new string Wrapper. */
  string view_hope_name = "hope";
  int view_hope_sfp = sfp++;
  string view_hope_str = "I sure hope these tests pass.";
  createStringView( root, view_hope_name, view_hope_sfp, view_hope_str );

  /* Create a new string array Wrapper. */
  string view_restart_name = "restart";
  int view_restart_sfp = sfp++;
  string_array view_restart_arr( 6 );
  view_restart_arr[0] = "Man ";
  view_restart_arr[1] = "this ";
  view_restart_arr[2] = "restart ";
  view_restart_arr[3] = "stuff ";
  view_restart_arr[4] = "better ";
  view_restart_arr[5] = "work.";
  createStringArrayView( root, view_restart_name, view_restart_sfp, view_restart_arr );

  /* Create a new group. */
  Group * strings_group = root->RegisterGroup( "strings" );
  strings_group->resize( group_size + 1 );

  /* Create a new string Wrapper. */
  string view_hello_name = "hello";
  int view_hello_sfp = sfp++;
  string view_hello_str = "Hello, how are you doing on this fine day?";
  createStringView( strings_group, view_hello_name, view_hello_sfp, view_hello_str );

  /* Create a new string Wrapper. */
  string view_goodbye_name = "goodbye";
  int view_goodbye_sfp = sfp++;
  string view_goodbye_str =
    "I hate this weather so I'm heading inside. Goodbye.";
  createStringView( strings_group,
                    view_goodbye_name,
                    view_goodbye_sfp,
                    view_goodbye_str );

  /* Create a new group. */
  Group * real64_group = root->RegisterGroup( "real64" );
  real64_group->resize( group_size + 2 );

  /* Create a new real64_array Wrapper. */
  string view_real641_name = "real641";
  int view_real641_sfp = sfp++;
  int view_real641_size = 1000;
  real64_array view_real641_data( view_real641_size );
  for( int i = 0; i < view_real641_size; i++ )
  {
    view_real641_data[i] = i * i / ( i + 5.0 );
  }
  createArrayView( real64_group,
                   view_real641_name,
                   view_real641_sfp,
                   view_real641_data );

  /* Create a new real64_array Wrapper. */
  string view_real642_name = "real642";
  int view_real642_sfp = sfp++;
  int view_real642_size = 1000;
  real64_array view_real642_data( view_real642_size );
  for( int i = 0; i < view_real642_size; i++ )
  {
    view_real642_data[i] = i * i / ( 5.0 + 5.0 * i + i * i );
  }
  createArrayView( real64_group,
                   view_real642_name,
                   view_real642_sfp,
                   view_real642_data );

  /* Create a new group. */
  Group * mixed_group = real64_group->RegisterGroup( "mixed" );
  mixed_group->resize( group_size + 3 );

  /* Create a new localIndex_array Wrapper. */
  string view_localIndex_name = "localIndex";
  int view_localIndex_sfp = sfp++;
  int view_localIndex_size = 953;
  localIndex_array view_localIndex_data( view_localIndex_size );
  for( localIndex i = 0; i < view_localIndex_size; i++ )
  {
    view_localIndex_data[i] = i * i - 100 * i + 3;
  }
  createArrayView( mixed_group,
                   view_localIndex_name,
                   view_localIndex_sfp,
                   view_localIndex_data );

  /* Create a new real32_array Wrapper. */
  string view_real32_name = "real32";
  int view_real32_sfp = sfp++;
  int view_real32_size = 782;
  real32_array view_real32_data( view_real32_size );
  for( int i = 0; i < view_real32_size; i++ )
  {
    view_real32_data[i] = ( i * i - 100.0f * i + 3.0f ) / ( i + 3.0f );
  }
  createArrayView( mixed_group, view_real32_name, view_real32_sfp, view_real32_data );

  /* Create a new string Wrapper. */
  string view_what_name = "what";
  int view_what_sfp = sfp++;
  string view_what_str = "What are you talking about? Who doesn't like storms?";
  createStringView( mixed_group, view_what_name, view_what_sfp, view_what_str );

  /* Create a new real64 Wrapper. */
  string view_pi_name = "pi";
  int view_pi_sfp = sfp++;
  real64 view_pi_value = 3.14159;
  createScalarView( mixed_group, view_pi_name, view_pi_sfp, view_pi_value );

  /* Create a new SortedArray<int> Wrapper. */
  string view_setlocalIndex_name = "view_setlocalIndex";
  localIndex view_setlocalIndex_sfp = sfp++;
  SortedArray< localIndex > view_setlocalIndex_set;
  view_setlocalIndex_set.insert( 4 );
  view_setlocalIndex_set.insert( 3 );
  view_setlocalIndex_set.insert( 10 );
  view_setlocalIndex_set.insert( 0 );
  view_setlocalIndex_set.insert( 1000 );
  view_setlocalIndex_set.insert( -1000 );
  createSetView( mixed_group,
                 view_setlocalIndex_name,
                 view_setlocalIndex_sfp,
                 view_setlocalIndex_set );

  /* Create a new SortedArray<string> Wrapper. */
  string view_setString_name = "view_setString";
  int view_setString_sfp = sfp++;
  SortedArray< string > view_setString_set;
  view_setString_set.insert( "zaa" );
  view_setString_set.insert( "aab" );
  view_setString_set.insert( "QD" );
  view_setString_set.insert( "az" );
  view_setString_set.insert( "g" );
  view_setString_set.insert( "this better be sorted" );
  createSetView( mixed_group,
                 view_setString_name,
                 view_setString_sfp,
                 view_setString_set );

  string view_real642d_name = "view_real642d";
  int view_real642d_sfp = sfp++;
  localIndex dim0 = 10;
  localIndex dim1 = 20;
  array2d< real64 > view_real642d_arr( dim0, dim1 );
  for( localIndex i = 0; i < dim0; i++ )
  {
    for( localIndex j = 0; j < dim1; j++ )
    {
      view_real642d_arr[i][j] = i * i / 3.0 + j;
    }
  }
  createArray2dView( mixed_group,
                     view_real642d_name,
                     view_real642d_sfp,
                     view_real642d_arr );

  string view_r1t2d_name = "view_r1t2d";
  int view_r1t2d_sfp = sfp++;
  dim0 = 10;
  dim1 = 20;
  array2d< R1Tensor > view_r1t2d_arr( dim0, dim1 );
  for( localIndex i = 0; i < dim0; i++ )
  {
    for( localIndex j = 0; j < dim1; j++ )
    {
      for( localIndex k = 0; k < 3; k++ )
      {
        view_r1t2d_arr[i][j] = i * i / 3.0 + j + i * j * k / 7.0;
      }
    }
  }
  createArray2dView( mixed_group, view_r1t2d_name, view_r1t2d_sfp, view_r1t2d_arr );

  /* Save the conduit tree */
  root->prepareToWrite();
  writeTree( path );
  root->finishWriting();

  /* Delete geos tree and reset conduit tree. */
  delete root;
  rootConduitNode.reset();

  /* Restore the conduit tree */
  loadTree( path );
  root = new Group( std::string( "data" ), nullptr );

  /* Create dual GEOS tree. Groups automatically register with the associated conduit::Node. */
  Wrapper< globalIndex_array > * view_globalIndex_new =
    root->registerWrapper< globalIndex_array >( view_globalIndex_name );
  Wrapper< string > * view_hope_new = root->registerWrapper< string >( view_hope_name );
  Wrapper< string_array > * view_restart_new =
    root->registerWrapper< string_array >( view_restart_name );

  Group * strings_group_new = root->RegisterGroup( "strings" );
  Wrapper< string > * view_hello_new =
    strings_group_new->registerWrapper< string >( view_hello_name );
  Wrapper< string > * view_goodbye_new =
    strings_group_new->registerWrapper< string >( view_goodbye_name );

  Group * real64_group_new = root->RegisterGroup( "real64" );
  Wrapper< real64_array > * view_real641_new =
    real64_group_new->registerWrapper< real64_array >( view_real641_name );
  Wrapper< real64_array > * view_real642_new =
    real64_group_new->registerWrapper< real64_array >( view_real642_name );

  Group * mixed_group_new = real64_group_new->RegisterGroup( "mixed" );
  Wrapper< localIndex_array > * view_localIndex_new =
    mixed_group_new->registerWrapper< localIndex_array >( view_localIndex_name );
  Wrapper< real32_array > * view_real32_new =
    mixed_group_new->registerWrapper< real32_array >( view_real32_name );
  Wrapper< string > * view_what_new =
    mixed_group_new->registerWrapper< string >( view_what_name );
  Wrapper< real64 > * view_pi_new =
    mixed_group_new->registerWrapper< real64 >( view_pi_name );
  Wrapper< SortedArray< localIndex > > * view_setlocalIndex_new =
    mixed_group_new->registerWrapper< SortedArray< localIndex > >(
      view_setlocalIndex_name );
  Wrapper< SortedArray< string > > * view_setString_new =
    mixed_group_new->registerWrapper< SortedArray< string > >( view_setString_name );
  Wrapper< array2d< real64 > > * view_real642d_new =
    mixed_group_new->registerWrapper< array2d< real64 > >( view_real642d_name );
  Wrapper< array2d< R1Tensor > > * view_r1t2d_new =
    mixed_group_new->registerWrapper< array2d< R1Tensor > >( view_r1t2d_name );

  /* Load the data */
  root->loadFromConduit();

  /* Group sizes should have carried over. */
  EXPECT_EQ( root->size(), group_size );
  EXPECT_EQ( strings_group_new->size(), group_size + 1 );
  EXPECT_EQ( real64_group_new->size(), group_size + 2 );
  EXPECT_EQ( mixed_group_new->size(), group_size + 3 );

  /* Check that Wrapper values were restored. */
  checkArrayView( view_globalIndex_new, view_globalIndex_sfp, view_globalIndex_data );
  checkStringView( view_hope_new, view_hope_sfp, view_hope_str );
  checkStringArrayView( view_restart_new, view_restart_sfp, view_restart_arr );
  checkStringView( view_hello_new, view_hello_sfp, view_hello_str );
  checkStringView( view_goodbye_new, view_goodbye_sfp, view_goodbye_str );
  checkArrayView( view_real641_new, view_real641_sfp, view_real641_data );
  checkArrayView( view_real642_new, view_real642_sfp, view_real642_data );
  checkArrayView( view_localIndex_new, view_localIndex_sfp, view_localIndex_data );
  checkArrayView( view_real32_new, view_real32_sfp, view_real32_data );
  checkStringView( view_what_new, view_what_sfp, view_what_str );
  checkScalarView( view_pi_new, view_pi_sfp, view_pi_value );
  checkSetView( view_setlocalIndex_new,
                view_setlocalIndex_sfp,
                view_setlocalIndex_set );
  checkSetView( view_setString_new, view_setString_sfp, view_setString_set );
  checkArray2dView( view_real642d_new, view_real642d_sfp, view_real642d_arr );
  checkArray2dView( view_r1t2d_new, view_r1t2d_sfp, view_r1t2d_arr );

  delete root;
}

} /* end namespace dataRepository */
} /* end namespace geosx */

int
main( int argc, char * argv[] )
{
  testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
