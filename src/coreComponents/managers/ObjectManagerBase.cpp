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

/**
 * @file ObjectManagerBase.cpp
 */

#include "ObjectManagerBase.hpp"

#include "common/TimingMacros.hpp"
#include "mesh/ExtrinsicMeshData.hpp"
#include "mpiCommunications/MpiWrapper.hpp"

namespace geosx
{
using namespace dataRepository;

ObjectManagerBase::ObjectManagerBase( std::string const & name, Group * const parent ) :
  Group( name, parent ),
  m_sets( groupKeyStruct::setsString, this ),
  m_neighborGroup( groupKeyStruct::neighborDataString, this ),
  m_localToGlobalMap(),
  m_globalToLocalMap(),
  m_isExternal(),
  m_ghostRank(),
  m_neighborData()
{
  RegisterGroup( groupKeyStruct::setsString, &m_sets );
  RegisterGroup( groupKeyStruct::neighborDataString, &m_neighborGroup );

  registerWrapper( viewKeyStruct::localToGlobalMapString, &m_localToGlobalMap )
    ->setApplyDefaultValue( -1 )
    ->setDescription(
      "Array that contains a map from localIndex to globalIndex." );

  registerWrapper( viewKeyStruct::globalToLocalMapString, &m_globalToLocalMap );

  registerWrapper( viewKeyStruct::isExternalString, &m_isExternal );

  registerWrapper( viewKeyStruct::ghostRankString, &m_ghostRank )
    ->setApplyDefaultValue( -2 )
    ->setPlotLevel( PlotLevel::LEVEL_0 );

  registerWrapper< array1d< integer > >( viewKeyStruct::domainBoundaryIndicatorString );

  m_sets.registerWrapper< SortedArray< localIndex > >(
    this->m_ObjectManagerBaseViewKeys.externalSet );
}

ObjectManagerBase::~ObjectManagerBase()
{}

ObjectManagerBase::CatalogInterface::CatalogType &
ObjectManagerBase::GetCatalog()
{
  static ObjectManagerBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void
ObjectManagerBase::CreateSet( const std::string & newSetName )
{
  m_sets.registerWrapper< SortedArray< localIndex > >( newSetName );
}

void
ObjectManagerBase::ConstructSetFromSetAndMap(
  SortedArrayView< localIndex const > const & inputSet,
  const array2d< localIndex > & map,
  const std::string & setName )
{
  SortedArray< localIndex > & newset =
    m_sets.getReference< SortedArray< localIndex > >( setName );
  newset.clear();

  localIndex const numObjects = size();
  GEOSX_ERROR_IF( map.size( 0 ) != numObjects,
                  "Size mismatch. " << map.size( 0 ) << " != " << numObjects );

  if( setName == "all" )
  {
    newset.reserve( numObjects );

    for( localIndex ka = 0; ka < numObjects; ++ka )
    {
      newset.insert( ka );
    }
  }
  else
  {
    localIndex const mapSize = map.size( 1 );
    for( localIndex ka = 0; ka < numObjects; ++ka )
    {
      if( std::all_of( &map( ka, 0 ), &map( ka, 0 ) + mapSize, [&]( localIndex const i ) {
            return inputSet.contains( i );
          } ) )
      {
        newset.insert( ka );
      }
    }
  }
}

void
ObjectManagerBase::ConstructSetFromSetAndMap(
  SortedArrayView< localIndex const > const & inputSet,
  const array1d< localIndex_array > & map,
  const std::string & setName )
{
  SortedArray< localIndex > & newset =
    m_sets.getReference< SortedArray< localIndex > >( setName );
  newset.clear();

  localIndex const numObjects = size();
  GEOSX_ERROR_IF( map.size() != numObjects,
                  "Size mismatch. " << map.size() << " != " << numObjects );

  if( setName == "all" )
  {
    newset.reserve( numObjects );

    for( localIndex ka = 0; ka < numObjects; ++ka )
    {
      newset.insert( ka );
    }
  }
  else
  {
    for( localIndex ka = 0; ka < numObjects; ++ka )
    {
      if( std::all_of( map[ka].begin(), map[ka].end(), [&]( localIndex const i ) {
            return inputSet.contains( i );
          } ) )
      {
        newset.insert( ka );
      }
    }
  }
}

void
ObjectManagerBase::ConstructSetFromSetAndMap(
  SortedArrayView< localIndex const > const & inputSet,
  ArrayOfArraysView< localIndex const > const & map,
  const std::string & setName )
{
  SortedArray< localIndex > & newset =
    m_sets.getReference< SortedArray< localIndex > >( setName );
  newset.clear();

  localIndex const numObjects = size();
  GEOSX_ERROR_IF( map.size() != numObjects,
                  "Size mismatch. " << map.size() << " != " << numObjects );

  if( setName == "all" )
  {
    newset.reserve( numObjects );

    for( localIndex ka = 0; ka < numObjects; ++ka )
    {
      newset.insert( ka );
    }
  }
  else
  {
    for( localIndex ka = 0; ka < numObjects; ++ka )
    {
      localIndex const * const values = map[ka];
      localIndex const numValues = map.sizeOfArray( ka );
      if( std::all_of( values, values + numValues, [&]( localIndex const i ) {
            return inputSet.contains( i );
          } ) )
      {
        newset.insert( ka );
      }
    }
  }
}

void
ObjectManagerBase::ConstructLocalListOfBoundaryObjects(
  localIndex_array & objectList ) const
{
  arrayView1d< integer const > const & isDomainBoundary =
    this->getReference< integer_array >(
      m_ObjectManagerBaseViewKeys.domainBoundaryIndicator );
  for( localIndex k = 0; k < size(); ++k )
  {
    if( isDomainBoundary[k] == 1 )
    {
      objectList.emplace_back( k );
    }
  }
}

void
ObjectManagerBase::ConstructGlobalListOfBoundaryObjects(
  globalIndex_array & objectList ) const
{
  arrayView1d< integer const > const & isDomainBoundary =
    this->getReference< integer_array >(
      m_ObjectManagerBaseViewKeys.domainBoundaryIndicator );
  for( localIndex k = 0; k < size(); ++k )
  {
    if( isDomainBoundary[k] == 1 )
    {
      objectList.emplace_back( this->m_localToGlobalMap[k] );
    }
  }
  std::sort( objectList.begin(), objectList.end() );
}

void
ObjectManagerBase::ConstructGlobalToLocalMap()
{
  GEOSX_MARK_FUNCTION;

  m_globalToLocalMap.clear();
  localIndex const N = size();
  m_globalToLocalMap.reserve( N );
  for( localIndex k = 0; k < N; ++k )
  {
    updateGlobalToLocalMap( k );
  }
}

localIndex
ObjectManagerBase::PackSize( string_array const & wrapperNames,
                             arrayView1d< localIndex const > const & packList,
                             integer const recursive,
                             bool on_device ) const
{
  localIndex packedSize = 0;
  buffer_unit_type * junk;
  packedSize +=
    this->PackPrivate< false >( junk, wrapperNames, packList, recursive, on_device );

  return packedSize;
}

localIndex
ObjectManagerBase::Pack( buffer_unit_type *& buffer,
                         string_array const & wrapperNames,
                         arrayView1d< localIndex const > const & packList,
                         integer const recursive,
                         bool on_device ) const
{
  localIndex packedSize = 0;

  packedSize +=
    this->PackPrivate< true >( buffer, wrapperNames, packList, recursive, on_device );

  return packedSize;
}

template< bool DOPACK >
localIndex
ObjectManagerBase::PackPrivate( buffer_unit_type *& buffer,
                                string_array const & wrapperNames,
                                arrayView1d< localIndex const > const & packList,
                                integer const recursive,
                                bool on_device ) const
{
  localIndex packedSize = 0;
  packedSize += bufferOps::Pack< DOPACK >( buffer, this->getName() );

  int const rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  packedSize += bufferOps::Pack< DOPACK >( buffer, rank );

  localIndex const numPackedIndices = packList.size();
  packedSize += bufferOps::Pack< DOPACK >( buffer, numPackedIndices );
  if( numPackedIndices > 0 )
  {
    packedSize += bufferOps::Pack< DOPACK >( buffer, string( "Wrappers" ) );

    string_array wrapperNamesForPacking;
    if( wrapperNames.size() == 0 )
    {
      SortedArray< localIndex > exclusionList;
      ViewPackingExclusionList( exclusionList );
      wrapperNamesForPacking.resize( this->wrappers().size() );
      localIndex count = 0;
      for( localIndex k = 0; k < this->wrappers().size(); ++k )
      {
        if( exclusionList.count( k ) == 0 )
        {
          wrapperNamesForPacking[count++] = wrappers().values()[k].first;
        }
      }
      wrapperNamesForPacking.resize( count );
    }
    else
    {
      wrapperNamesForPacking = wrapperNames;
    }

    packedSize += bufferOps::Pack< DOPACK >( buffer, wrapperNamesForPacking.size() );
    for( auto const & wrapperName : wrapperNamesForPacking )
    {
      dataRepository::WrapperBase const * const wrapper =
        this->getWrapperBase( wrapperName );
      if( wrapper != nullptr )
      {
        packedSize += bufferOps::Pack< DOPACK >( buffer, wrapperName );
        if( DOPACK )
        {
          packedSize += wrapper->PackByIndex( buffer, packList, on_device );
        }
        else
        {
          packedSize += wrapper->PackByIndexSize( packList, on_device );
        }
      }
      else
      {
        packedSize += bufferOps::Pack< DOPACK >( buffer, string( "nullptr" ) );
      }
    }
  }

  if( recursive > 0 )
  {
    packedSize += bufferOps::Pack< DOPACK >( buffer, string( "SubGroups" ) );
    packedSize += bufferOps::Pack< DOPACK >( buffer, this->GetSubGroups().size() );
    for( auto const & keyGroupPair : this->GetSubGroups() )
    {
      packedSize += bufferOps::Pack< DOPACK >( buffer, keyGroupPair.first );
      packedSize += keyGroupPair.second->Pack( buffer,
                                               wrapperNames,
                                               packList,
                                               recursive,
                                               on_device );
    }
  }

  packedSize += bufferOps::Pack< DOPACK >( buffer, this->getName() );

  return packedSize;
}

localIndex
ObjectManagerBase::Unpack( buffer_unit_type const *& buffer,
                           arrayView1d< localIndex > & packList,
                           integer const recursive,
                           bool on_device )
{
  localIndex unpackedSize = 0;
  string groupName;
  unpackedSize += bufferOps::Unpack( buffer, groupName );
  GEOSX_ERROR_IF_NE( groupName, this->getName() );

  int sendingRank;
  unpackedSize += bufferOps::Unpack( buffer, sendingRank );

  localIndex numUnpackedIndices;
  unpackedSize += bufferOps::Unpack( buffer, numUnpackedIndices );
  if( numUnpackedIndices > 0 )
  {
    string wrappersLabel;
    unpackedSize += bufferOps::Unpack( buffer, wrappersLabel );
    GEOSX_ERROR_IF_NE( wrappersLabel, "Wrappers" );

    localIndex numWrappers;
    unpackedSize += bufferOps::Unpack( buffer, numWrappers );
    for( localIndex a = 0; a < numWrappers; ++a )
    {
      string wrapperName;
      unpackedSize += bufferOps::Unpack( buffer, wrapperName );
      if( wrapperName != "nullptr" )
      {
        WrapperBase * const wrapper = this->getWrapperBase( wrapperName );
        unpackedSize += wrapper->UnpackByIndex( buffer, packList, on_device );
      }
    }
  }

  if( recursive > 0 )
  {
    string subGroups;
    unpackedSize += bufferOps::Unpack( buffer, subGroups );
    GEOSX_ERROR_IF_NE( subGroups, "SubGroups" );

    decltype( this->GetSubGroups().size() ) numSubGroups;
    unpackedSize += bufferOps::Unpack( buffer, numSubGroups );
    GEOSX_ERROR_IF_NE( numSubGroups, this->GetSubGroups().size() );

    for( localIndex i = 0; i < this->numSubGroups(); ++i )
    {
      string subGroupName;
      unpackedSize += bufferOps::Unpack( buffer, subGroupName );
      unpackedSize +=
        this->GetGroup( subGroupName )->Unpack( buffer, packList, recursive, on_device );
    }
  }

  unpackedSize += bufferOps::Unpack( buffer, groupName );
  GEOSX_ERROR_IF_NE( groupName, this->getName() );

  return unpackedSize;
}

template< bool DOPACK >
localIndex
ObjectManagerBase::PackParentChildMapsPrivate(
  buffer_unit_type *& buffer,
  arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  if( this->hasExtrinsicData< extrinsicMeshData::ParentIndex >() )
  {
    arrayView1d< localIndex const > const & parentIndex =
      this->getExtrinsicData< extrinsicMeshData::ParentIndex >();
    packedSize +=
      bufferOps::Pack< DOPACK >( buffer,
                                 string( extrinsicMeshData::ParentIndex::key ) );
    packedSize += bufferOps::Pack< DOPACK >( buffer,
                                             parentIndex,
                                             packList,
                                             this->m_localToGlobalMap,
                                             this->m_localToGlobalMap );
  }

  if( this->hasExtrinsicData< extrinsicMeshData::ChildIndex >() )
  {
    arrayView1d< localIndex const > const & childIndex =
      this->getExtrinsicData< extrinsicMeshData::ChildIndex >();
    packedSize +=
      bufferOps::Pack< DOPACK >( buffer, string( extrinsicMeshData::ChildIndex::key ) );
    packedSize += bufferOps::Pack< DOPACK >( buffer,
                                             childIndex,
                                             packList,
                                             this->m_localToGlobalMap,
                                             this->m_localToGlobalMap );
  }

  return packedSize;
}

template localIndex
ObjectManagerBase::PackParentChildMapsPrivate< true >(
  buffer_unit_type *& buffer,
  arrayView1d< localIndex const > const & packList ) const;
template localIndex
ObjectManagerBase::PackParentChildMapsPrivate< false >(
  buffer_unit_type *& buffer,
  arrayView1d< localIndex const > const & packList ) const;

localIndex
ObjectManagerBase::UnpackParentChildMaps( buffer_unit_type const *& buffer,
                                          localIndex_array & packList )
{
  localIndex unpackedSize = 0;

  if( this->hasExtrinsicData< extrinsicMeshData::ParentIndex >() )
  {
    arrayView1d< localIndex > const & parentIndex =
      this->getExtrinsicData< extrinsicMeshData::ParentIndex >();
    string shouldBeParentIndexString;
    unpackedSize += bufferOps::Unpack( buffer, shouldBeParentIndexString );
    GEOSX_ERROR_IF(
      shouldBeParentIndexString != extrinsicMeshData::ParentIndex::key,
      "value read from buffer is:" << shouldBeParentIndexString << ". It should be "
                                   << extrinsicMeshData::ParentIndex::key );
    unpackedSize += bufferOps::Unpack( buffer,
                                       parentIndex,
                                       packList,
                                       this->m_globalToLocalMap,
                                       this->m_globalToLocalMap );
  }

  if( this->hasExtrinsicData< extrinsicMeshData::ChildIndex >() )
  {
    arrayView1d< localIndex > const & childIndex =
      this->getExtrinsicData< extrinsicMeshData::ChildIndex >();
    string shouldBeChildIndexString;
    unpackedSize += bufferOps::Unpack( buffer, shouldBeChildIndexString );
    GEOSX_ERROR_IF(
      shouldBeChildIndexString != extrinsicMeshData::ChildIndex::key,
      "value read from buffer is:" << shouldBeChildIndexString << ". It should be "
                                   << extrinsicMeshData::ChildIndex::key );
    unpackedSize += bufferOps::Unpack( buffer,
                                       childIndex,
                                       packList,
                                       this->m_globalToLocalMap,
                                       this->m_globalToLocalMap );
  }

  return unpackedSize;
}

template< bool DOPACK >
localIndex
ObjectManagerBase::PackSets(
  buffer_unit_type *& buffer,
  arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;
  packedSize += bufferOps::Pack< DOPACK >( buffer, m_sets.getName() );

  packedSize += bufferOps::Pack< DOPACK >( buffer, m_sets.wrappers().size() );
  for( auto const & wrapperIter : m_sets.wrappers() )
  {
    string const & setName = wrapperIter.first;
    SortedArrayView< localIndex const > const & currentSet =
      m_sets.getReference< SortedArray< localIndex > >( setName );
    packedSize += bufferOps::Pack< DOPACK >( buffer, setName );
    packedSize +=
      bufferOps::Pack< DOPACK >( buffer,
                                 currentSet,
                                 packList,
                                 SortedArray< globalIndex >().toViewConst(),
                                 m_localToGlobalMap );
  }
  return packedSize;
}
//template<>
//localIndex ObjectManagerBase::PackSets<true>( buffer_unit_type * &,
//                                              arrayView1d<localIndex const> const & );
//template<>
//localIndex ObjectManagerBase::PackSets<false>( buffer_unit_type * &,
//                                               arrayView1d<localIndex const> const & );

localIndex
ObjectManagerBase::UnpackSets( buffer_unit_type const *& buffer )
{
  localIndex unpackedSize = 0;
  string name;
  unpackedSize += bufferOps::Unpack( buffer, name );
  GEOSX_ERROR_IF( name != m_sets.getName(),
                  "ObjectManagerBase::UnpackSets(): group names do not match" );

  localIndex numUnpackedSets;
  unpackedSize += bufferOps::Unpack( buffer, numUnpackedSets );
  for( localIndex a = 0; a < numUnpackedSets; ++a )
  {
    string setName;
    unpackedSize += bufferOps::Unpack( buffer, setName );
    SortedArray< localIndex > & targetSet =
      m_sets.getReference< SortedArray< localIndex > >( setName );

    SortedArray< globalIndex > junk;
    unpackedSize +=
      bufferOps::Unpack( buffer, targetSet, junk, this->m_globalToLocalMap, false );
  }

  return unpackedSize;
}

localIndex
ObjectManagerBase::PackGlobalMapsSize(
  arrayView1d< localIndex > const & packList,
  integer const recursive ) const
{
  buffer_unit_type * junk = nullptr;
  return PackGlobalMapsPrivate< false >( junk, packList, recursive );
}

localIndex
ObjectManagerBase::PackGlobalMaps( buffer_unit_type *& buffer,
                                   arrayView1d< localIndex > const & packList,
                                   integer const recursive ) const
{
  return PackGlobalMapsPrivate< true >( buffer, packList, recursive );
}

template< bool DOPACK >
localIndex
ObjectManagerBase::PackGlobalMapsPrivate(
  buffer_unit_type *& buffer,
  arrayView1d< localIndex const > const & packList,
  integer const recursive ) const
{
  localIndex packedSize = bufferOps::Pack< DOPACK >( buffer, this->getName() );

  // this doesn't link without the string()...no idea why.
  packedSize +=
    bufferOps::Pack< DOPACK >( buffer,
                               string( viewKeyStruct::localToGlobalMapString ) );

  int const rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  packedSize += bufferOps::Pack< DOPACK >( buffer, rank );

  localIndex const numPackedIndices = packList.size();
  packedSize += bufferOps::Pack< DOPACK >( buffer, numPackedIndices );

  if( numPackedIndices > 0 )
  {
    globalIndex_array globalIndices;
    globalIndices.resize( numPackedIndices );
    for( localIndex a = 0; a < numPackedIndices; ++a )
    {
      globalIndices[a] = this->m_localToGlobalMap[packList[a]];
    }
    packedSize += bufferOps::Pack< DOPACK >( buffer, globalIndices );
  }

  // FIXME is this the responsibility of this instance to do this?
  if( this->hasExtrinsicData< extrinsicMeshData::ParentIndex >() )
  {
    arrayView1d< localIndex const > const & parentIndex =
      this->getExtrinsicData< extrinsicMeshData::ParentIndex >();
    packedSize +=
      bufferOps::Pack< DOPACK >( buffer,
                                 string( extrinsicMeshData::ParentIndex::key ) );
    packedSize += bufferOps::Pack< DOPACK >( buffer,
                                             parentIndex,
                                             packList,
                                             this->m_localToGlobalMap,
                                             this->m_localToGlobalMap );
  }

  if( recursive > 0 )
  {
    packedSize += bufferOps::Pack< DOPACK >( buffer, string( "SubGroups" ) );
    packedSize += bufferOps::Pack< DOPACK >( buffer, this->GetSubGroups().size() );
    for( auto const & keyGroupPair : this->GetSubGroups() )
    {
      packedSize += bufferOps::Pack< DOPACK >( buffer, keyGroupPair.first );
      ObjectManagerBase const * const subObjectManager =
        Group::group_cast< ObjectManagerBase const * >( keyGroupPair.second );
      if( subObjectManager )
      {
        packedSize += subObjectManager->PackGlobalMapsPrivate< DOPACK >( buffer,
                                                                         packList,
                                                                         recursive );
      }
    }
  }

  packedSize += PackSets< DOPACK >( buffer, packList );

  return packedSize;
}

localIndex
ObjectManagerBase::UnpackGlobalMaps( buffer_unit_type const *& buffer,
                                     localIndex_array & packList,
                                     integer const recursive )
{
  GEOSX_MARK_FUNCTION;

  localIndex unpackedSize = 0;
  string groupName;
  unpackedSize += bufferOps::Unpack( buffer, groupName );
  GEOSX_ERROR_IF( groupName != this->getName(),
                  "ObjectManagerBase::Unpack(): group names do not match" );

  string localToGlobalString;
  unpackedSize += bufferOps::Unpack( buffer, localToGlobalString );
  GEOSX_ERROR_IF( localToGlobalString != viewKeyStruct::localToGlobalMapString,
                  "ObjectManagerBase::Unpack(): label incorrect" );

  int const rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  int sendingRank;
  unpackedSize += bufferOps::Unpack( buffer, sendingRank );

  localIndex numUnpackedIndices;
  unpackedSize += bufferOps::Unpack( buffer, numUnpackedIndices );

  if( numUnpackedIndices > 0 )
  {
    localIndex_array unpackedLocalIndices;
    unpackedLocalIndices.resize( numUnpackedIndices );

    globalIndex_array globalIndices;
    unpackedSize += bufferOps::Unpack( buffer, globalIndices );
    localIndex numNewIndices = 0;
    globalIndex_array newGlobalIndices;
    newGlobalIndices.reserve( numUnpackedIndices );
    localIndex const oldSize = this->size();
    for( localIndex a = 0; a < numUnpackedIndices; ++a )
    {
      // check to see if the object already exists by checking for the global
      // index in m_globalToLocalMap. If it doesn't, then add the object
      unordered_map< globalIndex, localIndex >::iterator iterG2L =
        m_globalToLocalMap.find( globalIndices[a] );
      if( iterG2L == m_globalToLocalMap.end() )
      {
        // object does not exist on this domain
        const localIndex newLocalIndex = oldSize + numNewIndices;

        // add the global index of the new object to the globalToLocal map
        m_globalToLocalMap[globalIndices[a]] = newLocalIndex;

        unpackedLocalIndices( a ) = newLocalIndex;

        newGlobalIndices.emplace_back( globalIndices[a] );

        ++numNewIndices;

        GEOSX_ERROR_IF( packList.size() != 0,
                        "ObjectManagerBase::Unpack(): packList specified, "
                        "but a new globalIndex is unpacked" );
      }
      else
      {
        // object already exists on this domain
        // get the local index of the node
        localIndex b = iterG2L->second;
        unpackedLocalIndices( a ) = b;
        if( ( sendingRank < rank && m_ghostRank[b] <= -1 ) ||
            ( sendingRank < m_ghostRank[b] ) )
        {
          m_ghostRank[b] = sendingRank;
        }
      }
    }

    // figure out new size of object container, and resize it
    const localIndex newSize = oldSize + numNewIndices;
    this->resize( newSize );

    // add the new indices to the maps.
    for( int a = 0; a < numNewIndices; ++a )
    {
      localIndex const b = oldSize + a;
      m_localToGlobalMap[b] = newGlobalIndices( a );
      m_ghostRank[b] = sendingRank;
    }

    packList = unpackedLocalIndices;
  }

  if( this->hasExtrinsicData< extrinsicMeshData::ParentIndex >() )
  {
    arrayView1d< localIndex > const & parentIndex =
      this->getExtrinsicData< extrinsicMeshData::ParentIndex >();
    string parentIndicesString;
    unpackedSize += bufferOps::Unpack( buffer, parentIndicesString );
    GEOSX_ERROR_IF( parentIndicesString != extrinsicMeshData::ParentIndex::key,
                    "ObjectManagerBase::Unpack(): label incorrect" );
    unpackedSize += bufferOps::Unpack( buffer,
                                       parentIndex,
                                       packList,
                                       this->m_globalToLocalMap,
                                       this->m_globalToLocalMap );
  }

  if( recursive > 0 )
  {
    string subGroups;
    unpackedSize += bufferOps::Unpack( buffer, subGroups );
    GEOSX_ERROR_IF( subGroups != "SubGroups",
                    "Group::Unpack(): group names do not match" );

    decltype( this->GetSubGroups().size() ) numSubGroups;
    unpackedSize += bufferOps::Unpack( buffer, numSubGroups );
    GEOSX_ERROR_IF( numSubGroups != this->GetSubGroups().size(),
                    "Group::Unpack(): incorrect number of subGroups" );

    for( auto const & index : this->GetSubGroups() )
    {
      GEOSX_UNUSED_VAR( index );
      string subGroupName;
      unpackedSize += bufferOps::Unpack( buffer, subGroupName );
      unpackedSize += this->GetGroup< ObjectManagerBase >( subGroupName )
                        ->UnpackGlobalMaps( buffer, packList, recursive );
    }
  }

  unpackedSize += UnpackSets( buffer );

  return unpackedSize;
}

void
ObjectManagerBase::ViewPackingExclusionList(
  SortedArray< localIndex > & exclusionList ) const
{
  exclusionList.insert(
    this->getWrapperIndex( viewKeyStruct::localToGlobalMapString ) );
  exclusionList.insert(
    this->getWrapperIndex( viewKeyStruct::globalToLocalMapString ) );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::ghostRankString ) );
  exclusionList.insert( this->getWrapperIndex( extrinsicMeshData::ParentIndex::key ) );
  exclusionList.insert( this->getWrapperIndex( extrinsicMeshData::ChildIndex::key ) );
}

localIndex
ObjectManagerBase::GetNumberOfGhosts() const
{
  localIndex rval = 0;
  for( localIndex i = 0; i < size(); ++i )
  {
    if( m_ghostRank[i] > -1 )
    {
      ++rval;
    }
  }
  return rval;
  //  return std::count_if( m_ghostRank.begin(), m_ghostRank.end(), [](integer i)->localIndex {return i>-1;} );
}

localIndex
ObjectManagerBase::GetNumberOfLocalIndices() const
{
  localIndex rval = 0;
  for( localIndex i = 0; i < size(); ++i )
  {
    if( m_ghostRank[i] <= -1 )
    {
      ++rval;
    }
  }
  return rval;
  //return std::count_if( m_ghostRank.begin(), m_ghostRank.end(), [](integer i)->localIndex {return i==-1;} );
}

void
ObjectManagerBase::SetReceiveLists()
{
  for( std::pair< int const, NeighborData > & pair : m_neighborData )
  {
    pair.second.ghostsToReceive().clear();
  }

  for( localIndex a = 0; a < size(); ++a )
  {
    if( m_ghostRank[a] > -1 )
    {
      getNeighborData( m_ghostRank[a] ).ghostsToReceive().emplace_back( a );
    }
  }
}

integer
ObjectManagerBase::SplitObject( localIndex const indexToSplit,
                                int const GEOSX_UNUSED_PARAM( rank ),
                                localIndex & newIndex )
{
  // if the object index has a zero sized childIndices entry, then this object can be split into two
  // new objects

  if( size() + 1 > capacity() )
  {
    reserve( static_cast< localIndex >( size() * m_overAllocationFactor ) );
  }

  // the new indices are tacked on to the end of the arrays
  newIndex = size();
  this->resize( newIndex + 1 );

  // copy the fields
  CopyObject( indexToSplit, newIndex );

  if( this->hasExtrinsicData< extrinsicMeshData::ParentIndex >() )
  {
    arrayView1d< localIndex > const & parentIndex =
      this->getExtrinsicData< extrinsicMeshData::ParentIndex >();
    parentIndex[newIndex] = indexToSplit;
  }

  if( this->hasExtrinsicData< extrinsicMeshData::ChildIndex >() )
  {
    arrayView1d< localIndex > const & childIndex =
      this->getExtrinsicData< extrinsicMeshData::ChildIndex >();
    childIndex[indexToSplit] = newIndex;
  }

  m_localToGlobalMap[newIndex] = -1;

  if( m_isExternal[indexToSplit] == 1 )
  {
    m_isExternal[indexToSplit] = 1;
    m_isExternal[newIndex] = 1;
  }
  else
  {
    m_isExternal[indexToSplit] = 2;
    m_isExternal[newIndex] = 2;
  }

  return 1;
}

void
ObjectManagerBase::inheritGhostRankFromParent( std::set< localIndex > const & indices )
{
  arrayView1d< localIndex const > const & parentIndex =
    this->getExtrinsicData< extrinsicMeshData::ParentIndex >();

  for( auto const a : indices )
  {
    m_ghostRank[a] = m_ghostRank[parentIndex[a]];
  }
}

void
ObjectManagerBase::CopyObject( const localIndex source,
                               const localIndex destination )
{
  for( auto & wrapper : wrappers() )
  {
    wrapper.second->copy( source, destination );
  }

  for( localIndex i = 0; i < m_sets.wrappers().size(); ++i )
  {
    SortedArray< localIndex > & targetSet =
      m_sets.getReference< SortedArray< localIndex > >( i );

#if !defined( __CUDA_ARCH__ )
    targetSet.move( LvArray::MemorySpace::CPU, true );
#endif

    if( targetSet.count( source ) > 0 )
    {
      targetSet.insert( destination );
    }
  }
}

void
ObjectManagerBase::SetMaxGlobalIndex()
{
  MpiWrapper::allReduce( &m_localMaxGlobalIndex,
                         &m_maxGlobalIndex,
                         1,
                         MPI_MAX,
                         MPI_COMM_GEOSX );
}

void
ObjectManagerBase::CleanUpMap( std::set< localIndex > const & targetIndices,
                               array1d< SortedArray< localIndex > > & upmap,
                               arrayView2d< localIndex const > const & downmap )
{
  for( auto const & targetIndex : targetIndices )
  {
    SortedArray< localIndex > eraseList;
    for( auto const & compositeIndex : upmap[targetIndex] )
    {
      bool hasTargetIndex = false;
      for( localIndex a = 0; a < downmap.size( 1 ); ++a )
      {
        localIndex const compositeLocalIndex = downmap[compositeIndex][a];
        if( compositeLocalIndex == targetIndex )
        {
          hasTargetIndex = true;
        }
      }
      if( !hasTargetIndex )
      {
        eraseList.insert( compositeIndex );
      }
    }
    for( auto const & val : eraseList )
    {
      upmap[targetIndex].remove( val );
    }
  }
}

void
ObjectManagerBase::CleanUpMap( std::set< localIndex > const & targetIndices,
                               ArrayOfSetsView< localIndex > const & upmap,
                               arrayView2d< localIndex const > const & downmap )
{
  std::vector< localIndex > eraseList;
  for( localIndex const targetIndex : targetIndices )
  {
    eraseList.clear();
    localIndex pos = 0;
    for( auto const & compositeIndex : upmap[targetIndex] )
    {
      bool hasTargetIndex = false;
      for( localIndex a = 0; a < downmap.size( 1 ); ++a )
      {
        localIndex const compositeLocalIndex = downmap[compositeIndex][a];
        if( compositeLocalIndex == targetIndex )
        {
          hasTargetIndex = true;
        }
      }

      if( !hasTargetIndex )
      {
        eraseList.emplace_back( pos );
      }

      ++pos;
    }

    localIndex const numUniqueIndices =
      LvArray::sortedArrayManipulation::makeSortedUnique( eraseList.begin(),
                                                          eraseList.end() );
    upmap.removeFromSet( targetIndex,
                         eraseList.begin(),
                         eraseList.begin() + numUniqueIndices );
  }
}

void
ObjectManagerBase::CleanUpMap(
  std::set< localIndex > const & targetIndices,
  array1d< SortedArray< localIndex > > & upmap,
  arrayView1d< arrayView1d< localIndex const > const > const & downmap )
{
  for( auto const & targetIndex : targetIndices )
  {
    SortedArray< localIndex > eraseList;
    for( auto const & compositeIndex : upmap[targetIndex] )
    {
      bool hasTargetIndex = false;
      for( localIndex a = 0; a < downmap[compositeIndex].size(); ++a )
      {
        localIndex const compositeLocalIndex = downmap[compositeIndex][a];
        if( compositeLocalIndex == targetIndex )
        {
          hasTargetIndex = true;
        }
      }
      if( !hasTargetIndex )
      {
        eraseList.insert( compositeIndex );
      }
    }
    for( auto const & val : eraseList )
    {
      upmap[targetIndex].remove( val );
    }
  }
}

void
ObjectManagerBase::CleanUpMap(
  std::set< localIndex > const & targetIndices,
  ArrayOfSetsView< localIndex > const & upmap,
  arrayView1d< arrayView1d< localIndex const > const > const & downmap )
{
  std::vector< localIndex > eraseList;
  for( localIndex const targetIndex : targetIndices )
  {
    eraseList.clear();
    localIndex pos = 0;
    for( localIndex const compositeIndex : upmap[targetIndex] )
    {
      bool hasTargetIndex = false;
      for( localIndex a = 0; a < downmap[compositeIndex].size(); ++a )
      {
        localIndex const compositeLocalIndex = downmap[compositeIndex][a];
        if( compositeLocalIndex == targetIndex )
        {
          hasTargetIndex = true;
        }
      }

      if( !hasTargetIndex )
      {
        eraseList.emplace_back( pos );
      }

      ++pos;
    }

    localIndex const numUniqueIndices =
      LvArray::sortedArrayManipulation::makeSortedUnique( eraseList.begin(),
                                                          eraseList.end() );
    upmap.removeFromSet( targetIndex,
                         eraseList.begin(),
                         eraseList.begin() + numUniqueIndices );
  }
}

void
ObjectManagerBase::CleanUpMap( std::set< localIndex > const & targetIndices,
                               ArrayOfSetsView< localIndex > const & upmap,
                               ArrayOfArraysView< localIndex const > const & downmap )
{
  std::vector< localIndex > eraseList;
  for( localIndex const targetIndex : targetIndices )
  {
    eraseList.clear();
    for( localIndex const compositeIndex : upmap[targetIndex] )
    {
      bool hasTargetIndex = false;
      for( localIndex const compositeLocalIndex : downmap[compositeIndex] )
      {
        if( compositeLocalIndex == targetIndex )
        {
          hasTargetIndex = true;
        }
      }

      if( !hasTargetIndex )
      {
        eraseList.emplace_back( compositeIndex );
      }
    }

    localIndex const numUniqueIndices =
      LvArray::sortedArrayManipulation::makeSortedUnique( eraseList.begin(),
                                                          eraseList.end() );
    upmap.removeFromSet( targetIndex,
                         eraseList.begin(),
                         eraseList.begin() + numUniqueIndices );
  }
}

void
ObjectManagerBase::enforceStateFieldConsistencyPostTopologyChange(
  std::set< localIndex > const & targetIndices )
{
  arrayView1d< localIndex const > const & childFaceIndices =
    getExtrinsicData< extrinsicMeshData::ChildIndex >();

  for( localIndex const targetIndex : targetIndices )
  {
    localIndex const childIndex = childFaceIndices[targetIndex];
    if( childIndex != -1 )
    {
      this->m_isExternal[targetIndex] = m_isExternal[childIndex];
    }
  }
}

void
ObjectManagerBase::moveSets( LvArray::MemorySpace const targetSpace )
{
  m_sets.forWrappers< SortedArray< localIndex > >( [&]( auto & wrapper ) {
    SortedArray< localIndex > & set = wrapper.reference();
    set.move( targetSpace );
  } );
}

} /* namespace geosx */
