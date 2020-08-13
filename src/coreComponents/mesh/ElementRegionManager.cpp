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

#include <map>
#include <vector>

#include "ElementRegionManager.hpp"

#include "common/TimingMacros.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "FaceElementRegion.hpp"
#include "FaceManager.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "CellBlockManager.hpp"
#include "meshUtilities/MeshManager.hpp"

namespace geosx
{
using namespace dataRepository;

ElementRegionManager::ElementRegionManager( string const & name, Group * const parent ) :
  ObjectManagerBase( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
  this->RegisterGroup< Group >(
    ElementRegionManager::groupKeyStruct::elementRegionsGroup );
}

ElementRegionManager::~ElementRegionManager()
{
  // TODO Auto-generated destructor stub
}

localIndex
ElementRegionManager::numCellBlocks() const
{
  localIndex numCellBlocks = 0;
  this->forElementSubRegions< ElementSubRegionBase >(
    [&]( ElementSubRegionBase const & ) { numCellBlocks += 1; } );
  return numCellBlocks;
}

void
ElementRegionManager::resize(
  integer_array const & numElements,
  string_array const & regionNames,
  string_array const & GEOSX_UNUSED_PARAM( elementTypes ) )
{
  localIndex const n_regions =
    LvArray::integerConversion< localIndex >( regionNames.size() );
  for( localIndex reg = 0; reg < n_regions; ++reg )
  {
    ElementRegionBase * elemRegion = this->GetRegion( reg );
    elemRegion->resize( numElements[reg] );
  }
}

void
ElementRegionManager::SetMaxGlobalIndex()
{
  forElementSubRegions< ElementSubRegionBase >(
    [this]( ElementSubRegionBase const & subRegion ) {
      m_localMaxGlobalIndex =
        std::max( m_localMaxGlobalIndex, subRegion.maxGlobalIndex() );
    } );

  MpiWrapper::allReduce( &m_localMaxGlobalIndex,
                         &m_maxGlobalIndex,
                         1,
                         MPI_MAX,
                         MPI_COMM_GEOSX );
}

Group *
ElementRegionManager::CreateChild( string const & childKey,
                                   string const & childName )
{
  GEOSX_ERROR_IF(
    !( CatalogInterface::hasKeyName( childKey ) ),
    "KeyName (" << childKey << ") not found in ObjectManager::Catalog" );
  GEOSX_LOG_RANK_0( "Adding Object " << childKey << " named " << childName
                                     << " from ObjectManager::Catalog." );
  Group * const elementRegions =
    this->GetGroup( ElementRegionManager::groupKeyStruct::elementRegionsGroup );
  return elementRegions->RegisterGroup(
    childName,
    CatalogInterface::Factory( childKey, childName, elementRegions ) );
}

void
ElementRegionManager::ExpandObjectCatalogs()
{
  ObjectManagerBase::CatalogInterface::CatalogType const & catalog =
    ObjectManagerBase::GetCatalog();
  for( ObjectManagerBase::CatalogInterface::CatalogType::const_iterator iter =
         catalog.begin();
       iter != catalog.end();
       ++iter )
  {
    string const key = iter->first;
    if( key.find( "ElementRegion" ) != std::string::npos )
    {
      this->CreateChild( key, key );
    }
  }
}

void
ElementRegionManager::SetSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
                                           xmlWrapper::xmlNode schemaParent,
                                           integer documentationType )
{
  xmlWrapper::xmlNode targetChoiceNode = schemaParent.child( "xsd:choice" );
  if( targetChoiceNode.empty() )
  {
    targetChoiceNode = schemaParent.prepend_child( "xsd:choice" );
    targetChoiceNode.append_attribute( "minOccurs" ) = "0";
    targetChoiceNode.append_attribute( "maxOccurs" ) = "unbounded";
  }

  std::set< string > names;
  this->forElementRegions( [&]( ElementRegionBase & elementRegion ) {
    names.insert( elementRegion.getName() );
  } );

  for( string const & name : names )
  {
    ElementRegionBase * const elementRegion = GetRegion( name );
    schemaUtilities::SchemaConstruction( elementRegion,
                                         schemaRoot,
                                         targetChoiceNode,
                                         documentationType );
  }
}

void
ElementRegionManager::GenerateMesh( Group * const cellBlockManager )
{
  this->forElementRegions< CellElementRegion >( [&]( CellElementRegion & elemRegion ) {
    elemRegion.GenerateMesh( cellBlockManager->GetGroup( keys::cellBlocks ) );
  } );
}

void
ElementRegionManager::GenerateCellToEdgeMaps( FaceManager const * const faceManager )
{
  /*
   * Create cell to edges map
   * I use the existing maps from cells to faces and from faces to edges.
   */
  localIndex faceIndex, edgeIndex;
  int count = 0;
  bool isUnique = true;

  this->forElementSubRegions< CellElementSubRegion >(
    [&]( CellElementSubRegion & subRegion ) {
      FixedOneToManyRelation & cellToEdges = subRegion.edgeList();
      FixedOneToManyRelation const & cellToFaces = subRegion.faceList();
      InterObjectRelation< ArrayOfArrays< localIndex > > const & faceToEdges =
        faceManager->edgeList();

      //loop over the cells
      for( localIndex kc = 0; kc < subRegion.size(); kc++ )
      {
        // loop over the faces
        // std::cout << "Element: " << kc << std::endl;
        count = 0;
        for( localIndex kf = 0; kf < subRegion.numFacesPerElement(); kf++ )
        {
          // loop over edges of each face
          faceIndex = cellToFaces[kc][kf];
          for( localIndex ke = 0; ke < faceToEdges.sizeOfArray( faceIndex ); ke++ )
          {
            isUnique = true;
            edgeIndex = faceToEdges[faceIndex][ke];

            //loop over edges that have already been added to the element.
            for( localIndex kec = 0; kec < count + 1; kec++ )
            {
              // make sure that the edge has not been counted yet
              if( cellToEdges( kc, kec ) == edgeIndex )
              {
                isUnique = false;
                break;
              }
            }
            if( isUnique )
            {
              // std::cout << edgeIndex << std::endl;
              cellToEdges( kc, count ) = edgeIndex;
              count++;
            }

          }  // end edge loop
        }    // end face loop
      }      // end cell loop
    } );
}

void
ElementRegionManager::GenerateAggregates( FaceManager const * const faceManager,
                                          NodeManager const * const nodeManager )
{
  this->forElementRegions< CellElementRegion >( [&]( CellElementRegion & elemRegion ) {
    elemRegion.GenerateAggregates( faceManager, nodeManager );
  } );
}

void
ElementRegionManager::GenerateWells( MeshManager * const meshManager,
                                     MeshLevel * const meshLevel )
{
  NodeManager * const nodeManager = meshLevel->getNodeManager();

  // get the offsets to construct local-to-global maps for well nodes and elements
  nodeManager->SetMaxGlobalIndex();
  globalIndex const nodeOffsetGlobal = nodeManager->maxGlobalIndex() + 1;
  localIndex const elemOffsetLocal = this->getNumberOfElements();
  globalIndex const elemOffsetGlobal = MpiWrapper::Sum( elemOffsetLocal );

  globalIndex wellElemCount = 0;
  globalIndex wellNodeCount = 0;

  // construct the wells one by one
  forElementRegions< WellElementRegion >( [&]( WellElementRegion & wellRegion ) {
    // get the global well geometry from the well generator
    string const generatorName = wellRegion.GetWellGeneratorName();
    InternalWellGenerator const * const wellGeometry =
      meshManager->GetGroup< InternalWellGenerator >( generatorName );

    GEOSX_ERROR_IF( wellGeometry == nullptr,
                    "InternalWellGenerator " << generatorName
                                             << " not found in well "
                                             << wellRegion.getName() );

    // generate the local data (well elements, nodes, perforations) on this well
    // note: each MPI rank knows the global info on the entire well (constructed earlier in InternalWellGenerator)
    // so we only need node and element offsets to construct the local-to-global maps in each wellElemSubRegion
    wellRegion.GenerateWell( *meshLevel,
                             *wellGeometry,
                             nodeOffsetGlobal + wellNodeCount,
                             elemOffsetGlobal + wellElemCount );

    // increment counters with global number of nodes and elements
    wellElemCount += wellGeometry->GetNumElements();
    wellNodeCount += wellGeometry->GetNumNodes();

    string const subRegionName = wellRegion.GetSubRegionName();
    WellElementSubRegion * const subRegion =
      wellRegion.GetGroup( ElementRegionBase::viewKeyStruct::elementSubRegions )
        ->GetGroup< WellElementSubRegion >( subRegionName );

    GEOSX_ERROR_IF( subRegion == nullptr,
                    "Subregion " << subRegionName << " not found in well "
                                 << wellRegion.getName() );

    globalIndex const numWellElemsGlobal = MpiWrapper::Sum( subRegion->size() );

    GEOSX_ERROR_IF( numWellElemsGlobal != wellGeometry->GetNumElements(),
                    "Invalid partitioning in well " << subRegionName );
  } );

  // communicate to rebuild global node info since we modified global ordering
  nodeManager->SetMaxGlobalIndex();
}

int
ElementRegionManager::PackSize(
  string_array const & wrapperNames,
  ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackPrivate< false >( junk, wrapperNames, packList );
}

int
ElementRegionManager::Pack(
  buffer_unit_type *& buffer,
  string_array const & wrapperNames,
  ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  return PackPrivate< true >( buffer, wrapperNames, packList );
}

template< bool DOPACK >
int
ElementRegionManager::PackPrivate(
  buffer_unit_type *& buffer,
  string_array const & wrapperNames,
  ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  int packedSize = 0;

  //  packedSize += Group::Pack( buffer, wrapperNames, {}, 0, 0);

  packedSize += bufferOps::Pack< DOPACK >( buffer, this->getName() );
  packedSize += bufferOps::Pack< DOPACK >( buffer, numRegions() );

  for( typename dataRepository::indexType kReg = 0; kReg < numRegions(); ++kReg )
  {
    ElementRegionBase const * const elemRegion = GetRegion( kReg );
    packedSize += bufferOps::Pack< DOPACK >( buffer, elemRegion->getName() );

    packedSize += bufferOps::Pack< DOPACK >( buffer, elemRegion->numSubRegions() );

    elemRegion->forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase const & subRegion ) {
        packedSize += bufferOps::Pack< DOPACK >( buffer, subRegion.getName() );

        arrayView1d< localIndex > const elemList = packList[kReg][esr];
        if( DOPACK )
        {
          packedSize += subRegion.Pack( buffer, wrapperNames, elemList, 0 );
        }
        else
        {
          packedSize += subRegion.PackSize( wrapperNames, elemList, 0 );
        }
      } );
  }

  return packedSize;
}
//template int ElementRegionManager::PackPrivate<true>( buffer_unit_type * &,
//                                                      string_array const &,
//                                                      ElementViewAccessor<arrayView1d<localIndex>> const & ) const;
//template int ElementRegionManager::PackPrivate<false>( buffer_unit_type * &,
//                                                      string_array const &,
//                                                      ElementViewAccessor<arrayView1d<localIndex>> const & ) const;

int
ElementRegionManager::Unpack( buffer_unit_type const *& buffer,
                              ElementViewAccessor< arrayView1d< localIndex > > & packList )
{
  return UnpackPrivate( buffer, packList );
}

int
ElementRegionManager::Unpack(
  buffer_unit_type const *& buffer,
  ElementReferenceAccessor< array1d< localIndex > > & packList )
{
  return UnpackPrivate( buffer, packList );
}

template< typename T >
int
ElementRegionManager::UnpackPrivate( buffer_unit_type const *& buffer,
                                     T & packList )
{
  int unpackedSize = 0;

  string name;
  unpackedSize += bufferOps::Unpack( buffer, name );

  GEOSX_ERROR_IF( name != this->getName(),
                  "Unpacked name (" << name << ") does not equal object name ("
                                    << this->getName() );

  localIndex numRegionsRead;
  unpackedSize += bufferOps::Unpack( buffer, numRegionsRead );

  for( localIndex kReg = 0; kReg < numRegionsRead; ++kReg )
  {
    string regionName;
    unpackedSize += bufferOps::Unpack( buffer, regionName );

    ElementRegionBase * const elemRegion = GetRegion( regionName );

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
    elemRegion->forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase & subRegion ) {
        string subRegionName;
        unpackedSize += bufferOps::Unpack( buffer, subRegionName );

        /// THIS IS WRONG??
        arrayView1d< localIndex > & elemList = packList[kReg][esr];

        unpackedSize += subRegion.Unpack( buffer, elemList, 0 );
      } );
  }

  return unpackedSize;
}

int
ElementRegionManager::PackGlobalMapsSize(
  ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackGlobalMapsPrivate< false >( junk, packList );
}

int
ElementRegionManager::PackGlobalMaps(
  buffer_unit_type *& buffer,
  ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  return PackGlobalMapsPrivate< true >( buffer, packList );
}
template< bool DOPACK >
int
ElementRegionManager::PackGlobalMapsPrivate(
  buffer_unit_type *& buffer,
  ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  int packedSize = 0;

  packedSize += bufferOps::Pack< DOPACK >( buffer, numRegions() );

  for( typename dataRepository::indexType kReg = 0; kReg < numRegions(); ++kReg )
  {
    ElementRegionBase const * const elemRegion = GetRegion( kReg );
    packedSize += bufferOps::Pack< DOPACK >( buffer, elemRegion->getName() );

    packedSize += bufferOps::Pack< DOPACK >( buffer, elemRegion->numSubRegions() );
    elemRegion->forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase const & subRegion ) {
        packedSize += bufferOps::Pack< DOPACK >( buffer, subRegion.getName() );

        arrayView1d< localIndex > const & elemList = packList[kReg][esr];
        if( DOPACK )
        {
          packedSize += subRegion.PackGlobalMaps( buffer, elemList, 0 );
        }
        else
        {
          packedSize += subRegion.PackGlobalMapsSize( elemList, 0 );
        }
      } );
  }

  return packedSize;
}

int
ElementRegionManager::UnpackGlobalMaps(
  buffer_unit_type const *& buffer,
  ElementViewAccessor< ReferenceWrapper< localIndex_array > > & packList )
{
  int unpackedSize = 0;

  localIndex numRegionsRead;
  unpackedSize += bufferOps::Unpack( buffer, numRegionsRead );

  packList.resize( numRegionsRead );
  for( localIndex kReg = 0; kReg < numRegionsRead; ++kReg )
  {
    string regionName;
    unpackedSize += bufferOps::Unpack( buffer, regionName );

    ElementRegionBase * const elemRegion = GetRegion( regionName );

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
    packList[kReg].resize( numSubRegionsRead );
    elemRegion->forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase & subRegion ) {
        string subRegionName;
        unpackedSize += bufferOps::Unpack( buffer, subRegionName );

        /// THIS IS WRONG
        localIndex_array & elemList = packList[kReg][esr].get();

        unpackedSize += subRegion.UnpackGlobalMaps( buffer, elemList, 0 );
      } );
  }

  return unpackedSize;
}

int
ElementRegionManager::PackUpDownMapsSize(
  ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate< false >( junk, packList );
}
int
ElementRegionManager::PackUpDownMapsSize(
  ElementReferenceAccessor< array1d< localIndex > > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate< false >( junk, packList );
}

int
ElementRegionManager::PackUpDownMaps(
  buffer_unit_type *& buffer,
  ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  return PackUpDownMapsPrivate< true >( buffer, packList );
}
int
ElementRegionManager::PackUpDownMaps(
  buffer_unit_type *& buffer,
  ElementReferenceAccessor< array1d< localIndex > > const & packList ) const
{
  return PackUpDownMapsPrivate< true >( buffer, packList );
}

template< bool DOPACK, typename T >
int
ElementRegionManager::PackUpDownMapsPrivate( buffer_unit_type *& buffer,
                                             T const & packList ) const
{
  int packedSize = 0;

  packedSize += bufferOps::Pack< DOPACK >( buffer, numRegions() );

  for( typename dataRepository::indexType kReg = 0; kReg < numRegions(); ++kReg )
  {
    ElementRegionBase const * const elemRegion = GetRegion( kReg );
    packedSize += bufferOps::Pack< DOPACK >( buffer, elemRegion->getName() );

    packedSize += bufferOps::Pack< DOPACK >( buffer, elemRegion->numSubRegions() );
    elemRegion->forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase const & subRegion ) {
        packedSize += bufferOps::Pack< DOPACK >( buffer, subRegion.getName() );

        arrayView1d< localIndex > const & elemList = packList[kReg][esr];
        if( DOPACK )
        {
          packedSize += subRegion.PackUpDownMaps( buffer, elemList );
        }
        else
        {
          packedSize += subRegion.PackUpDownMapsSize( elemList );
        }
      } );
  }

  return packedSize;
}
//template int
//ElementRegionManager::
//PackUpDownMapsPrivate<true>( buffer_unit_type * & buffer,
//                             ElementViewAccessor<arrayView1d<localIndex>> const & packList ) const;
//template int
//ElementRegionManager::
//PackUpDownMapsPrivate<false>( buffer_unit_type * & buffer,
//                             ElementViewAccessor<arrayView1d<localIndex>> const & packList ) const;

int
ElementRegionManager::UnpackUpDownMaps(
  buffer_unit_type const *& buffer,
  ElementReferenceAccessor< localIndex_array > & packList,
  bool const overwriteMap )
{
  int unpackedSize = 0;

  localIndex numRegionsRead;
  unpackedSize += bufferOps::Unpack( buffer, numRegionsRead );

  for( localIndex kReg = 0; kReg < numRegionsRead; ++kReg )
  {
    string regionName;
    unpackedSize += bufferOps::Unpack( buffer, regionName );

    ElementRegionBase * const elemRegion = GetRegion( regionName );

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
    elemRegion->forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const kSubReg, ElementSubRegionBase & subRegion ) {
        string subRegionName;
        unpackedSize += bufferOps::Unpack( buffer, subRegionName );

        /// THIS IS WRONG
        localIndex_array & elemList = packList[kReg][kSubReg];
        unpackedSize +=
          subRegion.UnpackUpDownMaps( buffer, elemList, false, overwriteMap );
      } );
  }

  return unpackedSize;
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase,
                        ElementRegionManager,
                        string const &,
                        Group * const )
}  // namespace geosx
