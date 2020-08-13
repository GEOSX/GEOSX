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

#include "CellElementSubRegion.hpp"

#include "constitutive/ConstitutiveManager.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

CellElementSubRegion::CellElementSubRegion( string const & name, Group * const parent ) :
  CellBlock( name, parent )
{
  registerWrapper( viewKeyStruct::constitutiveGroupingString,
                   &m_constitutiveGrouping )
    ->setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::constitutivePointVolumeFraction,
                   &m_constitutivePointVolumeFraction );

  registerWrapper( viewKeyStruct::dNdXString, &m_dNdX )
    ->setSizedFromParent( 1 )
    ->reference()
    .resizeDimension< 3 >( 3 );

  registerWrapper( viewKeyStruct::detJString, &m_detJ )
    ->setSizedFromParent( 1 )
    ->reference();
}

CellElementSubRegion::~CellElementSubRegion()
{
  // TODO Auto-generated destructor stub
}

void
CellElementSubRegion::CopyFromCellBlock( CellBlock * source )
{
  this->SetElementType( source->GetElementTypeString() );
  this->setNumNodesPerElement( source->numNodesPerElement() );
  this->setNumFacesPerElement( source->numFacesPerElement() );
  this->resize( source->size() );
  this->nodeList() = source->nodeList();

  arrayView1d< globalIndex const > const & sourceLocalToGlobal =
    source->localToGlobalMap();
  this->m_localToGlobalMap.resize( sourceLocalToGlobal.size() );
  for( localIndex i = 0; i < localToGlobalMap().size(); ++i )
  {
    this->m_localToGlobalMap[i] = sourceLocalToGlobal[i];
  }

  this->ConstructGlobalToLocalMap();
  source->forExternalProperties( [&]( dataRepository::WrapperBase * const wrapper ) {
    std::type_index typeIndex = std::type_index( wrapper->get_typeid() );
    rtTypes::ApplyArrayTypeLambda2(
      rtTypes::typeID( typeIndex ),
      true,
      [&]( auto type, auto GEOSX_UNUSED_PARAM( baseType ) ) {
        using fieldType = decltype( type );
        dataRepository::Wrapper< fieldType > & field =
          dataRepository::Wrapper< fieldType >::cast( *wrapper );
        const fieldType & fieldref = field.reference();
        this->registerWrapper(
          wrapper->getName(),
          &const_cast< fieldType & >( fieldref ) );  //TODO remove const_cast
        //      auto const & origFieldRef = field.reference();
        //      fieldType & fieldRef = this->registerWrapper<fieldType>( wrapper->getName() )->reference();
        //      fieldRef.resize( origFieldRef.size() );
      } );
  } );
}

void
CellElementSubRegion::ConstructSubRegionFromFaceSet(
  FaceManager const * const faceManager,
  string const & setName )
{
  SortedArrayView< localIndex const > const & targetSet =
    faceManager->sets().getReference< SortedArray< localIndex > >( setName );
  m_toFacesRelation.resize( 0, 2 );
  this->resize( targetSet.size() );
}

void
CellElementSubRegion::ViewPackingExclusionList(
  SortedArray< localIndex > & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList( exclusionList );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::nodeListString ) );
  //  exclusionList.insert(this->getWrapperIndex(this->viewKeys.edgeListString));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::faceListString ) );
}

localIndex
CellElementSubRegion::PackUpDownMapsSize(
  arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate< false >( junk, packList );
}

localIndex
CellElementSubRegion::PackUpDownMaps(
  buffer_unit_type *& buffer,
  arrayView1d< localIndex const > const & packList ) const
{
  return PackUpDownMapsPrivate< true >( buffer, packList );
}

template< bool DOPACK >
localIndex
CellElementSubRegion::PackUpDownMapsPrivate(
  buffer_unit_type *& buffer,
  arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           nodeList().Base().toViewConst(),
                                           m_unmappedGlobalIndicesInNodelist,
                                           packList,
                                           this->localToGlobalMap(),
                                           nodeList().RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           faceList().Base().toViewConst(),
                                           m_unmappedGlobalIndicesInFacelist,
                                           packList,
                                           this->localToGlobalMap(),
                                           faceList().RelatedObjectLocalToGlobal() );

  return packedSize;
}

localIndex
CellElementSubRegion::UnpackUpDownMaps(
  buffer_unit_type const *& buffer,
  localIndex_array & packList,
  bool const GEOSX_UNUSED_PARAM( overwriteUpMaps ),
  bool const GEOSX_UNUSED_PARAM( overwriteDownMaps ) )
{
  localIndex unPackedSize = 0;
  unPackedSize += bufferOps::Unpack( buffer,
                                     nodeList().Base().toView(),
                                     packList,
                                     m_unmappedGlobalIndicesInNodelist,
                                     this->globalToLocalMap(),
                                     nodeList().RelatedObjectGlobalToLocal() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     faceList().Base(),
                                     packList,
                                     m_unmappedGlobalIndicesInFacelist,
                                     this->globalToLocalMap(),
                                     faceList().RelatedObjectGlobalToLocal() );

  return unPackedSize;
}

void
CellElementSubRegion::FixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::FixUpDownMaps( nodeList(),
                                    m_unmappedGlobalIndicesInNodelist,
                                    clearIfUnmapped );

  ObjectManagerBase::FixUpDownMaps( faceList(),
                                    m_unmappedGlobalIndicesInFacelist,
                                    clearIfUnmapped );
}

} /* namespace geosx */
