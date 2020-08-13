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
 * @file ElementRegionManager.hpp
 */

#ifndef GEOSX_MESH_ELEMENTREGIONMANAGER_HPP
#define GEOSX_MESH_ELEMENTREGIONMANAGER_HPP

#include "CellBlock.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "CellElementRegion.hpp"
#include "CellElementSubRegion.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "dataRepository/ReferenceWrapper.hpp"
#include "FaceElementRegion.hpp"
#include "EmbeddedSurfaceRegion.hpp"
#include "fileIO/schema/schemaUtilities.hpp"
#include "WellElementRegion.hpp"
#include "EmbeddedSurfaceRegion.hpp"

namespace geosx
{
class MeshManager;

/**
 * @class ElementRegionManager
 * @brief The ElementRegionManager class provides an interface to ObjectManagerBase in order to manage ElementRegion
 * data
 */
class ElementRegionManager : public ObjectManagerBase
{
public:
  /**
   * Limit on max number of nodes for each element
   */
  constexpr static int maxNumNodesPerElem = 8;

  /**
   * @brief The ElementViewAccessor at the ElementRegionManager level is a 2D array of  VIEWTYPE.
   * @tparam VIEWTYPE data type
   */
  template< typename VIEWTYPE >
  using ElementViewAccessor = array1d< array1d< VIEWTYPE > >;

  /**
   * @brief The ElementViewAccessor at the ElementRegionManager level is a 2D array of ReferenceWrapper around VIEWTYPE.
   * @tparam VIEWTYPE data type
   */
  template< typename VIEWTYPE >
  using ElementReferenceAccessor = array1d< array1d< ReferenceWrapper< VIEWTYPE > > >;

  /**
   * @brief The MaterialViewAccessor at the ElementRegionManager level is a 3D array of VIEWTYPE.
   * @tparam VIEWTYPE data type
   * var[elementRegionIndex][elementSubRegionIndex][materialIndexInRegion]
   */
  template< typename VIEWTYPE >
  using MaterialViewAccessor = array1d< array1d< array1d< VIEWTYPE > > >;

  /**
   * @brief The ConstitutiveRelationAccessor at the ElementRegionManager level is a 3D array of CONSTITUTIVE_TYPE
   * @tparam CONSTITUTIVE_TYPE constitutive type
   */
  template< typename CONSTITUTIVE_TYPE >
  using ConstitutiveRelationAccessor =
    array1d< array1d< array1d< CONSTITUTIVE_TYPE * > > >;

  /**
   * @brief The function is to return the name of the ElementRegionManager in the object catalog
   * @return string that contains the catalog name used to register/lookup this class in  the object catalog
   */
  static const string
  CatalogName()
  {
    return "ZoneManager";
  }

  /**
   * @brief Virtual access to CatalogName()
   * @return string that contains the catalog name used to register/lookup this class in the object catalog
   */
  virtual const string
  getCatalogName() const override final
  {
    return ElementRegionManager::CatalogName();
  }

  /**
   * @brief Constructor.
   * @param [in] name the name of this ObjectManager
   * @param [in] parent the parent Group
   */
  ElementRegionManager( string const & name, Group * const parent );

  /**
   * @brief Destructor
   */
  virtual ~ElementRegionManager() override;

  /**
   * @brief Get the number of elements in all ElementSubRegions of type T.
   * @return number of elements
   */
  template< typename T = ElementSubRegionBase >
  localIndex
  getNumberOfElements() const
  {
    localIndex numElem = 0;
    this->forElementSubRegions< T >( [&]( ElementSubRegionBase const & cellBlock ) {
      numElem += cellBlock.size();
    } );
    return numElem;
  }

  //  void Initialize(  ){}

  /**
   * @brief Generate the mesh.
   * @param [in] cellBlockManager pointer to the CellBlockManager
   */
  void
  GenerateMesh( Group * const cellBlockManager );

  /**
   * @brief Generate the cell-to-edge map
   * @param [in] faceManager pointer to the FaceManager
   */
  void
  GenerateCellToEdgeMaps( FaceManager const * const faceManager );

  /**
   * @brief Generate the aggregates.
   * @param [in] faceManager pointer to the FaceManager
   * @param [in] nodeManager pointer to the NodeManager
   */
  void
  GenerateAggregates( FaceManager const * const faceManager,
                      NodeManager const * const nodeManager );

  /**
   * @brief Generate the wells.
   * @param [in] meshManager pointer to meshManager
   * @param [in] meshLevel pointer to meshLevel
   */
  void
  GenerateWells( MeshManager * const meshManager, MeshLevel * const meshLevel );

  /**
   * @brief Create a new ElementRegion object as a child of this group.
   * @param childKey catalog key of the new ElementRegion derived type to create
   * @param childName name of the new ElementRegion object
   * @return pointer to the created ElementRegion object
   */
  virtual Group *
  CreateChild( string const & childKey,
               string const & childName ) override;
  //  virtual void ReadXMLsub( xmlWrapper::xmlNode const & targetNode ) override;

  /**
   * @brief Expand any catalogs in the data structure
   */
  virtual void
  ExpandObjectCatalogs() override;

  /**
   * @brief Inform the schema generator of any deviations between the xml and GEOS data structures.
   * @param schemaRoot        XML node corresponding to the root
   * @param schemaParent      XML node for the parent node
   * @param documentationType type of XML schema generated
   */
  virtual void
  SetSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
                       xmlWrapper::xmlNode schemaParent,
                       integer documentationType ) override;

  using Group::resize;

  /**
   * @brief Set the number of elements for a set of element regions.
   * @param numElements list of the new element numbers
   * @param regionNames list of the element region names
   * @param elementTypes list of the element types
   */
  void
  resize( integer_array const & numElements,
          string_array const & regionNames,
          string_array const & elementTypes );

  /**
   * @brief Set the maximum local and global index.
   */
  void
  SetMaxGlobalIndex();

  /**
   * @brief Get a collection of element regions
   * @return reference to immutable subGroupMap
   */
  subGroupMap const &
  GetRegions() const
  {
    return this->GetGroup( groupKeyStruct::elementRegionsGroup )->GetSubGroups();
  }

  /**
   * @brief Get a collection of element regions.
   * @return reference to mutable subGroupMap
   */
  subGroupMap &
  GetRegions()
  {
    return this->GetGroup( groupKeyStruct::elementRegionsGroup )->GetSubGroups();
  }

  /**
   * @brief Get a element region.
   * @param regionName name of element region
   * @return pointer to const ElementRegionBase
   */
  template< typename T = ElementRegionBase >
  T const *
  GetRegion( string const & regionName ) const
  {
    return this->GetGroup( groupKeyStruct::elementRegionsGroup )
      ->GetGroup< T >( regionName );
  }

  /**
   * @brief Get a element region.
   * @param regionName name of element region
   * @return pointer to ElementRegionBase
   */
  template< typename T = ElementRegionBase >
  T *
  GetRegion( string const & regionName )
  {
    return this->GetGroup( groupKeyStruct::elementRegionsGroup )
      ->GetGroup< T >( regionName );
  }

  /**
   * @brief This is a const function to get a element region.
   * @param index index of element region
   * @return pointer to const ElementRegionBase
   */
  template< typename T = ElementRegionBase >
  T const *
  GetRegion( localIndex const index ) const
  {
    return this->GetGroup( groupKeyStruct::elementRegionsGroup )->GetGroup< T >( index );
  }

  /**
   * @brief This is a function to get a element region.
   * @param index index of element region
   * @return pointer to ElementRegionBase
   */
  template< typename T = ElementRegionBase >
  T *
  GetRegion( localIndex const index )
  {
    return this->GetGroup( groupKeyStruct::elementRegionsGroup )->GetGroup< T >( index );
  }

  /**
   * @brief Get number of the regions.
   * @return number of the regions
   */
  localIndex
  numRegions() const
  {
    return this->GetGroup( groupKeyStruct::elementRegionsGroup )->GetSubGroups().size();
  }

  /**
   * @brief Get number of the cell blocks.
   * @return number of the cell blocks
   */
  localIndex
  numCellBlocks() const;

  /**
   * @brief This function is used to launch kernel function over all the element regions with region type =
   * ElementRegionBase.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename REGIONTYPE = ElementRegionBase, typename... REGIONTYPES, typename LAMBDA >
  void
  forElementRegions( LAMBDA && lambda )
  {
    Group * const elementRegions =
      this->GetGroup( groupKeyStruct::elementRegionsGroup );
    elementRegions->forSubGroups< REGIONTYPE, REGIONTYPES... >(
      std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over all the element regions with region type =
   * ElementRegionBase.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename REGIONTYPE = ElementRegionBase, typename... REGIONTYPES, typename LAMBDA >
  void
  forElementRegions( LAMBDA && lambda ) const
  {
    Group const * const elementRegions =
      this->GetGroup( groupKeyStruct::elementRegionsGroup );
    elementRegions->forSubGroups< REGIONTYPE, REGIONTYPES... >(
      std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the target element regions with region type =
   * ElementRegionBase.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename REGIONTYPE = ElementRegionBase,
            typename... REGIONTYPES,
            typename LOOKUP_CONTAINER,
            typename LAMBDA >
  void
  forElementRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    Group * const elementRegions =
      this->GetGroup( groupKeyStruct::elementRegionsGroup );
    elementRegions->forSubGroups< REGIONTYPE, REGIONTYPES... >(
      targetRegions,
      std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over the target element regions with region type =
   * ElementRegionBase.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename REGIONTYPE = ElementRegionBase,
            typename... REGIONTYPES,
            typename LOOKUP_CONTAINER,
            typename LAMBDA >
  void
  forElementRegions( LOOKUP_CONTAINER const & targetRegions,
                     LAMBDA && lambda ) const
  {
    Group const * const elementRegions =
      this->GetGroup( groupKeyStruct::elementRegionsGroup );
    elementRegions->forSubGroups< REGIONTYPE, REGIONTYPES... >(
      targetRegions,
      std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over all the types of element regions.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void
  forElementRegionsComplete( LAMBDA lambda ) const
  {
    forElementRegionsComplete< CellElementRegion,
                               FaceElementRegion,
                               EmbeddedSurfaceRegion,
                               WellElementRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over all the types of element regions.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void
  forElementRegionsComplete( LAMBDA lambda )
  {
    forElementRegionsComplete< CellElementRegion,
                               FaceElementRegion,
                               EmbeddedSurfaceRegion,
                               WellElementRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over all the element regions that can be casted to one of
   * the specified region types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename REGIONTYPE, typename... REGIONTYPES, typename LAMBDA >
  void
  forElementRegionsComplete( LAMBDA lambda )
  {
    for( localIndex er = 0; er < this->numRegions(); ++er )
    {
      ElementRegionBase & elementRegion = *this->GetRegion( er );

      Group::applyLambdaToContainer< REGIONTYPE, REGIONTYPES... >(
        elementRegion,
        [&]( auto & castedRegion ) { lambda( er, castedRegion ); } );
    }
  }

  /**
   * @brief This const function is used to launch kernel function over all the element regions that can be casted to one
   * of the specified region types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename REGIONTYPE, typename... REGIONTYPES, typename LAMBDA >
  void
  forElementRegionsComplete( LAMBDA lambda ) const
  {
    for( localIndex er = 0; er < this->numRegions(); ++er )
    {
      ElementRegionBase const & elementRegion = *this->GetRegion( er );

      Group::applyLambdaToContainer< REGIONTYPE, REGIONTYPES... >(
        elementRegion,
        [&]( auto const & castedRegion ) { lambda( er, castedRegion ); } );
    }
  }

  /**
   * @brief This const function is used to launch kernel function over the specified target element regions.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void
  forElementRegionsComplete( LOOKUP_CONTAINER const & targetRegions,
                             LAMBDA lambda ) const
  {
    forElementRegionsComplete< CellElementRegion,
                               FaceElementRegion,
                               EmbeddedSurfaceRegion,
                               WellElementRegion >( targetRegions,
                                                    std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the specified target element regions.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void
  forElementRegionsComplete( LOOKUP_CONTAINER const & targetRegions,
                             LAMBDA lambda )
  {
    forElementRegionsComplete< CellElementRegion,
                               FaceElementRegion,
                               EmbeddedSurfaceRegion,
                               WellElementRegion >( targetRegions,
                                                    std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the specified target element regions with region type =
   * specified element region types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename REGIONTYPE, typename... REGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void
  forElementRegionsComplete( LOOKUP_CONTAINER const & targetRegions,
                             LAMBDA lambda )
  {
    forElementRegions< REGIONTYPE, REGIONTYPES... >(
      targetRegions,
      [&]( localIndex const targetIndex, auto & elementRegion ) {
        lambda( targetIndex, elementRegion.getIndexInParent(), elementRegion );
      } );
  }

  /**
   * @brief This const function is used to launch kernel function over the specified target element regions with region
   * type = specified element region types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename REGIONTYPE, typename... REGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void
  forElementRegionsComplete( LOOKUP_CONTAINER const & targetRegions,
                             LAMBDA lambda ) const
  {
    forElementRegions< REGIONTYPE, REGIONTYPES... >(
      targetRegions,
      [&]( localIndex const targetIndex, auto const & elementRegion ) {
        lambda( targetIndex, elementRegion.getIndexInParent(), elementRegion );
      } );
  }

  /**
   * @brief This function is used to launch kernel function over the element subregions of all the subregion types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void
  forElementSubRegions( LAMBDA && lambda )
  {
    forElementSubRegions< CellElementSubRegion,
                          FaceElementSubRegion,
                          EmbeddedSurfaceSubRegion,
                          WellElementSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over the element subregions of all the subregion
   * types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void
  forElementSubRegions( LAMBDA && lambda ) const
  {
    forElementSubRegions< CellElementSubRegion,
                          FaceElementSubRegion,
                          EmbeddedSurfaceSubRegion,
                          WellElementSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the specified target element subregions.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void
  forElementSubRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forElementSubRegions< CellElementSubRegion,
                          FaceElementSubRegion,
                          EmbeddedSurfaceSubRegion,
                          WellElementSubRegion >( targetRegions,
                                                  std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over the specified target element subregions.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void
  forElementSubRegions( LOOKUP_CONTAINER const & targetRegions,
                        LAMBDA && lambda ) const
  {
    forElementSubRegions< CellElementSubRegion,
                          FaceElementSubRegion,
                          EmbeddedSurfaceSubRegion,
                          WellElementSubRegion >( targetRegions,
                                                  std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the element subregions of the specified subregion
   * types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE, typename... SUBREGIONTYPES, typename LAMBDA >
  void
  forElementSubRegions( LAMBDA && lambda )
  {
    forElementSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >(
      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const,
                                                   localIndex const,
                                                   ElementRegionBase &,
                                                   auto & subRegion ) {
        lambda( subRegion );
      } );
  }

  /**
   * @brief This const function is used to launch kernel function over the element subregions of the specified subregion
   * types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE, typename... SUBREGIONTYPES, typename LAMBDA >
  void
  forElementSubRegions( LAMBDA && lambda ) const
  {
    forElementSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >(
      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const,
                                                   localIndex const,
                                                   ElementRegionBase const &,
                                                   auto const & subRegion ) {
        lambda( subRegion );
      } );
  }

  /**
   * @brief This function is used to launch kernel function over the specified target element subregions with the
   * specified subregion types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE,
            typename... SUBREGIONTYPES,
            typename LOOKUP_CONTAINER,
            typename LAMBDA >
  void
  forElementSubRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forElementSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >(
      targetRegions,
      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const targetIndex,
                                                   localIndex const,
                                                   localIndex const,
                                                   ElementRegionBase &,
                                                   auto & subRegion ) {
        lambda( targetIndex, subRegion );
      } );
  }

  /**
   * @brief This const function is used to launch kernel function over the specified target element subregions with the
   * specified subregion types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE,
            typename... SUBREGIONTYPES,
            typename LOOKUP_CONTAINER,
            typename LAMBDA >
  void
  forElementSubRegions( LOOKUP_CONTAINER const & targetRegions,
                        LAMBDA && lambda ) const
  {
    forElementSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >(
      targetRegions,
      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const targetIndex,
                                                   localIndex const,
                                                   localIndex const,
                                                   ElementRegionBase const &,
                                                   auto const & subRegion ) {
        lambda( targetIndex, subRegion );
      } );
  }

  /**
   * @brief This const function is used to launch kernel function over the element subregions of all subregion types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void
  forElementSubRegionsComplete( LAMBDA && lambda ) const
  {
    forElementSubRegionsComplete< CellElementSubRegion,
                                  FaceElementSubRegion,
                                  EmbeddedSurfaceSubRegion,
                                  WellElementSubRegion >(
      std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the element subregions of all subregion types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void
  forElementSubRegionsComplete( LAMBDA && lambda )
  {
    forElementSubRegionsComplete< CellElementSubRegion,
                                  FaceElementSubRegion,
                                  EmbeddedSurfaceSubRegion,
                                  WellElementSubRegion >(
      std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the specified target element subregions
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void
  forElementSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions,
                                LAMBDA && lambda )
  {
    forElementSubRegionsComplete< CellElementSubRegion,
                                  FaceElementSubRegion,
                                  EmbeddedSurfaceSubRegion,
                                  WellElementSubRegion >(
      targetRegions,
      std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over the specified target element subregions
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void
  forElementSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions,
                                LAMBDA && lambda ) const
  {
    forElementSubRegionsComplete< CellElementSubRegion,
                                  FaceElementSubRegion,
                                  EmbeddedSurfaceSubRegion,
                                  WellElementSubRegion >(
      targetRegions,
      std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over all the element subregions that can be casted to one of
   * the specified subregion types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE, typename... SUBREGIONTYPES, typename LAMBDA >
  void
  forElementSubRegionsComplete( LAMBDA && lambda )
  {
    for( localIndex er = 0; er < this->numRegions(); ++er )
    {
      ElementRegionBase & elementRegion = *this->GetRegion( er );

      for( localIndex esr = 0; esr < elementRegion.numSubRegions(); ++esr )
      {
        ElementSubRegionBase & subRegion = *elementRegion.GetSubRegion( esr );

        Group::applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >(
          subRegion,
          [&]( auto & castedSubRegion ) {
            lambda( er, esr, elementRegion, castedSubRegion );
          } );
      }
    }
  }

  /**
   * @brief This const function is used to launch kernel function over all the element subregions that can be casted to
   * one of the specified subregion types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE, typename... SUBREGIONTYPES, typename LAMBDA >
  void
  forElementSubRegionsComplete( LAMBDA && lambda ) const
  {
    for( localIndex er = 0; er < this->numRegions(); ++er )
    {
      ElementRegionBase const & elementRegion = *this->GetRegion( er );

      for( localIndex esr = 0; esr < elementRegion.numSubRegions(); ++esr )
      {
        ElementSubRegionBase const & subRegion = *elementRegion.GetSubRegion( esr );

        Group::applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >(
          subRegion,
          [&]( auto const & castedSubRegion ) {
            lambda( er, esr, elementRegion, castedSubRegion );
          } );
      }
    }
  }

  /**
   * @brief This function is used to launch kernel function over the specified target element subregions that can be
   * casted to one of the specified subregion types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE,
            typename... SUBREGIONTYPES,
            typename LOOKUP_CONTAINER,
            typename LAMBDA >
  void
  forElementSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions,
                                LAMBDA && lambda )
  {
    forElementRegions(
      targetRegions,
      [&]( localIndex const targetIndex, ElementRegionBase & elementRegion ) {
        localIndex const er = elementRegion.getIndexInParent();

        for( localIndex esr = 0; esr < elementRegion.numSubRegions(); ++esr )
        {
          ElementSubRegionBase & subRegion = *elementRegion.GetSubRegion( esr );

          Group::applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >(
            subRegion,
            [&]( auto & castedSubRegion ) {
              lambda( targetIndex, er, esr, elementRegion, castedSubRegion );
            } );
        }
      } );
  }

  /**
   * @brief This const function is used to launch kernel function over the specified target element subregions that can
   * be casted to one of the specified subregion types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE,
            typename... SUBREGIONTYPES,
            typename LOOKUP_CONTAINER,
            typename LAMBDA >
  void
  forElementSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions,
                                LAMBDA && lambda ) const
  {
    forElementRegions(
      targetRegions,
      [&]( localIndex const targetIndex, ElementRegionBase const & elementRegion ) {
        localIndex const er = elementRegion.getIndexInParent();

        for( localIndex esr = 0; esr < elementRegion.numSubRegions(); ++esr )
        {
          ElementSubRegionBase const & subRegion =
            *elementRegion.GetSubRegion( esr );

          Group::applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >(
            subRegion,
            [&]( auto const & castedSubRegion ) {
              lambda( targetIndex, er, esr, elementRegion, castedSubRegion );
            } );
        }
      } );
  }

  /**
   * @brief This is a const function to construct a ElementViewAccessor to access the data registered on the mesh.
   * @tparam VIEWTYPE data type
   * @param name view name of the data
   * @param neighborName neighbor data name
   * @return ElementViewAccessor that contains VIEWTYPE data
   */
  template< typename VIEWTYPE, typename LHS = VIEWTYPE >
  ElementViewAccessor< LHS >
  ConstructViewAccessor(
    string const & name,
    string const & neighborName = string() ) const;

  /**
   * @brief This is a function to construct a ElementViewAccessor to access the data registered on the mesh.
   * @tparam VIEWTYPE data type
   * @param name view name of the data
   * @param neighborName neighbor data name
   * @return ElementViewAccessor that contains VIEWTYPE data
   */
  template< typename VIEWTYPE, typename LHS = VIEWTYPE >
  ElementViewAccessor< LHS >
  ConstructViewAccessor(
    string const & name,
    string const & neighborName = string() );

  /**
   * @brief This is a function to construct a ElementViewAccessor to access array data registered on the mesh.
   * @tparam T data type
   * @tparam NDIM number of array dimensions
   * @param name view name of the data
   * @param neighborName neighbor data name
   * @return ElementViewAccessor that contains ArrayView<T const, NDIM> of data
   */
  template< typename T, int NDIM >
  ElementViewAccessor< ArrayView< T const, NDIM > >
  ConstructArrayViewAccessor(
    string const & name,
    string const & neighborName = string() ) const;

  /**
   * @brief This is a const function to construct a ElementViewAccessor to access the data registered on the mesh.
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @param neighborName neighbor data name
   * @return ElementViewAccessor that contains pointers to wrapped VIEWTYPE data
   */
  template< typename VIEWTYPE >
  ElementViewAccessor< ReferenceWrapper< VIEWTYPE > >
  ConstructReferenceAccessor(
    string const & viewName,
    string const & neighborName = string() ) const;

  /**
   * @brief This is a function to construct a ElementViewAccessor to access the data registered on the mesh.
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @param neighborName neighbor data name
   * @return ElementViewAccessor that contains pointers to wrapped VIEWTYPE data
   */
  template< typename VIEWTYPE >
  ElementViewAccessor< ReferenceWrapper< VIEWTYPE > >
  ConstructReferenceAccessor(
    string const & viewName,
    string const & neighborName = string() );

  /**
   * @brief This is a const function to construct a MaterialViewAccessor to access the material data.
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @param cm pointer to ConstitutiveManager
   * @return MaterialViewAccessor that contains VIEWTYPE data
   */
  template< typename VIEWTYPE, typename LHS = VIEWTYPE >
  MaterialViewAccessor< LHS >
  ConstructFullMaterialViewAccessor(
    string const & viewName,
    constitutive::ConstitutiveManager const * const cm ) const;

  /**
   * @brief This is a function to construct a MaterialViewAccessor to access the material data.
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @param cm pointer to ConstitutiveManager
   * @return MaterialViewAccessor that contains VIEWTYPE data
   */
  template< typename VIEWTYPE, typename LHS = VIEWTYPE >
  MaterialViewAccessor< LHS >
  ConstructFullMaterialViewAccessor(
    string const & viewName,
    constitutive::ConstitutiveManager const * const cm );

  /**
   * @brief This is a const function to construct a MaterialViewAccessor to access the material data for specified
   * regions/materials.
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @param regionNames list of region names
   * @param materialNames list of corresponding material names
   * @param allowMissingViews flag to indicate whether it is allowed to miss the specified material data in material
   * list
   * @return ElementViewAccessor that contains VIEWTYPE data
   */
  template< typename VIEWTYPE, typename LHS = VIEWTYPE >
  ElementViewAccessor< LHS >
  ConstructMaterialViewAccessor(
    string const & viewName,
    arrayView1d< string const > const & regionNames,
    arrayView1d< string const > const & materialNames,
    bool const allowMissingViews = false ) const;

  /**
   * @brief This is a function to construct a MaterialViewAccessor to access the material data for specified
   * regions/materials.
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @param regionNames list of region names
   * @param materialNames list of corresponding material names
   * @param allowMissingViews flag to indicate whether it is allowed to miss the specified material data in material
   * list
   * @return ElementViewAccessor that contains VIEWTYPE data
   */
  template< typename VIEWTYPE, typename LHS = VIEWTYPE >
  ElementViewAccessor< LHS >
  ConstructMaterialViewAccessor(
    string const & viewName,
    arrayView1d< string const > const & regionNames,
    arrayView1d< string const > const & materialNames,
    bool const allowMissingViews = false );

  /**
   * @brief Construct a view accessor for material data, assuming array as storage type
   * @tparam T underlying data type
   * @tparam NDIM number of array dimensions
   * @param viewName view name of the data
   * @param regionNames list of region names
   * @param materialNames list of corresponding material names
   * @param allowMissingViews flag to indicate whether it is allowed to miss the specified material data in material list
   * @return MaterialViewAccessor that contains the data views
   */
  template< typename T, int NDIM >
  ElementViewAccessor< ArrayView< T const, NDIM > >
  ConstructMaterialArrayViewAccessor(
    string const & viewName,
    arrayView1d< string const > const & regionNames,
    arrayView1d< string const > const & materialNames,
    bool const allowMissingViews = false ) const;

  /**
   * @brief Construct a ConstitutiveRelationAccessor.
   * @tparam CONSTITUTIVE_TYPE constitutive type
   * @param cm pointer to ConstitutiveManager
   * @return ConstitutiveRelationAccessor
   */
  template< typename CONSTITUTIVE_TYPE >
  ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE >
  ConstructFullConstitutiveAccessor(
    constitutive::ConstitutiveManager const * const cm ) const;

  /**
   * @brief Construct a ConstitutiveRelationAccessor.
   * @tparam CONSTITUTIVE_TYPE constitutive type
   * @param cm pointer to ConstitutiveManager
   * @return ConstitutiveRelationAccessor
   */
  template< typename CONSTITUTIVE_TYPE >
  ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE >
  ConstructFullConstitutiveAccessor(
    constitutive::ConstitutiveManager const * const cm );

  using Group::Pack;
  using Group::PackSize;
  using ObjectManagerBase::PackGlobalMaps;
  using ObjectManagerBase::PackGlobalMapsSize;
  using ObjectManagerBase::PackUpDownMaps;
  using ObjectManagerBase::PackUpDownMapsSize;
  using ObjectManagerBase::UnpackGlobalMaps;
  using ObjectManagerBase::UnpackUpDownMaps;

  /**
   * @brief Get the buffer size needed to pack a list of wrappers.
   * @param wrapperNames list of wrapper names
   * @param packList list of indices to pack
   * @return the size of the buffer required to pack the wrappers
   */
  int
  PackSize( string_array const & wrapperNames,
            ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Pack a list of wrappers to a buffer.
   * @param buffer pointer to the buffer to be packed
   * @param wrapperNames list of wrapper names
   * @param packList list of indices to pack
   * @return the size of data packed to the buffer
   */
  int
  Pack( buffer_unit_type *& buffer,
        string_array const & wrapperNames,
        ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /// @copydoc dataRepository::Group::Unpack
  using ObjectManagerBase::Unpack;

  /**
   * @brief Unpack a buffer.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to unpack
   * @return the size of data unpacked
   */
  int
  Unpack( buffer_unit_type const *& buffer,
          ElementViewAccessor< arrayView1d< localIndex > > & packList );

  /**
   * @brief Unpack a buffer.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to unpack
   * @return the size of data unpacked.
   */
  int
  Unpack( buffer_unit_type const *& buffer,
          ElementReferenceAccessor< array1d< localIndex > > & packList );

  /**
   * @brief Get the size of the buffer to be packed.
   * @param packList list of indices to pack
   * @return the size of the data packed
   */
  int
  PackGlobalMapsSize(
    ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Pack a buffer.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of the data packed
   */
  int
  PackGlobalMaps(
    buffer_unit_type *& buffer,
    ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Unpack a buffer.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to pack
   * @return the size of the data unpacked
   */
  int
  UnpackGlobalMaps(
    buffer_unit_type const *& buffer,
    ElementViewAccessor< ReferenceWrapper< localIndex_array > > & packList );

  /**
   * @brief Get the buffer size needed to pack element-to-node and element-to-face maps.
   * @param packList list of indices to pack
   * @return the size of data packed.
   */
  int
  PackUpDownMapsSize(
    ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Get the buffer size needed to pack element-to-node and element-to-face maps.
   * @param packList list of indices to pack
   * @return the size of data packed.
   */
  int
  PackUpDownMapsSize(
    ElementReferenceAccessor< array1d< localIndex > > const & packList ) const;

  /**
   * @brief Pack element-to-node and element-to-face maps.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of data packed.
   */
  int
  PackUpDownMaps(
    buffer_unit_type *& buffer,
    ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Pack element-to-node and element-to-face maps.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of data packed.
   */
  int
  PackUpDownMaps(
    buffer_unit_type *& buffer,
    ElementReferenceAccessor< array1d< localIndex > > const & packList ) const;

  /**
   * @brief Unpack element-to-node and element-to-face maps.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to pack
   * @param overwriteMap flag to indicate whether to overwrite the local map
   * @return the size of data packed.
   */
  int
  UnpackUpDownMaps( buffer_unit_type const *& buffer,
                    ElementReferenceAccessor< localIndex_array > & packList,
                    bool const overwriteMap );

  /**
   * @brief Group key associated with elementRegionsGroup
     struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
   */
  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    /// element regions group string key
    static constexpr auto elementRegionsGroup = "elementRegionsGroup";
  } m_ElementRegionManagerKeys;  ///< Element region manager keys

private:
  /**
   * @brief Pack a list of wrappers or get the buffer size needed to pack.
   * @param buffer pointer to the buffer to be packed
   * @param wrapperNames list of wrapper names
   * @param packList list of indices to pack
   * @return the size of the buffer required to pack the wrappers
   */
  template< bool DOPACK >
  int
  PackPrivate(
    buffer_unit_type *& buffer,
    string_array const & wrapperNames,
    ElementViewAccessor< arrayView1d< localIndex > > const & viewAccessor ) const;

  /**
   * @brief Pack a buffer or get the buffer size.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of the data packed
   */
  template< bool DOPACK >
  int
  PackGlobalMapsPrivate(
    buffer_unit_type *& buffer,
    ElementViewAccessor< arrayView1d< localIndex > > const & viewAccessor ) const;

  /**
   * @brief Pack element-to-node and element-to-face maps to a buffer or get the buffer size.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of the data packed
   */
  template< bool DOPACK, typename T >
  int
  PackUpDownMapsPrivate( buffer_unit_type *& buffer, T const & packList ) const;
  /**
   * @brief Unpack element-to-node and element-to-face maps.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to pack
   * @return the size of the data unpacked
   */
  template< typename T >
  int
  UnpackPrivate( buffer_unit_type const *& buffer, T & packList );

  /**
   * @brief Copy constructor.
   */
  ElementRegionManager( const ElementRegionManager & );

  /**
   * @brief Copy assignment operator.
   * @return reference to this object
   */
  ElementRegionManager &
  operator=( const ElementRegionManager & );
};

template< typename VIEWTYPE, typename LHS >
ElementRegionManager::ElementViewAccessor< LHS >
ElementRegionManager::ConstructViewAccessor( string const & viewName,
                                             string const & neighborName ) const
{
  ElementViewAccessor< LHS > viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg = 0; kReg < numRegions(); ++kReg )
  {
    ElementRegionBase const * const elemRegion = GetRegion( kReg );
    viewAccessor[kReg].resize( elemRegion->numSubRegions() );

    for( typename dataRepository::indexType kSubReg = 0;
         kSubReg < elemRegion->numSubRegions();
         ++kSubReg )
    {
      Group const * group = elemRegion->GetSubRegion( kSubReg );

      if( !neighborName.empty() )
      {
        group =
          group->GetGroup( ObjectManagerBase::groupKeyStruct::neighborDataString )
            ->GetGroup( neighborName );
      }

      if( group->hasWrapper( viewName ) &&
          group->getWrapperBase( viewName )->get_typeid() == typeid( VIEWTYPE ) )
      {
        viewAccessor[kReg][kSubReg] = group->getReference< VIEWTYPE >( viewName );
      }
    }
  }
  return viewAccessor;
}

template< typename VIEWTYPE, typename LHS >
ElementRegionManager::ElementViewAccessor< LHS >
ElementRegionManager::ConstructViewAccessor( string const & viewName,
                                             string const & neighborName )
{
  ElementViewAccessor< LHS > viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg = 0; kReg < numRegions(); ++kReg )
  {
    ElementRegionBase * const elemRegion = GetRegion( kReg );
    viewAccessor[kReg].resize( elemRegion->numSubRegions() );

    for( typename dataRepository::indexType kSubReg = 0;
         kSubReg < elemRegion->numSubRegions();
         ++kSubReg )
    {
      Group * group = elemRegion->GetSubRegion( kSubReg );

      if( !neighborName.empty() )
      {
        group =
          group->GetGroup( ObjectManagerBase::groupKeyStruct::neighborDataString )
            ->GetGroup( neighborName );
      }

      if( group->hasWrapper( viewName ) &&
          group->getWrapperBase( viewName )->get_typeid() == typeid( VIEWTYPE ) )
      {
        viewAccessor[kReg][kSubReg] = group->getReference< VIEWTYPE >( viewName );
      }
    }
  }
  return viewAccessor;
}

template< typename T, int NDIM >
ElementRegionManager::ElementViewAccessor< ArrayView< T const, NDIM > >
ElementRegionManager::ConstructArrayViewAccessor( string const & name,
                                                  string const & neighborName ) const
{
  return ConstructViewAccessor< Array< T, NDIM >, ArrayView< T const, NDIM > >(
    name,
    neighborName );
}

template< typename VIEWTYPE >
ElementRegionManager::ElementViewAccessor< ReferenceWrapper< VIEWTYPE > >
ElementRegionManager::ConstructReferenceAccessor( string const & viewName,
                                                  string const & neighborName ) const
{
  ElementViewAccessor< ReferenceWrapper< VIEWTYPE > > viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg = 0; kReg < numRegions(); ++kReg )
  {
    ElementRegionBase const * const elemRegion = GetRegion( kReg );
    viewAccessor[kReg].resize( elemRegion->numSubRegions() );

    for( typename dataRepository::indexType kSubReg = 0;
         kSubReg < elemRegion->numSubRegions();
         ++kSubReg )
    {
      Group const * group = elemRegion->GetSubRegion( kSubReg );

      if( !neighborName.empty() )
      {
        group =
          group->GetGroup( ObjectManagerBase::groupKeyStruct::neighborDataString )
            ->GetGroup( neighborName );
      }

      if( group->hasWrapper( viewName ) )
      {
        viewAccessor[kReg][kSubReg].set( group->getReference< VIEWTYPE >( viewName ) );
      }
    }
  }
  return viewAccessor;
}

template< typename VIEWTYPE >
ElementRegionManager::ElementViewAccessor< ReferenceWrapper< VIEWTYPE > >
ElementRegionManager::ConstructReferenceAccessor( string const & viewName,
                                                  string const & neighborName )
{
  ElementViewAccessor< ReferenceWrapper< VIEWTYPE > > viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg = 0; kReg < numRegions(); ++kReg )
  {
    ElementRegionBase * const elemRegion = GetRegion( kReg );
    viewAccessor[kReg].resize( elemRegion->numSubRegions() );

    for( typename dataRepository::indexType kSubReg = 0;
         kSubReg < elemRegion->numSubRegions();
         ++kSubReg )
    {
      Group * group = elemRegion->GetSubRegion( kSubReg );

      if( !neighborName.empty() )
      {
        group =
          group->GetGroup( ObjectManagerBase::groupKeyStruct::neighborDataString )
            ->GetGroup( neighborName );
      }

      if( group->hasWrapper( viewName ) )
      {
        viewAccessor[kReg][kSubReg].set( group->getReference< VIEWTYPE >( viewName ) );
      }
    }
  }
  return viewAccessor;
}

template< typename VIEWTYPE, typename LHS >
ElementRegionManager::MaterialViewAccessor< LHS >
ElementRegionManager::ConstructFullMaterialViewAccessor(
  string const & viewName,
  constitutive::ConstitutiveManager const * const cm ) const
{
  MaterialViewAccessor< LHS > accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg = 0; kReg < numRegions(); ++kReg )
  {
    ElementRegionBase const * const elemRegion = GetRegion( kReg );
    accessor[kReg].resize( elemRegion->numSubRegions() );

    for( localIndex kSubReg = 0; kSubReg < elemRegion->numSubRegions(); ++kSubReg )
    {
      ElementSubRegionBase const * const subRegion =
        elemRegion->GetSubRegion( kSubReg );
      dataRepository::Group const * const constitutiveGroup =
        subRegion->GetConstitutiveModels();

      accessor[kReg][kSubReg].resize( cm->numSubGroups() );

      for( localIndex matIndex = 0; matIndex < cm->numSubGroups(); ++matIndex )
      {
        string constitutiveName = cm->GetGroup( matIndex )->getName();
        dataRepository::Group const * const constitutiveRelation =
          constitutiveGroup->GetGroup( constitutiveName );
        if( constitutiveRelation != nullptr )
        {
          dataRepository::Wrapper< VIEWTYPE > const * const wrapper =
            constitutiveRelation->getWrapper< VIEWTYPE >( viewName );

          if( wrapper != nullptr )
          {
            accessor[kReg][kSubReg][matIndex] = wrapper->reference();
          }
        }
      }
    }
  }
  return accessor;
}

template< typename VIEWTYPE, typename LHS >
ElementRegionManager::MaterialViewAccessor< LHS >
ElementRegionManager::ConstructFullMaterialViewAccessor(
  string const & viewName,
  constitutive::ConstitutiveManager const * const cm )
{
  MaterialViewAccessor< LHS > accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg = 0; kReg < numRegions(); ++kReg )
  {
    ElementRegionBase * const elemRegion = GetRegion( kReg );
    accessor[kReg].resize( elemRegion->numSubRegions() );

    for( localIndex kSubReg = 0; kSubReg < elemRegion->numSubRegions(); ++kSubReg )
    {
      ElementSubRegionBase * const subRegion = elemRegion->GetSubRegion( kSubReg );
      dataRepository::Group * const constitutiveGroup =
        subRegion->GetConstitutiveModels();

      accessor[kReg][kSubReg].resize( cm->numSubGroups() );

      for( localIndex matIndex = 0; matIndex < cm->numSubGroups(); ++matIndex )
      {
        string constitutiveName = cm->GetGroup( matIndex )->getName();
        dataRepository::Group * const constitutiveRelation =
          constitutiveGroup->GetGroup( constitutiveName );
        if( constitutiveRelation != nullptr )
        {
          dataRepository::Wrapper< VIEWTYPE > * const wrapper =
            constitutiveRelation->getWrapper< VIEWTYPE >( viewName );

          if( wrapper != nullptr )
          {
            accessor[kReg][kSubReg][matIndex] = wrapper->reference();
          }
        }
      }
    }
  }
  return accessor;
}

template< typename VIEWTYPE, typename LHS >
ElementRegionManager::ElementViewAccessor< LHS >
ElementRegionManager::ConstructMaterialViewAccessor(
  string const & viewName,
  arrayView1d< string const > const & regionNames,
  arrayView1d< string const > const & materialNames,
  bool const allowMissingViews ) const
{
  GEOSX_ASSERT_EQ( regionNames.size(), materialNames.size() );
  ElementViewAccessor< LHS > accessor;

  // Resize the accessor to all regions and subregions
  accessor.resize( numRegions() );
  for( localIndex kReg = 0; kReg < numRegions(); ++kReg )
  {
    accessor[kReg].resize( GetRegion( kReg )->numSubRegions() );
  }

  subGroupMap const & regionMap = GetRegions();

  // Loop only over regions named and populate according to given material names
  for( localIndex k = 0; k < regionNames.size(); ++k )
  {
    localIndex const er = regionMap.getIndex( regionNames[k] );
    GEOSX_ERROR_IF_EQ_MSG( er,
                           subGroupMap::KeyIndex::invalid_index,
                           "Region not found: " << regionNames[k] );
    ElementRegionBase const & region = *GetRegion( er );

    region.forElementSubRegionsIndex(
      [&]( localIndex const esr, ElementSubRegionBase const & subRegion ) {
        dataRepository::Group const & constitutiveGroup =
          *subRegion.GetConstitutiveModels();
        dataRepository::Group const * const constitutiveRelation =
          constitutiveGroup.GetGroup( materialNames[k] );
        GEOSX_ERROR_IF( constitutiveRelation == nullptr,
                        "Material " << materialNames[k] << " not found in "
                                    << regionNames[k] << '/'
                                    << subRegion.getName() );
        dataRepository::Wrapper< VIEWTYPE > const * const wrapper =
          constitutiveRelation->getWrapper< VIEWTYPE >( viewName );
        GEOSX_ERROR_IF(
          !allowMissingViews && wrapper == nullptr,
          "Material " << materialNames[k] << " does not contain " << viewName );
        if( wrapper != nullptr )
        {
          accessor[er][esr] = wrapper->reference();
        }
      } );
  }
  return accessor;
}

template< typename VIEWTYPE, typename LHS >
ElementRegionManager::ElementViewAccessor< LHS >
ElementRegionManager::ConstructMaterialViewAccessor(
  string const & viewName,
  arrayView1d< string const > const & regionNames,
  arrayView1d< string const > const & materialNames,
  bool const allowMissingViews )
{
  GEOSX_ASSERT_EQ( regionNames.size(), materialNames.size() );
  ElementViewAccessor< LHS > accessor;

  // Resize the accessor to all regions and subregions
  accessor.resize( numRegions() );
  for( localIndex kReg = 0; kReg < numRegions(); ++kReg )
  {
    accessor[kReg].resize( GetRegion( kReg )->numSubRegions() );
  }

  subGroupMap const & regionMap = GetRegions();

  // Loop only over regions named and populate according to given material names
  for( localIndex k = 0; k < regionNames.size(); ++k )
  {
    localIndex const er = regionMap.getIndex( regionNames[k] );
    GEOSX_ERROR_IF_EQ_MSG( er,
                           subGroupMap::KeyIndex::invalid_index,
                           "Region not found: " << regionNames[k] );
    ElementRegionBase & region = *GetRegion( er );

    region.forElementSubRegionsIndex( [&]( localIndex const esr,
                                           ElementSubRegionBase & subRegion ) {
      dataRepository::Group & constitutiveGroup =
        *subRegion.GetConstitutiveModels();
      dataRepository::Group * const constitutiveRelation =
        constitutiveGroup.GetGroup( materialNames[k] );
      GEOSX_ERROR_IF( constitutiveRelation == nullptr,
                      "Material " << materialNames[k] << " not found in "
                                  << regionNames[k] << '/' << subRegion.getName() );
      dataRepository::Wrapper< VIEWTYPE > * const wrapper =
        constitutiveRelation->getWrapper< VIEWTYPE >( viewName );
      GEOSX_ERROR_IF(
        !allowMissingViews && wrapper == nullptr,
        "Material " << materialNames[k] << " does not contain " << viewName );
      if( wrapper != nullptr )
      {
        accessor[er][esr] = wrapper->reference();
      }
    } );
  }
  return accessor;
}

template< typename T, int NDIM >
ElementRegionManager::ElementViewAccessor< ArrayView< T const, NDIM > >
ElementRegionManager::ConstructMaterialArrayViewAccessor(
  string const & viewName,
  arrayView1d< string const > const & regionNames,
  arrayView1d< string const > const & materialNames,
  bool const allowMissingViews ) const
{
  return ConstructMaterialViewAccessor< Array< T, NDIM >, ArrayView< T const, NDIM > >(
    viewName,
    regionNames,
    materialNames,
    allowMissingViews );
}

template< typename CONSTITUTIVE_TYPE >
ElementRegionManager::ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE >
ElementRegionManager::ConstructFullConstitutiveAccessor(
  constitutive::ConstitutiveManager const * const cm ) const
{
  ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE > accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg = 0; kReg < numRegions(); ++kReg )
  {
    ElementRegionBase const * const elemRegion = GetRegion( kReg );
    accessor[kReg].resize( elemRegion->numSubRegions() );

    for( localIndex kSubReg = 0; kSubReg < elemRegion->numSubRegions(); ++kSubReg )
    {
      ElementSubRegionBase const * const subRegion =
        elemRegion->GetSubRegion( kSubReg );
      dataRepository::Group const * const constitutiveGroup =
        subRegion->GetConstitutiveModels();
      accessor[kReg][kSubReg].resize( cm->numSubGroups() );

      for( localIndex matIndex = 0; matIndex < cm->numSubGroups(); ++matIndex )
      {
        string const constitutiveName = cm->GetGroup( matIndex )->getName();

        CONSTITUTIVE_TYPE * const constitutiveRelation =
          constitutiveGroup->GetGroup< CONSTITUTIVE_TYPE >( constitutiveName );
        if( constitutiveRelation != nullptr )
        {
          accessor[kReg][kSubReg][matIndex] = constitutiveRelation;
        }
      }
    }
  }
  return accessor;
}

template< typename CONSTITUTIVE_TYPE >
ElementRegionManager::ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE >
ElementRegionManager::ConstructFullConstitutiveAccessor(
  constitutive::ConstitutiveManager const * const cm )
{
  ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE > accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg = 0; kReg < numRegions(); ++kReg )
  {
    ElementRegionBase * const elemRegion = GetRegion( kReg );
    accessor[kReg].resize( elemRegion->numSubRegions() );

    for( localIndex kSubReg = 0; kSubReg < elemRegion->numSubRegions(); ++kSubReg )
    {
      ElementSubRegionBase * const subRegion = elemRegion->GetSubRegion( kSubReg );
      dataRepository::Group * const constitutiveGroup =
        subRegion->GetConstitutiveModels();
      accessor[kReg][kSubReg].resize( cm->numSubGroups() );

      for( localIndex matIndex = 0; matIndex < cm->numSubGroups(); ++matIndex )
      {
        string const constitutiveName = cm->GetGroup( matIndex )->getName();

        CONSTITUTIVE_TYPE * const constitutiveRelation =
          constitutiveGroup->GetGroup< CONSTITUTIVE_TYPE >( constitutiveName );
        if( constitutiveRelation != nullptr )
        {
          accessor[kReg][kSubReg][matIndex] = constitutiveRelation;
        }
      }
    }
  }
  return accessor;
}

}  // namespace geosx
#endif /* GEOSX_MESH_ELEMENTREGIONMANAGER_HPP */
