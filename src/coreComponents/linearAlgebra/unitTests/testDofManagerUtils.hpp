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
 * @file testDofManagerUtils.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UNITTESTS_TESTDOFMANAGERUTILS_HPP
#define GEOSX_LINEARALGEBRA_UNITTESTS_TESTDOFMANAGERUTILS_HPP

#include <gtest/gtest.h>

namespace geosx
{
namespace testing
{
/**
 * @brief Set up a problem from an xml input buffer
 * @param problemManager the target problem manager
 * @param xmlInput       the XML input string
 */
void
setupProblemFromXML( ProblemManager * const problemManager,
                     char const * const xmlInput )
{
  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult =
    xmlDocument.load_buffer( xmlInput, strlen( xmlInput ) );
  if( !xmlResult )
  {
    GEOSX_LOG_RANK_0( "XML parsed with errors!" );
    GEOSX_LOG_RANK_0( "Error description: " << xmlResult.description() );
    GEOSX_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }

  int mpiSize = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
  dataRepository::Group * commandLine =
    problemManager->GetGroup< dataRepository::Group >(
      problemManager->groupKeys.commandLine );
  commandLine
    ->registerWrapper< integer >( problemManager->viewKeys.xPartitionsOverride.Key() )
    ->setApplyDefaultValue( mpiSize );

  xmlWrapper::xmlNode xmlProblemNode = xmlDocument.child( "Problem" );
  problemManager->InitializePythonInterpreter();
  problemManager->ProcessInputFileRecursive( xmlProblemNode );

  // Open mesh levels
  DomainPartition * domain = problemManager->getDomainPartition();
  MeshManager * meshManager =
    problemManager->GetGroup< MeshManager >( problemManager->groupKeys.meshManager );
  meshManager->GenerateMeshLevels( domain );

  ElementRegionManager * elementManager =
    domain->getMeshBody( 0 )->getMeshLevel( 0 )->getElemManager();
  xmlWrapper::xmlNode topLevelNode =
    xmlProblemNode.child( elementManager->getName().c_str() );
  elementManager->ProcessInputFileRecursive( topLevelNode );
  elementManager->PostProcessInputRecursive();

  problemManager->ProblemSetup();
}

/**
 * @brief Get a region list from an input list (empty region lists to mean all regions in the mesh).
 * @param mesh  pointer to the mesh
 * @param input input list of region names (may be empty)
 * @return the list of region names (same as input unless input is empty)
 *
 * Mainly used to allow empty region lists to mean all regions.
 */
string_array
getRegions( MeshLevel const * const mesh,
            std::vector< string > const & input )
{
  string_array regions;
  if( !input.empty() )
  {
    regions.insert( 0, input.begin(), input.end() );
  }
  else
  {
    mesh->getElemManager()->forElementRegions( [&]( ElementRegionBase const & region ) {
      regions.emplace_back( region.getName() );
    } );
  }
  return regions;
}

namespace internal
{
template< DofManager::Location LOC >
struct testMeshHelper
{};

template<>
struct testMeshHelper< DofManager::Location::Node >
{
  static auto constexpr managerKey =
    MeshLevel::groupStructKeys::nodeManagerString;
  static auto constexpr elemMapKey =
    ElementSubRegionBase::viewKeyStruct::nodeListString;
  template< typename SUBREGION >
  using ElemToObjMap = typename SUBREGION::NodeMapType;
};

template<>
struct testMeshHelper< DofManager::Location::Face >
{
  static auto constexpr managerKey =
    MeshLevel::groupStructKeys::faceManagerString;
  static auto constexpr elemMapKey =
    ElementSubRegionBase::viewKeyStruct::faceListString;
  template< typename SUBREGION >
  using ElemToObjMap = typename SUBREGION::FaceMapType;
};

template< int USD >
localIndex
size1( arrayView2d< localIndex const, USD > const & map,
       localIndex const GEOSX_UNUSED_PARAM( i0 ) )
{
  return map.size( 1 );
}

localIndex
size1( arrayView1d< arrayView1d< localIndex const > const > const & map,
       localIndex const i0 )
{
  return map[i0].size();
}

localIndex
size1( ArrayOfArraysView< localIndex const > const & map,
       localIndex const i0 )
{
  return map.sizeOfArray( i0 );
}

template< DofManager::Location LOC >
struct forLocalObjectsImpl
{
  template< typename LAMBDA >
  static void
  f( MeshLevel const * const mesh,
     string_array const & regions,
     LAMBDA lambda )
  {
    using helper = testMeshHelper< LOC >;
    ObjectManagerBase const * const manager =
      mesh->GetGroup< ObjectManagerBase >( helper::managerKey );

    arrayView1d< integer const > ghostRank =
      manager->getReference< array1d< integer > >(
        ObjectManagerBase::viewKeyStruct::ghostRankString );

    array1d< bool > visited( ghostRank.size() );

    mesh->getElemManager()->forElementSubRegions(
      regions,
      [&]( localIndex const, auto const & subRegion ) {
        using MapType = typename helper::template ElemToObjMap<
          std::remove_reference_t< decltype( subRegion ) > >;

        typename MapType::ViewTypeConst const & elemToObjMap =
          subRegion.template getReference< MapType >( helper::elemMapKey );

        for( localIndex k = 0; k < subRegion.size(); ++k )
        {
          for( localIndex a = 0; a < size1( elemToObjMap, k ); ++a )
          {
            localIndex const idx = elemToObjMap[k][a];
            if( ghostRank[idx] < 0 && !std::exchange( visited[idx], true ) )
            {
              lambda( idx );
            }
          }
        }
      } );
  }
};

template<>
struct forLocalObjectsImpl< DofManager::Location::Elem >
{
  template< typename LAMBDA >
  static void
  f( MeshLevel const * const mesh,
     string_array const & regions,
     LAMBDA lambda )
  {
    // make a list of regions
    ElementRegionManager const * const elemManager = mesh->getElemManager();

    elemManager->forElementSubRegionsComplete< ElementSubRegionBase >(
      regions,
      [&]( localIndex const,
           localIndex const er,
           localIndex const esr,
           ElementRegionBase const &,
           ElementSubRegionBase const & subRegion ) {
        arrayView1d< integer const > ghostRank =
          subRegion.getReference< array1d< integer > >(
            ObjectManagerBase::viewKeyStruct::ghostRankString );

        for( localIndex ei = 0; ei < subRegion.size(); ++ei )
        {
          if( ghostRank[ei] < 0 )
          {
            lambda( std::array< localIndex, 3 > { er, esr, ei } );
          }
        }
      } );
  }
};

}  // namespace internal

/**
 * @brief Apply a lambda to all locally owned mesh objects in the mesh.
 * @tparam LOC type of mesh location (Node, Element, etc.)
 * @tparam LAMBDA type of lambda
 * @param mesh    pointer to the mesh
 * @param regions list of input region names to loop over
 * @param lambda  the lambda to apply
 */
template< DofManager::Location LOC, typename LAMBDA >
void
forLocalObjects( MeshLevel const * const mesh,
                 array1d< string > const & regions,
                 LAMBDA && lambda )
{
  internal::forLocalObjectsImpl< LOC >::template f( mesh,
                                                    regions,
                                                    std::forward< LAMBDA >( lambda ) );
}

/**
 * @brief Count the number of local objects in the mesh.
 * @tparam LOC type of mesh location (Node, Element, etc.)
 * @param mesh    pointer to the mesh
 * @param regions list of input region names to loop over
 * @return the number of locally owned objects (e.g. nodes)
 */
template< DofManager::Location LOC >
localIndex
countLocalObjects( MeshLevel const * const mesh,
                   array1d< string > const & regions )
{
  localIndex numLocal = 0;
  forLocalObjects< LOC >( mesh, regions, [&]( auto const & ) { ++numLocal; } );
  return numLocal;
}

/**
 * @brief Create a TPFA-type sparsity pattern.
 * @param mesh        pointer to the mesh
 * @param dofIndexKey the DofManager key for the dof index array
 * @param regions     list of region names to include (if empty, all regions are used)
 * @param numComp     number of components per cell
 * @param sparsity    the matrix to be populated, must be properly sized.
 */
template< typename MATRIX >
void
makeSparsityTPFA( MeshLevel const * const mesh,
                  string const & dofIndexKey,
                  string_array const & regions,
                  localIndex const numComp,
                  MATRIX & sparsity )
{
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  FaceManager const * const faceManager = mesh->getFaceManager();

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofIndex =
    elemManager->ConstructViewAccessor< array1d< globalIndex >,
                                        arrayView1d< globalIndex const > >(
      dofIndexKey );

  // Make a set of target region indices to check face fluxes.
  SortedArray< localIndex > regionSet;
  elemManager->forElementRegions(
    regions,
    [&]( localIndex const, ElementRegionBase const & region ) {
      regionSet.insert( region.getIndexInParent() );
    } );

  // prepare data for assembly loop
  FaceManager::ElemMapType const & faceToElem = faceManager->toElementRelation();
  localIndex const numElem = faceToElem.size( 1 );

  array1d< globalIndex > localDofIndex( numElem * numComp );
  array2d< real64 > localValues( numElem * numComp, numElem * numComp );
  localValues.setValues< serialPolicy >( 1.0 );

  // Loop over faces and assemble TPFA-style "flux" contributions
  forLocalObjects< DofManager::Location::Face >( mesh, regions, [&]( localIndex const kf ) {
    localIndex const er0 = faceToElem.m_toElementRegion[kf][0];
    localIndex const er1 = faceToElem.m_toElementRegion[kf][1];

    if( er0 >= 0 && er1 >= 0 && regionSet.contains( er0 ) && regionSet.contains( er1 ) )
    {
      for( localIndex ke = 0; ke < numElem; ++ke )
      {
        localIndex const er = faceToElem.m_toElementRegion[kf][ke];
        localIndex const esr = faceToElem.m_toElementSubRegion[kf][ke];
        localIndex const ei = faceToElem.m_toElementIndex[kf][ke];

        for( localIndex c = 0; c < numComp; ++c )
        {
          localDofIndex[ke * numComp + c] = elemDofIndex[er][esr][ei] + c;
        }
      }
      sparsity.insert( localDofIndex, localDofIndex, localValues );
    }
  } );
}

/**
 * @brief Populate a FEM-type sparsity pattern.
 * @param mesh        pointer to the mesh
 * @param dofIndexKey the DofManager key for the dof index array
 * @param regions     list of region names to include (if empty, all regions are used)
 * @param numComp     number of components per node
 * @param sparsity    the matrix to be populated, must be properly sized.
 */
template< typename MATRIX >
void
makeSparsityFEM( MeshLevel const * const mesh,
                 string const & dofIndexKey,
                 string_array const & regions,
                 localIndex const numComp,
                 MATRIX & sparsity )
{
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  NodeManager const * const nodeManager = mesh->getNodeManager();

  arrayView1d< globalIndex const > nodeDofIndex =
    nodeManager->getReference< array1d< globalIndex > >( dofIndexKey );

  // perform assembly loop over elements
  elemManager->forElementSubRegions(
    regions,
    [&]( localIndex const, auto const & subRegion ) {
      using NodeMapType = typename TYPEOFREF( subRegion )::NodeMapType;
      typename NodeMapType::ViewTypeConst const & nodeMap =
        subRegion.template getReference< NodeMapType >(
          ElementSubRegionBase::viewKeyStruct::nodeListString );

      localIndex const numNode = subRegion.numNodesPerElement();
      array1d< globalIndex > localDofIndex( numNode * numComp );
      array2d< real64 > localValues( numNode * numComp, numNode * numComp );
      localValues.setValues< serialPolicy >( 1.0 );

      for( localIndex k = 0; k < subRegion.size(); ++k )
      {
        for( localIndex a = 0; a < numNode; ++a )
        {
          for( localIndex c = 0; c < numComp; ++c )
          {
            localDofIndex[a * numComp + c] = nodeDofIndex[nodeMap[k][a]] + c;
          }
        }

        sparsity.insert( localDofIndex, localDofIndex, localValues );
      }
    } );
}

/**
 * @brief Populate a FEM/FVM coupling sparsity.
 * @param mesh            pointer to the mesh
 * @param dofIndexKeyNode the DofManager key for the node-based dof index array
 * @param dofIndexKeyElem the DofManager key for the element-based dof index array
 * @param regions         list of region names to include (if empty, all regions are used)
 * @param numCompNode     number of components per node
 * @param numCompElem     number of components per element
 * @param sparsity        the matrix to be populated, must be properly sized.
 */
template< typename MATRIX >
void
makeSparsityFEM_FVM( MeshLevel const * const mesh,
                     string const & dofIndexKeyNode,
                     string const & dofIndexKeyElem,
                     string_array const & regions,
                     localIndex const numCompNode,
                     localIndex const numCompElem,
                     MATRIX & sparsity )
{
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  NodeManager const * const nodeManager = mesh->getNodeManager();

  arrayView1d< globalIndex const > nodeDofIndex =
    nodeManager->getReference< array1d< globalIndex > >( dofIndexKeyNode );

  // perform assembly loop over elements
  elemManager->forElementSubRegions(
    regions,
    [&]( localIndex const, auto const & subRegion ) {
      using NodeMapType = typename TYPEOFREF( subRegion )::NodeMapType;
      typename NodeMapType::ViewTypeConst const & nodeMap =
        subRegion.template getReference< NodeMapType >(
          ElementSubRegionBase::viewKeyStruct::nodeListString );

      arrayView1d< globalIndex const > elemDofIndex =
        subRegion.template getReference< array1d< globalIndex > >( dofIndexKeyElem );

      localIndex const numNode = subRegion.numNodesPerElement();

      array1d< globalIndex > localNodeDofIndex( numNode * numCompNode );
      array1d< globalIndex > localElemDofIndex( numCompElem );
      array2d< real64 > localValues1( numNode * numCompNode, numCompElem );
      localValues1.setValues< serialPolicy >( 1.0 );
      array2d< real64 > localValues2( numCompElem, numNode * numCompNode );
      localValues2.setValues< serialPolicy >( 1.0 );

      for( localIndex k = 0; k < subRegion.size(); ++k )
      {
        for( localIndex c = 0; c < numCompElem; ++c )
        {
          localElemDofIndex[c] = elemDofIndex[k] + c;
        }
        for( localIndex a = 0; a < numNode; ++a )
        {
          for( localIndex c = 0; c < numCompNode; ++c )
          {
            localNodeDofIndex[a * numCompNode + c] =
              nodeDofIndex[nodeMap[k][a]] + c;
          }
        }

        sparsity.insert( localNodeDofIndex, localElemDofIndex, localValues1 );
        sparsity.insert( localElemDofIndex, localNodeDofIndex, localValues2 );
      }
    } );
}

/**
 * @brief Create a mass matrix-type sparsity pattern (diagonal)
 * @param mesh        pointer to the mesh
 * @param dofIndexKey the DofManager key for the dof index array
 * @param regions     list of region names to include (if empty, all regions are used)
 * @param numComp     number of components per cell
 * @param sparsity    the matrix to be populated
 */
template< typename MATRIX >
void
makeSparsityMass( MeshLevel const * const mesh,
                  string const & dofIndexKey,
                  string_array const & regions,
                  localIndex const numComp,
                  MATRIX & sparsity )
{
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofIndex =
    elemManager->ConstructViewAccessor< array1d< globalIndex >,
                                        arrayView1d< globalIndex const > >(
      dofIndexKey );

  array1d< globalIndex > localDofIndex( numComp );
  array2d< real64 > localValues( numComp, numComp );
  localValues.setValues< serialPolicy >( 1.0 );

  forLocalObjects< DofManager::Location::Elem >( mesh, regions, [&]( auto const & idx ) {
    for( localIndex c = 0; c < numComp; ++c )
    {
      localDofIndex[c] = elemDofIndex[idx[0]][idx[1]][idx[2]] + c;
    }
    sparsity.insert( localDofIndex, localDofIndex, localValues );
  } );
}

/**
 * @brief Create a flux sparsity pattern (face-based dofs coupled via elements)
 * @param mesh        pointer to the mesh
 * @param dofIndexKey the DofManager key for the dof index array
 * @param regions     list of region names to include (if empty, all regions are used)
 * @param numComp     number of components per cell
 * @param sparsity    the matrix to be populated
 */
template< typename MATRIX >
void
makeSparsityFlux( MeshLevel const * const mesh,
                  string const & dofIndexKey,
                  string_array const & regions,
                  localIndex const numComp,
                  MATRIX & sparsity )
{
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  FaceManager const * const faceManager = mesh->getFaceManager();

  arrayView1d< globalIndex const > faceDofIndex =
    faceManager->getReference< array1d< globalIndex > >( dofIndexKey );

  // perform assembly loop over elements
  elemManager->forElementSubRegions(
    regions,
    [&]( localIndex const, auto const & subRegion ) {
      using FaceMapType = typename TYPEOFREF( subRegion )::FaceMapType;
      typename FaceMapType::ViewTypeConst const & faceMap =
        subRegion.template getReference< FaceMapType >(
          ElementSubRegionBase::viewKeyStruct::faceListString );

      localIndex const numFace = subRegion.numFacesPerElement();
      array1d< globalIndex > localDofIndex( numFace * numComp );
      array2d< real64 > localValues( numFace * numComp, numFace * numComp );
      localValues.setValues< serialPolicy >( 1.0 );

      for( localIndex k = 0; k < subRegion.size(); ++k )
      {
        for( localIndex a = 0; a < numFace; ++a )
        {
          for( localIndex c = 0; c < numComp; ++c )
          {
            localDofIndex[a * numComp + c] = faceDofIndex[faceMap[k][a]] + c;
          }
        }

        sparsity.insert( localDofIndex, localDofIndex, localValues );
      }
    } );
}

}  // namespace testing

}  // namespace geosx

#endif  //GEOSX_LINEARALGEBRA_UNITTESTS_TESTDOFMANAGERUTILS_HPP
