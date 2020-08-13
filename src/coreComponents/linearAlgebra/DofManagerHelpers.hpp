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
 * @file DofManagerHelpers.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_DOFMANAGERHELPERS_HPP
#define GEOSX_LINEARALGEBRA_DOFMANAGERHELPERS_HPP

namespace geosx
{
// unnamed namespace to avoid needless external linkage
namespace
{
/**
 * @brief A struct to abstract away some details of mesh access.
 * @tparam LOC type of mesh location
 */
template< DofManager::Location LOC >
struct MeshHelper
{};

template<>
struct MeshHelper< DofManager::Location::Node >
{
  using ManagerType = NodeManager;
  using LocalIndexType = localIndex;

  static LocalIndexType constexpr invalid_local_index { -1 };

  static constexpr auto managerGroupName =
    MeshLevel::groupStructKeys::nodeManagerString;
  static constexpr auto mapViewKey =
    ElementSubRegionBase::viewKeyStruct::nodeListString;
  static constexpr auto syncObjName = "node";

  template< typename MANAGER >
  using MapType = typename MANAGER::NodeMapType;
};

template<>
struct MeshHelper< DofManager::Location::Edge >
{
  using ManagerType = EdgeManager;
  using LocalIndexType = localIndex;

  static LocalIndexType constexpr invalid_local_index { -1 };

  static constexpr auto managerGroupName =
    MeshLevel::groupStructKeys::edgeManagerString;
  static constexpr auto mapViewKey =
    ElementSubRegionBase::viewKeyStruct::edgeListString;
  static constexpr auto syncObjName = "edge";

  template< typename MANAGER >
  using MapType = typename MANAGER::EdgeMapType;
};

template<>
struct MeshHelper< DofManager::Location::Face >
{
  using ManagerType = FaceManager;
  using LocalIndexType = localIndex;

  static LocalIndexType constexpr invalid_local_index { -1 };

  static constexpr auto managerGroupName =
    MeshLevel::groupStructKeys::faceManagerString;
  static constexpr auto mapViewKey =
    ElementSubRegionBase::viewKeyStruct::faceListString;
  static constexpr auto syncObjName = "face";

  template< typename MANAGER >
  using MapType = typename MANAGER::FaceMapType;
};

template<>
struct MeshHelper< DofManager::Location::Elem >
{
  using ManagerType = ElementSubRegionBase;
  using LocalIndexType = std::array< localIndex, 3 >;

  static LocalIndexType constexpr invalid_local_index { -1, -1, -1 };

  static constexpr auto managerGroupName =
    MeshLevel::groupStructKeys::elemManagerString;
  static constexpr auto syncObjName = "elems";

  template< typename MANAGER >
  using MapType = typename MANAGER::ElemMapType;
};

/**
 * @brief Compile-time limits on incidence between mesh objects.
 * @tparam LOC primary mesh object (location) type
 * @tparam LOC_ADJ adjacent mesh pbject type
 *
 * Typical usages:
 * - preallocate space in sparsity patterns
 * - choose size for stack arrays at compile time
 *
 * These limits may need to increase as we encounter more complex meshes.
 */
template< DofManager::Location LOC, DofManager::Location LOC_ADJ >
struct MeshIncidence
{};

/// Self-adjacency is always 1.
template< DofManager::Location LOC >
struct MeshIncidence< LOC, LOC >
{
  static localIndex constexpr max = 1;
};

/// Shortcut macro for declaring adjacency.
#define SET_MAX_MESH_INCIDENCE( LOC, LOC_ADJ, MAX )                                \
  template<>                                                                       \
  struct MeshIncidence< DofManager::Location::LOC, DofManager::Location::LOC_ADJ > \
  {                                                                                \
    static localIndex constexpr max = MAX;                                         \
  }

SET_MAX_MESH_INCIDENCE( Node, Edge, 10 );
SET_MAX_MESH_INCIDENCE( Node, Face, 10 );
SET_MAX_MESH_INCIDENCE( Node, Elem, 12 );

SET_MAX_MESH_INCIDENCE( Edge, Node, 2 );
SET_MAX_MESH_INCIDENCE( Edge, Face, 8 );
SET_MAX_MESH_INCIDENCE( Edge, Elem, 8 );

SET_MAX_MESH_INCIDENCE( Face, Node, 8 );
SET_MAX_MESH_INCIDENCE( Face, Edge, 4 );
SET_MAX_MESH_INCIDENCE( Face, Elem, 2 );

SET_MAX_MESH_INCIDENCE( Elem, Node, 8 );
SET_MAX_MESH_INCIDENCE( Elem, Edge, 12 );
SET_MAX_MESH_INCIDENCE( Elem, Face, 6 );

#undef SET_MAX_MESH_INCIDENCE

template< DofManager::Location LOC, typename MANAGER >
struct MapTypeHelper
{
  using type = typename MeshHelper< LOC >::template MapType< MANAGER >;
};

// return dummy type if target manager type does not declare a type alias to map to LOC objects
// this allows all switchyards to compile, but one shouldn't attempt to access a non-existent map
template< DofManager::Location LOC, typename MANAGER >
using MapType = typename MapTypeHelper< LOC, MANAGER >::type;

// some helper crust to extract underlying type from InterObjectRelation and the likes
template< typename T, bool >
struct BaseTypeHelper
{
  using type = T;
};

template< typename T >
struct BaseTypeHelper< T, true >
{
  using type = typename T::base_type;
};

HAS_ALIAS( base_type );

template< typename MAP >
using BaseType = typename BaseTypeHelper< MAP, HasAlias_base_type< MAP > >::type;

/**
 * @brief Helper struct that specializes access to various map types
 * @tparam MAP type of the map
 */
template< typename MAP >
struct MapHelperImpl
{};

template< typename T, typename PERMUTATION >
struct MapHelperImpl< array2d< T, PERMUTATION > >
{
  static localIndex
  size0( array2d< T, PERMUTATION > const & map )
  {
    return map.size( 0 );
  }

  static localIndex
  size1( array2d< T, PERMUTATION > const & map,
         localIndex const GEOSX_UNUSED_PARAM( i0 ) )
  {
    return map.size( 1 );
  }

  static T const &
  value( array2d< T, PERMUTATION > const & map,
         localIndex const i0,
         localIndex const i1 )
  {
    return map( i0, i1 );
  }

  static localIndex
  size0(
    arrayView2d< T const,
                 LvArray::typeManipulation::getStrideOneDimension( PERMUTATION {} ) > const & map )
  {
    return map.size( 0 );
  }

  static localIndex
  size1(
    arrayView2d< T const,
                 LvArray::typeManipulation::getStrideOneDimension( PERMUTATION {} ) > const & map,
    localIndex const GEOSX_UNUSED_PARAM( i0 ) )
  {
    return map.size( 1 );
  }

  static T const &
  value(
    arrayView2d< T const,
                 LvArray::typeManipulation::getStrideOneDimension( PERMUTATION {} ) > const & map,
    localIndex const i0,
    localIndex const i1 )
  {
    return map( i0, i1 );
  }
};

template< typename T >
struct MapHelperImpl< array1d< array1d< T > > >
{
  static localIndex
  size0( array1d< array1d< T > > const & map )
  {
    return map.size();
  }

  static localIndex
  size1( array1d< array1d< T > > const & map, localIndex const i0 )
  {
    return map[i0].size();
  }

  static T const &
  value( array1d< array1d< T > > const & map,
         localIndex const i0,
         localIndex const i1 )
  {
    return map[i0][i1];
  }

  static localIndex
  size0( arrayView1d< arrayView1d< T const > const > const & map )
  {
    return map.size();
  }

  static localIndex
  size1( arrayView1d< arrayView1d< T const > const > const & map,
         localIndex const i0 )
  {
    return map[i0].size();
  }

  static T const &
  value( arrayView1d< arrayView1d< T const > const > const & map,
         localIndex const i0,
         localIndex const i1 )
  {
    return map[i0][i1];
  }
};

template< typename T >
struct MapHelperImpl< ArrayOfArrays< T > >
{
  static localIndex
  size0( ArrayOfArraysView< T const > const & map )
  {
    return map.size();
  }

  static localIndex
  size1( ArrayOfArraysView< T const > const & map,
         localIndex const i0 )
  {
    return map.sizeOfArray( i0 );
  }

  static T const &
  value( ArrayOfArraysView< T const > const & map,
         localIndex const i0,
         localIndex const i1 )
  {
    return map( i0, i1 );
  }
};

template< typename T >
struct MapHelperImpl< array1d< SortedArray< T > > >
{
  static localIndex
  size0( array1d< SortedArray< T > > const & map )
  {
    return map.size();
  }

  static localIndex
  size1( array1d< SortedArray< T > > const & map, localIndex const i0 )
  {
    return map[i0].size();
  }

  static T const &
  value( array1d< SortedArray< T > > const & map,
         localIndex const i0,
         localIndex const i1 )
  {
    return map[i0][i1];
  }
};

template< typename T >
struct MapHelperImpl< ArrayOfSets< T > >
{
  static localIndex
  size0( ArrayOfSetsView< T const > const & map )
  {
    return map.size();
  }

  static localIndex
  size1( ArrayOfSetsView< T const > const & map, localIndex const i0 )
  {
    return map.sizeOfSet( i0 );
  }

  static T const &
  value( ArrayOfSetsView< T const > const & map,
         localIndex const i0,
         localIndex const i1 )
  {
    return map( i0, i1 );
  }
};

template< typename BASETYPE >
struct MapHelperImpl< ToElementRelation< BASETYPE > >
{
  using BaseHelper = MapHelperImpl< BASETYPE >;

  static localIndex
  size0( ToElementRelation< BASETYPE > const & map )
  {
    return BaseHelper::size0( map.m_toElementIndex );
  }

  static localIndex
  size1( ToElementRelation< BASETYPE > const & map,
         localIndex const i0 )
  {
    return BaseHelper::size1( map.m_toElementIndex, i0 );
  }

  static auto
  value( ToElementRelation< BASETYPE > const & map,
         localIndex const i0,
         localIndex const i1 )
  {
    return MeshHelper< DofManager::Location::Elem >::LocalIndexType {
      BaseHelper::value( map.m_toElementRegion, i0, i1 ),
      BaseHelper::value( map.m_toElementSubRegion, i0, i1 ),
      BaseHelper::value( map.m_toElementIndex, i0, i1 ) };
  }
};

/**
 * @brief Helper struct that specializes access to various map types
 * @tparam MAP type of the map
 *
 * @note We may need to strip off InterObjectRelation and get the underlying map type, hence the extra layer
 */
template< typename MAP >
using MapHelper = MapHelperImpl< BaseType< MAP > >;

/**
 * @brief Switchyard for Location (Node/Edge/Face/Elem)
 * @tparam LAMBDA generic functor type
 * @param loc location type
 * @param lambda functor to be called
 */
template< typename LAMBDA >
bool
LocationSwitch( DofManager::Location const loc, LAMBDA lambda )
{
  switch( loc )
  {
    case DofManager::Location::Node:
    {
      lambda(
        std::integral_constant< DofManager::Location, DofManager::Location::Node >() );
      return true;
    }
    case DofManager::Location::Edge:
    {
      lambda(
        std::integral_constant< DofManager::Location, DofManager::Location::Edge >() );
      return true;
    }
    case DofManager::Location::Face:
    {
      lambda(
        std::integral_constant< DofManager::Location, DofManager::Location::Face >() );
      return true;
    }
    case DofManager::Location::Elem:
    {
      lambda(
        std::integral_constant< DofManager::Location, DofManager::Location::Elem >() );
      return true;
    }
    default:
      return false;
  }
}

template< typename LAMBDA >
bool
LocationSwitch( DofManager::Location const loc1,
                DofManager::Location const loc2,
                LAMBDA lambda )
{
  bool ret2;
  bool const ret1 = LocationSwitch( loc1, [&]( auto const loc_type1 ) {
    ret2 = LocationSwitch( loc2, [&]( auto const loc_type2 ) {
      lambda( loc_type1, loc_type2 );
    } );
  } );
  return ret1 && ret2;
}

template< DofManager::Location LOC >
typename MeshHelper< LOC >::ManagerType const &
getObjectManager(
  MeshLevel const * const mesh )
{
  using ObjectManager = typename MeshHelper< LOC >::ManagerType;
  GEOSX_ASSERT( mesh != nullptr );
  ObjectManager const * manager =
    mesh->GetGroup< ObjectManager >( MeshHelper< LOC >::managerGroupName );
  GEOSX_ASSERT( manager != nullptr );
  return *manager;
}

template< DofManager::Location LOC >
typename MeshHelper< LOC >::ManagerType &
getObjectManager( MeshLevel * const mesh )
{
  using ObjectManager = typename MeshHelper< LOC >::ManagerType;
  return const_cast< ObjectManager & >(
    getObjectManager< LOC >( const_cast< MeshLevel const * >( mesh ) ) );
}

ObjectManagerBase const &
getObjectManager( DofManager::Location const loc,
                  MeshLevel const * const mesh )
{
  ObjectManagerBase const * manager = nullptr;
  LocationSwitch( loc, [&]( auto const LOC ) {
    manager = &getObjectManager< decltype( LOC )::value >( mesh );
  } );
  GEOSX_ASSERT( manager != nullptr );
  return *manager;
}

ObjectManagerBase &
getObjectManager( DofManager::Location const loc,
                  MeshLevel * const mesh )
{
  return const_cast< ObjectManagerBase & >(
    getObjectManager( loc, const_cast< MeshLevel const * >( mesh ) ) );
}

/**
 * @brief A helper struct implementing visitation of mesh objects and connected (adjacent) objects
 * @tparam LOC type of primary mesh object (location) to visit
 * @tparam CONN_LOC type of adjacent mesh object to visit from primary; can be the same as primary,
 *                  in which case each is visited only once and the same index is passed for both.
 * @tparam VISIT_GHOSTS whether to visit ghosted primary locations
 */
template< DofManager::Location LOC, DofManager::Location CONN_LOC, bool VISIT_GHOSTS >
struct MeshLoopHelper;

/**
 * @brief A specialization of MeshLoopHelper used when adjacent objects are NOT needed
 *        and the primary objects are NOT elements.
 * @tparam LOC type of primary mesh object (location) to visit
 * @tparam VISIT_GHOSTS whether to visit ghosted primary locations
 *
 * The algorithm will visit objects in order of increasing local index.
 */
template< DofManager::Location LOC, bool VISIT_GHOSTS >
struct MeshLoopHelper< LOC, LOC, VISIT_GHOSTS >
{
  template< typename... SUBREGIONTYPES, typename LAMBDA >
  static void
  visit( MeshLevel * const meshLevel,
         std::vector< std::string > const & regions,
         LAMBDA lambda )
  {
    // derive some useful type aliases
    using ObjectManagerLoc = typename MeshHelper< LOC >::ManagerType;

    // get access to location ghost rank (we don't want to visit ghosted locations
    ObjectManagerLoc const & objectManager = getObjectManager< LOC >( meshLevel );
    arrayView1d< integer const > const & ghostRank = objectManager.ghostRank();

    // create an array to track previously visited locations (to avoid multiple visits)
    array1d< bool > locationsVisited( objectManager.size() );
    locationsVisited.setValues< serialPolicy >( false );

    // create (and overallocate) an array to collect indicies to visit
    array1d< localIndex > locationsToVisit {};
    locationsToVisit.reserve( objectManager.size() );

    meshLevel->getElemManager()->forElementSubRegions< SUBREGIONTYPES... >(
      regions,
      [&]( localIndex const, auto const & subRegion ) {
        // derive some more useful, subregion-dependent type aliases
        using ElementSubRegionType = std::remove_reference_t< decltype( subRegion ) >;
        using ElemToLocMapType = MapType< LOC, ElementSubRegionType >;

        // get access to element-to-location map
        auto const & elemToLocMap =
          subRegion.template getReference< ElemToLocMapType >(
            MeshHelper< LOC >::mapViewKey );

        // loop over all elements (including ghosts, which may be necessary to access some locally owned locations)
        for( localIndex ei = 0; ei < subRegion.size(); ++ei )
        {
          // loop over all locations incident on an element
          for( localIndex a = 0;
               a < MapHelper< ElemToLocMapType >::size1( elemToLocMap, ei );
               ++a )
          {
            localIndex const locIdx =
              MapHelper< ElemToLocMapType >::value( elemToLocMap, ei, a );

            // check if we should visit this location
            if( ( VISIT_GHOSTS || ghostRank[locIdx] < 0 ) &&
                !std::exchange( locationsVisited[locIdx], true ) )
            {
              locationsToVisit.emplace_back( locIdx );
            }
          }
        }
      } );

    // optimize for the common case (e.g. all nodes of the mesh)
    if( locationsToVisit.size() == objectManager.size() )
    {
      for( localIndex locIdx = 0; locIdx < objectManager.size(); ++locIdx )
      {
        lambda( locIdx, locIdx, 0 );
      }
    }
    else
    {
      // I think this is faster than alternatives, but need to profile
      std::sort( locationsToVisit.begin(), locationsToVisit.end() );
      for( localIndex const locIdx : locationsToVisit )
      {
        lambda( locIdx, locIdx, 0 );
      }
    }
  }
};

/**
 * @brief A specialization of MeshLoopHelper used when NEITHER primary NOR adjacent objects are elements
 * @tparam LOC type of primary mesh object (location) to visit
 * @tparam CONN_LOC type of adjacent mesh object to visit from primary
 * @tparam VISIT_GHOSTS whether to visit ghosted primary locations
 */
template< DofManager::Location LOC, DofManager::Location CONN_LOC, bool VISIT_GHOSTS >
struct MeshLoopHelper
{
  template< typename... SUBREGIONTYPES, typename LAMBDA >
  static void
  visit( MeshLevel * const meshLevel,
         std::vector< std::string > const & regions,
         LAMBDA lambda )
  {
    // derive some useful type aliases
    using ObjectManagerLoc = typename MeshHelper< LOC >::ManagerType;
    using LocToConnMapType = MapType< CONN_LOC, ObjectManagerLoc >;

    ObjectManagerLoc const & objectManager = getObjectManager< LOC >( meshLevel );

    // get access to location-to-connected map
    auto const & locToConnMap =
      objectManager.template getReference< LocToConnMapType >(
        MeshHelper< CONN_LOC >::mapViewKey );

    // call the specialized version first, then add an extra loop over connected objects
    MeshLoopHelper< LOC, LOC, VISIT_GHOSTS >::template visit< SUBREGIONTYPES... >(
      meshLevel,
      regions,
      [&]( localIndex const locIdx, localIndex const, localIndex const ) {
        // loop over all connected locations
        for( localIndex b = 0;
             b < MapHelper< LocToConnMapType >::size1( locToConnMap, locIdx );
             ++b )
        {
          auto const connIdx =
            MapHelper< LocToConnMapType >::value( locToConnMap, locIdx, b );
          lambda( locIdx, connIdx, b );
        }
      } );
  }
};

/**
 * @brief A specialization of MeshLoopHelper used when primary objects are NOT elements, but adjacent are
 * @tparam LOC type of primary mesh object (location) to visit
 * @tparam VISIT_GHOSTS whether to visit ghosted primary locations
 *
 * @note We have to have this specialization because to-element maps are different from all other
 *       map types in that they are stored in data repository as three separate arrays.
 */
template< DofManager::Location LOC, bool VISIT_GHOSTS >
struct MeshLoopHelper< LOC, DofManager::Location::Elem, VISIT_GHOSTS >
{
  template< typename... SUBREGIONTYPES, typename LAMBDA >
  static void
  visit( MeshLevel * const meshLevel,
         std::vector< std::string > const & regions,
         LAMBDA lambda )
  {
    // derive some useful type aliases
    using ObjectManagerLoc = typename MeshHelper< LOC >::ManagerType;
    using ToElemMapType =
      BaseType< MapType< DofManager::Location::Elem, ObjectManagerLoc > >;

    // get mesh object manager for LOC to access maps
    ObjectManagerLoc const & objectManager = getObjectManager< LOC >( meshLevel );

    // access to location-to-element map
    auto const & elemRegionList =
      objectManager.template getReference< ToElemMapType >(
        ObjectManagerLoc::viewKeyStruct::elementRegionListString );
    auto const & elemSubRegionList =
      objectManager.template getReference< ToElemMapType >(
        ObjectManagerLoc::viewKeyStruct::elementSubRegionListString );
    auto const & elemIndexList =
      objectManager.template getReference< ToElemMapType >(
        ObjectManagerLoc::viewKeyStruct::elementListString );

    // call the specialized version first, then add an extra loop over connected elements
    MeshLoopHelper< LOC, LOC, VISIT_GHOSTS >::template visit< SUBREGIONTYPES... >(
      meshLevel,
      regions,
      [&]( localIndex const locIdx, localIndex const, localIndex const ) {
        // loop over all connected locations
        for( localIndex b = 0;
             b < MapHelper< ToElemMapType >::size1( elemIndexList, locIdx );
             ++b )
        {
          localIndex const er =
            MapHelper< ToElemMapType >::value( elemRegionList, locIdx, b );
          localIndex const esr =
            MapHelper< ToElemMapType >::value( elemSubRegionList, locIdx, b );
          localIndex const ei =
            MapHelper< ToElemMapType >::value( elemIndexList, locIdx, b );

          if( er >= 0 && esr >= 0 && ei >= 0 )
          {
            lambda(
              locIdx,
              MeshHelper< DofManager::Location::Elem >::LocalIndexType { er, esr, ei },
              b );
          }
        }
      } );
  }
};

/**
 * @brief A specialization of MeshLoopHelper used when primary objects are elements and adjacent are NOT
 * @tparam CONN_LOC type of adjacent mesh objects (locations)
 * @tparam VISIT_GHOSTS whether to visit ghosted primary locations
 *
 * @note We have to have this specialization because various subregion types might have different
 *       map types, so the correct map type must be derived within the subregion loop
 */
template< DofManager::Location CONN_LOC, bool VISIT_GHOSTS >
struct MeshLoopHelper< DofManager::Location::Elem, CONN_LOC, VISIT_GHOSTS >
{
  template< typename... SUBREGIONTYPES, typename LAMBDA >
  static void
  visit( MeshLevel * const meshLevel,
         std::vector< std::string > const & regions,
         LAMBDA lambda )
  {
    meshLevel->getElemManager()->forElementSubRegionsComplete< SUBREGIONTYPES... >(
      regions,
      [&]( localIndex const,
           localIndex const er,
           localIndex const esr,
           ElementRegionBase const &,
           auto const & subRegion ) {
        // derive subregion-dependent map type
        using ElemToConnMapType = MapType< CONN_LOC, TYPEOFREF( subRegion ) >;

        // get access to element-to-location map
        auto const & elemToConnMap =
          subRegion.template getReference< ElemToConnMapType >(
            MeshHelper< CONN_LOC >::mapViewKey );

        arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

        for( localIndex ei = 0; ei < subRegion.size(); ++ei )
        {
          if( VISIT_GHOSTS || elemGhostRank[ei] < 0 )
          {
            auto const elemIdx =
              MeshHelper< DofManager::Location::Elem >::LocalIndexType { er, esr, ei };

            for( localIndex a = 0;
                 a < MapHelper< ElemToConnMapType >::size1( elemToConnMap, ei );
                 ++a )
            {
              localIndex const locIdx =
                MapHelper< ElemToConnMapType >::value( elemToConnMap, ei, a );
              lambda( elemIdx, locIdx, a );
            }
          }
        }
      } );
  }
};

/**
 * @brief A specialization of MeshLoopHelper used when only loop over elements is required
 * @tparam VISIT_GHOSTS whether to visit ghosted elements
 */
template< bool VISIT_GHOSTS >
struct MeshLoopHelper< DofManager::Location::Elem, DofManager::Location::Elem, VISIT_GHOSTS >
{
  template< typename... SUBREGIONTYPES, typename LAMBDA >
  static void
  visit( MeshLevel * const meshLevel,
         std::vector< std::string > const & regions,
         LAMBDA && lambda )
  {
    meshLevel->getElemManager()->forElementSubRegionsComplete< SUBREGIONTYPES... >(
      regions,
      [&]( localIndex const,
           localIndex const er,
           localIndex const esr,
           ElementRegionBase const &,
           ElementSubRegionBase const & subRegion ) {
        arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
        localIndex const numElems = subRegion.size();

        for( localIndex ei = 0; ei < numElems; ++ei )
        {
          if( VISIT_GHOSTS || elemGhostRank[ei] < 0 )
          {
            auto const elemIdx =
              MeshHelper< DofManager::Location::Elem >::LocalIndexType { er, esr, ei };
            lambda( elemIdx, elemIdx, 0 );
          }
        }
      } );
  }
};

/**
 * @brief Visit mesh locations within active regions and their connected locations (e.g. nodes of an element)
 *
 * @tparam LOC type of location (Node/Edge/Face/Elem)
 * @tparam CONN_LOC type of connected location (Node/Edge/Face/Elem)
 * @tparam SUBREGIONTYPES variadic pack of subregion types to visit;
 *         if omitted, uses ElementRegionManager's default behavior
 * @tparam LAMBDA type of user-provided functor or lambda to call
 * @param mesh the mesh to loop over
 * @param regions list of region names to visit (assumed unique)
 * @param lambda functor or lambda to call
 *
 * Calls user-provided lambda with 3 parameters:
 * - mesh local index of primary location (one index for node/edge/face, reg/subreg/index triplet for elems)
 * - mesh local index of connected location
 * - index of connected location inside primary location's map (e.g. which node of an element)
 * If CONN_LOC == LOC, lambda is called only once per location, as lambda(loc, loc, 0)
 *
 * @note This function will exclude ghosted locations, but will not perform a similar check for
 *       connected locations - these may include ghosts, it is up to the user to filter as needed.
 *       Similarly, while primary loop is limited to @p regions, adjacent locations may not belong.
 */
template< DofManager::Location LOC,
          DofManager::Location CONN_LOC,
          bool VISIT_GHOSTS,
          typename... SUBREGIONTYPES,
          typename LAMBDA >
void
forMeshLocation( MeshLevel * const mesh,
                 std::vector< std::string > const & regions,
                 LAMBDA && lambda )
{
  MeshLoopHelper< LOC, CONN_LOC, VISIT_GHOSTS >::template visit< SUBREGIONTYPES... >(
    mesh,
    regions,
    std::forward< LAMBDA >( lambda ) );
}

/**
 * @brief A shortcut for forMeshLocation with CONN_LOC == LOC
 */
template< DofManager::Location LOC, bool VISIT_GHOSTS, typename... SUBREGIONTYPES, typename LAMBDA >
void
forMeshLocation( MeshLevel * const mesh,
                 std::vector< std::string > const & regions,
                 LAMBDA && lambda )
{
  forMeshLocation< LOC, LOC, VISIT_GHOSTS, SUBREGIONTYPES... >(
    mesh,
    regions,
    [&]( auto const locIdx, auto const, auto const ) { lambda( locIdx ); } );
}

/**
 * @brief Count objects of given type within a set of regions
 * @tparam LOC type of mesh objects (locations)
 * @tparam VISIT_GHOSTS whether ghosted objects should be counted
 * @tparam SUBREGIONTYPES types of subregions to include in counting
 * @param mesh the mesh level to use
 * @param regions a list of region names (assumed to be unique)
 * @return
 */
template< DofManager::Location LOC, bool VISIT_GHOSTS, typename... SUBREGIONTYPES >
localIndex
countMeshObjects( MeshLevel * const mesh,
                  std::vector< std::string > const & regions )
{
  localIndex count = 0;
  forMeshLocation< LOC, VISIT_GHOSTS, SUBREGIONTYPES... >(
    mesh,
    regions,
    [&]( auto const GEOSX_UNUSED_PARAM( connIdx ) ) { ++count; } );
  return count;
}

/**
 * @brief This template abstracts some differences between index arrays that live
 *        on elements vs other mesh objects: (de)registration, value access, etc.
 * @tparam LOC target location
 */
template< typename INDEX, DofManager::Location LOC >
struct IndexArrayHelper
{
  using IndexType = localIndex;
  using ArrayType = array1d< std::remove_const_t< INDEX > >;
  using ViewType = arrayView1d< INDEX >;
  using Accessor = ViewType const;
  using Mesh = add_const_if_t< MeshLevel, std::is_const< INDEX >::value >;

  template< typename... SUBREGIONTYPES >
  static void
  create( Mesh * const mesh,
          string const & key,
          string const & description,
          std::vector< std::string > const & GEOSX_UNUSED_PARAM( regions ) )
  {
    ObjectManagerBase & baseManager = getObjectManager< LOC >( mesh );
    baseManager.registerWrapper< ArrayType >( key )
      ->setApplyDefaultValue( -1 )
      ->setPlotLevel( dataRepository::PlotLevel::LEVEL_1 )
      ->setRestartFlags( dataRepository::RestartFlags::NO_WRITE )
      ->setDescription( description );
  }

  static Accessor
  get( Mesh * const mesh, string const & key )
  {
    return getObjectManager< LOC >( mesh ).template getReference< ArrayType >( key );
  }

  static inline INDEX
  value( Accessor & indexArray,
         typename MeshHelper< LOC >::LocalIndexType const i )
  {
    return indexArray[i];
  }

  static inline INDEX &
  reference( Accessor & indexArray,
             typename MeshHelper< LOC >::LocalIndexType const i )
  {
    return indexArray[i];
  }

  template< typename... SUBREGIONTYPES >
  static void
  remove( Mesh * const mesh,
          string const & key,
          std::vector< std::string > const & GEOSX_UNUSED_PARAM( regions ) )
  {
    getObjectManager< LOC >( mesh ).deregisterWrapper( key );
  }
};

/**
 * @brief Specialized version of IndexArrayHelper for elements.
 */
template< typename INDEX >
struct IndexArrayHelper< INDEX, DofManager::Location::Elem >
{
  using IndexType = INDEX;
  using ArrayType = array1d< std::remove_const_t< INDEX > >;
  using ViewType = arrayView1d< INDEX >;
  using Accessor = ElementRegionManager::ElementViewAccessor< ViewType > const;
  using Mesh = add_const_if_t< MeshLevel, std::is_const< INDEX >::value >;

  template< typename... SUBREGIONTYPES >
  static void
  create( Mesh * const mesh,
          string const & key,
          string const & description,
          std::vector< std::string > const & regions )
  {
    mesh->getElemManager()->template forElementSubRegions< SUBREGIONTYPES... >(
      regions,
      [&]( localIndex const, auto & subRegion ) {
        subRegion.template registerWrapper< ArrayType >( key )
          ->setApplyDefaultValue( -1 )
          ->setPlotLevel( dataRepository::PlotLevel::LEVEL_1 )
          ->setRestartFlags( dataRepository::RestartFlags::NO_WRITE )
          ->setDescription( description );
      } );
  }

  static Accessor
  get( Mesh * const mesh, string const & key )
  {
    return mesh->getElemManager()
      ->template ConstructViewAccessor< ArrayType, ViewType >( key );
  }

  static inline INDEX
  value(
    Accessor & indexArray,
    MeshHelper< DofManager::Location::Elem >::LocalIndexType const & e )
  {
    if( indexArray[std::get< 0 >( e )].empty() ||
        indexArray[std::get< 0 >( e )][std::get< 1 >( e )].empty() )
    {
      return -1;
    }
    return indexArray[std::get< 0 >( e )][std::get< 1 >( e )][std::get< 2 >( e )];
  }

  static inline INDEX &
  reference(
    Accessor & indexArray,
    MeshHelper< DofManager::Location::Elem >::LocalIndexType const & e )
  {
    return indexArray[std::get< 0 >( e )][std::get< 1 >( e )][std::get< 2 >( e )];
  }

  template< typename... SUBREGIONTYPES >
  static void
  remove( Mesh * const mesh,
          string const & key,
          std::vector< std::string > const & regions )
  {
    mesh->getElemManager()->template forElementSubRegions< SUBREGIONTYPES... >(
      regions,
      [&]( localIndex const, ElementSubRegionBase & subRegion ) {
        subRegion.deregisterWrapper( key );
      } );
  }
};

}  // namespace

}  // namespace geosx

#endif  //GEOSX_LINEARALGEBRA_DOFMANAGERHELPERS_HPP
