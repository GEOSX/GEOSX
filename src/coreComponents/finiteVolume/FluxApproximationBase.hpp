/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FluxApproximationBase.hpp
 *
 */

#ifndef GEOSX_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_
#define GEOSX_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_

#include "dataRepository/Group.hpp"
#include "finiteVolume/FluxStencil.hpp"
#include "CellElementStencilTPFA.hpp"
#include "FaceElementStencil.hpp"
#include "mesh/DomainPartition.hpp"

namespace geosx
{

/**
 * @struct CellDescriptor
 * @brief A structure containing a single cell (element) identifier triplet.
 */
struct CellDescriptor
{
  /// region index
  localIndex region;
  /// subregion index
  localIndex subRegion;
  /// cell index
  localIndex index;

  /**
   * @brief Comparison operator between two CellDescriptors.
   * @param[in] other the CellDescriptor to compare with
   * @return true if they represent the same mesh element
   */
  bool operator==( CellDescriptor const & other )
  {
    return( region==other.region && subRegion==other.subRegion && index==other.index );
  }
};

/**
 * @struct PointDescriptor
 * @brief A structure describing an arbitrary point participating in a stencil.
 *
 * Nodal and face center points are identified by local mesh index.
 * Cell center points are identified by a triplet <region,subregion,index>.
 *
 * The sad reality is, a boundary flux MPFA stencil may be comprised of a mix of
 * cell and face centroids, so we have to discriminate between them at runtime
 */
struct PointDescriptor
{
  /// Enum to classify the variable location
  enum class Tag { CELL, FACE, NODE };

  /// The tag
  Tag tag;

  /// union to characterize a PointDescriptor
  union
  {
    /// node index
    localIndex nodeIndex;
    /// face index
    localIndex faceIndex;
    /// CellDescriptor index
    CellDescriptor cellIndex;
  };
};

/**
 * @class FluxApproximationBase
 *
 * Base class for various flux approximation classes.
 * Stores the main and boundary stencils, construction is implemented in derived classes.
 * Main stencil is the one for cell-to-cell fluxes.
 * Boundary stencils are for Dirichlet boundary conditions
 */
class FluxApproximationBase : public dataRepository::Group
{
public:

  /// Alias for CatalogInterface, necessary declarations for factory instantiation of derived classes
  using CatalogInterface = dataRepository::CatalogInterface< FluxApproximationBase, string const &, Group * const >;
  /**
   * @brief Return the data type in the data repository.
   * @return the data type in the data repository
   */
  static typename CatalogInterface::CatalogType & getCatalog();

  FluxApproximationBase() = delete;

  /**
   * @brief Constructor.
   * @param name the name of the FluxApproximationBase in the data repository
   * @param parent the parent group of this group.
   */
  FluxApproximationBase( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Extract stencil stored under the mesh group.
   * @tparam TYPE type of Stencil to get
   * @param mesh the mesh level object
   * @param name name of the stencil object
   * @return reference to the stencil
   */
  template< typename TYPE >
  TYPE const & getStencil( MeshLevel const & mesh, string const & name ) const;

  /**
   * @copydoc getStencil(MeshLevel const &, string const &) const
   */
  template< typename TYPE >
  TYPE & getStencil( MeshLevel & mesh, string const & name ) const;

  /**
   * @brief Call a user-provided function for each stencil.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] mesh the mesh level containing the stencils
   * @param[in] lambda The LAMBDA function
   */
  template< typename LAMBDA >
  void forAllStencils( MeshLevel const & mesh, LAMBDA && lambda ) const;

  /**
   * @brief Call a user-provided function for the each stencil according to the provided TYPE.
   * @tparam TYPE The type to be passed to forWrappers
   * @tparam TYPES Other types to be passed to forWrappers
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] mesh the mesh level containing the stencils
   * @param[in] lambda The LAMBDA function
   */
  template< typename TYPE, typename ... TYPES, typename LAMBDA >
  void forStencils( MeshLevel const & mesh, LAMBDA && lambda ) const;

  /**
   * @brief Add a new fracture stencil.
   * @param[in,out] mesh the mesh on which to add the fracture stencil
   * @param[in] faceElementRegionName the face element region name
   * @param[in] initFlag if true initialize physical fields, like pressure
   */
  virtual void addToFractureStencil( MeshLevel & mesh,
                                     string const & faceElementRegionName,
                                     bool const initFlag ) const = 0;

  /**
   * @brief Add a new embedded fracture stencil.
   * @param[in,out] mesh the mesh on which to add the fracture stencil
   * @param[in] embeddedSurfaceRegionName the embedded surface element region name
   */
  virtual void addEDFracToFractureStencil( MeshLevel & mesh,
                                           string const & embeddedSurfaceRegionName ) const = 0;

  /**
   * @brief View keys.
   */
  struct viewKeyStruct
  {
    /// @return The key for fieldName
    static constexpr char const * fieldNameString() { return "fieldName"; }

    /// @return The key for coefficientName
    static constexpr char const * coeffNameString() { return "coefficientName"; }

    /// @return The key for coefficientName
    static constexpr char const * coefficientModelNamesString() { return "coefficientModelNames"; }

    /// @return The key for targetRegions
    static constexpr char const * targetRegionsString() { return "targetRegions"; }

    /// @return The key for areaRelTol
    static constexpr char const * areaRelativeToleranceString() { return "areaRelTol"; }

    /// @return The key for transMultiplier
    static constexpr char const * transMultiplierString() { return "TransMultiplier"; }


    // Keys below are for wrappers registered on MeshLevel, not the current object

    /// @return The key for cellStencil
    static constexpr char const * cellStencilString() { return "cellStencil"; }

    /// @return The key for fractureStencil
    static constexpr char const * fractureStencilString() { return "fractureStencil"; }
  };

  /**
   * @brief Group keys.
   */
  struct groupKeyStruct
  {
    /// @return Key under which the top-level group for all FV stencils will be registered on MeshLevel
    static constexpr auto stencilMeshGroupString() { return "finiteVolumeStencils"; }
  };

  /**
   * @brief Returns the target region name.
   * @return the target region name
   */
  string_array const & targetRegions() const { return m_targetRegions; }
  /**
   * @copydoc targetRegions() const
   */
  string_array & targetRegions()       { return m_targetRegions; }

protected:

  /// @copydoc geosx::dataRepository::Group::registerDataOnMesh
  virtual void registerDataOnMesh( Group & meshBodies ) override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  /**
   * @brief Register the wrapper for cell stencil on a mesh.
   * @param stencilGroup the group holding the stencil objects
   */
  virtual void registerCellStencil( Group & stencilGroup ) const = 0;

  /**
   * @brief Actual computation of the cell-to-cell stencil, to be overridden by implementations.
   * @param[in] mesh the mesh on which to perform the computation
   */
  virtual void computeCellStencil( MeshLevel & mesh ) const = 0;

  /**
   * @brief Register the wrapper for fracture stencil on a mesh.
   * @param stencilGroup the group holding the stencil objects
   */
  virtual void registerFractureStencil( Group & stencilGroup ) const = 0;

  /**
   * @brief Register the wrapper for boundary face stencil on a mesh.
   * @param stencilGroup the group holding the stencil objects
   * @param setName the face set name (used as the wrapper name)
   */
  virtual void registerBoundaryStencil( Group & stencilGroup,
                                        string const & setName ) const = 0;

  /**
   * @brief Allocate and populate a stencil to be used in boundary condition application
   * @param mesh the target mesh level
   * @param setName name of the face set, to be used as wrapper name for the produced stencil
   * @param faceSet set of face indices to use
   */
  virtual void computeBoundaryStencil( MeshLevel & mesh,
                                       string const & setName,
                                       SortedArrayView< localIndex const > const & faceSet ) const = 0;

  /// name of the primary solution field
  string m_fieldName;

  /// name of the coefficient field
  string m_coeffName;

  /// names of target regions to build the stencil for
  string_array m_coefficientModelNames;

  /// names of target regions to build the stencil for
  string_array m_targetRegions;

  /// relative tolerance
  real64 m_areaRelTol;

  /// length scale of the mesh body
  real64 m_lengthScale;

};

template< typename TYPE >
TYPE const & FluxApproximationBase::getStencil( MeshLevel const & mesh, string const & name ) const
{
  Group const & stencilGroup = mesh.getGroup( groupKeyStruct::stencilMeshGroupString() ).getGroup( getName() );
  return stencilGroup.getReference< TYPE >( name );
}

template< typename TYPE >
TYPE & FluxApproximationBase::getStencil( MeshLevel & mesh, string const & name ) const
{
  Group & stencilGroup = mesh.getGroup( groupKeyStruct::stencilMeshGroupString() ).getGroup( getName() );
  return stencilGroup.getReference< TYPE >( name );
}

template< typename LAMBDA >
void FluxApproximationBase::forAllStencils( MeshLevel const & mesh, LAMBDA && lambda ) const
{
  //TODO remove dependence on CellElementStencilTPFA and FaceElementStencil
  forStencils< CellElementStencilTPFA, FaceElementStencil >( mesh, std::forward< LAMBDA >( lambda ) );
}

template< typename TYPE, typename ... TYPES, typename LAMBDA >
void FluxApproximationBase::forStencils( MeshLevel const & mesh, LAMBDA && lambda ) const
{
  Group const & stencilGroup = mesh.getGroup( groupKeyStruct::stencilMeshGroupString() ).getGroup( getName() );
  stencilGroup.forWrappers< TYPE, TYPES... >( [&] ( auto const & wrapper )
  {
    lambda( wrapper.reference() );
  } );
}

} // namespace geosx

#endif //GEOSX_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_
