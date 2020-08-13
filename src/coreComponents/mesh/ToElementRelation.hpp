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
 * @file ToElementRelation.hpp
 */

#ifndef GEOSX_MESH_TOELEMENTRELATION_HPP_
#define GEOSX_MESH_TOELEMENTRELATION_HPP_

#include "InterObjectRelation.hpp"

namespace geosx
{
class ElementRegionManager;

/**
 * @brief A relationship to an element.
 * @tparam BASETYPE The underlying relation type to use to
 *                  store the relationsip information.
 */
template< typename BASETYPE >
class ToElementRelation
{
public:
  /// The type of the underlying relationship storage object.
  using base_type = BASETYPE;

  /// The default constructor
  ToElementRelation();
  /// The default destructor
  ~ToElementRelation();

  /**
   * @brief Resize the underlying relationship storage.
   * @tparam DIMS The types of each dimensions resize parameter.
   * @param newdims A parameter pack of appropriate size to resize each
   *                dimension of the relationship storage.
   */
  template< typename... DIMS >
  void
  resize( DIMS... newdims );

  /**
   * @brief Get the current size of the relationship storage.
   * @return The current size of the relationship storage.
   */
  localIndex
  size() const
  {
    return m_toElementRegion.size();
  }

  /**
   * @brief Get the size of a specific dimension of the relationship storage.
   * @param dim The dimension to get the storage size of.
   * @return The dimension size
   */
  localIndex
  size( int const dim ) const
  {
    return m_toElementRegion.size( dim );
  }

  /**
   * @brief Set the ElementRegionManager.
   * @param input The ElementRegionManager to set.
   */
  void
  setElementRegionManager( ElementRegionManager const * const input )
  {
    m_elemRegionManager = input;
  }

  /**
   * @brief Get the ElementRegionManager.
   * @return The current ElementRegionManager.
   */
  ElementRegionManager const *
  getElementRegionManager() const
  {
    return m_elemRegionManager;
  }

  /// The relationship between object indices and element regions.
  BASETYPE m_toElementRegion;
  /// The relationship between object indices and element subregions.
  BASETYPE m_toElementSubRegion;
  /// The relationship between object indices and element indices.
  BASETYPE m_toElementIndex;

  /// The current ElementRegionManager
  ElementRegionManager const * m_elemRegionManager;
};

/// @cond DO_NOT_DOCUMENT
template< typename BASETYPE >
ToElementRelation< BASETYPE >::ToElementRelation() :
  m_toElementRegion(),
  m_toElementSubRegion(),
  m_toElementIndex(),
  m_elemRegionManager( nullptr )
{}

template< typename BASETYPE >
ToElementRelation< BASETYPE >::~ToElementRelation()
{}

template< typename BASETYPE >
template< typename... DIMS >
void
ToElementRelation< BASETYPE >::resize( DIMS... newdims )
{
  m_toElementRegion.resize( newdims... );
  m_toElementSubRegion.resize( newdims... );
  m_toElementIndex.resize( newdims... );
}
/// @endcond

/// @brief A ToElementRelation where each object is related to the same number of elements.
typedef ToElementRelation< array2d< localIndex > > FixedToManyElementRelation;

/// @brief A ToElementRelation where each object is related to an arbitrary number of elements.
typedef ToElementRelation< ArrayOfArrays< localIndex > > OrderedVariableToManyElementRelation;

/**
 * @brief Remove an element relation from an object in the relation.
 * @param relation The relationship mapping to remove a single element relation from
 *                 a single object from.
 * @param firstIndex The object index to remove an element relation from.
 * @param er The element region to remove.
 * @param esr The element subregion to remove.
 * @param ei The element index to remove.
 */
void
erase( OrderedVariableToManyElementRelation & relation,
       localIndex const firstIndex,
       localIndex const er,
       localIndex const esr,
       localIndex const ei );

/**
 * @brief Insert an element relation for an object in the relation.
 * @param relation The relationship mapping to insert a single element relation for
 *                 a single object into.
 * @param firstIndex The object index to insert an element relation from.
 * @param er The element region to insert.
 * @param esr The element subregion to insert.
 * @param ei The element index to insert.
 */
void
insert( OrderedVariableToManyElementRelation & relation,
        localIndex const firstIndex,
        localIndex const er,
        localIndex const esr,
        localIndex const ei );

} /* namespace geosx */

#endif /* GEOSX_MESH_TOELEMENTRELATION_HPP_ */
