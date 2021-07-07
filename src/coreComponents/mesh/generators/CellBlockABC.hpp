/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2020-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_CELLBLOCKABC_HPP
#define GEOSX_CELLBLOCKABC_HPP

#include "dataRepository/Group.hpp"
#include "common/DataTypes.hpp"

#include <vector>

namespace geosx
{

/**
 * Abstract base class defining the information provided by any cell block implementation.
 * Mainly element type related information (number of faces per element...),
 * some element to nodes, faces and edges mappings, a local to global mapping and
 * some kind of accessors to external properties.
 *
 * The cells of CellBlockABC are all of the same kind.
 *
 * It's noteworthy that the CellBlockABC is immutable oriented.
 * The derived implementations need to have the modification/creation capabilities.
 */
class CellBlockABC : public dataRepository::Group
{
public:

  /**
   * @brief Constructor
   * @param name The name of this Group.
   * @param parent The parent Group.
   */
  CellBlockABC( string const & name,
                Group * const parent )
    :
    Group( name, parent )
  { }

  /**
   * @brief Get the type of element in this subregion.
   * @return a string specifying the type of element in this subregion
   *
   * See class FiniteElementBase for possible element type.
   */
  virtual string getElementTypeString() const = 0;

  /**
   * @brief Get the number of nodes per element.
   * @return number of nodes per element
   */
  virtual localIndex numNodesPerElement() const = 0;

  /**
   * @brief Get the number of faces per element.
   * @return number of faces per element
   */
  virtual localIndex numFacesPerElement() const = 0;

  /**
   * @brief Get the number of elements.
   * @return number of elements in the cell block
   */
  virtual localIndex numElements() const = 0;

  /**
   * @brief Get the element-to-node map.
   * @return a reference to the element-to-node map
   */
  virtual array2d< localIndex, cells::NODE_MAP_PERMUTATION > const & getElemToNode() const = 0;

  /**
   * @brief Get the element-to-edges map.
   * @return A view of the mapping.
   */
  virtual array2d< localIndex > const & getElemToEdges() const = 0;

  /**
   * @brief Get the element-to-faces map.
   * @return A view of the mapping.
   */
  virtual array2d< localIndex > const & getElemToFaces() const = 0;

  /**
   * @brief Get local to global map.
   * @return The mapping relationship as a array.
   */
  virtual arrayView1d< globalIndex const > localToGlobalMap() const = 0;

  /**
   * @brief Helper function to apply a lambda function over all the external properties of the subregion
   * @tparam LAMBDA the type of the lambda function
   * @param lambda lambda function that is applied to the wrappers of external properties
   */
  template< typename LAMBDA >
  void forExternalProperties( LAMBDA && lambda )
  {
    for( auto * wrapperBase: this->getExternalProperties() )
    {
      lambda( *wrapperBase );
    }
  } // TODO This is not abstract

private:
  /**
   * @brief Returns the external properties under the form of WrapperBase pointer to instances.
   * @return An iterable of pointers
   * @note I quite don't get why this cannot be const.
   *
   * In order not to expose the implementation details (external properties being stored as WrapperBase),
   * this abstract member function is made private.
   * Thus the list of pointers shall not be used anyhow by end-users that must use @see forExternalProperties.
   */
  virtual std::list< dataRepository::WrapperBase * > getExternalProperties() = 0;
};

}

#endif //GEOSX_CELLBLOCKABC_HPP
