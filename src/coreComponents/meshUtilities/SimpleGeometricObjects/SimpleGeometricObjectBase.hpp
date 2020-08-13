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

/*
 * @file SimpleGeometricObjectBase.hpp
 */

#ifndef SIMPLEGEOMETRICOBJECTS_H_
#define SIMPLEGEOMETRICOBJECTS_H_

//#include "common/Common.h"
#include "dataRepository/Group.hpp"
#include "codingUtilities/StringUtilities.hpp"
#include "dataRepository/ObjectCatalog.hpp"

class Function;

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const geometricObjects( "GeometricObjects" );
}
}  // namespace dataRepository

/**
 * @class SimpleGeometricObjectBase
 * @brief Base class for the geometric objects (box, plane, cylinder).
 */
class SimpleGeometricObjectBase : public dataRepository::Group
{
public:
  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  explicit SimpleGeometricObjectBase( std::string const & name,
                                      Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~SimpleGeometricObjectBase();

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static string
  CatalogName()
  {
    return "SimpleGeometricObjectBase";
  }

  /**
   * @brief Type alias for catalog interface used by this class.
   */
  using CatalogInterface =
    dataRepository::CatalogInterface< SimpleGeometricObjectBase,
                                      std::string const &,
                                      Group * const >;

  /**
   * @copydoc CatalogName()
   */
  static CatalogInterface::CatalogType &
  GetCatalog();

  ///@}

  /**
   * @brief Check if the input coordinates are in the object.
   * @param[in] coord the coordinates to test
   * @return true if the coordinates are in the object, false otherwise
   */
  virtual bool
  IsCoordInObject( const R1Tensor & coord ) const = 0;
};

}  // namespace geosx
#endif /* SIMPLEGEOMETRICOBJECTS_H_ */
