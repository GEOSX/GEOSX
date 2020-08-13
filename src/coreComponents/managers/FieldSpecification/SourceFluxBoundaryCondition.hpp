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
 * SourceFluxBoundaryCondition.hpp
 *
 */

#ifndef GEOSX_MANAGER_FIELDSPECIFICATION_SOURCEFLUXBOUNDARYCONDITION_HPP_
#define GEOSX_MANAGER_FIELDSPECIFICATION_SOURCEFLUXBOUNDARYCONDITION_HPP_

#include "FieldSpecificationBase.hpp"

namespace geosx
{
/**
 * @class SourceFluxBoundaryCondition
 * A class to manage Neumann boundary conditions
 */
class SourceFluxBoundaryCondition : public FieldSpecificationBase
{
public:
  /**
   * @brief constructor
   * @param name the name of the FieldSpecificationBase in the data repository
   * @param parent the parent group of this group.
   */
  SourceFluxBoundaryCondition( string const & name,
                               dataRepository::Group * const parent );

  /**
   * @brief destructor
   */
  SourceFluxBoundaryCondition() = delete;

  /**
   * @brief destructor
   */
  virtual ~SourceFluxBoundaryCondition() override;

  /**
   * @brief Called by Initialize() after to initializing sub-Groups.
   * @param[in] rootGroup A group that is passed in to the initialization functions
   *                      in order to facilitate the initialization.
   */
  virtual void
  InitializePreSubGroups( Group * const rootGroup ) override;

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string
  CatalogName()
  {
    return "SourceFlux";
  }

  virtual const string
  getCatalogName() const override
  {
    return SourceFluxBoundaryCondition::CatalogName();
  }
};

} /* namespace geosx */

#endif /* GEOSX_MANAGER_FIELDSPECIFICATION_SOURCEFLUXBOUNDARYCONDITION_HPP_ */
