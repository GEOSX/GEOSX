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
 * @file FiniteVolumeManager.hpp
 *
 */

#ifndef GEOSX_FINITEVOLUME_FINITEVOLUMEMANAGER_HPP_
#define GEOSX_FINITEVOLUME_FINITEVOLUMEMANAGER_HPP_

#include "dataRepository/Group.hpp"

namespace geosx
{
class DomainPartition;
class FluxApproximationBase;

/**
 * @class FiniteVolumeManager
 *
 * Class managing a finite volume discretization.
 */
class FiniteVolumeManager : public dataRepository::Group
{
public:
  /**
   * @brief Deleted default constructor.
   */
  FiniteVolumeManager() = delete;

  /**
   * @brief Constructor.
   * @param name the name of the FiniteVolumeManager in the data repository
   * @param parent the parent group of this group.
   */
  FiniteVolumeManager( string const & name, Group * const parent );

  /**
   * @brief Destructor.
   */
  virtual ~FiniteVolumeManager() override;

  virtual Group *
  CreateChild( string const & childKey,
               string const & childName ) override;

  virtual void
  ExpandObjectCatalogs() override;

  /**
   * @brief Return the FluxApproximation associated with the provided name.
   * @param[in] name the provided name
   * @return the FluxApproximation associated with the provided name
   */
  FluxApproximationBase const &
  getFluxApproximation( string const & name ) const;

  /**
   * @copydoc getFluxApproximation( string const & ) const
   */
  FluxApproximationBase &
  getFluxApproximation( string const & name );

private:
};

}  // namespace geosx

#endif  //GEOSX_FINITEVOLUME_FINITEVOLUMEMANAGER_HPP_
