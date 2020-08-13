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
 * @file SinglePhaseProppantBase.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEPROPPANTBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEPROPPANTBASE_HPP_

#include "SinglePhaseBase.hpp"

namespace geosx
{
class SinglePhaseProppantBase : public SinglePhaseBase
{
public:
  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseProppantBase( const std::string & name, Group * const parent );

  SinglePhaseProppantBase() = delete;

  /// deleted copy constructor
  SinglePhaseProppantBase( SinglePhaseProppantBase const & ) = delete;

  /// default move constructor
  SinglePhaseProppantBase( SinglePhaseProppantBase && ) = default;

  /// deleted assignment operator
  SinglePhaseProppantBase &
  operator=( SinglePhaseProppantBase const & ) = delete;

  /// deleted move operator
  SinglePhaseProppantBase &
  operator=( SinglePhaseProppantBase && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseProppantBase();

  virtual void
  UpdateFluidModel( Group & dataGroup,
                    localIndex const targetIndex ) const override;

protected:
  virtual void
  ValidateFluidModels( DomainPartition const & domain ) const override;

  virtual FluidPropViews
  getFluidProperties(
    constitutive::ConstitutiveBase const & fluid ) const override;

  virtual arrayView1d< real64 const > const &
  getPoreVolumeMult(
    ElementSubRegionBase const & subRegion ) const override;

private:
  virtual void
  ResetViewsPrivate( ElementRegionManager const & elemManager ) override;
};
}  // namespace geosx
#endif /* GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEPROPPANTBASE_HPP_ */
