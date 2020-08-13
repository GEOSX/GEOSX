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

#ifndef GEOSX_PHYSICSSOLVERS_PHYSICSSOLVERMANAGER_HPP_
#define GEOSX_PHYSICSSOLVERS_PHYSICSSOLVERMANAGER_HPP_

#include "dataRepository/Group.hpp"

namespace pugi
{
class xml_node;
}

namespace geosx
{
class SolverBase;

class PhysicsSolverManager : public dataRepository::Group
{
public:
  PhysicsSolverManager( std::string const & name, Group * const parent );

  virtual ~PhysicsSolverManager() override;

  virtual Group *
  CreateChild( string const & childKey,
               string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void
  ExpandObjectCatalogs() override;

  struct viewKeyStruct
  {
    constexpr static auto gravityVectorString = "gravityVector";
  } viewKeys;

  struct groupKeyStruct
  {
  } groupKeys;

  R1Tensor const &
  gravityVector() const
  {
    return m_gravityVector;
  }
  R1Tensor &
  gravityVector()
  {
    return m_gravityVector;
  }

private:
  PhysicsSolverManager() = delete;

  R1Tensor m_gravityVector;
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_PHYSICSSOLVERMANAGER_HPP_ */
