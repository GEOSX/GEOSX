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
 * @file Perforation.hpp
 */

#ifndef GEOSX_MESHUTILITIES_PERFORATION_HPP
#define GEOSX_MESHUTILITIES_PERFORATION_HPP

#include "dataRepository/Group.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
static constexpr auto perforation = "Perforation";
}
}  // namespace dataRepository

/**
 * @class Perforation
 *
 * This class describes a perforation with its location, well transmissibility  and corresponding well element
 */
class Perforation : public dataRepository::Group
{
public:
  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for Perforation Objects.
   * @param[in] name the name of this instantiation of Perforation in the repository
   * @param[in] parent the parent group of this instantiation of Perforation
   */
  explicit Perforation( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Default destructor.
   */
  ~Perforation() override;

  /**
   * @brief Deleted default constructor.
   */
  Perforation() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  Perforation( Perforation const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  Perforation( Perforation && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a perforation object
   */
  Perforation &
  operator=( Perforation const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a perforation object
   */
  Perforation &
  operator=( Perforation && ) = delete;

  ///@}

  /**
   * @name Getters
   */
  ///@{

  /**
   * @brief Get the linear distance between the well head and the perforation.
   * @return the distance between the well head and the perforation
   */
  real64 const &
  GetDistanceFromWellHead() const
  {
    return m_distanceFromHead;
  }

  /**
   * @brief Get the well Peaceman index at the perforation.
   * @return the well transmissibility
   */
  real64
  GetWellTransmissibility() const
  {
    return m_wellTransmissibility;
  }

  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the linear distance from well head
    static constexpr auto distanceFromHeadString = "distanceFromHead";
    /// String key for the well transmissibility at this perforation
    static constexpr auto wellTransmissibilityString = "transmissibility";
    /// ViewKey for the linear distance from well head
    dataRepository::ViewKey distanceFromHead = { distanceFromHeadString };
    /// ViewKey for the well transmissibility at this perforation
    dataRepository::ViewKey wellTransmissibility = { wellTransmissibilityString };
  }
  /// ViewKey struct for the Perforation class
  viewKeysPerforation;

protected:
  void
  PostProcessInput() override;

private:
  /// Linear distance from well head
  real64 m_distanceFromHead;

  /// Well transmissibility at this perforation
  real64 m_wellTransmissibility;
};

}  //namespace geosx

#endif  //GEOSX_MESHUTILITIES_PERFORATION_HPP
