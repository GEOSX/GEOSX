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
 * @file RelativePermeabilityBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP
#define GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{
namespace constitutive
{
class RelativePermeabilityBaseUpdate
{
public:
  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOSX_HOST_DEVICE
  localIndex
  numElems() const
  {
    return m_phaseRelPerm.size( 0 );
  }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOSX_HOST_DEVICE
  localIndex
  numGauss() const
  {
    return m_phaseRelPerm.size( 1 );
  }

  /**
   * @brief Get number of fluid phases.
   * @return number of phases
   */
  GEOSX_HOST_DEVICE
  localIndex
  numPhases() const
  {
    return m_phaseTypes.size();
  }

protected:
  RelativePermeabilityBaseUpdate( arrayView1d< integer const > const & phaseTypes,
                                  arrayView1d< integer const > const & phaseOrder,
                                  arrayView3d< real64 > const & phaseRelPerm,
                                  arrayView4d< real64 > const & dPhaseRelPerm_dPhaseVolFrac ) :
    m_phaseTypes( phaseTypes ),
    m_phaseOrder( phaseOrder ),
    m_phaseRelPerm( phaseRelPerm ),
    m_dPhaseRelPerm_dPhaseVolFrac( dPhaseRelPerm_dPhaseVolFrac )
  {}

  /// Default copy constructor
  RelativePermeabilityBaseUpdate( RelativePermeabilityBaseUpdate const & ) = default;

  /// Default move constructor
  RelativePermeabilityBaseUpdate( RelativePermeabilityBaseUpdate && ) = default;

  /// Deleted copy assignment operator
  RelativePermeabilityBaseUpdate &
  operator=(
    RelativePermeabilityBaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  RelativePermeabilityBaseUpdate &
  operator=( RelativePermeabilityBaseUpdate && ) = delete;

  arrayView1d< integer const > m_phaseTypes;
  arrayView1d< integer const > m_phaseOrder;

  arrayView3d< real64 > m_phaseRelPerm;
  arrayView4d< real64 > m_dPhaseRelPerm_dPhaseVolFrac;

private:
  GEOSX_HOST_DEVICE
  virtual void
  Compute(
    arraySlice1d< real64 const > const & phaseVolFraction,
    arraySlice1d< real64 > const & phaseRelPerm,
    arraySlice2d< real64 > const & dPhaseRelPerm_dPhaseVolFrac ) const = 0;

  GEOSX_HOST_DEVICE
  virtual void
  Update( localIndex const k,
          localIndex const q,
          arraySlice1d< real64 const > const & phaseVolFraction ) const = 0;
};

class RelativePermeabilityBase : public ConstitutiveBase
{
public:
  struct PhaseType
  {
    static constexpr integer OIL = 0;
    static constexpr integer GAS = 1;
    static constexpr integer WATER = 2;
    static constexpr integer MAX_NUM_PHASES = 3;
  };

  // order of the phase properties in the water-oil data
  struct WaterOilPairPhaseType
  {
    static constexpr integer WATER = 0;  // first water phase property
    static constexpr integer OIL = 1;    // second oil phase property
  };

  // order of the phase properties in the gas-oil data
  struct GasOilPairPhaseType
  {
    static constexpr integer GAS = 0;  // first gas phase property
    static constexpr integer OIL = 1;  // second oil phase property
  };

  RelativePermeabilityBase( std::string const & name,
                            dataRepository::Group * const parent );

  virtual ~RelativePermeabilityBase() override;

  void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  virtual void
  AllocateConstitutiveData(
    dataRepository::Group * const parent,
    localIndex const numConstitutivePointsPerParentIndex ) override;

  localIndex
  numFluidPhases() const
  {
    return m_phaseNames.size();
  }

  arrayView1d< string const > const &
  phaseNames() const
  {
    return m_phaseNames;
  }

  arrayView3d< real64 const > const &
  phaseRelPerm() const
  {
    return m_phaseRelPerm;
  }
  arrayView4d< real64 const > const &
  dPhaseRelPerm_dPhaseVolFraction() const
  {
    return m_dPhaseRelPerm_dPhaseVolFrac;
  }

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto phaseNamesString = "phaseNames";
    static constexpr auto phaseTypesString = "phaseTypes";
    static constexpr auto phaseOrderString = "phaseTypes";

    static constexpr auto phaseRelPermString = "phaseRelPerm";  // Kr
    static constexpr auto dPhaseRelPerm_dPhaseVolFractionString =
      "dPhaseRelPerm_dPhaseVolFraction";  // dKr_p/dS_p
  } viewKeysRelativePermeabilityBase;

protected:
  virtual void
  PostProcessInput() override;

  /**
   * @brief Function called internally to resize member arrays
   * @param size primary dimension (e.g. number of cells)
   * @param numPts secondary dimension (e.g. number of gauss points per cell)
   */
  void
  ResizeFields( localIndex const size, localIndex const numPts );

  // phase names read from input
  string_array m_phaseNames;

  // phase ordering info
  array1d< integer > m_phaseTypes;
  array1d< integer > m_phaseOrder;

  // output quantities
  array3d< real64 > m_phaseRelPerm;
  array4d< real64 > m_dPhaseRelPerm_dPhaseVolFrac;
};

}  // namespace constitutive

}  // namespace geosx

#endif  //GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP
