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
 * @file PerforationData.hpp
 */

#ifndef GEOSX_MESHUTILITIES_PERFORATIONDATA_HPP
#define GEOSX_MESHUTILITIES_PERFORATIONDATA_HPP

#include "dataRepository/Group.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "mesh/ToElementRelation.hpp"
#include "meshUtilities/InternalWellGenerator.hpp"

namespace geosx
{
class DomainPartition;
class MeshLevel;
class WellElementSubRegion;
class CellBlock;

/**
 * @class PerforationData
 * This class keeps track of all the local perforations on this rank
 */
class PerforationData : public ObjectManagerBase
{
public:
  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for PerforationData Objects.
   * @param[in] name name of this instantiation of PerforationData in the repository
   * @param[in] parent parent group of this instantiation of PerforationData
   */
  explicit PerforationData( string const & name,
                            dataRepository::Group * const parent );

  /**
   * @brief Default destructor.
   */
  ~PerforationData() override;

  /**
   * @brief Deleted default constructor.
   */
  PerforationData() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  PerforationData( PerforationData const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  PerforationData( PerforationData && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a PerforationData object
   */
  PerforationData &
  operator=( PerforationData const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a PerforationData object
   */
  PerforationData &
  operator=( PerforationData && ) = delete;

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
    return "PerforationData";
  }

  /**
   * @copydoc CatalogName()
   */
  virtual const string
  getCatalogName() const override
  {
    return CatalogName();
  }

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Set the global number of perforations used for well initialization.
   * @param[in] nPerfs global number of perforations (obtained from InternalWellGenerator)
   */
  void
  SetNumPerforationsGlobal( globalIndex nPerfs )
  {
    m_numPerforationsGlobal = nPerfs;
  }

  /**
   * @brief Get the global number of perforations (used for well initialization).
   * @return global number of perforations
   */
  globalIndex
  GetNumPerforationsGlobal() const
  {
    return m_numPerforationsGlobal;
  }

  /**
   * @brief Get perforation-to-mesh-element connectivity.
   * @return list of element region/subregion/index connected to each perforation
   */
  ToElementRelation< array1d< localIndex > > &
  GetMeshElements()
  {
    return m_toMeshElements;
  }

  /**
   * @brief Provide an immutable accessor to a const perforation-to-mesh-element connectivity.
   * @return list of element region/subregion/index connected to each perforation
   */
  ToElementRelation< array1d< localIndex > > const &
  GetMeshElements() const
  {
    return m_toMeshElements;
  }

  /**
   * @brief Get perforation-to-well-element connectivity.
   * @return list of well element index connected to each perforation
   */
  arrayView1d< localIndex > &
  GetWellElements()
  {
    return m_wellElementIndex;
  }

  /**
   * @brief Provide an immutable accessor to a const perforation-to-well-element connectivity.
   * @return list of well element index connected to each perforation
   */
  arrayView1d< localIndex const > const &
  GetWellElements() const
  {
    return m_wellElementIndex;
  }

  /**
   * @brief Get perforation locations.
   * @return list of perforation locations
   */
  arrayView1d< R1Tensor > &
  GetLocation()
  {
    return m_location;
  }

  /**
   * @brief Provide an immutable accessor to a const perforation location arrayView.
   * @return list of perforation locations
   */
  arrayView1d< R1Tensor const > const &
  GetLocation() const
  {
    return m_location;
  }

  /**
   * @brief Provide an immutable accessor to a const perforation well indices array.
   * @return list of perforation well indices
   */
  arrayView1d< real64 const > const &
  GetWellTransmissibility() const
  {
    return m_wellTransmissibility;
  }

  /**
   * @brief Get perforation well indices.
   * @return list of perforation well indices
   */
  arrayView1d< real64 > &
  GetWellTransmissibility()
  {
    return m_wellTransmissibility;
  }

  ///@}

  /**
   * @name Well transmissibility computation
   */
  ///@{

  /**
   * @brief Compute the well transmissibility for each local perforation on this well.
   * @param[in] mesh target mesh level
   * @param[in] wellElemSubRegion  subRegion corresponding to this well
   * @param[in] permeabilityKey key to access the permeability in the reservoir
   */
  void
  ComputeWellTransmissibility( MeshLevel const & mesh,
                               WellElementSubRegion const * const wellElemSubRegion,
                               string const & permeabilityKey );

  ///@}

  /**
   * @name Construction of the connectivity
   */
  ///@{

  /**
   * @brief Locate connected local mesh elements and resizes current object appropriately.
   * @param[in] mesh target mesh level
   * @param[in] wellGeometry  InternalWellGenerator containing the global well topology
   */
  void
  ConnectToMeshElements( MeshLevel const & mesh,
                         InternalWellGenerator const & wellGeometry );

  /**
   * @brief Connect each perforation to a local wellbore element.
   * @param[in] wellGeometry InternalWellGenerator containing the global well topology
   * @param[in] globalToLocalWellElementMap  global-to-local map of wellbore elements
   * @param[in] elemOffsetGlobal the offset of the first global well element ( = offset of last global mesh elem + 1 )
   */
  void
  ConnectToWellElements(
    InternalWellGenerator const & wellGeometry,
    unordered_map< globalIndex, localIndex > const & globalToLocalWellElementMap,
    globalIndex elemOffsetGlobal );

  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {
    /// String key for the global number of perforations
    static constexpr auto numPerforationsGlobalString = "numPerforationsGlobal";
    /// String key for the reservoir element region index
    static constexpr auto reservoirElementRegionString =
      "reservoirElementRegion";
    /// String key for the reservoir element subregion index
    static constexpr auto reservoirElementSubregionString =
      "reservoirElementSubregion";
    /// String key for the reservoir element index
    static constexpr auto reservoirElementIndexString = "reservoirElementIndex";
    /// String key for the well element index
    static constexpr auto wellElementIndexString = "wellElementIndex";
    /// String key for the perforation location
    static constexpr auto locationString = "location";
    /// String key for the well transmissibility
    static constexpr auto wellTransmissibilityString = "wellTransmissibility";

    /// ViewKey for the global number of perforations
    dataRepository::ViewKey numPerforationsGlobal = { numPerforationsGlobalString };
    /// ViewKey for the reservoir element region index
    dataRepository::ViewKey reservoirElementRegion = {
      reservoirElementRegionString };
    /// ViewKey for the reservoir element subregion index
    dataRepository::ViewKey reservoirElementSubregion = {
      reservoirElementSubregionString };
    /// ViewKey for the reservoir element index
    dataRepository::ViewKey reservoirElementIndex = { reservoirElementIndexString };
    /// ViewKey for the well element index
    dataRepository::ViewKey wellElementIndex = { wellElementIndexString };
    /// ViewKey for the well location
    dataRepository::ViewKey location = { locationString };
    /// ViewKey for the well transmissibility
    dataRepository::ViewKey wellTransmissibility = { wellTransmissibilityString };

  }
  /// ViewKey struct for the PerforationData class
  viewKeysPerforationData;

  /**
   * @brief struct to serve as a container for group strings and keys
   * @struct groupKeyStruct
   */
  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
  }
  /// groupKey struct for the PerforationData class
  groupKeysPerforationData;

private:
  /**
   * @name Helpers for transmissibility computation
   */
  ///@{

  /**
   * @brief Compute the approximate dimensions of the reservoir element containing a perforation.
   *        This is done by computing a bounding box containing the element.
   * @param[in] mesh target mesh level
   * @param[in] er  index of the element region containing the reservoir element
   * @param[in] esr  index of the element subRegion containing the reservoir element
   * @param[in] ei  index of the reservoir element
   * @param[out] dx dimension of the element in the x-direction
   * @param[out] dy dimension of the element in the y-direction
   * @param[out] dz dimension of the element in the z-direction
   */
  void
  GetReservoirElementDimensions( MeshLevel const & mesh,
                                 localIndex const er,
                                 localIndex const esr,
                                 localIndex const ei,
                                 real64 & dx,
                                 real64 & dy,
                                 real64 & dz ) const;

  /**
   * @brief Check if the well is along the x-, y-, or z- directions.
   * @param[in] vecWellElemCenterToPerf vector connecting the well element center to the perforation
   * @param[in] dx dimension of the element in the x-direction
   * @param[in] dy dimension of the element in the y-direction
   * @param[in] dz dimension of the element in the z-direction
   * @param[in] perm absolute permeability in the reservoir element
   * @param[out] d1 dimension of the element in the first direction
   * @param[out] d2 dimension of the element in the second direction
   * @param[out] h dimension of the element in the third direction
   * @param[out] k1 absolute permeability in the reservoir element (first direction)
   * @param[out] k2 absolute permeability in the reservoir element (second direction)
   */
  void
  DecideWellDirection( R1Tensor const & vecWellElemCenterToPerf,
                       real64 const & dx,
                       real64 const & dy,
                       real64 const & dz,
                       R1Tensor const & perm,
                       real64 & d1,
                       real64 & d2,
                       real64 & h,
                       real64 & k1,
                       real64 & k2 ) const;

  ///@}

  /// @cond DO_NOT_DOCUMENT
  void
  DebugLocalPerforations() const;
  /// @endcond

  /// Global number of perforations
  globalIndex m_numPerforationsGlobal;

  /// Indices of the mesh elements connected to perforations
  ToElementRelation< array1d< localIndex > > m_toMeshElements;

  /// Indices of the well elements to which perforations are attached
  array1d< localIndex > m_wellElementIndex;

  /// Location of the perforations
  array1d< R1Tensor > m_location;

  /// Well transmissibility at the perforations
  array1d< real64 > m_wellTransmissibility;
};

}  //namespace geosx

#endif  //GEOSX_MESHUTILITIES_PERFORATIONDATA_HPP
