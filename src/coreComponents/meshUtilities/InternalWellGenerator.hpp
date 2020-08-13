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
 * @file InternalWellGenerator.hpp
 *
 */

#ifndef GEOSX_MESHUTILITIES_INTERNALWELLGENERATOR_HPP_
#define GEOSX_MESHUTILITIES_INTERNALWELLGENERATOR_HPP_

#include "dataRepository/Group.hpp"
#include "meshUtilities/MeshGeneratorBase.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const nodeCoords = "polylineNodeCoords";
string const segmentConn = "polylineSegmentConn";
string const nElems = "numElementsPerSegment";
string const radius = "radius";
string const wellRegionName = "wellRegionName";
string const wellControlsName = "wellControlsName";
string const meshBodyName = "meshName";
}  // namespace keys
}  // namespace dataRepository

/**
 * @class InternalWellGenerator
 *
 * This class processes the data of a single well from the XML and generates the well geometry
 */
class InternalWellGenerator : public MeshGeneratorBase
{
public:
  /**
   * @brief Struct to define the top and bottom node of a segment.
   */
  struct NodeLocation
  {
    static constexpr integer TOP = 0;    /**< Top of the well */
    static constexpr integer BOTTOM = 1; /**< Bottom of the well */
  };

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  InternalWellGenerator( const std::string & name, Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~InternalWellGenerator() override;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this type in the catalog
   */
  static string
  CatalogName()
  {
    return "InternalWell";
  }

  ///@}

  /**
   * @name Overriding functions defined in MeshGeneratorBase and above
   */
  ///@{

  /// not implemented
  virtual void
  GenerateElementRegions( DomainPartition & ) override
  {}

  /**
   * @brief Creates a new sub-Group using the ObjectCatalog functionality.
   * @param[in] childKey The name of the new object type's key in the
   *                     ObjectCatalog.
   * @param[in] childName The name of the new object in the collection of
   *                      sub-Groups.
   * @return A pointer to the new Group created by this function.
   */
  virtual Group *
  CreateChild( string const & childKey,
               string const & childName ) override;

  /**
   * @brief Expand any catalogs in the data structure.
   */
  virtual void
  ExpandObjectCatalogs() override;

  /**
   * @brief Main function of the class that generates the well geometry
   * @param[in] domain the domain object
   */
  virtual void
  GenerateMesh( DomainPartition * const domain ) override;

  /// not implemented
  virtual void
  GetElemToNodesRelationInBox( std::string const &,
                               int const *,
                               int const &,
                               int *,
                               int const ) override
  {}

  /// not implemented
  virtual void
  RemapMesh( dataRepository::Group * const ) override
  {}

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  // getters for element data

  /**
   * @brief Get the global number of well elements.
   * @return the global number of elements
   */
  globalIndex
  GetNumElements() const
  {
    return m_numElems;
  }

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  arrayView1d< R1Tensor const > const &
  GetElemCoords() const
  {
    return m_elemCenterCoords;
  }

  /**
   * @brief Get the global indices mapping an element to the next.
   * @return list providing the global index of the next element for each element
   */
  arrayView1d< globalIndex const > const &
  GetNextElemIndex() const
  {
    return m_nextElemId;
  }

  /**
   * @brief Get the global indices mapping an element to the previous ones.
   * @return list providing the global indices of the previous elements for each element
   */
  arrayView1d< arrayView1d< globalIndex const > const > const &
  GetPrevElemIndices() const
  {
    return m_prevElemId.toViewConst();
  }

  /**
   * @brief Get the global indices of the well nodes nodes connected to each element.
   * @return list providing the global index of the well nodes for each well element
   */
  arrayView2d< globalIndex const > const &
  GetElemToNodesMap() const
  {
    return m_elemToNodesMap;
  }

  /**
   * @brief Get the volume of the well elements.
   * @return list of volumes of the well elements
   */
  arrayView1d< real64 const > const &
  GetElemVolume() const
  {
    return m_elemVolume;
  }

  /**
   * @brief Get the radius in the well.
   * @return the radius in the well
   */
  real64
  GetElementRadius() const
  {
    return m_radius;
  }

  // getters for node data

  /**
   * @brief Get the global number of well nodes.
   * @return the global number of nodes
   */
  globalIndex
  GetNumNodes() const
  {
    return m_numNodes;
  }

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  arrayView1d< R1Tensor const > const &
  GetNodeCoords() const
  {
    return m_nodeCoords;
  }

  // getters for perforation data

  /**
   * @brief Get the global number of perforations on this well.
   * @return the global number of elements
   */
  globalIndex
  GetNumPerforations() const
  {
    return m_numPerforations;
  }

  /**
   * @brief Get the locations of the perforations.
   * @return list of locations of all the perforations on the well
   */
  arrayView1d< R1Tensor const > const &
  GetPerfCoords() const
  {
    return m_perfCoords;
  }

  /**
   * @brief Get the well transmissibility at the perforations.
   * @return list of well transmissibility at all the perforations on the well
   */
  arrayView1d< real64 const > const &
  GetPerfTransmissibility() const
  {
    return m_perfTransmissibility;
  }

  /**
   * @brief Get the global indices of the well elements connected to each perforation.
   * @return list providing the global index of the connected well element for each perforation
   */
  arrayView1d< globalIndex const > const &
  GetPerfElemIndex() const
  {
    return m_perfElemId;
  }

  ///@}

protected:
  /**
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  void
  PostProcessInput() override final;

private:
  /**
   * @name Helper functions to construct the geometry of the well
   */
  ///@{

  /**
   * @brief Map each polyline node to the polyline segment(s) it is connected to.
   */
  void
  ConstructPolylineNodeToSegmentMap();

  /**
   * @brief Find the head node of the well (i.e., top node of the polyline).
   */
  void
  FindPolylineHeadNodeIndex();

  /**
   * @brief Discretize the polyline by placing well elements.
   */
  void
  DiscretizePolyline();

  /**
   * @brief Map each perforation to a well element.
   */
  void
  ConnectPerforationsToWellElements();

  /**
   * @brief Merge perforations on the elements with multiple perforations.
   */
  void
  MergePerforations();

  /**
   * @brief At a given node, find the next segment going in the direction of the bottom of the well.
   * @param[in] topSegId index of the top segment
   * @param[in] currentNodeId index of the current node
   */
  globalIndex
  GetNextSegmentIndex( globalIndex topSegId,
                       globalIndex currentNodeId ) const;

  ///@}

  /// @cond DO_NOT_DOCUMENT
  void
  DebugWellGeometry() const;
  /// @endcond

  // XML Input

  /// Coordinates of the polyline
  array2d< real64 > m_inputPolyNodeCoords;

  /// Connectivity between the polyline nodes
  array2d< globalIndex > m_segmentToPolyNodeMap;

  /// Number of well elements per polyline interval
  int m_numElemsPerSegment;

  /// Radius area of the well (assumed to be valid for the entire well)
  real64 m_radius;

  /// Name of the corresponding well region
  string m_wellRegionName;

  /// Name of the constraints associated with this well
  string m_wellControlsName;

  /// Name of the mesh body associated with this well
  string m_meshBodyName;

  // Geometry of the well (later passed to the WellElementSubRegion)

  // well element data

  /// Global number of well elements
  globalIndex m_numElems;

  /// Physical location of the center of the well element
  array1d< R1Tensor > m_elemCenterCoords;

  /// Global index of the next well element
  array1d< globalIndex > m_nextElemId;

  /// Global indices of the prev well elements (maybe need multiple prevs for branching)
  array1d< array1d< globalIndex > > m_prevElemId;

  /// Connectivity between elements and nodes
  array2d< globalIndex > m_elemToNodesMap;

  /// Volume of well elements
  array1d< real64 > m_elemVolume;

  // well node data

  /// Number of nodes per well element
  globalIndex const m_numNodesPerElem;

  /// Global number of well nodes
  globalIndex m_numNodes;

  /// Physical location of the nodes
  array1d< R1Tensor > m_nodeCoords;

  // perforation data

  /// Global number of perforations
  globalIndex m_numPerforations;

  /// Absolute physical location of the perforation
  array1d< R1Tensor > m_perfCoords;

  /// Well Peaceman index at the perforation
  array1d< real64 > m_perfTransmissibility;

  /// Global index of the well element
  array1d< globalIndex > m_perfElemId;

  // Auxiliary data

  // Number of physical dimensions
  const int m_nDims;

  /// Coordinates of the polyline nodes in R1Tensor format
  array1d< R1Tensor > m_polyNodeCoords;

  /// Map from the polyline nodes to the polyline nodes
  array1d< SortedArray< globalIndex > > m_polyNodeToSegmentMap;

  /// Index of the node at the well head
  globalIndex m_polylineHeadNodeId;

  /// Physical location of the polyline node wrt to well head
  array1d< real64 > m_nodeDistFromHead;

  // Perforation data

  /// List of perforation names
  string_array m_perforationList;

  /// Physical location of the perforation wrt to well head
  array1d< real64 > m_perfDistFromHead;
};
}  // namespace geosx

#endif /* GEOSX_MESHUTILITIES_INTERNALWELLGENERATOR_HPP_ */
