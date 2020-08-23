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
 * @file TwoPointFluxApproximationWithGraph.hpp
 */

#ifndef GEOSX_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATIONWITHGRAPH_HPP_
#define GEOSX_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATIONWITHGRAPH_HPP_

#include "finiteVolume/TwoPointFluxApproximation.hpp"
#include "mesh/GraphFromText.hpp"

namespace geosx
{

/**
 * @class TwoPointFluxApproximationWithGraph
 *
 * Provides management of the interior stencil points when using a two-point flux approximation. It uses a graph to build stencil instead of a mesh
 */
class TwoPointFluxApproximationWithGraph : public TwoPointFluxApproximation
{
public:

  /**
   * @brief Static Factory Catalog Functions.
   * @return the catalog name
   */
  static std::string CatalogName() { return "TwoPointFluxApproximationWithGraph"; }

  TwoPointFluxApproximationWithGraph() = delete;

  /**
   * @brief Constructor.
   * @param name the name of the TwoPointFluxApproximationWithGraph in the data repository
   * @param parent the parent group of this group.
   */
  TwoPointFluxApproximationWithGraph( std::string const & name, dataRepository::Group * const parent );

  void recoverGraph() const;
  

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
     /// The key for fieldName
    static constexpr auto fieldNameString             = "fieldName";
    /// The key for coefficientName
    static constexpr auto coeffNameString             = "coefficientName";
    /// The key for targetRegions
    static constexpr auto targetRegionsString         = "targetRegions";
    /// The key for areaRelTol
    static constexpr auto areaRelativeToleranceString = "areaRelTol";

    // Keys below are for wrappers registered on MeshLevel, not the current object

    /// The key for cellStencil
    static constexpr auto cellStencilString           = "cellStencil";
    /// The key for fractureStencil
    static constexpr auto fractureStencilString       = "fractureStencil";

    static constexpr auto graphString = "graph";
     } viewKeys;
  /// @endcond

  
protected:

  virtual void computeCellStencil(MeshLevel & mesh) const override;

private:
  string m_graphString;
  GraphFromText * m_graph;
};

}


#endif //GEOSX_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATIONWITHGRAPH_HPP_
