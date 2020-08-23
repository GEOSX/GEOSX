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
 * @file Edge.hpp
 */

#ifndef GEOSX_MESH_EDGE_HPP_
#define GEOSX_MESH_EDGE_HPP_
#include "mesh/Vertice.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"


namespace geosx
{

/**
 * @class Edge
 *
 * Represents connection between two points for the GraphBase.
 */
class Edge
{
public:
  /**
   * @brief Constructor for Edge object
   * @param [in] index of the edge
   * @param [in] neighbour1 one of the two Vertices forming the connection
   * @param [in] neighbour2 the other Vertice
   * 
   */   
  Edge( int const index, Vertice* neighbour1, Vertice* neighbour2);

  /// Destructor
  virtual ~Edge();

  localIndex getIndice() const { return ind; }

  Vertice* getN1() const { return n1; }

  Vertice* getN2() const { return n2; }

  
private:
  localIndex ind;
  Vertice* n1;
  Vertice* n2;
};

} /* namespace geosx */

#endif /* GEOSX_MESH_EDGE_HPP_ */
