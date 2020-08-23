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
 * @file Edge.cpp
 */
#include <string>
#include <iostream>
#include "common/Path.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertice.hpp"

namespace geosx
{


Edge::Edge( int const index, Vertice* neighbour1, Vertice* neighbour2):
  ind(index),n1(neighbour1),n2(neighbour2)
  {    
    std::cout<<"Constructing edge "<<ind<<" between vertices "<<n1->getIndice()<<" and "<<n2->getIndice()<<"\n";
  }


Edge::~Edge()
{}

} /* namespace geosx */
