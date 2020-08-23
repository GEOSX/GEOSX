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
 * @file Vertice.cpp
 */
#include <string>
#include <iostream>
#include "common/Path.hpp"
#include "Vertice.hpp"

namespace geosx
{



Vertice::Vertice( const int index):
  ind(index)
  {    
    std::cout<<"Constructing vertice "<<ind<<"\n";
  }


Vertice::~Vertice()
{}

} /* namespace geosx */
