/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "PartitionBase.hpp"

namespace geosx
{

PartitionBase::PartitionBase( const unsigned int numPartitions,
                              const unsigned int thisPartition )
  :
  m_size( numPartitions ),
  m_rank( thisPartition ),
  m_numColors( 1 )
{}

PartitionBase::~PartitionBase()
{}

}
