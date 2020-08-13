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
 * @file PartitionBase.cpp
 */

#include "mpiCommunications/PartitionBase.hpp"

#include "managers/DomainPartition.hpp"
#include "mesh/NodeManager.hpp"

#include <limits.h>

namespace geosx
{
using namespace dataRepository;

PartitionBase::PartitionBase() :
  m_domain( nullptr ),
  m_hasLocalGhosts( false )
{
  //maxComm
}

PartitionBase::PartitionBase( const unsigned int numPartitions,
                              const unsigned int thisPartiton ) :
  m_size( numPartitions ),
  m_rank( thisPartiton ),
  m_color( 0 ),
  m_numColors( 1 ),
  m_domain( nullptr ),
  m_t1( 0.0 ),
  m_t2( 0.0 ),
  m_t3( 0.0 ),
  m_t4( 0.0 ),
  m_hasLocalGhosts( false )
{
  //maxComm
}

PartitionBase::~PartitionBase()
{}

void
PartitionBase::SetDomain( DomainPartition * domain )
{
  // set the const pointer "m_domain" by casting away the the const
  // on the address of the pointer, and modifying what the address of
  // the pointer.
  DomainPartition ** temp = const_cast< DomainPartition ** >( &m_domain );
  *temp = domain;
}

}  // namespace geosx
