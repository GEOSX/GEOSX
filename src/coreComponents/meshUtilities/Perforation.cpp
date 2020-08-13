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
 * @file Perforation.cpp
 */

#include "Perforation.hpp"
#include "dataRepository/InputFlags.hpp"

namespace geosx
{
using namespace dataRepository;

Perforation::Perforation( string const & name, Group * const parent ) :
  Group( name, parent ),
  m_distanceFromHead( 0 ),
  m_wellTransmissibility( 0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::distanceFromHeadString, &m_distanceFromHead )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Linear distance from well head to the perforation" );

  registerWrapper( viewKeyStruct::wellTransmissibilityString,
                   &m_wellTransmissibility )
    ->setApplyDefaultValue( -1.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Perforation transmissibility" );
}

void
Perforation::PostProcessInput()
{
  GEOSX_ERROR_IF( m_distanceFromHead <= 0,
                  "Invalid distance well head to perforation " << getName() );
}

Perforation::~Perforation()
{}

}  //namespace geosx
