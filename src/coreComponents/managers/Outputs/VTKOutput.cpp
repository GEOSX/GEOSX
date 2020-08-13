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
 * @file VTKOutput.cpp
 */

#include "VTKOutput.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{
using namespace dataRepository;

VTKOutput::VTKOutput( std::string const & name, Group * const parent ) :
  OutputBase( name, parent ),
  m_plotFileRoot(),
  m_writeFaceMesh(),
  m_plotLevel(),
  m_writer( name )
{
  registerWrapper( viewKeysStruct::plotFileRoot, &m_plotFileRoot )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "" );

  registerWrapper( viewKeysStruct::writeFEMFaces, &m_writeFaceMesh )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "" );

  registerWrapper( viewKeysStruct::plotLevel, &m_plotLevel )
    ->setApplyDefaultValue( 1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "" );

  registerWrapper( viewKeysStruct::binaryString, &m_writeBinaryData )
    ->setApplyDefaultValue( 1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Output the data in binary format" );
}

VTKOutput::~VTKOutput()
{}

void
VTKOutput::Execute( real64 const time_n,
                    real64 const GEOSX_UNUSED_PARAM( dt ),
                    integer const cycleNumber,
                    integer const GEOSX_UNUSED_PARAM( eventCounter ),
                    real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                    Group * domain )
{
  DomainPartition * domainPartition = Group::group_cast< DomainPartition * >( domain );
  if( m_writeBinaryData )
  {
    m_writer.SetOutputMode( vtk::VTKOutputMode::BINARY );
  }
  else
  {
    m_writer.SetOutputMode( vtk::VTKOutputMode::ASCII );
  }
  m_writer.SetPlotLevel( m_plotLevel );
  m_writer.Write( time_n, cycleNumber, *domainPartition );
}

REGISTER_CATALOG_ENTRY( OutputBase, VTKOutput, std::string const &, Group * const )
} /* namespace geosx */
