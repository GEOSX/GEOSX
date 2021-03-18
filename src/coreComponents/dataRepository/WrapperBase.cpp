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

/** @file */

#include "WrapperBase.hpp"

#include "Group.hpp"
#include "RestartFlags.hpp"


namespace geosx
{
namespace dataRepository
{


WrapperBase::WrapperBase( string const & name,
                          Group & parent ):
  m_name( name ),
  m_parent( &parent ),
  m_sizedFromParent( 1 ),
  m_restart_flags( RestartFlags::WRITE_AND_READ ),
  m_plotLevel( PlotLevel::NOPLOT ),
  m_inputFlag( InputFlags::INVALID ),
  m_successfulReadFromInput( false ),
  m_description(),
  m_registeringObjects(),
  m_conduitNode( parent.getConduitNode()[ name ] )
{}


WrapperBase::~WrapperBase()
{}

void WrapperBase::resize()
{
  resize( m_parent->size());
}

void WrapperBase::copyWrapperAttributes( WrapperBase const & source )
{
  m_sizedFromParent = source.m_sizedFromParent;
  m_restart_flags = source.m_restart_flags;
  m_plotLevel  = source.m_plotLevel;
  m_inputFlag = source.m_inputFlag;
  m_description = source.m_description;
}

string WrapperBase::getPath() const
{
  // In the Conduit node heirarchy everything begins with 'Problem', we should change it so that
  // the ProblemManager actually uses the root Conduit Node but that will require a full rebaseline.
  string const noProblem = m_conduitNode.path().substr( sizeof( "Problem" ) -1 );
  return noProblem.empty() ? "/" : noProblem;
}

string WrapperBase::dumpInputOptions( bool const outputHeader ) const
{
  string rval;
  char temp[1000] = {0};
  if( outputHeader )
  {
    sprintf( temp, "  |         name         |  opt/req  | Description \n" );
    rval.append( temp );
    sprintf( temp, "  |----------------------|-----------|-----------------------------------------\n" );
    rval.append( temp );
  }

  if( getInputFlag() == InputFlags::OPTIONAL ||
      getInputFlag() == InputFlags::REQUIRED )
  {
    sprintf( temp, "  | %20s | %9s | %s \n",
             getName().c_str(),
             InputFlagToString( getInputFlag()).c_str(),
             getDescription().c_str() );
    rval.append( temp );
  }

  return rval;
}

#if defined(USE_TOTALVIEW_OUTPUT)
int WrapperBase::setTotalviewDisplay() const
{
  //std::cout<<"exectuing WrapperBase::setTotalviewDisplay()"<<std::endl;
//  TV_ttf_add_row("TYPE", TV_ttf_type_ascii_string, type.c_str() );
  TV_ttf_add_row( "m_name", LvArray::system::demangle< string >().c_str(), &m_name );
  TV_ttf_add_row( "m_parent", LvArray::system::demangle< Group >().c_str(), m_parent );
  TV_ttf_add_row( "m_sizedFromParent", "int", &m_sizedFromParent );
  TV_ttf_add_row( "m_restart_flags", LvArray::system::demangle< RestartFlags >().c_str(), &m_restart_flags );
  TV_ttf_add_row( "m_plotLevel", LvArray::system::demangle< PlotLevel >().c_str(), &m_plotLevel );
  TV_ttf_add_row( "m_inputFlag", LvArray::system::demangle< InputFlags >().c_str(), &m_inputFlag );
  TV_ttf_add_row( "m_description", LvArray::system::demangle< string >().c_str(), &m_description );
  size_t junk = m_registeringObjects.size();
  TV_ttf_add_row( "m_registeringObjects",
                  totalview::format< string, size_t >( 1, &junk ).c_str(),
                  m_registeringObjects.data() );

  return 0;
}
#endif


}
} /* namespace geosx */

#if defined(USE_TOTALVIEW_OUTPUT)
/**
 * @brief Global function correlated with WrapperBase to be called by Totalview when displaying
 *        a WrapperBase as a VieWrapper<T>
 * @param wrapper A pointer to the wrapper that will be displayed.
 * @return 0
 */
int TV_ttf_display_type( const geosx::dataRepository::WrapperBase * wrapper )
{
  if( wrapper!=nullptr )
  {
    //std::cout<<"displaying WrapperBase "<<wrapper->getName()<<" as "<<wrapper->totalviewTypeName()<<std::endl;
// keep this and try to make it work later on.
//    rval = TV_ttf_add_row( "casted_this", wrapper->totalviewTypeName().c_str(), wrapper );
    wrapper->setTotalviewDisplay();
  }
  return 0;
}
#endif
