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
 * @file GraphFromText.cpp
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "common/Path.hpp"
#include "GraphFromText.hpp"
#include "managers/Functions/FunctionManager.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertice.hpp"
#include "common/Logger.hpp"

namespace geosx
{

using namespace dataRepository;


GraphFromText::GraphFromText( const std::string & name,
                              Group * const parent ):
  GraphBase( name, parent ),m_file("")
  {
     setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

    // This enables logLevel filtering
    enableLogLevelInput();

    registerWrapper( viewKeyStruct::fileString, &m_file )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the graph file." );
    registerWrapper( viewKeyStruct::meshString, &m_meshString )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Name of the mesh to link to the graph." );

    
  }


GraphFromText::~GraphFromText()
{}

void GraphFromText::GetTargetReferences()
{
  if( !m_meshString.empty())
  {
    Group * tmp = this->GetGroupByPath( m_meshString );
    m_mesh = tmp;
    GEOSX_ERROR_IF( m_mesh == nullptr, "The event " << m_meshString << " does not exist or it is not executable." );
  }
}


void GraphFromText::GenerateGraph()
{
  GetTargetReferences();
  std::ifstream infile(m_file);
  std::string line;
  std::vector<Vertice*> vertices;
  //std::vector<Edge*> edges;
  while (std::getline(infile, line))
  {
    std::istringstream iss(line);
    int a, b, c;
    if (!(iss >> a >> b >> c)) { break; } // error
    std::cout<<a<<" "<<b<<" "<<c<<"\n" ;
    int place_b=-1;
    int place_c=-1;
    long unsigned int i=0;
    while (i<vertices.size() && (place_b==-1 || place_c==-1))
    {
      if (vertices[i]->getIndice()==b)
      {
        place_b=i;
      }
      if (vertices[i]->getIndice()==c)
      {
        place_c=i;
      }
      ++i;
    }
    Vertice* v1;
    Vertice* v2;
    if(place_b==-1)
    {
      v1=new Vertice(b);
      vertices.push_back(v1);
    }
    else
    {
      v1=vertices[place_b];
    }
    if(place_c==-1)
    {
      v2=new Vertice(c);
      vertices.push_back(v2);
    }
    else
    {
      v2=vertices[place_c];
    }
    Edge* e1= new Edge(a, v1, v2);
    m_edges.push_back(e1);

  }
}

REGISTER_CATALOG_ENTRY( GraphBase, GraphFromText, std::string const &, Group * const )
} /* namespace geosx */
