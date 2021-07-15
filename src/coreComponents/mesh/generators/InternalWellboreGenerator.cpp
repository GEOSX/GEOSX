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

/**
 * @file InternalWellboreGenerator.cpp
 */

#include "InternalWellboreGenerator.hpp"

#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"

namespace geosx
{
using namespace dataRepository;

InternalWellboreGenerator::InternalWellboreGenerator( string const & name,
                                                      Group * const parent ):
  InternalMeshGenerator( name, parent ),
  m_trajectory(),
  m_cartesianOuterBoundary(),
  m_isFullAnnulus( false ),
  m_autoSpaceRadialElems(),
  m_radialCoords( this->m_setCoords[0] )
{

  getWrapper< array1d< real64 > >( viewKeyStruct::xCoordsString() ).
    setInputFlag( InputFlags::FALSE );

  getWrapper< array1d< real64 > >( viewKeyStruct::yCoordsString() ).
    setInputFlag( InputFlags::FALSE );

  getWrapper< array1d< integer > >( viewKeyStruct::xElemsString() ).
    setInputFlag( InputFlags::FALSE );

  getWrapper< array1d< integer > >( viewKeyStruct::yElemsString() ).
    setInputFlag( InputFlags::FALSE );


  registerWrapper( viewKeyStruct::radiusString(), &( m_vertices[0] ) ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Wellbore radius" );

  registerWrapper( viewKeyStruct::thetaString(), &( m_vertices[1] ) ).
    setApplyDefaultValue( 360.0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Tangent angle defining geometry size: 90 for quarter, 180 for half and 360 for full wellbore geometry" );

  registerWrapper( viewKeyStruct::rElemsString(), &( m_nElems[0] ) ).
    setApplyDefaultValue( 10 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Number of elements in the radial direction" );

  registerWrapper( viewKeyStruct::tElemsString(), &( m_nElems[1] ) ).
    setApplyDefaultValue( 40 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Number of elements in the tangent direction" );

  // TODO to enable the use of radial bias
  registerWrapper( viewKeyStruct::rBiasString(), &( m_nElemBias[0] ) ).
    setApplyDefaultValue( -0.8 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Bias of element sizes in the radial direction" );

  registerWrapper( viewKeyStruct::trajectoryString(), &m_trajectory ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Coordinates defining the wellbore trajectory" );

  registerWrapper( viewKeyStruct::cartesianOuterBoundaryString(), &m_cartesianOuterBoundary ).
    setApplyDefaultValue( 1000000 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Enforce a Cartesian aligned outer boundary on the outer "
                    "block starting with the radial block specified in this value" );

  registerWrapper( viewKeyStruct::autoSpaceRadialElemsString(), &m_autoSpaceRadialElems ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Automatically set number and spacing of elements in the radial direction. "
                    "This overrides the values of nr" );

}

void InternalWellboreGenerator::postProcessInput()
{

  GEOSX_ERROR_IF( m_nElems[1].size() > 1,
                  "Only one block in the theta direction is currently supported. "
                  "This is specified by the nt keyword in InternalWellbore" );

  GEOSX_ERROR_IF( m_nElems[2].size() > 1,
                  "Only one block in the z direction is currently supported. "
                  "This is specified by the nz keyword in InternalWellbore" );



  GEOSX_ERROR_IF( m_trajectory.size( 0 ) != 2 || m_trajectory.size( 1 ) != 3,
                  "Input for trajectory should be specified in the form of "
                  "{ { xbottom, ybottom, zbottom }, { xtop, ytop, ztop } }." );

  // Project trajectory to bottom and top of the wellbore
  real64 trajectoryVector[3] = {0};
  for( int i=0; i<3; ++i )
  {
    trajectoryVector[i] = m_trajectory[1][i] - m_trajectory[0][i];
  }
  LvArray::tensorOps::normalize< 3 >( trajectoryVector );
  real64 const scaleb = ( m_vertices[2][0] - m_trajectory[0][2] ) / trajectoryVector[2];
  real64 const scalet = ( m_vertices[2][1] - m_trajectory[1][2] ) / trajectoryVector[2];
  for( int i=0; i<3; ++i )
  {
    m_trajectory[0][i] = m_trajectory[0][i] + scaleb * trajectoryVector[i];
    m_trajectory[1][i] = m_trajectory[1][i] + scalet * trajectoryVector[i];
  }

  arrayView1d< real64 const > const theta = m_vertices[1];
  real64 const dTheta = theta.back() - theta[0];

  // enable full annulus corrections if the mesh is 360 degrees
  if( dTheta >= 360.0 - 1e-99 )
  {
    m_isFullAnnulus = true;
  }

  // automatically set radial coordinates
  if( m_autoSpaceRadialElems.size() > 0 )
  {
    localIndex const numRadialBlocks = m_vertices[0].size() - 1;

    // loop over blocks in the radial direction (i-direction).
    for( localIndex iBlock=0; iBlock<numRadialBlocks; ++iBlock )
    {
      real64 const rInner = m_vertices[0][iBlock];
      real64 const rOuter = m_vertices[0][iBlock+1];

      // if we are on the first block, always insert the inner radius as the
      // fixed inner boundary.
      if( iBlock==0 )
      {
        m_radialCoords.emplace_back( rInner );
      }

      // Are we going to auto-size this block??
      if( m_autoSpaceRadialElems[iBlock] == 1 )
      {
        // We have to set the starting index for this block so that we skip any
        // values that are fixed from the last block when we scale the radial
        // coordinates.
        localIndex const startingIndex = m_radialCoords.size() - 1;

        // keep a count of actual number of radial elements...we will resize
        // the number of elements later.
        localIndex actualNumberOfRadialElements = 0;

        // This is the factor so that all the coordinates end up s.t. the inner
        // and outer radial coordinates are preserved.
        real64 scalingFactor = 0;

        // Loop over an excessive number of elements in the radial direction.
        // This bound needs to be more than we will end up with.
        for( localIndex i=0; i<(m_nElems[0][iBlock]+1)*100; ++i )
        {
          // approximate the theta direction size at i, and approximate the
          // radius at (i+1).
          real64 const t_i =  m_radialCoords.back() * ( 2 * M_PI * dTheta / 360 ) / m_nElems[1][0];
          real64 const r_ip1_0 = ( m_radialCoords.back() +  t_i );

          // approximate the theta direction size at the approximated radius.
          real64 const tElemSize_ip1_0 =  r_ip1_0 * ( 2 * M_PI * dTheta / 360 ) / m_nElems[1][0];

          // set the next radius a some combination of the inner and outer radius for this "element".
          constexpr real64 c = 0.5;
          real64 const r_ip1 = m_radialCoords.back() + ( 1.0 - c ) * t_i + c * tElemSize_ip1_0;

          // if the radius of the next layer is bigger than rOuter, we figure
          // out where to cut off the layer.
          if( r_ip1 > rOuter )
          {
            real64 const overshoot = r_ip1 - rOuter;
            real64 const undershoot = rOuter - m_radialCoords.back();

            if( overshoot < undershoot )
            {
              // use and append the overshot value
              m_radialCoords.emplace_back( r_ip1 );
              actualNumberOfRadialElements = i+1;
              scalingFactor =  ( rOuter - rInner ) / ( r_ip1 - rInner );
            }
            else
            {
              // use the undershot value...throw away the overshot value.
              actualNumberOfRadialElements = i;
              scalingFactor = ( rOuter - rInner ) / ( m_radialCoords.back() - rInner );
            }
            break;
          }
          else
          {
            m_radialCoords.emplace_back( r_ip1 );
          }
        }

        // set the number of actual radial elements specified by the auto
        // spacing
        m_nElems[0][iBlock] = actualNumberOfRadialElements;

        // scale the coordinates to ensure they align with the inner and outer
        // radius of the block.
        for( localIndex i=startingIndex; i<m_radialCoords.size(); ++i )
        {
          m_radialCoords[i] = ( m_radialCoords[i] - rInner ) * scalingFactor + rInner;
        }
      }
      else
      {
        // even/fixed spacing option
        real64 min = m_vertices[0][iBlock];
        real64 max = m_vertices[0][iBlock+1];
        real64 const h = (max-min) / m_nElems[0][iBlock];
        localIndex const numNodes = m_nElems[0][iBlock] + 1;
        for( localIndex i=1; i<numNodes; ++i )
        {
          m_radialCoords.emplace_back( min +  i * h );
        }
      }
    }
  }

  InternalMeshGenerator::postProcessInput();
}

void InternalWellboreGenerator::reduceNumNodesForPeriodicBoundary( SpatialPartition & partition,
                                                                   integer ( & numNodesInDir )[3] )
{
  if( m_isFullAnnulus )
  {
    if( partition.m_Partitions[1]==0 )
    {
      numNodesInDir[1] -= 1;
    }
    else if( partition.m_Partitions[1]>2 )
    {
      partition.m_Periodic[1] = 1;
    }
  }

}

void InternalWellboreGenerator::
  setNodeGlobalIndicesOnPeriodicBoundary( SpatialPartition & partition,
                                          int ( & globalIJK )[3] )
{

  GEOSX_UNUSED_VAR( partition );
  if( m_isFullAnnulus )
  {
    // TODO this only works for single theta partition.
    if( globalIJK[1] == m_nElems[1].back() )
    {
      globalIJK[1] = 0;
    }
  }
}

void InternalWellboreGenerator::setConnectivityForPeriodicBoundaries( int ( & globalIJK )[3],
                                                                      integer const ( &numNodesInDir )[3],
                                                                      int const ( &firstElemIndexInPartition )[3],
                                                                      localIndex ( & nodeOfBox )[8] )
{
  if( m_isFullAnnulus )
  {
    setConnectivityForPeriodicBoundary( 1,
                                        globalIJK,
                                        numNodesInDir,
                                        firstElemIndexInPartition,
                                        nodeOfBox );
  }
}

void InternalWellboreGenerator::coordinateTransformation( NodeManager & nodeManager )
{
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  Group & nodeSets = nodeManager.sets();
  SortedArray< localIndex > & xnegNodes = nodeSets.getReference< SortedArray< localIndex > >( string( "xneg" ) );
  SortedArray< localIndex > & xposNodes = nodeSets.getReference< SortedArray< localIndex > >( string( "xpos" ) );
  SortedArray< localIndex > & ynegNodes = nodeSets.getReference< SortedArray< localIndex > >( string( "yneg" ) );
  SortedArray< localIndex > & yposNodes = nodeSets.getReference< SortedArray< localIndex > >( string( "ypos" ) );

  SortedArray< localIndex > & rnegNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "rneg" ) ).reference();
  SortedArray< localIndex > & rposNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "rpos" ) ).reference();
  SortedArray< localIndex > & tnegNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "tneg" ) ).reference();
  SortedArray< localIndex > & tposNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "tpos" ) ).reference();

  real64 const cartesianMappingInnerRadius = m_cartesianOuterBoundary < m_vertices[0].size() ?
                                             m_vertices[0][m_cartesianOuterBoundary] :
                                             1e99;

  // Map to radial mesh
  for( localIndex a = 0; a<nodeManager.size(); ++a )
  {
    real64 meshTheta = X[a][1] * M_PI / 180.0;
    int meshAxis = static_cast< int >( round( meshTheta * 2.0 / M_PI ) );
    real64 meshPhi = fabs( meshTheta - meshAxis * M_PI / 2.0 );
    real64 meshRout = m_max[0] / cos( meshPhi );
    real64 meshRact;

    if( X[a][0] > cartesianMappingInnerRadius )
    {
      real64 const cartesianScaling = ( meshRout - cartesianMappingInnerRadius ) / ( m_max[0] - cartesianMappingInnerRadius );
      meshRact = cartesianScaling * ( X[a][0] - cartesianMappingInnerRadius ) + cartesianMappingInnerRadius;
    }
    else
    {
      meshRact = X[a][0];
    }

    // Wellbore nodesets
    if( isEqual( X[a][0], m_min[0], m_coordinatePrecision ) )
    {
      rnegNodes.insert( a );
    }

    if( isEqual( X[a][0], m_max[0], m_coordinatePrecision ) )
    {
      rposNodes.insert( a );
    }

    if( isEqual( X[a][1], m_min[1], m_coordinatePrecision ) )
    {
      tnegNodes.insert( a );
    }
    if( isEqual( X[a][1], m_max[1], m_coordinatePrecision ) )
    {
      tposNodes.insert( a );
    }

    X[a][0] = meshRact * cos( meshTheta );
    X[a][1] = meshRact * sin( meshTheta );

    // Add mapped values to nodesets
    if( m_cartesianOuterBoundary<m_vertices[0].size() )
    {
      if( isEqual( X[a][0], -1 * m_max[0], m_coordinatePrecision ) )
      {
        xnegNodes.insert( a );
      }
      if( isEqual( X[a][0], m_max[0], m_coordinatePrecision ) )
      {
        xposNodes.insert( a );
      }
      if( isEqual( X[a][1], -1 * m_max[0], m_coordinatePrecision ) )
      {
        ynegNodes.insert( a );
      }
      if( isEqual( X[a][1], m_max[0], m_coordinatePrecision ) )
      {
        yposNodes.insert( a );
      }
    }
  }

  // Map to inclined wellbore
  {
    for( int localNodeIndex=0; localNodeIndex<nodeManager.size(); ++localNodeIndex )
    {
      // Get Cartesian coordinates of a reference centered vertical wellbore
      real64 & xCoord = X( localNodeIndex, 0 );
      real64 & yCoord = X( localNodeIndex, 1 );
      real64 const & zCoord = X( localNodeIndex, 2 );

      // Compute cylindrical coordinates of a reference centered vertical wellbore
      real64 rCoord = sqrt( xCoord * xCoord + yCoord * yCoord );

      {
        real64 tCoord;

        if( rCoord < m_coordinatePrecision )
        {
          tCoord = 0.0;
        }
        else if( yCoord >= 0.0 )
        {
          tCoord = acos( xCoord/rCoord );
        }
        else
        {
          tCoord = 2.0 * M_PI - acos( xCoord/rCoord );
        }

        tCoord *= 180.0 / M_PI;

        // Radial distance of the outer square boundary of a reference centered vertical wellbore
        real64 meshTheta = tCoord * M_PI / 180.0;
        int meshAxis = static_cast< int >( round( meshTheta * 2.0 / M_PI ) );
        real64 meshPhi = fabs( meshTheta - meshAxis * M_PI / 2.0 );
        real64 meshRout = m_cartesianOuterBoundary < m_vertices[0].size() ? m_max[0] / cos( meshPhi ) : m_max[0];

        // Wellbore trajectory
        real64 xTopCenter = m_trajectory[0][0];
        real64 yTopCenter = m_trajectory[0][1];
        real64 zTop = m_min[2];

        real64 xBottomCenter = m_trajectory[1][0];
        real64 yBottomCenter = m_trajectory[1][1];
        real64 zBottom = m_max[2];

        real64 dx = xBottomCenter - xTopCenter;
        real64 dy = yBottomCenter - yTopCenter;
        real64 dz = zBottom - zTop;
        real64 dr = sqrt( dx*dx + dy*dy );
        real64 dl = sqrt( dr*dr + dz*dz );

        // Azimuth of the wellbore from x-axis
        real64 theta0;

        if( dr < m_coordinatePrecision )
        {
          theta0 = 0.0;
        }
        else if( dy>=0.0 )
        {
          theta0 = acos( dx/dr );
        }
        else
        {
          theta0 = 2.0 * M_PI - acos( dx/dr );
        }

        // The horizontal section of an inclined wellbore is an ellipse
        // The principle directions of this ellipse are defined by dTheta = 0, and PI/2
        real64 dTheta = meshTheta - theta0;
        real64 tanDTheta = tan( dTheta );

        // Transform radial coordinate regarding the elliptical shape of the wellbore section in the horizontal plane
        // This transformation ensures that the outer square boundary is unchanged
        // TODO create a function in ComputationalGeometry class for this pure geometrical transformation
        real64 transformCoeff = sqrt ( ( 1.0 + tanDTheta * tanDTheta )/( dz*dz/dl/dl + tanDTheta * tanDTheta ) );

        // rCoordRelative = 1 on the borehole, = 0 on the outer boundary
        // A power parameter round( dl/dz ) >= 1 is added
        // to control the decay of rCoordRelative to 0 at the farfield zone.
        // This consideration is userful for highly deviated wellbore problem
        // of which a large transformation is required for the near borehole zone.
        real64 rCoordRelative = pow( ( meshRout - rCoord ) / ( meshRout - m_min[0] ), round( dl/dz ) );

        // Transform the radial coordinate such that:
        // rCoordTransform = rCoord * transformCoeff on the borehole
        // rCoordTransform = rCoord on the boundary (the boundary is not transformed)
        real64 rCoordTransform = rCoord * ( rCoordRelative * ( transformCoeff - 1.0 ) + 1.0 );

        // Compute transformed cartesian coordinates
        xCoord = rCoordTransform * cos( meshTheta );
        yCoord = rCoordTransform * sin( meshTheta );

        // Moving the coordinate in the horizontal plane with respect to the center of the wellbore section
        // This transformation ensures that the ourter square boundary is unchanged
        real64 zRatio = ( zCoord - zTop ) / ( zBottom -zTop );
        real64 xCenter = xTopCenter + (xBottomCenter - xTopCenter) * zRatio;
        real64 yCenter = yTopCenter + (yBottomCenter - yTopCenter) * zRatio;

        xCoord += xCenter * ( meshRout - rCoord ) / ( meshRout - m_min[0] * transformCoeff );
        yCoord += yCenter * ( meshRout - rCoord ) / ( meshRout - m_min[0] * transformCoeff );
      }
    }
  }
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, InternalWellboreGenerator, string const &, Group * const )
} /* namespace geosx */
