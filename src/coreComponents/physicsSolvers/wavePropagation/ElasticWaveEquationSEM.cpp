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
 * @file ElasticWaveEquationSEM.cpp
 */

#include "ElasticWaveEquationSEM.hpp"

#include "dataRepository/KeyNames.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geosx
{

using namespace dataRepository;

ElasticWaveEquationSEM::ElasticWaveEquationSEM( const std::string & name,
                                                  Group * const parent ):
  SolverBase( name,
              parent )
{

  registerWrapper( viewKeyStruct::sourceCoordinatesString(), &m_sourceCoordinates ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Coordinates (x,y,z) of the sources" );

  registerWrapper( viewKeyStruct::sourceNodeIdsString(), &m_sourceNodeIds ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Indices of the nodes (in the right order) for each source point" );

  registerWrapper( viewKeyStruct::sourceConstantsString(), &m_sourceConstants_x ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Constant part of the source for the nodes listed in m_sourceNodeIds" );

  registerWrapper( viewKeyStruct::sourceConstantsString(), &m_sourceConstants_y ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Constant part of the source for the nodes listed in m_sourceNodeIds" );

  registerWrapper( viewKeyStruct::sourceConstantsString(), &m_sourceConstants_z ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Constant part of the source for the nodes listed in m_sourceNodeIds" );


  registerWrapper( viewKeyStruct::sourceIsLocalString(), &m_sourceIsLocal ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Flag that indicates whether the source is local to this MPI rank" );


  registerWrapper( viewKeyStruct::timeSourceFrequencyString(), &m_timeSourceFrequency ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Central frequency for the time source" );


  registerWrapper( viewKeyStruct::receiverCoordinatesString(), &m_receiverCoordinates ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Coordinates (x,y,z) of the receivers" );

  registerWrapper( viewKeyStruct::receiverNodeIdsString(), &m_receiverNodeIds ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Indices of the nodes (in the right order) for each receiver point" );

  registerWrapper( viewKeyStruct::sourceConstantsString(), &m_sourceConstants ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Constant part of the receiver for the nodes listed in m_receiverNodeIds" );

  registerWrapper( viewKeyStruct::receiverIsLocalString(), &m_receiverIsLocal ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Flag that indicates whether the receiver is local to this MPI rank" );

    registerWrapper( viewKeyStruct::rickerOrderString(), &m_rickerOrder ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 2 ).
    setDescription( "Flag that indicates the order of the Ricker to be used o, 1 or 2. Order 2 by default" );

     registerWrapper( viewKeyStruct::outputSismoTraceString(), &m_outputSismoTrace ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag that indicates if we write the sismo trace in a file .txt, 0 no output, 1 otherwise" );

  registerWrapper( viewKeyStruct::displacementNp1AtReceiversString(), &m_displacementNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep" );


}

ElasticWaveEquationSEM::~ElasticWaveEquationSEM()
{
  // TODO Auto-generated destructor stub
}

void ElasticWaveEquationSEM::postProcessInput()
{

  GEOSX_ERROR_IF( m_sourceCoordinates.size( 1 ) != 3,
                  "Invalid number of physical coordinates for the sources" );

  GEOSX_ERROR_IF( m_receiverCoordinates.size( 1 ) != 3,
                  "Invalid number of physical coordinates for the receivers" );

  localIndex const numNodesPerElem = 8;

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );

  m_sourceNodeIds.resizeDimension< 0, 1 >( numSourcesGlobal, numNodesPerElem );
  m_sourceConstants_x.resizeDimension< 0, 1 >( numSourcesGlobal, numNodesPerElem );
  m_sourceConstants_y.resizeDimension< 0, 1 >( numSourcesGlobal, numNodesPerElem );
  m_sourceConstants_z.resizeDimension< 0, 1 >( numSourcesGlobal, numNodesPerElem );
  m_sourceIsLocal.resizeDimension< 0 >( numSourcesGlobal );

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_receiverNodeIds.resizeDimension< 0, 1 >( numReceiversGlobal, numNodesPerElem );
  m_receiverConstants.resizeDimension< 0, 1 >( numReceiversGlobal, numNodesPerElem );
  m_receiverIsLocal.resizeDimension< 0 >( numReceiversGlobal );

  m_displacementNp1AtReceivers.resizeDimension< 0, 1 >( numReceiversGlobal, 3 );

}

void ElasticWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{

  meshBodies.forSubGroups< MeshBody >( [&]( MeshBody & meshBody )
  {

    MeshLevel & meshLevel =  meshBody.getMeshLevel( 0 );

    NodeManager & nodeManager = meshLevel.getNodeManager();

    nodeManager.registerExtrinsicData< extrinsicMeshData::Displacementx_nm1,
                                       extrinsicMeshData::Displacementy_nm1,
                                       extrinsicMeshData::Displacementz_nm1,
                                       extrinsicMeshData::Displacementx_n,
                                       extrinsicMeshData::Displacementy_n,
                                       extrinsicMeshData::Displacementz_n,
                                       extrinsicMeshData::Displacementx_np1,
                                       extrinsicMeshData::Displacementy_np1,
                                       extrinsicMeshData::Displacementz_np1,
                                       extrinsicMeshData::ForcingRHS_x,
                                       extrinsicMeshData::ForcingRHS_y,
                                       extrinsicMeshData::ForcingRHS_z,
                                       extrinsicMeshData::MassVector,
                                       extrinsicMeshData::DampingVector_x,
                                       extrinsicMeshData::DampingVector_y,
                                       extrinsicMeshData::DampingVector_z,
                                       extrinsicMeshData::StiffnessVector_x,
                                       extrinsicMeshData::StiffnessVector_y,
                                       extrinsicMeshData::StiffnessVector_z,
                                       extrinsicMeshData::FreeSurfaceNodeIndicator >( this->getName() );

    FaceManager & faceManager = meshLevel.getFaceManager();
    faceManager.registerExtrinsicData< extrinsicMeshData::FreeSurfaceFaceIndicator >( this->getName() );

    ElementRegionManager & elemManager = meshLevel.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerExtrinsicData< extrinsicMeshData::MediumVelocityVp >( this->getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::MediumVelocityVs >( this->getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::MediumDensity >( this->getName() );

      subRegion.registerExtrinsicData< extrinsicMeshData::LameCoefficientLambda >(this->getName());
      subRegion.registerExtrinsicData< extrinsicMeshData::LameCoefficientMu >(this->getName());
    } );

  } );
}


void ElasticWaveEquationSEM::precomputeSourceAndReceiverTerm( MeshLevel & mesh )
{
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

  arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();
  arrayView2d< localIndex > const sourceNodeIds = m_sourceNodeIds.toView();
  arrayView2d< real64 > const sourceConstants_x = m_sourceConstants_x.toView();
  arrayView2d< real64 > const sourceConstants_y = m_sourceConstants_y.toView();
  arrayView2d< real64 > const sourceConstants_z = m_sourceConstants_z.toView();
  arrayView1d< localIndex > const sourceIsLocal = m_sourceIsLocal.toView();
  sourceNodeIds.setValues< serialPolicy >( -1 );
  sourceConstants_x.setValues< serialPolicy >( -1 );
  sourceConstants_y.setValues< serialPolicy >( -1 );
  sourceConstants_z.setValues< serialPolicy >( -1 );
  sourceIsLocal.setValues< serialPolicy >( 0 );

  arrayView2d< real64 const > const receiverCoordinates = m_receiverCoordinates.toViewConst();
  arrayView2d< localIndex > const receiverNodeIds = m_receiverNodeIds.toView();
  arrayView2d< real64 > const receiverConstants = m_receiverConstants.toView();
  arrayView1d< localIndex > const receiverIsLocal = m_receiverIsLocal.toView();
  receiverNodeIds.setValues< serialPolicy >( -1 );
  receiverConstants.setValues< serialPolicy >( -1 );
  receiverIsLocal.setValues< serialPolicy >( 0 );

  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const,
                                       ElementRegionBase & elemRegion )
  {
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const,
                                                                       CellElementSubRegion & elementSubRegion )
    {

      GEOSX_THROW_IF( elementSubRegion.getElementTypeString() != "C3D8",
                      "Invalid type of element, the elastic solver is designed for hexahedral meshes only (C3D8) ",
                      InputError );

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
        localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();
        array1d< array1d< localIndex > > faceNodes( numFacesPerElem );

        forAll< serialPolicy >( elementSubRegion.size(), [=, &elementSubRegion] ( localIndex const k )
        {

          for( localIndex kf = 0; kf < numFacesPerElem; ++kf )
          {
            elementSubRegion.getFaceNodes( k, kf, faceNodes[kf] );
          }

          /// loop over all the source that haven't been found yet
          forAll< serialPolicy >( sourceCoordinates.size( 0 ), [=] ( localIndex const isrc )
          {
            if( sourceIsLocal[isrc] == 0 )
            {
              real64 const coords[3] = { sourceCoordinates[isrc][0],
                                         sourceCoordinates[isrc][1],
                                         sourceCoordinates[isrc][2] };

              real64 xLocal[numNodesPerElem][3];

              for( localIndex a=0; a< numNodesPerElem; ++a )
              {
                for( localIndex i=0; i<3; ++i )
                {
                  xLocal[a][i] = X( elemsToNodes( k, a ), i );
                }
              }

              real64 coordsOnRefElem[3]{};
              bool const sourceFound = computeCoordinatesOnReferenceElement< FE_TYPE >( coords, coordsOnRefElem, k, faceNodes, elemsToNodes, X );
              if( sourceFound )
              {
                sourceIsLocal[isrc] = 1;

                for( localIndex c=0; c<2; ++c )
                {
                  for( localIndex b=0; b<2; ++b )
                 {
                   for( localIndex a=0; a<2; ++a )
                   {
                     real64 const Grad[3] = { finiteElement::LagrangeBasis1::gradient( a, coordsOnRefElem[0])*
                                    finiteElement::LagrangeBasis1::value( b, coordsOnRefElem[1])*
                                    finiteElement::LagrangeBasis1::value( c, coordsOnRefElem[2]),
                                    finiteElement::LagrangeBasis1::value( a, coordsOnRefElem[0])*
                                    finiteElement::LagrangeBasis1::gradient( b, coordsOnRefElem[1])*
                                    finiteElement::LagrangeBasis1::value( c, coordsOnRefElem[2]),
                                    finiteElement::LagrangeBasis1::value( a, coordsOnRefElem[0])*
                                    finiteElement::LagrangeBasis1::value( b, coordsOnRefElem[1])*
                                    finiteElement::LagrangeBasis1::gradient( c, coordsOnRefElem[2])};

                     localIndex const nodeIndex = finiteElement::LagrangeBasis1::TensorProduct3D::linearIndex( a, b, c );

                     // real64 const detJ = finiteElement.template getGradN< FE_TYPE >( k, nodeIndex, xLocal, gradN );

                     real64 invJ[3][3]={{0}};
                     FE_TYPE::invJacobianTransformation( nodeIndex, xLocal, invJ );
                     sourceNodeIds[isrc][nodeIndex] = elemsToNodes[k][nodeIndex];
                     sourceConstants_x[isrc][nodeIndex] = Grad[0] * invJ[0][0] + Grad[1] * invJ[0][1] + Grad[2] * invJ[0][2];
                     sourceConstants_y[isrc][nodeIndex] = Grad[0] * invJ[1][0] + Grad[1] * invJ[1][1] + Grad[2] * invJ[1][2];
                     sourceConstants_z[isrc][nodeIndex] = Grad[0] * invJ[2][0] + Grad[1] * invJ[2][1] + Grad[2] * invJ[2][2];

                    }
                  }
                }
              }
            }
          } ); // End loop over all source


          /// loop over all the receiver that haven't been found yet
          //
          forAll< serialPolicy >( receiverCoordinates.size( 0 ), [=] ( localIndex const ircv )
          {
            if( receiverIsLocal[ircv] == 0 )
            {
              real64 const coords[3] = { receiverCoordinates[ircv][0],
                                         receiverCoordinates[ircv][1],
                                         receiverCoordinates[ircv][2] };

              real64 coordsOnRefElem[3]{};
              bool const receiverFound = computeCoordinatesOnReferenceElement< FE_TYPE >( coords, coordsOnRefElem, k, faceNodes, elemsToNodes, X );
              if( receiverFound )
              {
                receiverIsLocal[ircv] = 1;

                real64 Ntest[8];
                finiteElement::LagrangeBasis1::TensorProduct3D::value( coordsOnRefElem, Ntest );

                for( localIndex a=0; a< numNodesPerElem; ++a )
                {
                  receiverNodeIds[ircv][a] = elemsToNodes[k][a];
                  receiverConstants[ircv][a] = Ntest[a];
                }
              }
            }
          } ); // End loop over reciever

        } );// End loop over elements
      } );
    } );
  } );
}


void ElasticWaveEquationSEM::addSourceToRightHandSide( real64 const & time_n, arrayView1d< real64 > const rhs_x, arrayView1d< real64 > const rhs_y, arrayView1d< real64 > const rhs_z )
{
  arrayView2d< localIndex const > const sourceNodeIds = m_sourceNodeIds.toViewConst();
  arrayView2d< real64 const > const sourceConstants_x   = m_sourceConstants_x.toViewConst();
  arrayView2d< real64 const > const sourceConstants_y   = m_sourceConstants_y.toViewConst();
  arrayView2d< real64 const > const sourceConstants_z   = m_sourceConstants_z.toViewConst();
  arrayView1d< localIndex const > const sourceIsLocal = m_sourceIsLocal.toViewConst();

  real64 const fi = evaluateRicker( time_n, this->m_timeSourceFrequency, this->m_rickerOrder );


  forAll< serialPolicy >( m_sourceConstants_x.size( 0 ), [=] ( localIndex const isrc )
  {
    if( sourceIsLocal[isrc] == 1 )
    {
      for( localIndex inode = 0; inode < m_sourceConstants_x.size( 1 ); ++inode )
      {
        rhs_x[sourceNodeIds[isrc][inode]] = sourceConstants_x[isrc][inode] * fi;
        rhs_y[sourceNodeIds[isrc][inode]] = sourceConstants_y[isrc][inode] * fi;
        rhs_z[sourceNodeIds[isrc][inode]] = sourceConstants_z[isrc][inode] * fi;
      }
    }
  } );
}


void ElasticWaveEquationSEM::computeSismoTrace( localIndex const isismo, arrayView1d< real64 > const ux_np1, arrayView1d< real64 > const uy_np1,arrayView1d< real64 > const uz_np1 )
{
  arrayView2d< localIndex const > const receiverNodeIds = m_receiverNodeIds.toViewConst();
  arrayView2d< real64 const > const receiverConstants   = m_receiverConstants.toViewConst();
  arrayView1d< localIndex const > const receiverIsLocal = m_receiverIsLocal.toViewConst();

  arrayView2d< real64 > const u_rcvs   = m_displacementNp1AtReceivers.toView();

  // for( localIndex ircv = 0; ircv < receiverConstants.size( 0 ); ++ircv )
  forAll< serialPolicy >( receiverConstants.size( 0 ), [=] ( localIndex const ircv )
  {
    if( receiverIsLocal[ircv] == 1 )
    {
      u_rcvs[ircv][0] = 0.0;
      u_rcvs[ircv][1] = 0.0;
      u_rcvs[ircv][2] = 0.0;

      for( localIndex inode = 0; inode < receiverConstants.size( 1 ); ++inode )
      {
        u_rcvs[ircv][0] += ux_np1[receiverNodeIds[ircv][inode]]*receiverConstants[ircv][inode];
        u_rcvs[ircv][1] += uy_np1[receiverNodeIds[ircv][inode]]*receiverConstants[ircv][inode];
        u_rcvs[ircv][2] += uz_np1[receiverNodeIds[ircv][inode]]*receiverConstants[ircv][inode];
      }
    }
  } );

  forAll< serialPolicy >( receiverConstants.size( 0 ), [=] ( localIndex const ircv )
  {
    if( this->m_outputSismoTrace == 1 )
    {
      if( receiverIsLocal[ircv] == 1 )
      {
        // Note: this "manual" output to file is temporary
        //       It should be removed as soon as we can use TimeHistory to output data not registered on the mesh
        // TODO: remove the (sprintf+saveSismo) and replace with TimeHistory
        char filename[50];
        sprintf( filename, "sismoTraceUxReceiver%0ld.txt", ircv );
        this->saveSismo( isismo, u_rcvs[ircv][0], filename );
        sprintf( filename, "sismoTraceUyReceiver%0ld.txt", ircv );
        this->saveSismo( isismo, u_rcvs[ircv][1], filename );
        sprintf( filename, "sismoTraceUzReceiver%0ld.txt", ircv );
        this->saveSismo( isismo, u_rcvs[ircv][2], filename );

      }
    }
  } );
}


void ElasticWaveEquationSEM::saveSismo( localIndex isismo, real64 val_displacement, char *filename )
{
  std::ofstream f( filename, std::ios::app );
  f<< isismo << " " << val_displacement << std::endl;
  f.close();
}


void ElasticWaveEquationSEM::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const * const
  feDiscretization = feDiscretizationManager.getGroupPointer< FiniteElementDiscretization >( m_discretizationName );
  GEOSX_THROW_IF( feDiscretization == nullptr,
                  getName() << ": FE discretization not found: " << m_discretizationName,
                  InputError );
}


void ElasticWaveEquationSEM::initializePostInitialConditionsPreSubGroups()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  real64 const time = 0.0;
  applyFreeSurfaceBC( time, domain );
  applyABC( time, domain );
  precomputeSourceAndReceiverTerm( mesh );

  NodeManager & nodeManager = mesh.getNodeManager();
  FaceManager & faceManager = mesh.getFaceManager();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

  /// Get table containing all the face normals

  arrayView1d< real64 > const mass = nodeManager.getExtrinsicData< extrinsicMeshData::MassVector >();

  /// get array of indicators: 1 if face is on the free surface; 0 otherwise
  arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceFaceIndicator >();


  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const,
                                       ElementRegionBase & elemRegion )
  {
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const,
                                                                       CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      arrayView1d< real64 > const rho = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumDensity>();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        /// Loop over elements
        forAll< serialPolicy >( elemsToNodes.size( 0 ), [=] ( localIndex const k )
        {

          constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
          constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

          real64 xLocal[numNodesPerElem][3];
          for( localIndex a=0; a< numNodesPerElem; ++a )
          {
            for( localIndex i=0; i<3; ++i )
            {
              xLocal[a][i] = X( elemsToNodes( k, a ), i );
            }
          }

          real64 N[numNodesPerElem];
          real64 gradN[ numNodesPerElem ][ 3 ];


          for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
          {
            FE_TYPE::calcN( q, N );
            real64 const detJ = finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

            for( localIndex a=0; a< numNodesPerElem; ++a )
            {
              mass[elemsToNodes[k][a]] += rho[k] * detJ * N[a];
            }
          }

        } ); // end loop over element
      } );
    } );
  } );

}


void ElasticWaveEquationSEM::applyFreeSurfaceBC( real64 const time, DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getNodeManager();

  arrayView1d< real64 > const ux_nm1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementx_nm1 >();
  arrayView1d< real64 > const uy_nm1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementy_nm1 >();
  arrayView1d< real64 > const uz_nm1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementz_nm1 >();
  arrayView1d< real64 > const ux_n = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementx_n >();
  arrayView1d< real64 > const uy_n = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementy_n >();
  arrayView1d< real64 > const uz_n = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementz_n >();
  arrayView1d< real64 > const ux_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementx_np1 >();
  arrayView1d< real64 > const uy_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementy_np1 >();
  arrayView1d< real64 > const uz_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementz_np1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// set array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceFaceIndicator >();

  /// set array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceNodeIndicator >();

  fsManager.apply( time,
                   domain,
                   "faceManager",
                   string( "FreeSurface" ),
                   [&]( FieldSpecificationBase const & bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group &,
                        string const & )
  {
    string const & functionName = bc.getFunctionName();

    if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
    {
      real64 const value = bc.getScale();

      freeSurfaceFaceIndicator.setValues< serialPolicy >( 0 );
      freeSurfaceNodeIndicator.setValues< serialPolicy >( 0 );

      for( localIndex i = 0; i < targetSet.size(); ++i )
      {
        localIndex const kf = targetSet[ i ];
        freeSurfaceFaceIndicator[kf] = 1;

        localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
        for( localIndex a=0; a < numNodes; ++a )
        {
          localIndex const dof = faceToNodeMap( kf, a );
          freeSurfaceNodeIndicator[dof] = 1;

          ux_np1[dof] = value;
          uy_np1[dof] = value;
          uz_np1[dof] = value;
          ux_n[dof]   = value;
          uy_n[dof]   = value;
          uz_n[dof]   = value;
          ux_nm1[dof] = value;
          uy_nm1[dof] = value;
          uz_nm1[dof] = value;
        }
      }
    }
    else
    {
      GEOSX_ERROR( "This option is not supported yet" );
    }
  } );
}

void ElasticWaveEquationSEM::applyABC( real64 const time, DomainPartition & domain )
{

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NodeManager & nodeManager = mesh.getNodeManager();
  FaceManager & faceManager = mesh.getFaceManager();

  /// get the array of indicators: 1 if the face is on the boundary; 0 otherwise
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

  /// Get table containing all the face normals
  arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();

  FaceManager::ElemMapType const & faceToElem = faceManager.toElementRelation();
  arrayView2d< localIndex const > const & faceToElemIndex = faceToElem.m_toElementIndex;


  ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

  /// damping matrix to be computed for each dof in the boundary of the mesh
  arrayView1d< real64 > const damping_x = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector_x >();
  arrayView1d< real64 > const damping_y = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector_y >();
  arrayView1d< real64 > const damping_z = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector_z >();

  damping_x.setValues< serialPolicy >( 0.0 );
  damping_y.setValues< serialPolicy >( 0.0 );
  damping_z.setValues< serialPolicy >( 0.0 );

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  fsManager.apply( time,
                   domain,
                   "faceManager",
                   string( "ABC" ),
                   [&]( FieldSpecificationBase const & bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group &,
                        string const & )
  {
    string const & functionName = bc.getFunctionName();

    if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
    {

      forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const,
                                       ElementRegionBase & elemRegion )
      {

        elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const,
                                                                       CellElementSubRegion & elementSubRegion )
        {

          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

          arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();

          arrayView1d< real64 > const rho = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumDensity>();
          arrayView1d< real64 > const vp = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocityVp >();
          arrayView1d< real64 > const vs = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocityVs >();


          finiteElement::FiniteElementBase const &
          fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

          finiteElement::dispatch3D( fe,
                                     [&]
                                     ( auto const finiteElement )
          {
            using FE_TYPE = TYPEOFREF( finiteElement );

            //localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();
            localIndex const numNodesPerFace = 4;

            constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
            constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;


            for( localIndex i = 0; i < targetSet.size(); ++i )
            {

              localIndex const kf = targetSet[ i ];

              localIndex const k = faceToElemIndex[kf][0];

              real64 xLocal[numNodesPerElem][3];
              for( localIndex a=0; a< numNodesPerElem; ++a )
              {
                for( localIndex b=0; b<3; ++b )
                {
                  xLocal[a][b] = X( elemsToNodes( k, a ), b );
                }
              }

              real64 N[numNodesPerElem];
              real64 gradN[ numNodesPerElem ][ 3 ];

              for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
              {
                FE_TYPE::calcN( q, N );
                real64 const detJ = finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

                real64 invJ[3][3]={{0}};
                FE_TYPE::invJacobianTransformation( q, xLocal, invJ );

                for( localIndex a=0; a < numNodesPerFace; ++a )
                {
                  /// compute ds=||detJ*invJ*normalFace_{kfe}||
                  real64 tmp[3]={0};
                  real64 ds = 0.0;
                  for( localIndex b=0; b<3; ++b )
                  {
                    for( localIndex j = 0; j < 3; ++j )
                    {
                      tmp[b] += invJ[j][b]*faceNormal[kf][j];
                    }
                    ds +=tmp[b]*tmp[b];
                  }
                  ds = std::sqrt( ds );

                  localIndex numNodeGl = facesToNodes[kf][a];


                  // Damping in x=xpos direction
                  if ( faceNormal[kf][0] > 0.0 )
                  {
                    real64 const alpha_x = rho[k] * vp[k];
                    real64 const alpha_y = rho[k] * vs[k];
                    real64 const alpha_z = rho[k] * vs[k];
                    damping_x[numNodeGl] += alpha_x*detJ*ds*N[a];
                    damping_y[numNodeGl] += alpha_y*detJ*ds*N[a];
                    damping_z[numNodeGl] += alpha_z*detJ*ds*N[a];

                  }
                  // Damping in x=xneg direction
                  if ( faceNormal[kf][0] < 0.0 )
                  {
                    real64 const alpha_x = rho[k] * vp[k];
                    real64 const alpha_y = rho[k] * vs[k];
                    real64 const alpha_z = rho[k] * vs[k];
                    damping_x[numNodeGl] += alpha_x*detJ*ds*N[a];
                    damping_y[numNodeGl] += alpha_y*detJ*ds*N[a];
                    damping_z[numNodeGl] += alpha_z*detJ*ds*N[a];

                  }
                  //Damping in y=ypos direction
                  if ( faceNormal[kf][1] > 0.0 )
                  {
                    real64 const alpha_x = rho[k] * vs[k];
                    real64 const alpha_y = rho[k] * vp[k];
                    real64 const alpha_z = rho[k] * vs[k];
                    damping_x[numNodeGl] += alpha_x*detJ*ds*N[a];
                    damping_y[numNodeGl] += alpha_y*detJ*ds*N[a];
                    damping_z[numNodeGl] += alpha_z*detJ*ds*N[a];

                  }
                  //Damping in y=yneg direction
                  if ( faceNormal[kf][1] < 0.0 )
                  {
                    real64 const alpha_x = rho[k] * vs[k];
                    real64 const alpha_y = rho[k] * vp[k];
                    real64 const alpha_z = rho[k] * vs[k];
                    damping_x[numNodeGl] += alpha_x*detJ*ds*N[a];
                    damping_y[numNodeGl] += alpha_y*detJ*ds*N[a];
                    damping_z[numNodeGl] += alpha_z*detJ*ds*N[a];

                  }
                  //Damping in z=zpos direction
                  if ( faceNormal[kf][2] > 0.0 )
                  {
                    real64 const alpha_x = rho[k] * vs[k];
                    real64 const alpha_y = rho[k] * vs[k];
                    real64 const alpha_z = rho[k] * vp[k];
                    damping_x[numNodeGl] += alpha_x*detJ*ds*N[a];
                    damping_y[numNodeGl] += alpha_y*detJ*ds*N[a];
                    damping_z[numNodeGl] += alpha_z*detJ*ds*N[a];

                  }
                  //Damping in z=zneg direction
                  if ( faceNormal[kf][2] < 0.0 )
                  {
                    real64 const alpha_x = rho[k] * vs[k];
                    real64 const alpha_y = rho[k] * vs[k];
                    real64 const alpha_z = rho[k] * vp[k];
                    damping_x[numNodeGl] += alpha_x*detJ*ds*N[a];
                    damping_y[numNodeGl] += alpha_y*detJ*ds*N[a];
                    damping_z[numNodeGl] += alpha_z*detJ*ds*N[a];

                  }
                }
              }
            }
          } );
        } );
      } );
    }
  } );
}


real64 ElasticWaveEquationSEM::solverStep( real64 const & time_n,
                                            real64 const & dt,
                                            integer const cycleNumber,
                                            DomainPartition & domain )
{
  return explicitStep( time_n, dt, cycleNumber, domain );
}


real64 ElasticWaveEquationSEM::evaluateRicker( real64 const & time_n, real64 const & f0, localIndex order )
{
  real64 const o_tpeak = 1.0/f0;
  real64 pulse = 0.0;
  if((time_n <= -0.9*o_tpeak) || (time_n >= 2.9*o_tpeak))
    return pulse;

  constexpr real64 pi = M_PI;
  real64 const lam = (f0*pi)*(f0*pi);

  switch( order )
  {
    case 2:
    {
      pulse = 2.0*lam*(2.0*lam*(time_n-o_tpeak)*(time_n-o_tpeak)-1.0)*exp( -lam*(time_n-o_tpeak)*(time_n-o_tpeak));
    }
    break;
    case 1:
    {
      pulse = -2.0*lam*(time_n-o_tpeak)*exp( -lam*(time_n-o_tpeak)*(time_n-o_tpeak));
    }
    break;
    case 0:
    {
      pulse = -(time_n-o_tpeak)*exp( -lam*(time_n-o_tpeak)*(time_n-o_tpeak) );
    }
    break;
    default:
      GEOSX_ERROR( "This option is not supported yet, rickerOrder must be 0, 1 or 2" );
  }

  return pulse;
}




real64 ElasticWaveEquationSEM::explicitStep( real64 const & time_n,
                                              real64 const & dt,
                                              integer const cycleNumber,
                                              DomainPartition & domain )
{

  GEOSX_UNUSED_VAR( time_n, dt, cycleNumber );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  if( time_n < 0.5*dt )
  {
    applyFreeSurfaceBC( time_n, domain );
    postProcessInput();
    precomputeSourceAndReceiverTerm( mesh );
  }

  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real64 const > const mass = nodeManager.getExtrinsicData< extrinsicMeshData::MassVector >();
  arrayView1d< real64 const > const damping_x = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector_x >();
  arrayView1d< real64 const > const damping_y = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector_y >();
  arrayView1d< real64 const > const damping_z = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector_z >();

  arrayView1d< real64 > const ux_nm1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementx_nm1 >();
  arrayView1d< real64 > const uy_nm1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementy_nm1 >();
  arrayView1d< real64 > const uz_nm1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementz_nm1 >();
  arrayView1d< real64 > const ux_n = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementx_n >();
  arrayView1d< real64 > const uy_n = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementy_n >();
  arrayView1d< real64 > const uz_n = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementz_n >();
  arrayView1d< real64 > const ux_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementx_np1 >();
  arrayView1d< real64 > const uy_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementy_np1 >();
  arrayView1d< real64 > const uz_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementz_np1 >();

  /// get array of indicators: 1 if node on free surface; 0 otherwise
  arrayView1d< localIndex const > const freeSurfaceNodeIndicator = nodeManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceNodeIndicator >();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

  arrayView1d< real64 > const stiffnessVector_x = nodeManager.getExtrinsicData< extrinsicMeshData::StiffnessVector_x >();
  arrayView1d< real64 > const stiffnessVector_y = nodeManager.getExtrinsicData< extrinsicMeshData::StiffnessVector_y >();
  arrayView1d< real64 > const stiffnessVector_z = nodeManager.getExtrinsicData< extrinsicMeshData::StiffnessVector_z >();

  arrayView1d< real64 > const rhs_x = nodeManager.getExtrinsicData< extrinsicMeshData::ForcingRHS_x >();
  arrayView1d< real64 > const rhs_y = nodeManager.getExtrinsicData< extrinsicMeshData::ForcingRHS_y >();
  arrayView1d< real64 > const rhs_z = nodeManager.getExtrinsicData< extrinsicMeshData::ForcingRHS_z >();

  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const,
                                       ElementRegionBase & elemRegion )
  {
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const,
                                                                       CellElementSubRegion & elementSubRegion )
    {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      arrayView1d< real64 > const vp = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocityVp >();
      arrayView1d< real64 > const vs = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocityVs >();
      arrayView1d< real64 > const rho = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumDensity >();

      arrayView1d< real64 > const lambda = elementSubRegion.getExtrinsicData< extrinsicMeshData::LameCoefficientLambda >();
      arrayView1d< real64 > const mu = elementSubRegion.getExtrinsicData< extrinsicMeshData::LameCoefficientMu >();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        forAll< serialPolicy >( elemsToNodes.size( 0 ), [=] ( localIndex const k )
        {

          constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
          constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

          // Computation of the Lamé coefficients lambda and mu
          mu[k] = rho[k] * vs[k] * vs[k];
          lambda[k] = rho[k] * vp[k] * vp[k] - 2.0*mu[k];

          real64 xLocal[numNodesPerElem][3] = {{0.0}};

          for( localIndex i=0; i<numNodesPerElem; ++i )
          {
            for( localIndex a=0; a<3; ++a )
            {
              xLocal[i][a] = X[elemsToNodes[k][i]][a];
            }

          }

          real64 gradN[ numNodesPerElem ][ 3 ];

          for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
          {
            real64 const detJ = finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

	          for( localIndex i=0; i<numNodesPerElem; ++i )
	          {
		          for(localIndex j=0; j<numNodesPerElem; ++j)
		          {

                real64 const Rxx_ij = detJ* ((lambda[k]+2.0*mu[k])*gradN[j][0]*gradN[i][0] + mu[k] * gradN[j][1]*gradN[i][1] + mu[k] * gradN[j][2]*gradN[i][2]);
                real64 const Ryy_ij = detJ* ((lambda[k]+2.0*mu[k])*gradN[j][1]*gradN[i][1] + mu[k] * gradN[j][0]*gradN[i][0] + mu[k] * gradN[j][2]*gradN[i][2]);
                real64 const Rzz_ij = detJ*((lambda[k]+2.0*mu[k])*gradN[j][2]*gradN[i][2] + mu[k] * gradN[j][1]*gradN[i][1] + mu[k] * gradN[j][0]*gradN[i][0]);
                real64 const Rxy_ij =  detJ*(mu[k] * gradN[j][1]*gradN[i][0] + lambda[k] * gradN[j][0]*gradN[i][1]);
                real64 const Ryx_ij =  detJ*(mu[k] * gradN[j][0]*gradN[i][1] + lambda[k] * gradN[j][1]*gradN[i][0]);
                real64 const Rxz_ij =  detJ*(mu[k] * gradN[j][2]*gradN[i][0] + lambda[k] * gradN[j][0]*gradN[i][2]);
                real64 const Rzx_ij =  detJ*(mu[k] * gradN[j][0]*gradN[i][2] + lambda[k] * gradN[j][2]*gradN[i][0]);
                real64 const Ryz_ij =  detJ* (mu[k] * gradN[j][2]*gradN[i][1] + lambda[k] * gradN[j][1]*gradN[i][2]);
                real64 const Rzy_ij =  detJ*(mu[k] * gradN[j][1]*gradN[i][2] + lambda[k] * gradN[j][2]*gradN[i][1]);


                stiffnessVector_x[elemsToNodes[k][i]] +=  (Rxx_ij * ux_n[elemsToNodes[k][j]] + Rxy_ij*uy_n[elemsToNodes[k][j]] + Rxz_ij*uz_n[elemsToNodes[k][j]]);
                stiffnessVector_y[elemsToNodes[k][i]] += (Ryx_ij * ux_n[elemsToNodes[k][j]] + Ryy_ij*uy_n[elemsToNodes[k][j]] + Ryz_ij*uz_n[elemsToNodes[k][j]]);
                stiffnessVector_z[elemsToNodes[k][i]] += ( Rzx_ij * ux_n[elemsToNodes[k][j]] + Rzy_ij*uy_n[elemsToNodes[k][j]] + Rzz_ij*uz_n[elemsToNodes[k][j]]);

		          }
	          }
          }

        } );

      } );
    } );
  } );

  addSourceToRightHandSide( time_n, rhs_x, rhs_y, rhs_z );

  /// Calculate your time integrators
  real64 dt2 = dt*dt;
  forAll< serialPolicy >( nodeManager.size(), [=] ( localIndex const a )
  {
    if( freeSurfaceNodeIndicator[a]!=1 )
    {

       ux_np1[a] = ux_n[a];
       ux_np1[a] *= 2.0*mass[a];
       ux_np1[a] -= (mass[a]-0.5*dt*damping_x[a])*ux_nm1[a];
       ux_np1[a] += dt2*(rhs_x[a]-stiffnessVector_x[a]);
       ux_np1[a] /= mass[a]+0.5*dt*damping_x[a];
       uy_np1[a] = uy_n[a];
       uy_np1[a] *= 2.0*mass[a];
       uy_np1[a] -= (mass[a]-0.5*dt*damping_y[a])*uy_nm1[a];
       uy_np1[a] += dt2*(rhs_y[a]-stiffnessVector_y[a]);
       uy_np1[a] /= mass[a]+0.5*dt*damping_y[a];
       uz_np1[a] = uz_n[a];
       uz_np1[a] *= 2.0*mass[a];
       uz_np1[a] -= (mass[a]-0.5*dt*damping_z[a])*uz_nm1[a];
       uz_np1[a] += dt2*(rhs_z[a]-stiffnessVector_z[a]);
       uz_np1[a] /= mass[a]+0.5*dt*damping_z[a];
    }
  } );

  /// Synchronize pressure fields
  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( "displacementx_np1" );
  fieldNames["node"].emplace_back( "displacementy_np1" );
  fieldNames["node"].emplace_back( "displacementz_np1" );

  CommunicationTools & syncFields = CommunicationTools::getInstance();
  syncFields.synchronizeFields( fieldNames,
                                domain.getMeshBody( 0 ).getMeshLevel( 0 ),
                                domain.getNeighbors(),
                                true );

  for( localIndex a=0; a<nodeManager.size(); ++a )
  {
    ux_nm1[a] = ux_n[a];
    uy_nm1[a] = uy_n[a];
    uz_nm1[a] = uz_n[a];
    ux_n[a] = ux_np1[a];
    uy_n[a] = uy_np1[a];
    uz_n[a] = uz_np1[a];

    stiffnessVector_x[a] = 0.0;
    stiffnessVector_y[a] = 0.0;
    stiffnessVector_z[a] = 0.0;
    rhs_x[a] = 0.0;
    rhs_y[a] = 0.0;
    rhs_z[a] = 0.0;
  }

  computeSismoTrace( cycleNumber, ux_np1, uy_np1, uz_np1 );


  return dt;

}

REGISTER_CATALOG_ENTRY( SolverBase, ElasticWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geosx */
