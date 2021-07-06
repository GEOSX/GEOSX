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
 * @file AcousticWaveEquationSEM.cpp
 */

#include "AcousticWaveEquationSEM.hpp"

#include "dataRepository/KeyNames.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geosx
{

using namespace dataRepository;

AcousticWaveEquationSEM::AcousticWaveEquationSEM( const std::string & name,
                                                  Group * const parent ):
  SolverBase( name,
              parent )
{
  GEOSX_MARK_FUNCTION;

  registerWrapper( viewKeyStruct::sourceCoordinatesString(), &m_sourceCoordinates ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Coordinates (x,y,z) of the sources" );

  registerWrapper( viewKeyStruct::sourceNodeIdsString(), &m_sourceNodeIds ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Indices of the nodes (in the right order) for each source point" );

  registerWrapper( viewKeyStruct::sourceConstantsString(), &m_sourceConstants ).
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

  registerWrapper( viewKeyStruct::pressureNp1AtReceiversString(), &m_pressureNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Pressure value at each receiver for each timestep" );

  registerWrapper( viewKeyStruct::rickerOrderString(), &m_rickerOrder ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 2 ).
    setDescription( "Flag that indicates the order of the Ricker to be used o, 1 or 2. Order 2 by default" );

  registerWrapper( viewKeyStruct::outputSismoTraceString(), &m_outputSismoTrace ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag that indicates if we write the sismo trace in a file .txt, 0 no output, 1 otherwise" );

  registerWrapper( viewKeyStruct::dtSismoTraceString(), &m_dtSismoTrace ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Time step size to compute the sismo trace" );

  registerWrapper( viewKeyStruct::nSampleSismoTraceString(), &m_nSampleSismoTrace ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( " Number of sismo trace to be computed" );

  registerWrapper( viewKeyStruct::indexSismoTraceString(), &m_indexSismoTrace ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( " Index of the sismo trace " );
}

AcousticWaveEquationSEM::~AcousticWaveEquationSEM()
{
  // TODO Auto-generated destructor stub
}

void AcousticWaveEquationSEM::postProcessInput()
{
  GEOSX_MARK_FUNCTION;

  GEOSX_THROW_IF( m_sourceCoordinates.size( 1 ) != 3,
                  "Invalid number of physical coordinates for the sources",
                  InputError );

  GEOSX_THROW_IF( m_receiverCoordinates.size( 1 ) != 3,
                  "Invalid number of physical coordinates for the receivers",
                  InputError );

  localIndex const numNodesPerElem = 8;

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );

  m_sourceNodeIds.resizeDimension< 0, 1 >( numSourcesGlobal, numNodesPerElem );
  m_sourceConstants.resizeDimension< 0, 1 >( numSourcesGlobal, numNodesPerElem );
  m_sourceIsLocal.resizeDimension< 0 >( numSourcesGlobal );

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_receiverNodeIds.resizeDimension< 0, 1 >( numReceiversGlobal, numNodesPerElem );
  m_receiverConstants.resizeDimension< 0, 1 >( numReceiversGlobal, numNodesPerElem );
  m_receiverIsLocal.resizeDimension< 0 >( numReceiversGlobal );

  m_pressureNp1AtReceivers.resizeDimension< 0,1 >( numReceiversGlobal,m_nSampleSismoTrace );

}

void AcousticWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{
  GEOSX_MARK_FUNCTION;

  meshBodies.forSubGroups< MeshBody >( [&]( MeshBody & meshBody )
  {

    MeshLevel & meshLevel =  meshBody.getMeshLevel( 0 );

    NodeManager & nodeManager = meshLevel.getNodeManager();

    nodeManager.registerExtrinsicData< extrinsicMeshData::Pressure_nm1,
                                       extrinsicMeshData::Pressure_n,
                                       extrinsicMeshData::Pressure_np1,
                                       extrinsicMeshData::ForcingRHS,
                                       extrinsicMeshData::MassVector,
                                       extrinsicMeshData::DampingVector,
                                       extrinsicMeshData::StiffnessVector,
                                       extrinsicMeshData::FreeSurfaceNodeIndicator >( this->getName() );

    FaceManager & faceManager = meshLevel.getFaceManager();
    faceManager.registerExtrinsicData< extrinsicMeshData::FreeSurfaceFaceIndicator >( this->getName() );

    ElementRegionManager & elemManager = meshLevel.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerExtrinsicData< extrinsicMeshData::MediumVelocity >( this->getName() );
    } );

  } );
}


void AcousticWaveEquationSEM::precomputeSourceAndReceiverTerm( MeshLevel & mesh )
{
  GEOSX_MARK_FUNCTION;

  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

  arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();
  arrayView2d< localIndex > const sourceNodeIds = m_sourceNodeIds.toView();
  arrayView2d< real64 > const sourceConstants = m_sourceConstants.toView();
  arrayView1d< localIndex > const sourceIsLocal = m_sourceIsLocal.toView();
  sourceNodeIds.setValues< OMP_EXEC_POLICY >( -1 );
  sourceConstants.setValues< OMP_EXEC_POLICY >( -1 );
  sourceIsLocal.zero();

  arrayView2d< real64 const > const receiverCoordinates = m_receiverCoordinates.toViewConst();
  arrayView2d< localIndex > const receiverNodeIds = m_receiverNodeIds.toView();
  arrayView2d< real64 > const receiverConstants = m_receiverConstants.toView();
  arrayView1d< localIndex > const receiverIsLocal = m_receiverIsLocal.toView();
  receiverNodeIds.setValues< OMP_EXEC_POLICY >( -1 );
  receiverConstants.setValues< OMP_EXEC_POLICY >( -1 );
  receiverIsLocal.zero();

  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const,
                                       ElementRegionBase & elemRegion )
  {
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const,
                                                                       CellElementSubRegion & elementSubRegion )
    {

      GEOSX_THROW_IF( elementSubRegion.getElementTypeString() != "C3D8",
                      "Invalid type of element, the acoustic solver is designed for hexahedral meshes only (C3D8) ",
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

        forAll< OMP_EXEC_POLICY >( elementSubRegion.size(), [=, &elementSubRegion] ( localIndex const k )
        {
	  array1d< array1d< localIndex > > faceNodes( numFacesPerElem );

          for( localIndex kf = 0; kf < numFacesPerElem; ++kf )
          {
            elementSubRegion.getFaceNodes( k, kf, faceNodes[kf] );
          }

          /// loop over all the source that haven't been found yet
          for( localIndex isrc=0; isrc < sourceCoordinates.size( 0 ); isrc++ )
          {
            if( sourceIsLocal[isrc] == 0 )
            {
              real64 const coords[3] = { sourceCoordinates[isrc][0],
                                         sourceCoordinates[isrc][1],
                                         sourceCoordinates[isrc][2] };

              real64 coordsOnRefElem[3]{};
              bool const sourceFound = computeCoordinatesOnReferenceElement< FE_TYPE >( coords, coordsOnRefElem, k, faceNodes, elemsToNodes, X );
              if( sourceFound )
              {
                sourceIsLocal[isrc] = 1;
                real64 Ntest[8];
                finiteElement::LagrangeBasis1::TensorProduct3D::value( coordsOnRefElem, Ntest );


                for( localIndex a=0; a< numNodesPerElem; ++a )
                {
                  sourceNodeIds[isrc][a] = elemsToNodes[k][a];
                  sourceConstants[isrc][a] = Ntest[a];
                }
              }
            }
          } // End loop over all source


          /// loop over all the receiver that haven't been found yet
          for( localIndex ircv=0; ircv < receiverCoordinates.size( 0 ); ircv++ )
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
          } // End loop over receiver

        } ); // End loop over elements
      } );
    } );
  } );
}


void AcousticWaveEquationSEM::addSourceToRightHandSide( real64 const & time_n, arrayView1d< real64 > const rhs )
{
  GEOSX_MARK_FUNCTION;

  arrayView2d< localIndex const > const sourceNodeIds = m_sourceNodeIds.toViewConst();
  arrayView2d< real64 const > const sourceConstants   = m_sourceConstants.toViewConst();
  arrayView1d< localIndex const > const sourceIsLocal = m_sourceIsLocal.toViewConst();

  real64 const fi = evaluateRicker( time_n, this->m_timeSourceFrequency, this->m_rickerOrder );

  forAll< EXEC_POLICY >( sourceConstants.size( 0 ), [=] GEOSX_HOST_DEVICE ( localIndex const isrc )
  {
    if( sourceIsLocal[isrc] == 1 )
    {
      for( localIndex inode = 0; inode < sourceConstants.size( 1 ); ++inode )
      {
        rhs[sourceNodeIds[isrc][inode]] = sourceConstants[isrc][inode] * fi;
      }
    }
  } );
}


void AcousticWaveEquationSEM::computeSismoTrace(real64 const time_n, real64 const dt, localIndex const iSismo, arrayView1d< real64 > const pressure_np1, arrayView1d< real64 > const pressure_n )
{
  GEOSX_MARK_FUNCTION;

  real64 const time_sismo = m_dtSismoTrace*iSismo;
  real64 const time_np1 = time_n+dt;
  arrayView2d< localIndex const > const receiverNodeIds = m_receiverNodeIds.toViewConst();
  arrayView2d< real64 const > const receiverConstants   = m_receiverConstants.toViewConst();
  arrayView1d< localIndex const > const receiverIsLocal = m_receiverIsLocal.toViewConst();

  arrayView2d< real64 > const p_rcvs   = m_pressureNp1AtReceivers.toView();

  forAll< serialPolicy >( receiverConstants.size( 0 ), [=] ( localIndex const ircv )
  {
    if( receiverIsLocal[ircv] == 1 )
    {
      p_rcvs[ircv][iSismo] = 0.0;
      real64 ptmp_np1 = 0.0;
      real64 ptmp_n = 0.0;
      for( localIndex inode = 0; inode < receiverConstants.size( 1 ); ++inode )
      {
	ptmp_np1 += pressure_np1[receiverNodeIds[ircv][inode]]*receiverConstants[ircv][inode];
	ptmp_n += pressure_n[receiverNodeIds[ircv][inode]]*receiverConstants[ircv][inode];
      }
      p_rcvs[ircv][iSismo] = ((time_np1 - time_sismo)*ptmp_n+(time_sismo-time_n)*ptmp_np1)/dt;
    }
  } );

}

/// Use for now until we get the same functionality in TimeHistory
void AcousticWaveEquationSEM::saveSismo( arrayView2d< real64 const > const pressure_receivers )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< localIndex const > const receiverIsLocal = m_receiverIsLocal.toViewConst();

  forAll< serialPolicy >( pressure_receivers.size( 0 ), [=] ( localIndex const ircv )
  {
    if( receiverIsLocal[ircv] == 1 )
    {
      char filename[50];
      sprintf( filename, "sismoTraceReceiver%0d.txt", static_cast< int >( ircv ) );
      std::ofstream f( filename, std::ios::app );

      for( localIndex iSismo=0; iSismo < pressure_receivers.size( 1 ); iSismo++ )
      {
	//real64 timeSismo = iSismo*m_dtSismoTrace;
	f<< iSismo << " " << pressure_receivers[ircv][iSismo] << std::endl;
      }

      f.close();
    }
  } );
}


void AcousticWaveEquationSEM::initializePreSubGroups()
{
  GEOSX_MARK_FUNCTION;

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


void AcousticWaveEquationSEM::initializePostInitialConditionsPreSubGroups()
{
  GEOSX_MARK_FUNCTION;

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  real64 const time = 0.0;
  applyFreeSurfaceBC( time, domain );
  precomputeSourceAndReceiverTerm( mesh );

  NodeManager & nodeManager = mesh.getNodeManager();
  FaceManager & faceManager = mesh.getFaceManager();

  /// get the array of indicators: 1 if the face is on the boundary; 0 otherwise
  arrayView1d< integer > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

  /// Get table containing all the face normals
  arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();
  ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

  arrayView1d< real64 > const mass = nodeManager.getExtrinsicData< extrinsicMeshData::MassVector >();

  /// damping matrix to be computed for each dof in the boundary of the mesh
  arrayView1d< real64 > const damping = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector >();

  damping.zero();

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

      arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();

      arrayView1d< real64 > const c = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocity >();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();
        localIndex const numNodesPerFace = 4;

        /// Loop over elements
        forAll< serialPolicy >( elemsToNodes.size( 0 ), [=] ( localIndex const k )
        {
          constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
          constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

          real64 const invC2 = 1.0 / ( c[k] * c[k] );
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
              mass[elemsToNodes[k][a]] +=  invC2 * detJ * N[a];
            }
          }

          real64 const alpha = 1.0/c[k];

          for( localIndex kfe=0; kfe< numFacesPerElem; ++kfe )
          {
            localIndex const numFaceGl = elemsToFaces[k][kfe];

            /// Face on the domain boundary and not on free surface
            if( facesDomainBoundaryIndicator[numFaceGl]==1 && freeSurfaceFaceIndicator[numFaceGl]!=1 )
            {
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
                  for( localIndex i=0; i<3; ++i )
                  {
                    for( localIndex j = 0; j < 3; ++j )
                    {
                      tmp[i] += invJ[j][i]*faceNormal[numFaceGl][j];
                    }
                    ds +=tmp[i]*tmp[i];
                  }
                  ds = std::sqrt( ds );

                  localIndex numNodeGl = facesToNodes[numFaceGl][a];
                  damping[numNodeGl] += alpha*detJ*ds*N[a];
                }
              }
            }
          }
        } ); // end loop over element
      } );
    } );
  } );

}


void AcousticWaveEquationSEM::applyFreeSurfaceBC( real64 const time, DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getNodeManager();

  arrayView1d< real64 > const p_nm1 = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_nm1 >();
  arrayView1d< real64 > const p_n = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_n >();
  arrayView1d< real64 > const p_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_np1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceFaceIndicator >();

  /// array of indicators: 1 if a node is on on free surface; 0 otherwise
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

      freeSurfaceFaceIndicator.zero();
      freeSurfaceNodeIndicator.zero();

      for( localIndex i = 0; i < targetSet.size(); ++i )
      {
        localIndex const kf = targetSet[ i ];
        freeSurfaceFaceIndicator[kf] = 1;

        localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
        for( localIndex a=0; a < numNodes; ++a )
        {
          localIndex const dof = faceToNodeMap( kf, a );
          freeSurfaceNodeIndicator[dof] = 1;

          p_np1[dof] = value;
          p_n[dof]   = value;
          p_nm1[dof] = value;
        }
      }
    }
    else
    {
      GEOSX_ERROR( "This option is not supported yet" );
    }
  } );
}



real64 AcousticWaveEquationSEM::solverStep( real64 const & time_n,
                                            real64 const & dt,
                                            integer const cycleNumber,
                                            DomainPartition & domain )
{
  return explicitStep( time_n, dt, cycleNumber, domain );
}

real64 AcousticWaveEquationSEM::evaluateRicker( real64 const & time_n, real64 const & f0, localIndex order )
{
  GEOSX_MARK_FUNCTION;

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
      pulse = -(time_n-o_tpeak)*exp( -2*lam*(time_n-o_tpeak)*(time_n-o_tpeak) );
    }
    break;
    default:
      GEOSX_ERROR( "This option is not supported yet, rickerOrder must be 0, 1 or 2" );
  }

  return pulse;
}

template< typename FE_TYPE >
struct StiffnessVectorKernel
{
public:
  StiffnessVectorKernel( FE_TYPE const & finiteElementSpace )
    :
    m_finiteElementSpace( finiteElementSpace )
  {}

  FE_TYPE const & m_finiteElementSpace;

  template< typename POLICY>
  void compute( arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
		arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
		arrayView1d< real64 const > const p_n,
		arrayView1d< real64 > const stiffnessVector )
  {
    constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
    constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

    forAll< POLICY >( elemsToNodes.size( 0 ), [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      real64 xLocal[numNodesPerElem][3];

      for( localIndex a=0; a< numNodesPerElem; ++a )
      {
	for( localIndex i=0; i<3; ++i )
	{
	  xLocal[a][i] = X( elemsToNodes( k, a ), i );
	}
      }

      real64 gradN[ numNodesPerElem ][ 3 ];
      for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
      {
	real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

	for( localIndex i=0; i<numNodesPerElem; ++i )
	{
	  for( localIndex j=0; j<numNodesPerElem; ++j )
	  {
	    real64 const Rh_ij = detJ * LvArray::tensorOps::AiBi< 3 >( gradN[ i ], gradN[ j ] );

	    stiffnessVector[elemsToNodes[k][i]] += Rh_ij*p_n[elemsToNodes[k][j]];
	  }
	}
      }
    } );
  }
};


real64 AcousticWaveEquationSEM::explicitStep( real64 const & time_n,
                                              real64 const & dt,
                                              integer const cycleNumber,
                                              DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( time_n, dt, cycleNumber );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  ///Use to reinit params if pygeox set time_n to zero
  if( time_n < 0.5*dt )
  {
    applyFreeSurfaceBC( time_n, domain );
    postProcessInput();
    precomputeSourceAndReceiverTerm( mesh );
  }

  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real64 const > const mass = nodeManager.getExtrinsicData< extrinsicMeshData::MassVector >();
  arrayView1d< real64 const > const damping = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector >();

  arrayView1d< real64 > const p_nm1 = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_nm1 >();
  arrayView1d< real64 > const p_n = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_n >();
  arrayView1d< real64 > const p_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_np1 >();

  /// get array of indicators: 1 if node on free surface; 0 otherwise
  arrayView1d< localIndex const > const freeSurfaceNodeIndicator = nodeManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceNodeIndicator >();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

  arrayView1d< real64 > const stiffnessVector = nodeManager.getExtrinsicData< extrinsicMeshData::StiffnessVector >();

  arrayView1d< real64 > const rhs = nodeManager.getExtrinsicData< extrinsicMeshData::ForcingRHS >();

  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const,
                                       ElementRegionBase & elemRegion )
  {
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const,
                                                                       CellElementSubRegion & elementSubRegion )
    {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );
	GEOSX_MARK_SCOPE (stiffness);
	StiffnessVectorKernel<FE_TYPE> kernel( finiteElement );
	kernel.template compute< EXEC_POLICY >( elemsToNodes,
						X,
						p_n,
						stiffnessVector );

      } );
    } );
  } );

  addSourceToRightHandSide( time_n, rhs );

  /// Calculate your time integrators
  real64 const dt2 = dt*dt;
  {
  GEOSX_MARK_SCOPE (updateP);
  forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    if( freeSurfaceNodeIndicator[a]!=1 )
    {
      p_np1[a] = p_n[a];
      p_np1[a] *= 2.0*mass[a];
      p_np1[a] -= (mass[a]-0.5*dt*damping[a])*p_nm1[a];
      p_np1[a] += dt2*(rhs[a]-stiffnessVector[a]);
      p_np1[a] /= mass[a]+0.5*dt*damping[a];
    }
  } );
  }

  /// Synchronize pressure fields
  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( "pressure_np1" );

  CommunicationTools & syncFields = CommunicationTools::getInstance();
  syncFields.synchronizeFields( fieldNames,
                                domain.getMeshBody( 0 ).getMeshLevel( 0 ),
                                domain.getNeighbors(),
                                true );

  forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOSX_HOST_DEVICE (localIndex a )
  {
    p_nm1[a]=p_n[a];
    p_n[a] = p_np1[a];

    stiffnessVector[a] = 0.0;
    rhs[a] = 0.0;
  });

  real64 checkSismo = m_dtSismoTrace*m_indexSismoTrace;
  real64 epsilonLoc = 1e-12;
  if( (time_n-epsilonLoc) <= checkSismo &&  checkSismo < (time_n + dt) )
  {
    computeSismoTrace( time_n, dt, m_indexSismoTrace, p_np1, p_n );
    m_indexSismoTrace ++;
  }

  // Note: this "manual" output to file is temporary
  //       It should be removed as soon as we can use TimeHistory to output data not registered on the mesh
  // TODO: remove the (sprintf+saveSismo) and replace with TimeHistory
  if( m_outputSismoTrace == 1 )
  {
    if( m_indexSismoTrace % m_nSampleSismoTrace == 0 )
    {
      arrayView2d< real64 const > const p_rcvs = m_pressureNp1AtReceivers.toView();
      this->saveSismo( p_rcvs );
    }
  }

  return dt;
}


REGISTER_CATALOG_ENTRY( SolverBase, AcousticWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geosx */
