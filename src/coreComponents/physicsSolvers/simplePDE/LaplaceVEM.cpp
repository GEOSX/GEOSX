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
 * @file LaplaceVEM.cpp
 */

// Source includes
#include "LaplaceVEM.hpp"
#include "virtualElement/ConformingVirtualElementOrder1.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{}
}

using namespace dataRepository;

/*----------------------------------------------------------------------------------
 * LaplaceVEM: Solving Laplace's partial differential equation with virtual elements
 * ---------------------------------------------------------------------------------
 *
 * What does this solver do?
 * --------------------------
 *
 * This solver finds a solution f(x,y,z) to the Laplace equation: div ( grad ( f )) = 0
 * This common elliptic PDE represents the solution of a steady-state heat transfer, for instance.
 *
 * Where can I find an example of what it does?
 * --------------------------------------------
 *
 * Integrated tests associated to this solver are found in the ./integratedTests/ folder
 * These tests consist of computing the steady-state temperature profile in a simple cube-shaped domain
 * with fixed temperatures applied on two opposite cube faces ("Dirichlet" boundary conditions: imposing a value).
 * Feel free to run these tests cases, check out the XML input files, and inspect the output.
 *
 * Implementation: before we start:
 * ---------------------------------
 * In this implementation, the solution function (called above f) is called m_fieldName.
 * The variable m_fieldName is a string that points to a data container (an array) that
 * holds the numerical values of the PDE solution for each location at which f is evaluated.
 *
 * Let's take a look at the implementation step by step.
 *
 * ---------------------------------------------------------------------------------
 */


/* CONSTRUCTOR
   First, let us inspect the constructor of a "LaplaceVEM" object.
   This constructor does three important things:
   1 - It constructs an instance of the LaplaceVEM class (here: using the SolverBase constructor and passing through the arguments).
   2 - It sets some default values for the LaplaceVEM-specific private variables (here: m_fieldName and m_timeIntegrationOption).
   3 - It creates and activates a "registerWrapper" for each private variable.
   This is where the private variables are declared either as REQUIRED or OPTIONAL.
   An error is thrown if a REQUIRED variable is not specified in the XML file,
   along with the description of this variable and possible enum values if relevant.
   The description that is set is used in auto-generated documentation and console error messages.
 */

LaplaceVEM::LaplaceVEM( const string & name,
                        Group * const parent ):
  LaplaceBaseH1( name, parent )
{}

// Destructor
LaplaceVEM::~LaplaceVEM()
{
  // TODO Auto-generated destructor stub
}

/* SETUP SYSTEM
   Setting up the system using the base class method
 */
void LaplaceVEM::setupSystem( DomainPartition & domain,
                              DofManager & dofManager,
                              CRSMatrix< real64, globalIndex > & localMatrix,
                              array1d< real64 > & localRhs,
                              array1d< real64 > & localSolution,
                              bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;
  SolverBase::setupSystem( domain, dofManager, localMatrix, localRhs, localSolution, setSparsity );
}

/*
   ASSEMBLE SYSTEM
   This is the most important method to assemble the matrices needed before sending them to our solver.
   For a system A.x = B (with x the unknown), here, we use:
   - A : "localMatrix" this represents a Compressed Row Storage (optimized for sparse) matrix of real64 values associated with their index,
   - B : "localRhs" this represents a vector (1d array) of real64 numbers specified at the equation's right-hand side.
   The "local" prefix indicates that we are working on a local problem here, and the parallelization is performed at a higher level.
   This assembly step collects all the information needed to create the matrices localMatrix and localRhs, and the computation of values
   is done in a specific Laplace kernel optimized for parallel performance. Here we:
   1 - identify and point to the mesh of this domain,
   2 - find the node manager of this mesh,
   3 - extract the indices of the nodes that will be solved for (ie. the degrees of freedom or "dof")
   4 - pass all this information to a Laplace-specific finite element computation kernel.
   The call to the kernel is a templated call designed for performance (we will not explain the kernel here).
   See the implementation in LaplaceFEMKernel.cpp.
 */
void LaplaceVEM::assembleSystem( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                 real64 const GEOSX_UNUSED_PARAM( dt ),
                                 DomainPartition & domain,
                                 DofManager const & dofManager,
                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs )
{
  using VEM = virtualElement::ConformingVirtualElementOrder1< m_maxCellNodes, m_maxFaceNodes >;

  MeshLevel & mesh = domain.getMeshBodies().getGroup< MeshBody >( 0 ).getMeshLevel( 0 );
  NodeManager & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();

  // Get geometric properties used to compute projectors.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > nodesCoords =
    nodeManager.referencePosition();
  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const faceToEdgeMap = faceManager.edgeList().toViewConst();
  arrayView2d< localIndex const > const edgeToNodeMap = edgeManager.nodeList().toViewConst();
  arrayView2d< real64 const > const faceCenters = faceManager.faceCenter();
  arrayView2d< real64 const > const faceNormals = faceManager.faceNormal();
  arrayView1d< real64 const > const faceAreas = faceManager.faceArea();
  string const dofKey = dofManager.getKey( m_fieldName );
  arrayView1d< globalIndex const > const & dofIndex =
    nodeManager.getReference< array1d< globalIndex > >( dofKey );

  real64 const diffusion = 1.0;
  globalIndex const rankOffset = dofManager.rankOffset();

  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const,
                                       ElementRegionBase & elementRegion )
  {
    elementRegion.forElementSubRegions< CellElementSubRegion >
      ( [&]( CellElementSubRegion const & elemSubRegion )
    {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > elemToNodeMap = elemSubRegion.nodeList().toViewConst();
      arrayView2d< localIndex const > const elementToFaceMap = elemSubRegion.faceList().toViewConst();
      arrayView2d< real64 const > elemCenters = elemSubRegion.getElementCenter();
      arrayView1d< real64 const > elemVolumes = elemSubRegion.getElementVolume();
      arrayView1d< integer const > const & elemGhostRank = elemSubRegion.ghostRank();
      localIndex const numCells = elemSubRegion.size();
      forAll< parallelDevicePolicy< 32 > >( numCells, [=] GEOSX_HOST_DEVICE
                                              ( localIndex const cellIndex )
      {
        VEM::BasisData basisData;

        if( elemGhostRank[cellIndex] < 0 )
        {
          real64 cellVolume = elemVolumes[cellIndex];
          real64 const cellCenter[3] { elemCenters( cellIndex, 0 ),
                                       elemCenters( cellIndex, 1 ),
                                       elemCenters( cellIndex, 2 ) };
          // VEM::computeProjectors( cellIndex,
          //                         nodesCoords,
          //                         elemToNodeMap,
          //                         elementToFaceMap,
          //                         faceToNodeMap,
          //                         faceToEdgeMap,
          //                         edgeToNodeMap,
          //                         faceCenters,
          //                         faceNormals,
          //                         faceAreas,
          //                         cellCenter,
          //                         cellVolume,
          //                         basisData );

          VEM virtualElement;
          virtualElement.processLocalGeometry( cellIndex,
                                               nodesCoords,
                                               elemToNodeMap,
                                               elementToFaceMap,
                                               faceToNodeMap,
                                               faceToEdgeMap,
                                               edgeToNodeMap,
                                               faceCenters,
                                               faceNormals,
                                               faceAreas,
                                               cellCenter,
                                               cellVolume );
          basisData = virtualElement.m_basisData;

          real64 derivativesIntMean[VEM::maxSupportPoints][3] { { 0.0 } };
          globalIndex elemDofIndex[VEM::maxSupportPoints] { 0 };
          real64 element_matrix[VEM::maxSupportPoints][VEM::maxSupportPoints] { { 0.0 } };
          localIndex const numSupportPoints = VEM::getNumSupportPoints( basisData );
          for( localIndex a = 0; a < numSupportPoints; ++a )
          {
            elemDofIndex[a] = dofIndex[ elemToNodeMap( cellIndex, a ) ];
            for( localIndex b = 0; b < numSupportPoints; ++b )
            {
              element_matrix[a][b] = VEM::calcStabilizationValue( a, b, basisData );
            }
            localRhs[a] = 0.0;
          }
          for( localIndex q = 0; q < VEM::numQuadraturePoints; ++q )
          {
            VEM::calcGradN( q, basisData, derivativesIntMean );
            for( localIndex a = 0; a < numSupportPoints; ++a )
            {
              for( localIndex b = 0; b < numSupportPoints; ++b )
              {
                element_matrix[a][b] += diffusion*VEM::transformedQuadratureWeight( q, basisData ) *
                                        (derivativesIntMean[a][0] * derivativesIntMean[b][0] +
                                         derivativesIntMean[a][1] * derivativesIntMean[b][1] +
                                         derivativesIntMean[a][2] * derivativesIntMean[b][2] );
              }
            }
          }
          for( localIndex a = 0; a < numSupportPoints; ++a )
          {
            localIndex const dof = LvArray::integerConversion< localIndex >
                                     ( elemDofIndex[a] - rankOffset );
            if( dof < 0 || dof >= localMatrix.numRows() )
              continue;
            localMatrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >
              ( dof, elemDofIndex, element_matrix[ a ], numSupportPoints );
          }
        }
      } );
    } );
  } );

}

REGISTER_CATALOG_ENTRY( SolverBase, LaplaceVEM, string const &, Group * const )

} /* namespace geosx */
