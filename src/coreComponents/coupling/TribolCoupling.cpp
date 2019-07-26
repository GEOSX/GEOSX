/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file TribolCoupling.cpp
 */

#if HAVE_TRIBOLCOUPLING

#include "TribolCoupling.hpp"
#include "dataRepository/KeyNames.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/MeshBody.hpp"

#define GLOBALID_IS_64BIT 1
#include "worlds/SlideWorldAdapter.h"

#include <algorithm>
#include "tribol/tribol.hpp"
#include "axom/mint/mesh/CellTypes.hpp"

static void GEOSXSlideWorldErrorHandler(const char* msg, int etype, int)
{
   switch (etype) {
      default: // undefined values default to VERR_WARN
      case VERR_WARN:
         std::cout<<"***** TRIBOL COUPLING WARNING " << std::endl << msg << std::endl ;
         break ;
      case VERR_FATAL:
         GEOS_ERROR(msg);
         break ;
   }
}

namespace geosx
{

vista::View *s_tribolProblem = nullptr ;
vista::View *s_tribolDomain = nullptr ;
vista::View *s_srcNodes = nullptr ;
vista::View *s_srcFaces = nullptr ;
vista::View *s_srcBricks = nullptr ;
vista::View *s_slideWorldSourceNodes = nullptr ;
vista::View *s_slideWorldSourceFaces = nullptr ;

void TribolCoupling::Initialize(dataRepository::ManagedGroup * eventManager, dataRepository::ManagedGroup * domain)
{
  vista::ErrSetMask(VERR_CONVERSION | VERR_TRUNCATION | VERR_NOTCREATED | VERR_INTERNAL | VERR_NILNAME | VERR_INVALID_INDEXSET) ;
  vista::ErrHandler(GEOSXSlideWorldErrorHandler) ;

  const int nodesPerFace = 4 ;
  const int nodesPerElem = 8 ;

  real64& currentTime = *(eventManager->getData<real64>("time"));
  real64& dt = *(eventManager->getData<real64>("dt"));
  integer& cycle = *(eventManager->getData<integer>("cycle"));

  // Currently we need the previous dt, but we assume fixed dt for now.
  // We probably want to pass in prevDt in the update anyway.
  s_tribolProblem = new vista::View("geosx", 0) ;
  s_tribolDomain = s_tribolProblem->viewCreate("DomainWorld", 0)->viewCreate("domain", 0, 0) ;

  // Assuming only a single element region
  MeshLevel * const meshLevel = domain->group_cast<DomainPartition*>()->getMeshBody(0)->getMeshLevel(0);
  CellBlockSubRegion * const subRegion = meshLevel->getElemManager()->GetRegion(0)->GetSubRegion(0);
  NodeManager * const nodeManager = meshLevel->getNodeManager();
  FaceManager * const faceManager = meshLevel->getFaceManager();

  int numFaces = faceManager->size() ;
  int numNodes = nodeManager->size() ;
  int numBricks = subRegion->size() ;
  const globalIndex_array faceMap = faceManager->m_localToGlobalMap ;
  globalIndex_array nodeMap = nodeManager->m_localToGlobalMap ;
  globalIndex_array brickMap = subRegion->m_localToGlobalMap ;

  const array2d<localIndex> & facesToElems = faceManager->elementList();
  const OrderedVariableOneToManyRelation & facesToNodes = faceManager->nodeList();
  const integer_array isExternalFace = faceManager->m_isExternal ;
  const integer_array faceGhostRank = faceManager->m_ghostRank ;

  globalID *extFaceMap ;
  int *slideNodeMap ;
  int *extFacesToNodes ;
  int *extFacesToElement ;
  vista::MemAlloc(numFaces*nodesPerFace, &slideNodeMap) ;
  vista::MemAlloc(numFaces, &extFaceMap) ;
  vista::MemAlloc(numFaces*nodesPerFace, &extFacesToNodes) ;
  vista::MemAlloc(numFaces*2, &extFacesToElement) ;

  int numSlideNodes = 0 ;
  int numExtFaces = 0 ;

  for (int i = 0 ; i < numFaces ; ++i) {
     // Only include external, non-ghost faces/nodes.
     if (isExternalFace[i] && faceGhostRank[i] < 0) {
        extFaceMap[numExtFaces] = globalID((GIDTYPE)faceMap[i]) ;
        for (int j = 0 ; j < nodesPerFace ; ++j) {
           int nodeIdx = facesToNodes[i][j] ;
           extFacesToNodes[numExtFaces*4+j] = nodeIdx ;
           slideNodeMap[numSlideNodes++] = nodeIdx ;
        }
        for (int j = 0 ; j < 2 ; ++j) {
           extFacesToElement[numExtFaces*2+j] = facesToElems[i][j] ;
        }

        ++numExtFaces ;
     }
  }
  std::sort(slideNodeMap, slideNodeMap + numSlideNodes) ;
  int * newEnd = std::unique(slideNodeMap, slideNodeMap + numSlideNodes) ;

  numSlideNodes = newEnd - slideNodeMap ;
  vista::MemRealloc(numSlideNodes, &slideNodeMap) ;
  vista::MemRealloc(numExtFaces, &extFaceMap) ;
  vista::MemRealloc(numExtFaces*nodesPerFace, &extFacesToNodes) ;
  vista::MemRealloc(numExtFaces*2, &extFacesToElement) ;

  s_srcFaces = s_tribolDomain->viewCreate("extFaces", 1, &numExtFaces, new vista::IndexSet(vista::VISTA_ACQUIRES, numExtFaces, extFaceMap)) ;
  s_srcNodes = s_tribolDomain->viewCreate("nodes", 1, &numNodes, new vista::IndexSet(vista::VISTA_COPIES, numNodes, (globalID*) nodeMap.data())) ;
  s_srcBricks = s_tribolDomain->viewCreate("bricks", 1, &numBricks, new vista::IndexSet(vista::VISTA_COPIES, numBricks, (globalID*) brickMap.data())) ;

  vista::View* localNodes = s_srcNodes->viewCreate("local", 1) ;

  // We mimic the ALE3D slide data structure for convenience here.
  // This could be optimized by creating a more generic SetSourceData call.
  s_slideWorldSourceNodes = s_srcNodes->attach(new vista::View("LSslaveNodes0", 1, &numSlideNodes, new vista::IndexSet(numSlideNodes, slideNodeMap, vista::VISTA_ACQUIRES))) ;

  s_slideWorldSourceFaces = s_srcFaces->attach(new vista::View("LSslaveFaces0", 1, &numExtFaces, new vista::IndexSet(numExtFaces))) ;
  // since we are including all extFaces, we can just copy the extFacesToNodes relation for our slide facesToDomainNodes relation.
  s_slideWorldSourceFaces->relationCreateFixed("facesToDomainNodes", nodesPerFace, extFacesToNodes, s_srcNodes, vista::VISTA_COPIES) ;

  s_srcFaces->relationCreateFixed("facesToNodes", nodesPerFace, extFacesToNodes, s_srcNodes, vista::VISTA_ACQUIRES) ;
  s_srcFaces->relationCreateFixed("facesToElement", 2, extFacesToElement, s_srcBricks, vista::VISTA_ACQUIRES) ;

  const array1d< R1Tensor > & X = nodeManager->referencePosition();
  
  real64 *x0 = s_srcNodes->fieldCreateReal("x0")->Real() ;
  real64 *y0 = s_srcNodes->fieldCreateReal("y0")->Real() ;
  real64 *z0 = s_srcNodes->fieldCreateReal("z0")->Real() ;

  // These are used directly on the srcNodes by the slide decomposition,
  // so they need to have the proper tree structure.
  s_srcNodes->fieldCreateReal("x") ;
  s_srcNodes->fieldCreateReal("y") ;
  s_srcNodes->fieldCreateReal("z") ;
  s_srcNodes->fieldCreateReal("xd") ;
  s_srcNodes->fieldCreateReal("yd") ;
  s_srcNodes->fieldCreateReal("zd") ;
  s_srcNodes->fieldCreateReal("local:xdd") ;
  s_srcNodes->fieldCreateReal("local:ydd") ;
  s_srcNodes->fieldCreateReal("local:zdd") ;

  const FixedOneToManyRelation& elemsToNodes = subRegion->nodeList() ;

  // We need initial positions for the nodes of all bricks adjacent to the
  // slide faces (this is needed to determine the orientation of the faces).

  for (int i = 0 ; i < numExtFaces ; ++i) {
     int faceElem = extFacesToElement[2*i] ;
     localIndex const * const elemNodes = elemsToNodes[faceElem] ;

     for (int j = 0 ; j < nodesPerElem ; ++j) {
        int nodeIdx = elemNodes[j] ;
        x0[nodeIdx] = X[nodeIdx][0] ;
        y0[nodeIdx] = X[nodeIdx][1] ;
        z0[nodeIdx] = X[nodeIdx][2] ;
     }
  }

  s_srcNodes->fieldCreateReal("fx") ;
  s_srcNodes->fieldCreateReal("fy") ;
  s_srcNodes->fieldCreateReal("fz") ;

  CopyPositionsToTribolSourceData(nodeManager) ;
  CopyForcesToTribolSourceData(nodeManager) ;

  int *bricksToNodes ;
  vista::MemAlloc(numBricks*nodesPerElem, &bricksToNodes) ;
  int bricksToNodesIdx = 0 ;

  for (int i = 0 ; i < numBricks ; ++i) {
     localIndex const * const elemNodes = elemsToNodes[i] ;

     for (int j = 0 ; j < nodesPerElem ; ++j) {
         bricksToNodes[bricksToNodesIdx++] = elemNodes[j] ;
     }
  }

  s_srcBricks->relationCreateFixed("bricksToNodes", nodesPerElem, bricksToNodes, s_srcNodes, vista::VISTA_ACQUIRES) ;

  SlideWorldAdapter::CreateWorld(MPI_COMM_GEOSX, MPI_COMM_WORLD,
                                 4000, // comm tag 
                                 3, // dimension
                                 0, // axisym
                                 1, // numSS
                                 nodesPerFace,
                                 nodesPerElem,
                                 &cycle,
                                 &dt,
                                 &dt, // prevDt
                                 1, // noParamInput
                                 true, // sort views
                                 s_tribolProblem) ;

  SlideWorldAdapter::SetSourceData(s_tribolProblem) ;
}

void TribolCoupling::ApplyTribolForces(dataRepository::ManagedGroup * domain,
                                       real64 const& time_n,
                                       real64 const& dt,
                                       const int cycleNumber)
{
   MeshLevel * const meshLevel = domain->group_cast<DomainPartition*>()->getMeshBody(0)->getMeshLevel(0);
   NodeManager * const nodeManager = meshLevel->getNodeManager();

   CopyPositionsToTribolSourceData(nodeManager) ;
   CopyVelocitiesToTribolSourceData(nodeManager) ;
   CopyAccelerationsToTribolSourceData(nodeManager) ;
   CopyForcesToTribolSourceData(nodeManager) ;

   int meshDim = 3 ;
   int nodesPerFace = 4 ;

   // set the cellType
   const int cellType = (meshDim == 2) ? (int)axom::mint::CellType::SEGMENT : (int)axom::mint::CellType::QUAD;

   SlideWorldAdapter::Setup(0) ;
   SlideWorldAdapter::TransferToWorld() ;

   vista::View *slideWorldData = SlideWorldAdapter::GetWorldData() ;

   // Just handling one slide for now.
   vista::View *slideWorldNodesSlave = slideWorldData->view("slidenodes0s") ;
   vista::View *slideWorldFacesSlave = slideWorldData->view("slidefaces0s") ;

   int numSlaveFaces = slideWorldFacesSlave->length() ;

   const vista::FixedSizeRelation *facesToNodesSlave = slideWorldFacesSlave->relationFixed("facesToNodes") ;
   const int* facesToNodesSlaveData = facesToNodesSlave->data() ;

   real64 *fxs = slideWorldNodesSlave->fieldReal("fx") ;
   real64 *fys = slideWorldNodesSlave->fieldReal("fy") ;
   real64 *fzs = slideWorldNodesSlave->fieldReal("fz") ;
   const real64 *xs = slideWorldNodesSlave->fieldReal("x") ;
   const real64 *ys = slideWorldNodesSlave->fieldReal("y") ;
   const real64 *zs = slideWorldNodesSlave->fieldReal("z") ;

   // initialize the contact library
   tribol::initialize(meshDim, MPI_COMM_WORLD);

   // set penalty stiffness, set to negative to account for flipped sign in the
   // contact nodal forces. This should be corrected in the normal that is
   // used in tribol penalty.
   tribol::setPenaltyStiffness(1.0);

   // register the current configuration mesh
   const int slaveMesh = 0 ;
   tribol::registerMesh(slaveMesh, numSlaveFaces, nodesPerFace,
                        facesToNodesSlaveData, cellType, xs, ys, zs);

   // register nodal response (i.e. contact forces)
   tribol::registerNodalResponse(slaveMesh, fxs, fys, fzs);

   // register faces in contact
   int maxNumPairs = numSlaveFaces*(numSlaveFaces-1) ;
   int *meshId1 = new int [maxNumPairs] ;
   int *meshId2 = new int [maxNumPairs] ;
   int *pairType = new int [maxNumPairs] ;
   int *pairIndex1 = new int [maxNumPairs] ;
   int *pairIndex2 = new int [maxNumPairs] ;

   int numPairs = 0 ;

   for (int i = 0 ; i < numSlaveFaces ; ++i) {
      // Ensure only faces with all nodes defined are included.
      // We need this because the faces and nodes are decomposed separately.
      bool validFace = true ;
      const int * faceNodesSlave = (*facesToNodesSlave)[i] ;
      for (int k = 0 ; k < nodesPerFace ; ++k) {
         if (faceNodesSlave[k] == -1) {
            validFace = false ;
            break ;
         }
      }

      if (validFace) {
         int startIdx = i + 1 ;
         int endIdx = numSlaveFaces ;
         int jMesh = slaveMesh ;
         const vista::FixedSizeRelation *facesToNodes = facesToNodesSlave ;

         for (int j = startIdx ; j < endIdx; ++j) {
            bool validFace2 = true ;
            const int * faceNodes2 = (*facesToNodes)[j] ;
            for (int k = 0 ; k < nodesPerFace ; ++k) {
               if (faceNodes2[k] == -1) {
                  validFace2 = false ;
                  break ;
               }
            }

            if (validFace2) {
               meshId1[numPairs] = slaveMesh ;
               meshId2[numPairs] = jMesh ;
               pairType[numPairs] = cellType ;
               pairIndex1[numPairs] = i ;
               pairIndex2[numPairs] = j ;
               ++numPairs ;
            }
         }
      }
   }

   //tribol::setInterfacePairs(numPairs, meshId1, pairType, pairIndex1, meshId2, pairType, pairIndex2) ;

   // register the coupling scheme
   tribol::registerCouplingScheme(0, slaveMesh, slaveMesh,
                                  tribol::SURFACE_TO_SURFACE,
                                  tribol::AUTO,
                                  tribol::COMMON_PLANE,
                                  tribol::FRICTIONLESS,
                                  tribol::PENALTY,
                                  tribol::BINNING_CARTERIAN_PRODUCT);

   // call the contact update routine. For the serial problem it doesn't matter if this is
   // in the domain loop or outside
   bool outputCycle = false ;
   int err = tribol::update(cycleNumber, time_n, dt, outputCycle);

   if (err) {
      GEOS_ERROR( "TRIBOL update error" );
   }

   tribol::finalize();

   delete[] meshId1 ;
   delete[] meshId2 ;
   delete[] pairType ;
   delete[] pairIndex1 ;
   delete[] pairIndex2 ;

   SlideWorldAdapter::TransferFromWorld() ;

   SlideWorldAdapter::Cleanup() ;

   CopyAccelerationsFromTribolSourceData(nodeManager) ;
}

void TribolCoupling::CopyPositionsToTribolSourceData(NodeManager const * const nodeManager)
{
  const array1d< R1Tensor > & X = nodeManager->referencePosition();
  const array1d< R1Tensor > & u = nodeManager->getReference< array1d< R1Tensor > >("TotalDisplacement");
  const int numSlideNodes = s_slideWorldSourceNodes->length() ;
  const int *slideMap = s_slideWorldSourceNodes->map() ;
  real64 *x = s_srcNodes->fieldReal("x") ;
  real64 *y = s_srcNodes->fieldReal("y") ;
  real64 *z = s_srcNodes->fieldReal("z") ;

  for (int i = 0 ; i < numSlideNodes ; ++i) {
     int nodeIdx = slideMap[i] ;
     const R1Tensor & Xi = X[nodeIdx] ;
     const R1Tensor & ui = u[nodeIdx] ;
     x[nodeIdx] = Xi[0] + ui[0] ;
     y[nodeIdx] = Xi[1] + ui[1] ;
     z[nodeIdx] = Xi[2] + ui[2] ;
  }
}

void TribolCoupling::CopyVelocitiesToTribolSourceData(NodeManager const * const nodeManager)
{
  const array1d< R1Tensor > & v = nodeManager->getReference< array1d< R1Tensor > >("Velocity");
  const int numSlideNodes = s_slideWorldSourceNodes->length() ;
  const int *slideMap = s_slideWorldSourceNodes->map() ;
  real64 *xd = s_srcNodes->fieldReal("xd") ;
  real64 *yd = s_srcNodes->fieldReal("yd") ;
  real64 *zd = s_srcNodes->fieldReal("zd") ;

  for (int i = 0 ; i < numSlideNodes ; ++i) {
     int nodeIdx = slideMap[i] ;
     const R1Tensor & vi = v[nodeIdx] ;
     xd[nodeIdx] = vi[0] ;
     yd[nodeIdx] = vi[1] ;
     zd[nodeIdx] = vi[2] ;
  }
}

void TribolCoupling::CopyAccelerationsToTribolSourceData(NodeManager const * const nodeManager)
{
  const array1d< R1Tensor > & a = nodeManager->getReference< array1d< R1Tensor > >("Acceleration");
  const int numSlideNodes = s_slideWorldSourceNodes->length() ;
  const int *slideMap = s_slideWorldSourceNodes->map() ;
  real64 *xdd = s_srcNodes->fieldReal("local:xdd") ;
  real64 *ydd = s_srcNodes->fieldReal("local:ydd") ;
  real64 *zdd = s_srcNodes->fieldReal("local:zdd") ;

  for (int i = 0 ; i < numSlideNodes ; ++i) {
     int nodeIdx = slideMap[i] ;
     const R1Tensor & ai = a[nodeIdx] ;
     xdd[nodeIdx] = ai[0] ;
     ydd[nodeIdx] = ai[1] ;
     zdd[nodeIdx] = ai[2] ;
  }
}

void TribolCoupling::CopyForcesToTribolSourceData(NodeManager const * const nodeManager)
{
  const array1d< R1Tensor > & a = nodeManager->getReference< array1d< R1Tensor > >( dataRepository::keys::Acceleration );
  const array1d< real64 > & nodalMass = nodeManager->getReference< array1d< real64 > >( dataRepository::keys::Mass );
  const int numSlideNodes = s_slideWorldSourceNodes->length() ;
  const int *slideMap = s_slideWorldSourceNodes->map() ;
  real64 *fx = s_srcNodes->fieldReal("fx") ;
  real64 *fy = s_srcNodes->fieldReal("fy") ;
  real64 *fz = s_srcNodes->fieldReal("fz") ;

  for (int i = 0 ; i < numSlideNodes ; ++i) {
     int nodeIdx = slideMap[i] ;
     const real64 nmass = nodalMass[nodeIdx] ;
     const R1Tensor & ai = a[nodeIdx] ;
     fx[nodeIdx] = ai[0]*nmass ;
     fy[nodeIdx] = ai[1]*nmass ;
     fz[nodeIdx] = ai[2]*nmass ;
  }
}

void TribolCoupling::CopyAccelerationsFromTribolSourceData(NodeManager * const nodeManager)
{
  array1d< R1Tensor > & a = nodeManager->getReference< array1d< R1Tensor > >( dataRepository::keys::Acceleration );
  const array1d< real64 > & nodalMass = nodeManager->getReference< array1d< real64 > >( dataRepository::keys::Mass );
  const int numSlideNodes = s_slideWorldSourceNodes->length() ;
  const int *slideMap = s_slideWorldSourceNodes->map() ;
  const real64 *fx = s_srcNodes->fieldReal("fx") ;
  const real64 *fy = s_srcNodes->fieldReal("fy") ;
  const real64 *fz = s_srcNodes->fieldReal("fz") ;

  for (int i = 0 ; i < numSlideNodes ; ++i) {
     int nodeIdx = slideMap[i] ;
     const real64 nmass = nodalMass[nodeIdx] ;
     R1Tensor & ai = a[nodeIdx] ;
     ai[0] = fx[nodeIdx]/nmass ;
     ai[1] = fy[nodeIdx]/nmass ;
     ai[2] = fz[nodeIdx]/nmass ;
  }
}

void TribolCoupling::SyncTermination(int* terminate)
{
   SlideWorldAdapter::SyncTermination(terminate) ;
}

void TribolCoupling::SyncTimestep(real64* newDt)
{
    // Currently does not change GEOS timestep.
    // Will also be incorrect if timestep is adjusted by maxTime.
    SlideWorldAdapter::GetRequiredTimestep(newDt) ;

    // Ignore GEOS requests
    int plotThisCycle = 0 ;
    int dumpThisCycle = 0 ;
    int advectCountThisCycle = 0 ;
    int slideDeleteThisCycle = 0 ;

    SlideWorldAdapter::SyncControl(&slideDeleteThisCycle, &advectCountThisCycle, &plotThisCycle, &dumpThisCycle) ;

    if (slideDeleteThisCycle) {
        SlideWorldAdapter::Cleanup();
        SlideWorldAdapter::SetSourceData(s_tribolProblem, true);
    }
}

void TribolCoupling::Cleanup()
{
   SlideWorldAdapter::DestroyWorld() ;
   delete s_tribolProblem ;
   s_tribolProblem = nullptr ;
   s_tribolDomain = nullptr ;
   s_srcNodes = nullptr ;
   s_srcFaces = nullptr ;
   s_srcBricks = nullptr ;
   s_slideWorldSourceNodes = nullptr ;
   s_slideWorldSourceFaces = nullptr ;
}

#ifdef GEOSX_USE_MPI
// @author Tony De Groot
void TribolCoupling::InitCommSubset(MPI_Comm const mpiComm,
                                    MPI_Comm *myComm,
                                    MPI_Comm *otherComm,
                                    int myCode)
{
   int numProcs, myProcNum ;
   MPI_Group mpiGroup, myGroup, otherGroup;

   // Each processor sends an identification message to processor 0
   // Processor 0 determines which code is on which processor

   MPI_Comm_rank(mpiComm, &myProcNum);
   MPI_Comm_size(mpiComm, &numProcs);
   MPI_Comm_group(mpiComm, &mpiGroup);

   // Create list of code identities
   int *myCodes = new int[numProcs] ;
   int *allCodes = new int[numProcs] ;

   MPI_Allgather(&myCode, 1, MPI_INT, allCodes, 1, MPI_INT, mpiComm);

   // Make a list of ALE3D procs from the list of all procs
   int numMyProcs = 0;
   bool otherGroupNeeded = false;
   for (int i=0; i<numProcs; ++i) {
      if (allCodes[i]==myCode) {
         myCodes[numMyProcs] = i;
         ++numMyProcs;
      }
      else {
         otherGroupNeeded = true;
      }
   }
   if (otherGroupNeeded) {
      // Create the communicators
      MPI_Group_excl(mpiGroup, numMyProcs, myCodes, &otherGroup);
      MPI_Comm_create(mpiComm, otherGroup, otherComm);
   }

   MPI_Group_incl(mpiGroup, numMyProcs, myCodes, &myGroup);
   MPI_Comm_create(mpiComm, myGroup, myComm);

   delete [] allCodes ;
   delete [] myCodes ;
}
#endif

} /* namespace geosx */

#endif
