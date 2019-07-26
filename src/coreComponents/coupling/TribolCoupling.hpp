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


#ifndef SRC_COMPONENTS_CORE_SRC_TRIBOLCOUPLING_HPP_
#define SRC_COMPONENTS_CORE_SRC_TRIBOLCOUPLING_HPP_

#if HAVE_TRIBOLCOUPLING

#include "dataRepository/ManagedGroup.hpp"
#include "mesh/NodeManager.hpp"

namespace geosx
{

/**
 * @class TribolCoupling
 *
 * A class for managing Tribol code coupling.
 */
class TribolCoupling
{
public:
   static void Initialize(dataRepository::ManagedGroup * eventManager, dataRepository::ManagedGroup * domain) ;

   static void ApplyTribolForces(dataRepository::ManagedGroup * domain,
                                 real64 const& time_n, real64 const& dt, const int cycleNumber) ;
   static void CopyPositionsToTribolSourceData(NodeManager const * const nodeManager) ;
   static void CopyVelocitiesToTribolSourceData(NodeManager const * const nodeManager) ;
   static void CopyAccelerationsToTribolSourceData(NodeManager const * const nodeManager) ;
   static void CopyForcesToTribolSourceData(NodeManager const * const nodeManager) ;
   static void CopyAccelerationsFromTribolSourceData(NodeManager * const nodeManager) ;

   static void SyncTermination(int* terminate) ;
   static void SyncTimestep(real64* newDt) ;
   static void Cleanup() ;

#ifdef GEOSX_USE_MPI
   static void InitCommSubset(MPI_Comm const mpiComm, MPI_Comm *myComm, MPI_Comm *otherComm, int myCode) ;
#endif
};

} /* namespace geosx */

#endif

#endif /* SRC_COMPONENTS_CORE_SRC_TRIBOLCOUPLING_HPP_ */
