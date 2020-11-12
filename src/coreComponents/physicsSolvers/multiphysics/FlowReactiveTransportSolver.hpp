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
 * @file FlowReactiveTransportSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_COUPLEDSOLVERS_FLOWREACTIVETRANSPORTSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_COUPLEDSOLVERS_FLOWREACTIVETRANSPORTSOLVER_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{

class ReactiveTransport;
class GeochemicalModel;
class FlowSolverBase;

class FlowReactiveTransportSolver : public SolverBase
{
public:
  FlowReactiveTransportSolver( const std::string & name,
                               Group * const parent );
  ~FlowReactiveTransportSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "FlowReactiveTransport"; }

  virtual void RegisterDataOnMesh( dataRepository::Group * const MeshBodies ) override final;

  virtual real64
  SolverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override;

  virtual void
  ResetStateToBeginningOfStep( DomainPartition & domain ) override;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto reactiveTransportSolverNameString = "reactiveTransportSolverName";
    constexpr static auto flowSolverNameString = "flowSolverName";
    constexpr static auto geochemicalModelNameString = "geochemicalModelName";

  } flowReactiveTransportSolverViewKeys;


  void PreStepUpdate( real64 const & time_n,
                      real64 const & dt,
                      int const cycleNumber,
                      DomainPartition & domain );

  void PostStepUpdate( real64 const & time_n,
                       real64 const & dt,
                       int const cycleNumber,
                       DomainPartition & domain );

protected:

  virtual void PostProcessInput() override final;

private:

  string m_flowSolverName;
  string m_reactiveTransportSolverName;
  string m_geochemicalModelName;

  FlowSolverBase * m_flowSolver;
  ReactiveTransport * m_reactiveTransportSolver;
  GeochemicalModel * m_geochemicalModel;

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_COUPLEDSOLVERS_FLOWREACTIVETRANSPORTSOLVER_HPP_ */
