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


#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

#include "common/DataTypes.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/multiphysics/ReservoirSolverBase.hpp"
#include "physicsSolvers/multiphysics/SinglePhaseReservoir.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVM.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWell.hpp"

using namespace geosx;
using namespace geosx::dataRepository;
using namespace geosx::constitutive;
using namespace geosx::testing;

CommandLineOptions g_commandLineOptions;

char const * xmlInput =
  "<Problem>\n"
  "  <Solvers gravityVector=\"0.0, 0.0, -9.81\">\n"
  "    <SinglePhaseReservoir name=\"reservoirSystem\"\n"
  "               flowSolverName=\"singlePhaseFlow\"\n"
  "               wellSolverName=\"singlePhaseWell\"\n"
  "               logLevel=\"1\"\n"
  "               targetRegions=\"{Region1,wellRegion1,wellRegion2}\">\n"
  "      <NonlinearSolverParameters newtonMaxIter=\"40\"/>\n"
  "      <LinearSolverParameters solverType=\"direct\"\n"
  "                              logLevel=\"2\"/>\n"
  "    </SinglePhaseReservoir>\n"
  "    <SinglePhaseFVM name=\"singlePhaseFlow\"\n"
  "                             logLevel=\"1\"\n"
  "                             discretization=\"singlePhaseTPFA\"\n"
  "                             fluidNames=\"{water}\"\n"
  "                             solidNames=\"{rock}\"\n"
  "                             permeabilityNames=\"{permeabilityModel}\"\n"
  "                             targetRegions=\"{Region1}\">\n"
  "    </SinglePhaseFVM>\n"
  "    <SinglePhaseWell name=\"singlePhaseWell\"\n"
  "                     logLevel=\"1\"\n"
  "                     fluidNames=\"{water}\"\n"
  "                     targetRegions=\"{wellRegion1,wellRegion2}\">\n"
  "        <WellControls name=\"wellControls1\"\n"
  "                      type=\"producer\"\n"
  "                      referenceElevation=\"2\"\n"
  "                      control=\"BHP\"\n"
  "                      targetBHP=\"5e5\"\n"
  "                      targetTotalRate=\"1e-3\"/>\n"
  "        <WellControls name=\"wellControls2\"\n"
  "                      type=\"injector\"\n"
  "                      referenceElevation=\"2\"\n"
  "                      control=\"totalVolRate\" \n"
  "                      targetBHP=\"2e7\"\n"
  "                      targetTotalRate=\"1e-4\"/>\n"
  "    </SinglePhaseWell>\n"
  "  </Solvers>\n"
  "  <Mesh>\n"
  "    <InternalMesh name=\"mesh1\"\n"
  "                  elementTypes=\"{C3D8}\" \n"
  "                  xCoords=\"{0, 5}\"\n"
  "                  yCoords=\"{0, 1}\"\n"
  "                  zCoords=\"{0, 1}\"\n"
  "                  nx=\"{3}\"\n"
  "                  ny=\"{1}\"\n"
  "                  nz=\"{1}\"\n"
  "                  cellBlockNames=\"{cb1}\"/>\n"
  "    <InternalWell name=\"well_producer1\"\n"
  "                  wellRegionName=\"wellRegion1\"\n"
  "                  wellControlsName=\"wellControls1\"\n"
  "                  meshName=\"mesh1\"\n"
  "                  polylineNodeCoords=\"{ {4.5, 0,  2  },\n"
  "                                         {4.5, 0,  0.5} }\"\n"
  "                  polylineSegmentConn=\"{ {0, 1} }\"\n"
  "                  radius=\"0.1\"\n"
  "                  numElementsPerSegment=\"1\">\n"
  "        <Perforation name=\"producer1_perf1\"\n"
  "                     distanceFromHead=\"1.45\"/> \n"
  "    </InternalWell> \n"
  "    <InternalWell name=\"well_injector1\"\n"
  "                  wellRegionName=\"wellRegion2\"\n"
  "                  wellControlsName=\"wellControls2\"\n"
  "                  meshName=\"mesh1\"\n"
  "                  polylineNodeCoords=\"{ {0.5, 0, 2  },\n"
  "                                         {0.5, 0, 0.5} }\"\n"
  "                  polylineSegmentConn=\"{ {0, 1} }\"\n"
  "                  radius=\"0.1\"\n"
  "                  numElementsPerSegment=\"1\">\n"
  "        <Perforation name=\"injector1_perf1\"\n"
  "                     distanceFromHead=\"1.45\"/> \n"
  "    </InternalWell>\n"
  "  </Mesh>\n"
  "  <NumericalMethods>\n"
  "    <FiniteVolume>\n"
  "      <TwoPointFluxApproximation name=\"singlePhaseTPFA\"\n"
  "                                 fieldName=\"pressure\"\n"
  "                                 coefficientName=\"permeability\"/>\n"
  "    </FiniteVolume>\n"
  "  </NumericalMethods>\n"
  "  <ElementRegions>\n"
  "    <CellElementRegion name=\"Region1\"\n"
  "                       cellBlocks=\"{cb1}\"\n"
  "                       materialList=\"{water, rock, permeabilityModel}\"/>\n"
  "    <WellElementRegion name=\"wellRegion1\"\n"
  "                       materialList=\"{water}\"/> \n"
  "    <WellElementRegion name=\"wellRegion2\"\n"
  "                       materialList=\"{water}\"/> \n"
  "  </ElementRegions>\n"
  "  <Constitutive>\n"
  "    <CompressibleSinglePhaseFluid name=\"water\"\n"
  "                                  defaultDensity=\"1000\"\n"
  "                                  defaultViscosity=\"0.001\"\n"
  "                                  referencePressure=\"0.0\"\n"
  "                                  referenceDensity=\"1000\"\n"
  "                                  compressibility=\"5e-10\"\n"
  "                                  referenceViscosity=\"0.001\"\n"
  "                                  viscosibility=\"0.0\"/>\n"
  "    <CompressibleRock name=\"rock\"\n"
  "                      referencePressure=\"0.0\"\n"
  "                      defaultReferencePorosity=\"0.05\"\n"
  "                      compressibility=\"1e-9\"/>\n"
  "    <ConstantPermeability name=\"permeabilityModel\"\n"
  "                          permeabilityComponents=\"{ 2e-16, 2e-16, 2e-16}\"/>\n"
  "  </Constitutive>\n"
  "  <FieldSpecifications>\n"
  "    <FieldSpecification name=\"initialPressure\"\n"
  "                        initialCondition=\"1\"\n"
  "                        setNames=\"{all}\"\n"
  "                        objectPath=\"ElementRegions/Region1/cb1\"\n"
  "                        fieldName=\"pressure\"\n"
  "                        scale=\"5e6\"/>\n"
  "  </FieldSpecifications>\n"
  "</Problem>";

template< typename LAMBDA >
void testNumericalJacobian( SinglePhaseReservoir & solver,
                            DomainPartition & domain,
                            real64 const perturbParameter,
                            real64 const relTol,
                            LAMBDA && assembleFunction )
{
  SinglePhaseWell & wellSolver = dynamicCast< SinglePhaseWell & >( *solver.getWellSolver() );
  SinglePhaseFVM< SinglePhaseBase > & flowSolver = dynamicCast< SinglePhaseFVM< SinglePhaseBase > & >( *solver.getFlowSolver() );

  CRSMatrix< real64, globalIndex > const & jacobian = solver.getLocalMatrix();
  array1d< real64 > const & residual = solver.getLocalRhs();
  DofManager const & dofManager = solver.getDofManager();

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elemManager = mesh.getElemManager();

  // assemble the analytical residual
  solver.resetStateToBeginningOfStep( domain );

  residual.zero();
  jacobian.zero();

  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );
  residual.move( LvArray::MemorySpace::host, false );

  // copy the analytical residual
  array1d< real64 > residualOrig( residual );

  // create the numerical jacobian
  jacobian.move( LvArray::MemorySpace::host );
  CRSMatrix< real64, globalIndex > jacobianFD( jacobian );
  jacobianFD.zero();

  string const & resDofKey  = dofManager.getKey( wellSolver.resElementDofName() );
  string const & wellDofKey = dofManager.getKey( wellSolver.wellElementDofName() );

  // at this point we start assembling the finite-difference block by block

  ////////////////////////////////////////////////
  // Step 1) Compute the terms in J_RR and J_WR //
  ////////////////////////////////////////////////

  for( localIndex er = 0; er < elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    elemRegion.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      // get the dof numbers and ghosting information
      arrayView1d< globalIndex const > const & dofNumber =
        subRegion.getReference< array1d< globalIndex > >( resDofKey );

      // get the primary variables on reservoir elements
      arrayView1d< real64 const > const & pres =
        subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::pressureString() );
      arrayView1d< real64 > const & dPres =
        subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::deltaPressureString() );
      pres.move( LvArray::MemorySpace::host, false );

      // a) compute all the derivatives wrt to the pressure in RESERVOIR elem ei
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        {
          solver.resetStateToBeginningOfStep( domain );

          // here is the perturbation in the pressure of the element
          real64 const dP = perturbParameter * (pres[ei] + perturbParameter);
          dPres.move( LvArray::MemorySpace::host, true );
          dPres[ei] = dP;

          // after perturbing, update the pressure-dependent quantities in the reservoir
          flowSolver.forTargetSubRegions( mesh, [&]( localIndex const targetIndex2,
                                                     ElementSubRegionBase & subRegion2 )
          {
            flowSolver.updateFluidState( subRegion2, targetIndex2 );
          } );
          wellSolver.forTargetSubRegions< WellElementSubRegion >( mesh, [&]( localIndex const targetIndex3,
                                                                             WellElementSubRegion & subRegion3 )
          {
            wellSolver.updateSubRegionState( subRegion3, targetIndex3 );
          } );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 dofNumber[ei],
                                 dP,
                                 jacobianFD.toViewConstSizes() );
        }
      }
    } );
  }

  /////////////////////////////////////////////////
  // Step 2) Compute the terms in J_RW and J_WW //
  /////////////////////////////////////////////////

  // loop over the wells
  wellSolver.forTargetSubRegions< WellElementSubRegion >( mesh, [&]( localIndex const targetIndex,
                                                                     WellElementSubRegion & subRegion )
  {
    // get the degrees of freedom and ghosting information
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    // get the primary variables on well elements
    arrayView1d< real64 const > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( SinglePhaseWell::viewKeyStruct::pressureString() );
    arrayView1d< real64 > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( SinglePhaseWell::viewKeyStruct::deltaPressureString() );
    wellElemPressure.move( LvArray::MemorySpace::host, false );

    arrayView1d< real64 const > const & connRate =
      subRegion.getReference< array1d< real64 > >( SinglePhaseWell::viewKeyStruct::connRateString() );
    arrayView1d< real64 > const & dConnRate =
      subRegion.getReference< array1d< real64 > >( SinglePhaseWell::viewKeyStruct::deltaConnRateString() );
    connRate.move( LvArray::MemorySpace::host, false );

    // a) compute all the derivatives wrt to the pressure in WELL elem iwelem
    for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
    {
      {
        solver.resetStateToBeginningOfStep( domain );

        // here is the perturbation in the pressure of the well element
        real64 const dP = perturbParameter * ( wellElemPressure[iwelem] + perturbParameter );
        dWellElemPressure.move( LvArray::MemorySpace::host, true );
        dWellElemPressure[iwelem] = dP;

        // after perturbing, update the pressure-dependent quantities in the well
        wellSolver.updateSubRegionState( subRegion, targetIndex );

        residual.zero();
        jacobian.zero();
        assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

        // consider mass balance eq lid in RESERVOIR elems and WELL elems
        //      this is computing J_RW and J_WW
        fillNumericalJacobian( residual.toViewConst(),
                               residualOrig.toViewConst(),
                               wellElemDofNumber[iwelem] + SinglePhaseWell::ColOffset::DPRES,
                               dP,
                               jacobianFD.toViewConstSizes() );
      }
    }

    // b) compute all the derivatives wrt to the connection in WELL elem iwelem
    for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
    {
      {
        solver.resetStateToBeginningOfStep( domain );

        // here is the perturbation in the pressure of the well element
        real64 const dRate = perturbParameter * ( connRate[iwelem] + perturbParameter );
        dConnRate.move( LvArray::MemorySpace::host, true );
        dConnRate[iwelem] = dRate;

        // after perturbing, update the rate-dependent quantities in the well (well controls)
        wellSolver.updateSubRegionState( subRegion, targetIndex );

        residual.zero();
        jacobian.zero();
        assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

        // consider mass balance eq lid in RESERVOIR elems and WELL elems
        //      this is computing J_RW and J_WW
        fillNumericalJacobian( residual.toViewConst(),
                               residualOrig.toViewConst(),
                               wellElemDofNumber[iwelem] + SinglePhaseWell::ColOffset::DRATE,
                               dRate,
                               jacobianFD.toViewConstSizes() );
      }
    }
  } );

  // assemble the analytical jacobian
  solver.resetStateToBeginningOfStep( domain );

  residual.zero();
  jacobian.zero();
  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

  compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol );
}

class SinglePhaseReservoirSolverTest : public ::testing::Test
{
public:

  SinglePhaseReservoirSolverTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< SinglePhaseReservoir >( "reservoirSystem" );

    DomainPartition & domain = state.getProblemManager().getDomainPartition();

    solver->setupSystem( domain,
                         solver->getDofManager(),
                         solver->getLocalMatrix(),
                         solver->getLocalRhs(),
                         solver->getLocalSolution() );

    solver->implicitStepSetup( time, dt, domain );
  }

  static real64 constexpr time = 0.0;
  static real64 constexpr dt = 1e4;
  static real64 constexpr eps = std::numeric_limits< real64 >::epsilon();

  GeosxState state;
  SinglePhaseReservoir * solver;
};

real64 constexpr SinglePhaseReservoirSolverTest::time;
real64 constexpr SinglePhaseReservoirSolverTest::dt;
real64 constexpr SinglePhaseReservoirSolverTest::eps;

TEST_F( SinglePhaseReservoirSolverTest, jacobianNumericalCheck_Perforation )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleCouplingTerms( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

TEST_F( SinglePhaseReservoirSolverTest, jacobianNumericalCheck_Flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->getWellSolver()->assembleFluxTerms( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

TEST_F( SinglePhaseReservoirSolverTest, jacobianNumericalCheck_PressureRel )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->getWellSolver()->formPressureRelations( domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
