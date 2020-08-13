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

#include "physicsSolvers/fluidFlow/unitTests/testCompFlowUtils.hpp"

#include "common/DataTypes.hpp"
#include "managers/initialization.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/DomainPartition.hpp"
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
  "                             targetRegions=\"{Region1}\">\n"
  "    </SinglePhaseFVM>\n"
  "    <SinglePhaseWell name=\"singlePhaseWell\"\n"
  "                     logLevel=\"1\"\n"
  "                     fluidNames=\"{water}\"\n"
  "                     targetRegions=\"{wellRegion1,wellRegion2}\">\n"
  "        <WellControls name=\"wellControls1\"\n"
  "                      type=\"producer\"\n"
  "                      control=\"BHP\"\n"
  "                      targetBHP=\"5e5\"\n"
  "                      targetRate=\"1e-3\"/>\n"
  "        <WellControls name=\"wellControls2\"\n"
  "                      type=\"injector\"\n"
  "                      control=\"liquidRate\" \n"
  "                      targetBHP=\"2e7\"\n"
  "                      targetRate=\"1e-4\"/>\n"
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
  "                       materialList=\"{water, rock}\"/>\n"
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
  "    <PoreVolumeCompressibleSolid name=\"rock\"\n"
  "                                 referencePressure=\"0.0\"\n"
  "                                 compressibility=\"1e-9\"/>\n"
  "  </Constitutive>\n"
  "  <FieldSpecifications>\n"
  "    <FieldSpecification name=\"permx\"\n"
  "                        component=\"0\"\n"
  "                        initialCondition=\"1\"  \n"
  "                        setNames=\"{all}\"\n"
  "                        objectPath=\"ElementRegions/Region1/cb1\"\n"
  "                        fieldName=\"permeability\"\n"
  "                        scale=\"2.0e-16\"/>\n"
  "    <FieldSpecification name=\"permy\"\n"
  "                        component=\"1\"\n"
  "                        initialCondition=\"1\"\n"
  "                        setNames=\"{all}\"\n"
  "                        objectPath=\"ElementRegions/Region1/cb1\"\n"
  "                        fieldName=\"permeability\"\n"
  "                        scale=\"2.0e-16\"/>\n"
  "    <FieldSpecification name=\"permz\"\n"
  "                        component=\"2\"\n"
  "                        initialCondition=\"1\"\n"
  "                        setNames=\"{all}\"\n"
  "                        objectPath=\"ElementRegions/Region1/cb1\"\n"
  "                        fieldName=\"permeability\"\n"
  "                        scale=\"2.0e-16\"/>\n"
  "    <FieldSpecification name=\"referencePorosity\"\n"
  "                        initialCondition=\"1\"\n"
  "                        setNames=\"{all}\"\n"
  "                        objectPath=\"ElementRegions/Region1/cb1\"\n"
  "                        fieldName=\"referencePorosity\"\n"
  "                        scale=\"0.05\"/>\n"
  "    <FieldSpecification name=\"initialPressure\"\n"
  "                        initialCondition=\"1\"\n"
  "                        setNames=\"{all}\"\n"
  "                        objectPath=\"ElementRegions/Region1/cb1\"\n"
  "                        fieldName=\"pressure\"\n"
  "                        scale=\"5e6\"/>\n"
  "  </FieldSpecifications>\n"
  "</Problem>";

template< typename LAMBDA >
void
testNumericalJacobian( SinglePhaseReservoir & solver,
                       DomainPartition & domain,
                       real64 const perturbParameter,
                       real64 const relTol,
                       LAMBDA && assembleFunction )
{
  SinglePhaseWell & wellSolver =
    *solver.GetWellSolver()->group_cast< SinglePhaseWell * >();
  SinglePhaseFVM< SinglePhaseBase > & flowSolver =
    *solver.GetFlowSolver()->group_cast< SinglePhaseFVM< SinglePhaseBase > * >();

  CRSMatrix< real64, globalIndex > const & jacobian = solver.getLocalMatrix();
  array1d< real64 > const & residual = solver.getLocalRhs();
  DofManager const & dofManager = solver.getDofManager();

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager & elemManager = *mesh.getElemManager();

  // assemble the analytical residual
  solver.ResetStateToBeginningOfStep( domain );

  residual.setValues< parallelDevicePolicy<> >( 0.0 );
  jacobian.setValues< parallelDevicePolicy<> >( 0.0 );

  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );
  residual.move( LvArray::MemorySpace::CPU, false );

  // copy the analytical residual
  array1d< real64 > residualOrig( residual );

  // create the numerical jacobian
  CRSMatrix< real64, globalIndex > jacobianFD( jacobian );
  jacobianFD.setValues< parallelDevicePolicy<> >( 0.0 );

  string const & resDofKey = dofManager.getKey( wellSolver.ResElementDofName() );
  string const & wellDofKey = dofManager.getKey( wellSolver.WellElementDofName() );

  // at this point we start assembling the finite-difference block by block

  ////////////////////////////////////////////////
  // Step 1) Compute the terms in J_RR and J_WR //
  ////////////////////////////////////////////////

  for( localIndex er = 0; er < elemManager.numRegions(); ++er )
  {
    ElementRegionBase * const elemRegion = elemManager.GetRegion( er );
    elemRegion->forElementSubRegions< CellElementSubRegion >(
      [&]( CellElementSubRegion & subRegion ) {
        // get the dof numbers and ghosting information
        arrayView1d< globalIndex const > const & dofNumber =
          subRegion.getReference< array1d< globalIndex > >( resDofKey );

        // get the primary variables on reservoir elements
        arrayView1d< real64 const > const & pres =
          subRegion.getReference< array1d< real64 > >(
            FlowSolverBase::viewKeyStruct::pressureString );
        arrayView1d< real64 > const & dPres =
          subRegion.getReference< array1d< real64 > >(
            FlowSolverBase::viewKeyStruct::deltaPressureString );
        pres.move( LvArray::MemorySpace::CPU, false );

        // a) compute all the derivatives wrt to the pressure in RESERVOIR elem ei
        for( localIndex ei = 0; ei < subRegion.size(); ++ei )
        {
          {
            solver.ResetStateToBeginningOfStep( domain );

            // here is the perturbation in the pressure of the element
            real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
            dPres.move( LvArray::MemorySpace::CPU, true );
            dPres[ei] = dP;

            // after perturbing, update the pressure-dependent quantities in the reservoir
            flowSolver.forTargetSubRegions(
              mesh,
              [&]( localIndex const targetIndex2, ElementSubRegionBase & subRegion2 ) {
                flowSolver.UpdateState( subRegion2, targetIndex2 );
              } );
            wellSolver.forTargetSubRegions< WellElementSubRegion >(
              mesh,
              [&]( localIndex const targetIndex3, WellElementSubRegion & subRegion3 ) {
                wellSolver.UpdateState( subRegion3, targetIndex3 );
              } );

            residual.setValues< parallelDevicePolicy<> >( 0.0 );
            jacobian.setValues< parallelDevicePolicy<> >( 0.0 );
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
  wellSolver.forTargetSubRegions< WellElementSubRegion >(
    mesh,
    [&]( localIndex const targetIndex, WellElementSubRegion & subRegion ) {
      // get the degrees of freedom and ghosting information
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );

      // get the primary variables on well elements
      arrayView1d< real64 const > const & wellElemPressure =
        subRegion.getReference< array1d< real64 > >(
          SinglePhaseWell::viewKeyStruct::pressureString );
      arrayView1d< real64 > const & dWellElemPressure =
        subRegion.getReference< array1d< real64 > >(
          SinglePhaseWell::viewKeyStruct::deltaPressureString );
      wellElemPressure.move( LvArray::MemorySpace::CPU, false );

      arrayView1d< real64 const > const & connRate =
        subRegion.getReference< array1d< real64 > >(
          SinglePhaseWell::viewKeyStruct::connRateString );
      arrayView1d< real64 > const & dConnRate =
        subRegion.getReference< array1d< real64 > >(
          SinglePhaseWell::viewKeyStruct::deltaConnRateString );
      connRate.move( LvArray::MemorySpace::CPU, false );

      // a) compute all the derivatives wrt to the pressure in WELL elem iwelem
      for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
      {
        {
          solver.ResetStateToBeginningOfStep( domain );

          // here is the perturbation in the pressure of the well element
          real64 const dP =
            perturbParameter * ( wellElemPressure[iwelem] + perturbParameter );
          dWellElemPressure.move( LvArray::MemorySpace::CPU, true );
          dWellElemPressure[iwelem] = dP;

          // after perturbing, update the pressure-dependent quantities in the well
          wellSolver.UpdateState( subRegion, targetIndex );

          residual.setValues< parallelDevicePolicy<> >( 0.0 );
          jacobian.setValues< parallelDevicePolicy<> >( 0.0 );
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          // consider mass balance eq lid in RESERVOIR elems and WELL elems
          //      this is computing J_RW and J_WW
          fillNumericalJacobian(
            residual.toViewConst(),
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
          solver.ResetStateToBeginningOfStep( domain );

          // here is the perturbation in the pressure of the well element
          real64 const dRate =
            perturbParameter * ( connRate[iwelem] + perturbParameter );
          dConnRate.move( LvArray::MemorySpace::CPU, true );
          dConnRate[iwelem] = dRate;

          residual.setValues< parallelDevicePolicy<> >( 0.0 );
          jacobian.setValues< parallelDevicePolicy<> >( 0.0 );
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          // consider mass balance eq lid in RESERVOIR elems and WELL elems
          //      this is computing J_RW and J_WW
          fillNumericalJacobian(
            residual.toViewConst(),
            residualOrig.toViewConst(),
            wellElemDofNumber[iwelem] + SinglePhaseWell::ColOffset::DRATE,
            dRate,
            jacobianFD.toViewConstSizes() );
        }
      }
    } );

  // assemble the analytical jacobian
  solver.ResetStateToBeginningOfStep( domain );

  residual.setValues< parallelDevicePolicy<> >( 0.0 );
  jacobian.setValues< parallelDevicePolicy<> >( 0.0 );
  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

  compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol );
}

class SinglePhaseReservoirSolverTest : public ::testing::Test
{
public:
  SinglePhaseReservoirSolverTest() :
    problemManager( std::make_unique< ProblemManager >( "Problem", nullptr ) )
  {}

protected:
  void
  SetUp() override
  {
    setupProblemFromXML( *problemManager, xmlInput );
    solver =
      problemManager->GetPhysicsSolverManager().GetGroup< SinglePhaseReservoir >(
        "reservoirSystem" );

    DomainPartition & domain = *problemManager->getDomainPartition();

    solver->SetupSystem( domain,
                         solver->getDofManager(),
                         solver->getLocalMatrix(),
                         solver->getLocalRhs(),
                         solver->getLocalSolution() );

    solver->ImplicitStepSetup( time, dt, domain );
  }

  static real64 constexpr time = 0.0;
  static real64 constexpr dt = 1e4;
  static real64 constexpr eps = std::numeric_limits< real64 >::epsilon();

  std::unique_ptr< ProblemManager > problemManager;
  SinglePhaseReservoir * solver;
};

real64 constexpr SinglePhaseReservoirSolverTest::time;
real64 constexpr SinglePhaseReservoirSolverTest::dt;
real64 constexpr SinglePhaseReservoirSolverTest::eps;

TEST_F( SinglePhaseReservoirSolverTest, jacobianNumericalCheck_Perforation )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1;  // 10% error margin

  DomainPartition & domain = *problemManager->getDomainPartition();

  testNumericalJacobian(
    *solver,
    domain,
    perturb,
    tol,
    [&]( CRSMatrixView< real64, globalIndex const > const & localMatrix,
         arrayView1d< real64 > const & localRhs ) {
      solver->AssembleCouplingTerms( time,
                                     dt,
                                     domain,
                                     solver->getDofManager(),
                                     localMatrix,
                                     localRhs );
    } );
}

TEST_F( SinglePhaseReservoirSolverTest, jacobianNumericalCheck_Flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1;  // 10% error margin

  DomainPartition & domain = *problemManager->getDomainPartition();

  testNumericalJacobian(
    *solver,
    domain,
    perturb,
    tol,
    [&]( CRSMatrixView< real64, globalIndex const > const & localMatrix,
         arrayView1d< real64 > const & localRhs ) {
      solver->GetWellSolver()->AssembleFluxTerms( time,
                                                  dt,
                                                  domain,
                                                  solver->getDofManager(),
                                                  localMatrix,
                                                  localRhs );
    } );
}

TEST_F( SinglePhaseReservoirSolverTest, jacobianNumericalCheck_PressureRel )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1;  // 10% error margin

  DomainPartition & domain = *problemManager->getDomainPartition();

  testNumericalJacobian(
    *solver,
    domain,
    perturb,
    tol,
    [&]( CRSMatrixView< real64, globalIndex const > const & localMatrix,
         arrayView1d< real64 > const & localRhs ) {
      solver->GetWellSolver()->FormPressureRelations( domain,
                                                      solver->getDofManager(),
                                                      localMatrix,
                                                      localRhs );
    } );
}

int
main( int argc, char ** argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
