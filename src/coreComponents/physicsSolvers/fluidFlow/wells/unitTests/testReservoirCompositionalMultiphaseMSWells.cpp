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
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/multiphysics/ReservoirSolverBase.hpp"
#include "physicsSolvers/multiphysics/CompositionalMultiphaseReservoir.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFlow.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"

using namespace geosx;
using namespace geosx::dataRepository;
using namespace geosx::constitutive;
using namespace geosx::testing;

char const * xmlInput =
  "<Problem>\n"
  "  <Solvers gravityVector=\"0.0, 0.0, -9.81\">\n"
  "    <CompositionalMultiphaseReservoir name=\"reservoirSystem\"\n"
  "               flowSolverName=\"compositionalMultiphaseFlow\"\n"
  "               wellSolverName=\"compositionalMultiphaseWell\"\n"
  "               logLevel=\"1\"\n"
  "               targetRegions=\"{Region1,wellRegion1,wellRegion2}\">\n"
  "      <NonlinearSolverParameters newtonMaxIter=\"40\"/>\n"
  "      <LinearSolverParameters solverType=\"direct\"\n"
  "                              logLevel=\"2\"/>\n"
  "    </CompositionalMultiphaseReservoir>\n"
  "    <CompositionalMultiphaseFlow name=\"compositionalMultiphaseFlow\"\n"
  "                                 logLevel=\"1\"\n"
  "                                 discretization=\"fluidTPFA\"\n"
  "                                 targetRegions=\"{Region1}\"\n"
  "                                 fluidNames=\"{fluid1}\"\n"
  "                                 solidNames=\"{rock}\"\n"
  "                                 relPermNames=\"{relperm}\"\n"
  "                                 temperature=\"297.15\"\n"
  "                                 useMass=\"0\">\n"
  "    </CompositionalMultiphaseFlow>\n"
  "    <CompositionalMultiphaseWell name=\"compositionalMultiphaseWell\"\n"
  "                                 logLevel=\"1\"\n"
  "                                 "
  "targetRegions=\"{wellRegion1,wellRegion2}\"\n"
  "                                 fluidNames=\"{fluid1}\"\n"
  "                                 relPermNames=\"{relperm}\"\n"
  "                                 wellTemperature=\"297.15\"\n"
  "                                 useMass=\"0\">\n"
  "        <WellControls name=\"wellControls1\"\n"
  "                      type=\"producer\"\n"
  "                      control=\"BHP\"\n"
  "                      targetBHP=\"4e6\"\n"
  "                      targetRate=\"1\"/>\n"
  "        <WellControls name=\"wellControls2\"\n"
  "                      type=\"injector\"\n"
  "                      control=\"liquidRate\" \n"
  "                      targetBHP=\"2e7\"\n"
  "                      targetRate=\"1e-5\" \n"
  "                      injectionStream=\"{0.1, 0.1, 0.1, 0.7}\"/>\n"
  "    </CompositionalMultiphaseWell>\n"
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
  "    </InternalWell>\n"
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
  "      <TwoPointFluxApproximation name=\"fluidTPFA\"\n"
  "                                 fieldName=\"pressure\"\n"
  "                                 coefficientName=\"permeability\"/>\n"
  "    </FiniteVolume>\n"
  "  </NumericalMethods>\n"
  "  <ElementRegions>\n"
  "    <CellElementRegion name=\"Region1\"\n"
  "                       cellBlocks=\"{cb1}\"\n"
  "                       materialList=\"{fluid1, rock, relperm}\"/>\n"
  "    <WellElementRegion name=\"wellRegion1\"\n"
  "                       materialList=\"{fluid1, relperm}\"/> \n"
  "    <WellElementRegion name=\"wellRegion2\"\n"
  "                       materialList=\"{fluid1, relperm}\"/> \n"
  "  </ElementRegions>\n"
  "  <Constitutive>\n"
  "    <CompositionalMultiphaseFluid name=\"fluid1\"\n"
  "                                  phaseNames=\"{oil, gas}\"\n"
  "                                  equationsOfState=\"{PR, PR}\"\n"
  "                                  componentNames=\"{N2, C10, C20, H2O}\"\n"
  "                                  componentCriticalPressure=\"{34e5, "
  "25.3e5, 14.6e5, 220.5e5}\"\n"
  "                                  componentCriticalTemperature=\"{126.2, "
  "622.0, 782.0, 647.0}\"\n"
  "                                  componentAcentricFactor=\"{0.04, 0.443, "
  "0.816, 0.344}\"\n"
  "                                  componentMolarWeight=\"{28e-3, 134e-3, "
  "275e-3, 18e-3}\"\n"
  "                                  componentVolumeShift=\"{0, 0, 0, 0}\"\n"
  "                                  componentBinaryCoeff=\"{ {0, 0, 0, 0},\n"
  "                                                          {0, 0, 0, 0},\n"
  "                                                          {0, 0, 0, 0},\n"
  "                                                          {0, 0, 0, 0} "
  "}\"/>\n"
  "    <PoreVolumeCompressibleSolid name=\"rock\"\n"
  "                                 referencePressure=\"0.0\"\n"
  "                                 compressibility=\"1e-9\"/>\n"
  "    <BrooksCoreyRelativePermeability name=\"relperm\"\n"
  "                                     phaseNames=\"{oil, gas}\"\n"
  "                                     phaseMinVolumeFraction=\"{0.1, "
  "0.15}\"\n"
  "                                     phaseRelPermExponent=\"{2.0, 2.0}\"\n"
  "                                     phaseRelPermMaxValue=\"{0.8, 0.9}\"/>\n"
  "  </Constitutive>\n"
  "  <FieldSpecifications>\n"
  "    <FieldSpecification name=\"permx\"\n"
  "               component=\"0\"\n"
  "               initialCondition=\"1\"  \n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region1/cb1\"\n"
  "               fieldName=\"permeability\"\n"
  "               scale=\"2.0e-16\"/>\n"
  "    <FieldSpecification name=\"permy\"\n"
  "               component=\"1\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region1/cb1\"\n"
  "               fieldName=\"permeability\"\n"
  "               scale=\"2.0e-16\"/>\n"
  "    <FieldSpecification name=\"permz\"\n"
  "               component=\"2\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region1/cb1\"\n"
  "               fieldName=\"permeability\"\n"
  "               scale=\"2.0e-16\"/>\n"
  "    <FieldSpecification name=\"referencePorosity\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region1/cb1\"\n"
  "               fieldName=\"referencePorosity\"\n"
  "               scale=\"0.05\"/>\n"
  "    <FieldSpecification name=\"initialPressure\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region1/cb1\"\n"
  "               fieldName=\"pressure\"\n"
  "               scale=\"5e6\"/>\n"
  "    <FieldSpecification name=\"initialComposition_N2\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region1/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"0\"\n"
  "               scale=\"0.099\"/>\n"
  "    <FieldSpecification name=\"initialComposition_C10\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region1/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"1\"\n"
  "               scale=\"0.3\"/>\n"
  "    <FieldSpecification name=\"initialComposition_C20\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region1/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"2\"\n"
  "               scale=\"0.6\"/>\n"
  "    <FieldSpecification name=\"initialComposition_H20\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region1/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"3\"\n"
  "               scale=\"0.001\"/>\n"
  "  </FieldSpecifications>\n"
  "</Problem>";

template< typename LAMBDA >
void
testNumericalJacobian( CompositionalMultiphaseReservoir & solver,
                       DomainPartition & domain,
                       real64 const perturbParameter,
                       real64 const relTol,
                       LAMBDA && assembleFunction )
{
  CompositionalMultiphaseWell & wellSolver =
    *solver.GetWellSolver()->group_cast< CompositionalMultiphaseWell * >();
  CompositionalMultiphaseFlow & flowSolver =
    *solver.GetFlowSolver()->group_cast< CompositionalMultiphaseFlow * >();

  localIndex const NC = flowSolver.numFluidComponents();

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

  string const resDofKey = dofManager.getKey( wellSolver.ResElementDofName() );
  string const wellDofKey = dofManager.getKey( wellSolver.WellElementDofName() );

  // at this point we start assembling the finite-difference block by block

  ////////////////////////////////////////////////
  // Step 1) Compute the terms in J_RR and J_WR //
  ////////////////////////////////////////////////

  for( localIndex er = 0; er < elemManager.numRegions(); ++er )
  {
    ElementRegionBase * const elemRegion = elemManager.GetRegion( er );
    elemRegion->forElementSubRegionsIndex< CellElementSubRegion >(
      [&]( localIndex const, CellElementSubRegion & subRegion ) {
        // get the degrees of freedom and ghosting information
        arrayView1d< globalIndex const > const & dofNumber =
          subRegion.getReference< array1d< globalIndex > >( resDofKey );

        // get the primary variables on the reservoir elements
        arrayView1d< real64 const > const & pres =
          subRegion.getReference< array1d< real64 > >(
            CompositionalMultiphaseFlow::viewKeyStruct::pressureString );
        arrayView1d< real64 > const & dPres =
          subRegion.getReference< array1d< real64 > >(
            CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString );
        pres.move( LvArray::MemorySpace::CPU, false );

        arrayView2d< real64 const > const & compDens =
          subRegion.getReference< array2d< real64 > >(
            CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );
        arrayView2d< real64 > const & dCompDens =
          subRegion.getReference< array2d< real64 > >(
            CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );
        compDens.move( LvArray::MemorySpace::CPU, false );

        // a) compute all the derivatives wrt to the pressure in RESERVOIR elem ei
        for( localIndex ei = 0; ei < subRegion.size(); ++ei )
        {
          real64 totalDensity = 0.0;
          for( localIndex ic = 0; ic < NC; ++ic )
          {
            totalDensity += compDens[ei][ic];
          }

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

          for( localIndex jc = 0; jc < NC; ++jc )
          {
            solver.ResetStateToBeginningOfStep( domain );

            real64 const dRho = perturbParameter * totalDensity;
            dCompDens.move( LvArray::MemorySpace::CPU, true );
            dCompDens[ei][jc] = dRho;

            flowSolver.forTargetSubRegions(
              mesh,
              [&]( localIndex const targetIndex2, ElementSubRegionBase & subRegion2 ) {
                flowSolver.UpdateState( subRegion2, targetIndex2 );
              } );

            residual.setValues< parallelDevicePolicy<> >( 0.0 );
            jacobian.setValues< parallelDevicePolicy<> >( 0.0 );
            assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

            fillNumericalJacobian( residual.toViewConst(),
                                   residualOrig.toViewConst(),
                                   dofNumber[ei] + jc + 1,
                                   dRho,
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
      // get the degrees of freedom, ghosting info and next well elem index
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );

      // get the primary variables on the well elements
      arrayView1d< real64 > const & wellElemPressure =
        subRegion.getReference< array1d< real64 > >(
          CompositionalMultiphaseWell::viewKeyStruct::pressureString );
      arrayView1d< real64 > const & dWellElemPressure =
        subRegion.getReference< array1d< real64 > >(
          CompositionalMultiphaseWell::viewKeyStruct::deltaPressureString );
      wellElemPressure.move( LvArray::MemorySpace::CPU, false );

      arrayView2d< real64 const > const & wellElemCompDens =
        subRegion.getReference< array2d< real64 > >(
          CompositionalMultiphaseWell::viewKeyStruct::globalCompDensityString );
      arrayView2d< real64 > const & dWellElemCompDens =
        subRegion.getReference< array2d< real64 > >(
          CompositionalMultiphaseWell::viewKeyStruct::deltaGlobalCompDensityString );
      wellElemCompDens.move( LvArray::MemorySpace::CPU, false );

      arrayView1d< real64 const > const & connRate =
        subRegion.getReference< array1d< real64 > >(
          CompositionalMultiphaseWell::viewKeyStruct::mixtureConnRateString );
      arrayView1d< real64 > const & dConnRate =
        subRegion.getReference< array1d< real64 > >(
          CompositionalMultiphaseWell::viewKeyStruct::deltaMixtureConnRateString );
      connRate.move( LvArray::MemorySpace::CPU, false );

      // a) compute all the derivatives wrt to the pressure in WELL elem iwelem
      for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
      {
        real64 wellElemTotalDensity = 0.0;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          wellElemTotalDensity += wellElemCompDens[iwelem][ic];
        }

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

          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 wellElemDofNumber[iwelem] +
                                   CompositionalMultiphaseWell::ColOffset::DPRES,
                                 dP,
                                 jacobianFD.toViewConstSizes() );
        }

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          solver.ResetStateToBeginningOfStep( domain );

          real64 const dRho = perturbParameter * wellElemTotalDensity;
          dWellElemCompDens.move( LvArray::MemorySpace::CPU, true );
          dWellElemCompDens[iwelem][jc] = dRho;

          wellSolver.UpdateStateAll( domain );

          residual.setValues< parallelDevicePolicy<> >( 0.0 );
          jacobian.setValues< parallelDevicePolicy<> >( 0.0 );
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          fillNumericalJacobian(
            residual.toViewConst(),
            residualOrig.toViewConst(),
            wellElemDofNumber[iwelem] +
              CompositionalMultiphaseWell::ColOffset::DCOMP + jc,
            dRho,
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

          fillNumericalJacobian(
            residual.toViewConst(),
            residualOrig.toViewConst(),
            wellElemDofNumber[iwelem] +
              CompositionalMultiphaseWell::ColOffset::DCOMP + NC,
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

class CompositionalMultiphaseReservoirSolverTest : public ::testing::Test
{
public:
  CompositionalMultiphaseReservoirSolverTest() :
    problemManager( std::make_unique< ProblemManager >( "Problem", nullptr ) )
  {}

protected:
  void
  SetUp() override
  {
    setupProblemFromXML( *problemManager, xmlInput );
    solver = problemManager->GetPhysicsSolverManager()
               .GetGroup< CompositionalMultiphaseReservoir >( "reservoirSystem" );

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
  CompositionalMultiphaseReservoir * solver;
};

real64 constexpr CompositionalMultiphaseReservoirSolverTest::time;
real64 constexpr CompositionalMultiphaseReservoirSolverTest::dt;
real64 constexpr CompositionalMultiphaseReservoirSolverTest::eps;

#if 0

TEST_F( CompositionalMultiphaseReservoirSolverTest, jacobianNumericalCheck_Perforation )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = *problemManager->getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->AssembleCouplingTerms( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

#endif

TEST_F( CompositionalMultiphaseReservoirSolverTest, jacobianNumericalCheck_Flux )
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

TEST_F( CompositionalMultiphaseReservoirSolverTest,
        jacobianNumericalCheck_VolumeBalance )
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
      solver->GetWellSolver()->AssembleVolumeBalanceTerms( time,
                                                           dt,
                                                           domain,
                                                           solver->getDofManager(),
                                                           localMatrix,
                                                           localRhs );
    } );
}

TEST_F( CompositionalMultiphaseReservoirSolverTest,
        jacobianNumericalCheck_PressureRel )
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
