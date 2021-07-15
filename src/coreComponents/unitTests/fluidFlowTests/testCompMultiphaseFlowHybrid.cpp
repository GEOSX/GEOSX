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

#include "constitutive/fluid/MultiFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "mainInterface/initialization.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

using namespace geosx;
using namespace geosx::dataRepository;
using namespace geosx::constitutive;
using namespace geosx::testing;

CommandLineOptions g_commandLineOptions;

char const * xmlInput =
  "<Problem>\n"
  "  <Solvers gravityVector=\"0.0, 0.0, -9.81\">\n"
  "    <CompositionalMultiphaseHybridFVM name=\"compflow\"\n"
  "                                 logLevel=\"0\"\n"
  "                                 discretization=\"fluidHM\"\n"
  "                                 targetRegions=\"{Region}\"\n"
  "                                 fluidNames=\"{fluid1}\"\n"
  "                                 solidNames=\"{rock}\"\n"
  "                                 permeabilityNames=\"{rockPerm}\"\n"
  "                                 relPermNames=\"{relperm}\"\n"
  "                                 temperature=\"297.15\"\n"
  "                                 useMass=\"1\">\n"
  "                                 \n"
  "      <NonlinearSolverParameters newtonTol=\"1.0e-6\"\n"
  "                                 newtonMaxIter=\"2\"/>\n"
  "      <LinearSolverParameters solverType=\"gmres\"\n"
  "                              krylovTol=\"1.0e-10\"/>\n"
  "    </CompositionalMultiphaseHybridFVM>\n"
  "  </Solvers>\n"
  "  <Mesh>\n"
  "    <InternalMesh name=\"mesh1\"\n"
  "                  elementTypes=\"{C3D8}\" \n"
  "                  xCoords=\"{0, 1}\"\n"
  "                  yCoords=\"{0, 1}\"\n"
  "                  zCoords=\"{0, 10}\"\n"
  "                  nx=\"{1}\"\n"
  "                  ny=\"{1}\"\n"
  "                  nz=\"{4}\"\n"
  "                  cellBlockNames=\"{cb1}\"/>\n"
  "  </Mesh>\n"
  "  <NumericalMethods>\n"
  "    <FiniteVolume>\n"
  "      <HybridMimeticDiscretization name=\"fluidHM\"\n"
  "                                   innerProductType=\"beiraoDaVeigaLipnikovManzini\"\n"
  "                                   coefficientName=\"permeability\"/>\n"
  "    </FiniteVolume>\n"
  "  </NumericalMethods>\n"
  "  <ElementRegions>\n"
  "    <CellElementRegion name=\"Region\" cellBlocks=\"{cb1}\" materialList=\"{fluid1, rock, relperm, rockPerm}\" />\n"
  "  </ElementRegions>\n"
  "  <Constitutive>\n"
  "    <CompositionalMultiphaseFluid name=\"fluid1\"\n"
  "                                  phaseNames=\"{oil, gas}\"\n"
  "                                  equationsOfState=\"{PR, PR}\"\n"
  "                                  componentNames=\"{N2, C10, C20, H2O}\"\n"
  "                                  componentCriticalPressure=\"{34e5, 25.3e5, 14.6e5, 220.5e5}\"\n"
  "                                  componentCriticalTemperature=\"{126.2, 622.0, 782.0, 647.0}\"\n"
  "                                  componentAcentricFactor=\"{0.04, 0.443, 0.816, 0.344}\"\n"
  "                                  componentMolarWeight=\"{28e-3, 134e-3, 275e-3, 18e-3}\"\n"
  "                                  componentVolumeShift=\"{0, 0, 0, 0}\"\n"
  "                                  componentBinaryCoeff=\"{ {0, 0, 0, 0},\n"
  "                                                          {0, 0, 0, 0},\n"
  "                                                          {0, 0, 0, 0},\n"
  "                                                          {0, 0, 0, 0} }\"/>\n"
  "    <PoreVolumeCompressibleSolid name=\"rock\"\n"
  "                                 referencePressure=\"0.0\"\n"
  "                                 compressibility=\"1e-9\"/>\n"
  "    <BrooksCoreyRelativePermeability name=\"relperm\"\n"
  "                                     phaseNames=\"{oil, gas}\"\n"
  "                                     phaseMinVolumeFraction=\"{0.1, 0.15}\"\n"
  "                                     phaseRelPermExponent=\"{2.0, 2.0}\"\n"
  "                                     phaseRelPermMaxValue=\"{0.8, 0.9}\"/>\n"
  "  <ConstantPermeability name=\"rockPerm\"\n"
  "                      permeabilityComponents=\"{2.0e-16, 2.0e-16, 2.0e-16}\"/> \n"
  "  </Constitutive>\n"
  "  <FieldSpecifications>\n"
  "    <FieldSpecification name=\"referencePorosity\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region/cb1\"\n"
  "               fieldName=\"referencePorosity\"\n"
  "               scale=\"0.05\"/>\n"
  "    <FieldSpecification name=\"initialPressure\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region/cb1\"\n"
  "               fieldName=\"pressure\"\n"
  "               functionName=\"initialPressureFunc\"\n"
  "               scale=\"5e6\"/>\n"
  "    <FieldSpecification name=\"initialFacePressure\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"faceManager\"\n"
  "               fieldName=\"facePressure\"\n"
  "               functionName=\"initialFacePressureFunc\"\n"
  "               scale=\"5e6\"/>\n"
  "    <FieldSpecification name=\"initialComposition_N2\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"0\"\n"
  "               scale=\"0.099\"/>\n"
  "    <FieldSpecification name=\"initialComposition_C10\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"1\"\n"
  "               scale=\"0.3\"/>\n"
  "    <FieldSpecification name=\"initialComposition_C20\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"2\"\n"
  "               scale=\"0.6\"/>\n"
  "    <FieldSpecification name=\"initialComposition_H20\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"3\"\n"
  "               scale=\"0.001\"/>\n"
  "  </FieldSpecifications>\n"
  "  <Functions>\n"
  "    <TableFunction name=\"initialPressureFunc\"\n"
  "                   inputVarNames=\"{elementCenter}\"\n"
  "                   coordinates=\"{0.0, 2.0, 4.0, 6.0, 8.0, 10.0 }\"\n"
  "                   values=\"{ 1.0, 0.5, 0.2, 3.0, 2.1, 1.0 }\"/>\n"
  "    <TableFunction name=\"initialFacePressureFunc\"\n"
  "                   inputVarNames=\"{faceCenter}\"\n"
  "                   coordinates=\"{0.0, 2.0, 4.0, 6.0, 8.0, 10.0 }\"\n"
  "                   values=\"{ 2.0, 0.1, 0.2, 2.1, 2.1, 0.1 }\"/>\n"
  "  </Functions>"
  "</Problem>";

template< typename LAMBDA >
void testNumericalJacobian( CompositionalMultiphaseHybridFVM & solver,
                            DomainPartition & domain,
                            real64 const perturbParameter,
                            real64 const relTol,
                            LAMBDA assembleFunction )
{
  localIndex const NC = solver.numFluidComponents();

  CRSMatrix< real64, globalIndex > const & jacobian = solver.getLocalMatrix();
  array1d< real64 > const & residual = solver.getLocalRhs();
  DofManager const & dofManager = solver.getDofManager();

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

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

  string const elemDofKey = dofManager.getKey( CompositionalMultiphaseHybridFVM::viewKeyStruct::elemDofFieldString() );

  solver.forTargetSubRegions( mesh, [&]( localIndex const,
                                         ElementSubRegionBase & subRegion )
  {
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
    arrayView1d< globalIndex const > const & elemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( elemDofKey );

    arrayView1d< real64 const > const & pres =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseHybridFVM::viewKeyStruct::pressureString() );
    pres.move( LvArray::MemorySpace::host, false );

    arrayView1d< real64 > const & dPres =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseHybridFVM::viewKeyStruct::deltaPressureString() );

    arrayView2d< real64 const, compflow::USD_COMP > const & compDens =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( CompositionalMultiphaseHybridFVM::viewKeyStruct::globalCompDensityString() );
    compDens.move( LvArray::MemorySpace::host, false );

    arrayView2d< real64, compflow::USD_COMP > const & dCompDens =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( CompositionalMultiphaseHybridFVM::viewKeyStruct::deltaGlobalCompDensityString() );

    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      if( elemGhostRank[ei] >= 0 )
      {
        continue;
      }

      real64 totalDensity = 0.0;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        totalDensity += compDens[ei][ic];
      }

      {
        solver.resetStateToBeginningOfStep( domain );

        real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
        dPres.move( LvArray::MemorySpace::host, true );
        dPres[ei] = dP;

        solver.forTargetSubRegions( mesh, [&]( localIndex const targetIndex2,
                                               ElementSubRegionBase & subRegion2 )
        {
          solver.updateState( subRegion2, targetIndex2 );
        } );

        residual.zero();
        jacobian.zero();
        assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

        fillNumericalJacobian( residual.toViewConst(),
                               residualOrig.toViewConst(),
                               elemDofNumber[ei],
                               dP,
                               jacobianFD.toViewConstSizes() );
      }

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        solver.resetStateToBeginningOfStep( domain );

        real64 const dRho = perturbParameter * totalDensity;
        dCompDens.move( LvArray::MemorySpace::host, true );
        dCompDens[ei][jc] = dRho;

        solver.forTargetSubRegions( mesh, [&]( localIndex const targetIndex2,
                                               ElementSubRegionBase & subRegion2 )
        {
          solver.updateState( subRegion2, targetIndex2 );
        } );

        residual.zero();
        jacobian.zero();
        assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

        fillNumericalJacobian( residual.toViewConst(),
                               residualOrig.toViewConst(),
                               elemDofNumber[ei] + jc + 1,
                               dRho,
                               jacobianFD.toViewConstSizes() );
      }
    }
  } );

  FaceManager & faceManager = mesh.getFaceManager();

  // get the face-based pressure
  arrayView1d< real64 > const & facePres =
    faceManager.getReference< array1d< real64 > >( CompositionalMultiphaseHybridFVM::viewKeyStruct::facePressureString() );
  facePres.move( LvArray::MemorySpace::host, false );
  arrayView1d< real64 > const & dFacePres =
    faceManager.getReference< array1d< real64 > >( CompositionalMultiphaseHybridFVM::viewKeyStruct::deltaFacePressureString() );

  string const faceDofKey = dofManager.getKey( CompositionalMultiphaseHybridFVM::viewKeyStruct::faceDofFieldString() );

  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager.getReference< array1d< globalIndex > >( faceDofKey );
  faceDofNumber.move( LvArray::MemorySpace::host );

  arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
  faceGhostRank.move( LvArray::MemorySpace::host );

  for( localIndex iface = 0; iface < faceManager.size(); ++iface )
  {
    if( faceGhostRank[iface] >= 0 )
    {
      continue;
    }

    solver.resetStateToBeginningOfStep( domain );

    real64 const dFP = perturbParameter * ( facePres[iface] + perturbParameter );
    dFacePres.move( LvArray::MemorySpace::host, true );
    dFacePres[iface] = dFP;

    residual.zero();
    jacobian.zero();
    assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

    fillNumericalJacobian( residual.toViewConst(),
                           residualOrig.toViewConst(),
                           faceDofNumber[iface],
                           dFP,
                           jacobianFD.toViewConstSizes() );
  }

  // assemble the analytical jacobian
  solver.resetStateToBeginningOfStep( domain );

  residual.zero();
  jacobian.zero();
  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

  compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol );
}

class CompositionalMultiphaseHybridFlowTest : public ::testing::Test
{
public:

  CompositionalMultiphaseHybridFlowTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< CompositionalMultiphaseHybridFVM >( "compflow" );

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
  CompositionalMultiphaseHybridFVM * solver;
};

real64 constexpr CompositionalMultiphaseHybridFlowTest::time;
real64 constexpr CompositionalMultiphaseHybridFlowTest::dt;
real64 constexpr CompositionalMultiphaseHybridFlowTest::eps;


TEST_F( CompositionalMultiphaseHybridFlowTest, jacobianNumericalCheck_flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 5e-3; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleFluxTerms( dt, domain, solver->getDofManager(), localMatrix, localRhs );
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
