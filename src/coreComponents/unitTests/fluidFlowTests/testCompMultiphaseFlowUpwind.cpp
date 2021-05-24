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

// Source includes
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mainInterface/initialization.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFlowUpwindHelperKernels.hpp"
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

#include "unitTests/fluidFlowTests/testFlowKernelHelpers.hpp"
// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::CompositionalMultiphaseFlowUpwindHelperKernels;
using UHelpers = geosx::CompositionalMultiphaseFlowUpwindHelperKernels::UpwindHelpers;

using namespace geosx::dataRepository;
using namespace geosx::constitutive;
using namespace geosx::testing;

CommandLineOptions g_commandLineOptions;

// Sphinx start after input XML

char const * xmlInput =
  "<Problem>\n"
  "  <Solvers gravityVector=\"-9.81, 0.0, 0.0\">\n"
  "    <CompositionalMultiphaseFVM name=\"compflow\"\n"
  "                                 logLevel=\"0\"\n"
  "                                 discretization=\"fluidTPFA\"\n"
  "                                 targetRegions=\"{Region2}\"\n"
  "                                 fluidNames=\"{fluid1}\"\n"
  "                                 solidNames=\"{rock}\"\n"
  "                                 relPermNames=\"{relperm}\"\n"
  "                                 capPressureNames=\"{cappressure}\"\n"
  "                                 temperature=\"297.15\"\n"
  "                                 useMass=\"1\">\n"
  "                                 \n"
  "      <NonlinearSolverParameters newtonTol=\"1.0e-6\"\n"
  "                                 newtonMaxIter=\"2\"/>\n"
  "      <LinearSolverParameters solverType=\"gmres\"\n"
  "                              krylovTol=\"1.0e-10\"/>\n"
  "    </CompositionalMultiphaseFVM>\n"
  "  </Solvers>\n"
  "  <Mesh>\n"
  "    <InternalMesh name=\"mesh1\"\n"
  "                  elementTypes=\"{C3D8}\" \n"
  "                  xCoords=\"{0, 3}\"\n"
  "                  yCoords=\"{0, 1}\"\n"
  "                  zCoords=\"{0, 1}\"\n"
  "                  nx=\"{13}\"\n"
  "                  ny=\"{1}\"\n"
  "                  nz=\"{1}\"\n"
  "                  cellBlockNames=\"{cb1}\"/>\n"
  "  </Mesh>\n"
  "  <NumericalMethods>\n"
  "    <FiniteVolume>\n"
  "      <TwoPointFluxApproximation name=\"fluidTPFA\"\n"
  "                                 fieldName=\"pressure\"\n"
  "                                 coefficientName=\"permeability\"/>\n"
  "    </FiniteVolume>\n"
  "  </NumericalMethods>\n"
  "  <ElementRegions>\n"
  "    <CellElementRegion name=\"Region2\" cellBlocks=\"{cb1}\" materialList=\"{fluid1, rock, relperm, cappressure}\" />\n"
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
  "    <BrooksCoreyCapillaryPressure name=\"cappressure\"\n"
  "                                  phaseNames=\"{oil, gas}\"\n"
  "                                  phaseMinVolumeFraction=\"{0.2, 0.05}\"\n"
  "                                  phaseCapPressureExponentInv=\"{4.25, 3.5}\"\n"
  "                                  phaseEntryPressure=\"{0., 1e8}\"\n"
  "                                  capPressureEpsilon=\"0.0\"/> \n"
  "  </Constitutive>\n"
  "  <FieldSpecifications>\n"
  "    <FieldSpecification name=\"permx\"\n"
  "               component=\"0\"\n"
  "               initialCondition=\"1\"  \n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"permeability\"\n"
  "               scale=\"2.0e-16\"/>\n"
  "    <FieldSpecification name=\"permy\"\n"
  "               component=\"1\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"permeability\"\n"
  "               scale=\"2.0e-16\"/>\n"
  "    <FieldSpecification name=\"permz\"\n"
  "               component=\"2\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"permeability\"\n"
  "               scale=\"2.0e-16\"/>\n"
  "    <FieldSpecification name=\"referencePorosity\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"referencePorosity\"\n"
  "               scale=\"0.05\"/>\n"
  "    <FieldSpecification name=\"initialPressure\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"pressure\"\n"
  "               functionName=\"initialPressureFunc\"\n"
  "               scale=\"5e6\"/>\n"
  "    <FieldSpecification name=\"initialComposition_N2\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"0\"\n"
  "               scale=\"0.099\"/>\n"
  "    <FieldSpecification name=\"initialComposition_C10\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"1\"\n"
  "               functionName=\"initialC10\"\n"
  "               scale=\"1\"/>\n"
  "    <FieldSpecification name=\"initialComposition_C20\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"2\"\n"
  "               functionName=\"initialC20\"\n"
  "               scale=\"1\"/>\n"
  "    <FieldSpecification name=\"initialComposition_H20\"\n"
  "               initialCondition=\"1\"\n"
  "               setNames=\"{all}\"\n"
  "               objectPath=\"ElementRegions/Region2/cb1\"\n"
  "               fieldName=\"globalCompFraction\"\n"
  "               component=\"3\"\n"
  "               scale=\"0.001\"/>\n"
  "  </FieldSpecifications>\n"
  "  <Functions>\n"
  "    <TableFunction name=\"initialPressureFunc\"\n"
  "                   inputVarNames=\"{elementCenter}\"\n"
  "                   coordinates=\"{0.0,3.0}\"\n"
  "                   values=\"{1.0, 0.5}\"/>\n"
  "    <TableFunction name=\"initialC10\"\n"
  "                   inputVarNames=\"{elementCenter}\"\n"
  "                   coordinates=\"{0.0, 3.0}\"\n"
  "                   values=\"{0.3, 0.6}\"/>\n"
  "    <TableFunction name=\"initialC20\"\n"
  "                   inputVarNames=\"{elementCenter}\"\n"
  "                   coordinates=\"{0.0, 3.0}\"\n"
  "                   values=\"{0.6, 0.3}\"/>\n"
  "  </Functions>"
  "</Problem>";


// Useful macro
#ifndef DEFINE_CAPB_FIELDS
#define DEFINE_CAP_FIELDS( active ) \
  int const capPressureFlag = active; \
  auto phaseCapPressureView = getElementAccessor< true, 3, real64 >( &cap, \
                                                                     CapillaryPressureBase::viewKeyStruct::phaseCapPressureString(), \
                                                                     MAX_STENCIL, \
                                                                     seri[iconn], \
                                                                     sesri[iconn], \
                                                                     sei[iconn], \
                                                                     1, \
                                                                     NP ); \
  auto dPhaseCapPressure_dPhaseVolFracView = getElementAccessor< true, 4, real64 >( &cap, \
                                                                                    CapillaryPressureBase::viewKeyStruct::dPhaseCapPressure_dPhaseVolFractionString(), \
                                                                                    MAX_STENCIL, \
                                                                                    seri[iconn], \
                                                                                    sesri[iconn], \
                                                                                    sei[iconn], \
                                                                                    1, \
                                                                                    NP, \
                                                                                    NP );
#endif


#ifndef DEFINE_FLUID_FIELDS
#define DEFINE_FLUID_FIELDS() \
  auto phaseMassDensView = getElementAccessor< true, 3, real64 >( &fluid, \
                                                                  MultiFluidBase::viewKeyStruct::phaseMassDensityString(), \
                                                                  MAX_STENCIL, \
                                                                  seri[iconn], \
                                                                  sesri[iconn], \
                                                                  sei[iconn], 1, \
                                                                  NP ); \
  auto dPhaseMassDens_dPView = getElementAccessor< true, 3, real64 >( &fluid, \
                                                                      MultiFluidBase::viewKeyStruct::dPhaseMassDensity_dPressureString(), \
                                                                      MAX_STENCIL, \
                                                                      seri[iconn], \
                                                                      sesri[iconn], \
                                                                      sei[iconn], \
                                                                      1, \
                                                                      NP ); \
  auto dPhaseMassDens_dCView = getElementAccessor< true, 4, real64 >( &fluid, \
                                                                      MultiFluidBase::viewKeyStruct::dPhaseMassDensity_dGlobalCompFractionString(), \
                                                                      MAX_STENCIL, \
                                                                      seri[iconn], \
                                                                      sesri[iconn], \
                                                                      sei[iconn], \
                                                                      1, NP, \
                                                                      NC );
#endif

#ifndef DEFINE_SUBR_FIELDS
#define DEFINE_SUBR_FIELDS( A ) \
  auto presView = getElementAccessor< true, 1, real64 >( &subRegion, \
                                                         CompositionalMultiphaseFVM::viewKeyStruct::pressureString(), \
                                                         MAX_STENCIL, \
                                                         seri[iconn], \
                                                         sesri[iconn], \
                                                         sei[iconn] ); \
  auto dPresView = getElementAccessor< true, 1, real64 >( &subRegion, \
                                                          CompositionalMultiphaseFVM::viewKeyStruct::deltaPressureString(), \
                                                          MAX_STENCIL, \
                                                          seri[iconn], \
                                                          sesri[iconn], \
                                                          sei[iconn] ); \
  Array< Array< Array< double, 1 >, 1 >, 1 > gravCoefView; \
  if( A ) \
  { \
    gravCoefView = getElementAccessor< true, 1, real64 >( &subRegion, \
                                                          CompositionalMultiphaseFVM::viewKeyStruct::gravityCoefString(), \
                                                          MAX_STENCIL, \
                                                          seri[iconn], \
                                                          sesri[iconn], \
                                                          sei[iconn] ); \
  } \
  else \
  { \
    real64 const temp[1] = {0.0}; \
    gravCoefView = AccessorHelper< true >::template makeElementAccessor< 1 >( &(temp[0]), \
                                                                              MAX_STENCIL, \
                                                                              seri[iconn], \
                                                                              sesri[iconn], \
                                                                              sei[iconn] ); \
  } \
  auto phaseMobView = getElementAccessor< true, 2, real64 >( &subRegion, \
                                                             CompositionalMultiphaseFVM::viewKeyStruct::phaseMobilityString(), \
                                                             MAX_STENCIL, \
                                                             seri[iconn], \
                                                             sesri[iconn], \
                                                             sei[iconn], NP ); \
  auto dPhaseMob_dPView = getElementAccessor< true, 2, real64 >( &subRegion, \
                                                                 CompositionalMultiphaseFVM::viewKeyStruct::dPhaseMobility_dPressureString(), \
                                                                 MAX_STENCIL, \
                                                                 seri[iconn], \
                                                                 sesri[iconn], \
                                                                 sei[iconn], \
                                                                 NP ); \
  auto dPhaseMob_dCView = getElementAccessor< true, 3, real64 >( &subRegion, \
                                                                 CompositionalMultiphaseFVM::viewKeyStruct::dPhaseMobility_dGlobalCompDensityString(), \
                                                                 MAX_STENCIL, \
                                                                 seri[iconn], \
                                                                 sesri[iconn], \
                                                                 sei[iconn], NP, \
                                                                 NC ); \
  auto \
    dCompFrac_dCompDensView = getElementAccessor< true, 3, real64 >( &subRegion, \
                                                                     CompositionalMultiphaseFVM::viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString(), \
                                                                     MAX_STENCIL, \
                                                                     seri[iconn], \
                                                                     sesri[iconn], \
                                                                     sei[iconn], \
                                                                     NC, \
                                                                     NC ); \
  auto dPhaseVolFrac_dPView = getElementAccessor< true, 2, real64 >( &subRegion, \
                                                                     CompositionalMultiphaseFVM::viewKeyStruct::dPhaseVolumeFraction_dPressureString(), \
                                                                     MAX_STENCIL, \
                                                                     seri[iconn], \
                                                                     sesri[iconn], \
                                                                     sei[iconn], \
                                                                     NP ); \
  auto dPhaseVolFrac_dCView = getElementAccessor< true, 3, real64 >( &subRegion, \
                                                                     CompositionalMultiphaseFVM::viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString(), \
                                                                     MAX_STENCIL, \
                                                                     seri[iconn], \
                                                                     sesri[iconn], \
                                                                     sei[iconn], \
                                                                     NP, \
                                                                     NC );
#endif

//// Sphinx end before input XML
template< bool FULL, localIndex N, typename T, typename ... DIMS >
auto getElementAccessor( Group const * const group,
                         std::string const & key,
                         localIndex const stencilSize,
                         arraySlice1d< localIndex const > const & stencilRegIndices,
                         arraySlice1d< localIndex const > const & stencilSubRegIndices,
                         arraySlice1d< localIndex const > const & stencilElemIndices,
                         DIMS... otherDims )
{

  ArrayView< T const, N > data;
  if( ElementSubRegionBase const * const subRegion = dynamicCast< ElementSubRegionBase const * const >( group ) )
  {
    data = subRegion->getReference< Array< T, N > >( key );
  }
  else if( MultiFluidBase const * const fluid = dynamicCast< MultiFluidBase const * const >( group ) )
  {
    data = fluid->getReference< Array< T, N > >( key );
  }
  else if( CapillaryPressureBase const * const cap = dynamicCast< CapillaryPressureBase const * const >( group ) )
  {
    data = cap->getReference< Array< T, N > >( key );
  }

  data.move( LvArray::MemorySpace::CPU, false );
  auto view = AccessorHelper< FULL >::template makeElementAccessor< N, T >( data.data() + stencilElemIndices[0],
                                                                            stencilSize,
                                                                            stencilRegIndices,
                                                                            stencilSubRegIndices,
                                                                            stencilElemIndices,
                                                                            otherDims ... );
  return view;

}

template< localIndex NC, localIndex NP >
void testCompositionalStandardUpwind( CompositionalMultiphaseFVM & solver,
                                      DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  solver.forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                         ElementSubRegionBase & subRegion )
  {
    SCOPED_TRACE( subRegion.getParent().getParent().getName() + "/" + subRegion.getName() );

    string const & fluidName = solver.fluidModelNames()[targetIndex];
    string const & capName = solver.capPresModelNames()[targetIndex];
    string const & discretizationName = solver.getDiscretization();

    MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
    CapillaryPressureBase & cap = subRegion.getConstitutiveModel< CapillaryPressureBase >( capName );

    // reset the solver state to zero out variable updates
    solver.resetStateToBeginningOfStep( domain );

    FluxApproximationBase const & fluxApprox =
      domain.getNumericalMethodManager().getFiniteVolumeManager().getFluxApproximation( discretizationName );
    fluxApprox.forStencils< CellElementStencilTPFA >( mesh, [&]( auto const & stencil )
    {

      auto const & weights = stencil.getWeights();
      auto const & seri = stencil.getElementRegionIndices();
      auto const & sesri = stencil.getElementSubRegionIndices();
      auto const & sei = stencil.getElementIndices();
      localIndex constexpr NUM_ELEMS = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
      localIndex constexpr MAX_STENCIL = CellElementStencilTPFA::MAX_STENCIL_SIZE;


      for( localIndex iconn = 0; iconn < stencil.size(); ++iconn )
      {
        DEFINE_SUBR_FIELDS( 1 )
        //include cap pressure
        DEFINE_CAP_FIELDS( 1 );
        // -- fluid fetched fields
        DEFINE_FLUID_FIELDS()

        real64 totFlux = 0;

        localIndex k_up[NP]{};
        real64 phaseFlux[NP]{};
        real64 dPhaseFlux_dP[NP][MAX_STENCIL]{};
        real64 dPhaseFlux_dC[NP][MAX_STENCIL][NC]{};

        real64 gravHead[NP]{};
        real64 dGravHead_dP[NP][NUM_ELEMS]{};
        real64 dGravHead_dC[NP][NUM_ELEMS][NC]{};
        real64 dProp_dC[NP][NC]{};

        real64 capHead[NP]{};
        real64 dCapHead_dP[NP][MAX_STENCIL]{};
        real64 dCapHead_dC[NP][MAX_STENCIL][NC]{};

        for( localIndex ip = 0; ip < NP; ++ip )
        {


          UHelpers::formPPUVelocity< NC, NUM_ELEMS, MAX_STENCIL >( NP,
                                                                   ip,
                                                                   MAX_STENCIL,
                                                                   seri[iconn],
                                                                   sesri[iconn],
                                                                   sei[iconn],
                                                                   weights[iconn],
                                                                   presView.toNestedViewConst(),
                                                                   dPresView.toNestedViewConst(),
                                                                   gravCoefView.toNestedViewConst(),
                                                                   phaseMobView.toNestedViewConst(),
                                                                   dPhaseMob_dPView.toNestedViewConst(),
                                                                   dPhaseMob_dCView.toNestedViewConst(),
                                                                   dPhaseVolFrac_dPView.toNestedViewConst(),
                                                                   dPhaseVolFrac_dCView.toNestedViewConst(),
                                                                   dCompFrac_dCompDensView.toNestedViewConst(),
                                                                   phaseMassDensView.toNestedViewConst(),
                                                                   dPhaseMassDens_dPView.toNestedViewConst(),
                                                                   dPhaseMassDens_dCView.toNestedViewConst(),
                                                                   phaseCapPressureView.toNestedViewConst(),
                                                                   dPhaseCapPressure_dPhaseVolFracView.toNestedViewConst(),
                                                                   capPressureFlag,
                                                                   k_up[ip],
                                                                   phaseFlux[ip],
                                                                   dPhaseFlux_dP[ip],
                                                                   dPhaseFlux_dC[ip] );

          totFlux += phaseFlux[ip];

          UHelpers::formPotential< NC, CompositionalMultiphaseFlowUpwindHelperKernels::term::Gravity, NUM_ELEMS, MAX_STENCIL >::compute(
            NP,
            ip,
            MAX_STENCIL,
            seri[iconn],
            sesri[iconn],
            sei[iconn],
            weights[iconn],
            totFlux,
            gravCoefView.toNestedViewConst(),
            dCompFrac_dCompDensView.toNestedViewConst(),
            phaseMassDensView.toNestedViewConst(),
            dPhaseMassDens_dPView.toNestedViewConst(),
            dPhaseMassDens_dCView.toNestedViewConst(),
            dPhaseVolFrac_dPView.toNestedViewConst(),
            dPhaseVolFrac_dCView.toNestedViewConst(),
            phaseCapPressureView.toNestedViewConst(),
            dPhaseCapPressure_dPhaseVolFracView.toNestedViewConst(),
            gravHead[ip],
            dGravHead_dP[ip],
            dGravHead_dC[ip],
            dProp_dC[ip]
            );

          UHelpers::formPotential< NC, CompositionalMultiphaseFlowUpwindHelperKernels::term::Capillary, NUM_ELEMS, MAX_STENCIL >::compute(
            NP,
            ip,
            MAX_STENCIL,
            seri[iconn],
            sesri[iconn],
            sei[iconn],
            weights[iconn],
            totFlux,
            gravCoefView.toNestedViewConst(),
            dCompFrac_dCompDensView.toNestedViewConst(),
            phaseMassDensView.toNestedViewConst(),
            dPhaseMassDens_dPView.toNestedViewConst(),
            dPhaseMassDens_dCView.toNestedViewConst(),
            dPhaseVolFrac_dPView.toNestedViewConst(),
            dPhaseVolFrac_dCView.toNestedViewConst(),
            phaseCapPressureView.toNestedViewConst(),
            dPhaseCapPressure_dPhaseVolFracView.toNestedViewConst(),
            capHead[ip],
            dCapHead_dP[ip],
            dCapHead_dC[ip],
            dProp_dC[ip]
            );
        }                                  //standar loop

        localIndex k_up_scheme[NP]{};
        real64 fflow[NP]{};
        real64 dFflow_dP[NP][MAX_STENCIL]{};
        real64 dFflow_dC[NP][MAX_STENCIL][NC]{};
        for( localIndex ip = 0; ip < NP; ++ip )
        {


          UHelpers::formFracFlow< NC, NUM_ELEMS, MAX_STENCIL,
                                  CompositionalMultiphaseFlowUpwindHelperKernels::term::Viscous,
                                  CompositionalMultiphaseFlowUpwindHelperKernels::PhasePotentialUpwind >( NP,
                                                                                                          ip,
                                                                                                          MAX_STENCIL,
                                                                                                          seri[iconn],
                                                                                                          sesri[iconn],
                                                                                                          sei[iconn],
                                                                                                          weights[iconn],
                                                                                                          totFlux,
                                                                                                          presView.toNestedViewConst(),
                                                                                                          dPresView.toNestedViewConst(),
                                                                                                          gravCoefView.toNestedViewConst(),
                                                                                                          dCompFrac_dCompDensView.toNestedViewConst(),
                                                                                                          phaseMassDensView.toNestedViewConst(),
                                                                                                          dPhaseMassDens_dPView.toNestedViewConst(),
                                                                                                          dPhaseMassDens_dCView.toNestedViewConst(),
                                                                                                          phaseMobView.toNestedViewConst(),
                                                                                                          dPhaseMob_dPView.toNestedViewConst(),
                                                                                                          dPhaseMob_dCView.toNestedViewConst(),
                                                                                                          dPhaseVolFrac_dPView.toNestedViewConst(),
                                                                                                          dPhaseVolFrac_dCView.toNestedViewConst(),
                                                                                                          phaseCapPressureView.toNestedViewConst(),
                                                                                                          dPhaseCapPressure_dPhaseVolFracView.toNestedViewConst(),
                                                                                                          capPressureFlag,
                                                                                                          k_up_scheme[ip],
                                                                                                          fflow[ip],
                                                                                                          dFflow_dP[ip],
                                                                                                          dFflow_dC[ip] );

          EXPECT_EQ( k_up_scheme[ip], k_up[ip] );
        }                                  //loop on scheme to test


        //test fluxes values
        for( localIndex ip = 0; ip < NP; ++ip )
        {
          real64 phaseFluxScheme{};
          phaseFluxScheme += fflow[ip] * totFlux;
          for( localIndex jp = 0; jp < NP; ++jp )
          {
            if( ip != jp )
            {
              phaseFluxScheme -=
                fflow[ip] * phaseMobView[seri[iconn][0]][sesri[iconn][0]][sei[iconn][k_up_scheme[jp]]][jp]
                * ( (gravHead[ip] + capHead[ip] ) - (gravHead[jp] + capHead[jp]) );

            }
          }

          real64 const relTol = 1e-12;
          real64 const absTol = 1e-18;
          checkRelativeError( phaseFlux[ip], phaseFluxScheme, relTol, absTol );

        }//compare flux construction


      }//inner most loop on connexion
    } );//second level lambda on stencil list
  } );//outer lambda on subregions

}//EOfunc

template< localIndex NC, localIndex NP, term T, bool isViscous = (T==term::Viscous), bool isGrav = (T==term::Gravity), bool isCapillary = (T==term::Capillary) >
void testCompositionalUpwindHUPU( CompositionalMultiphaseFVM & solver,
                                  DomainPartition & domain )
{

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  solver.forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                         ElementSubRegionBase & subRegion )
  {
    SCOPED_TRACE( subRegion.getParent().getParent().getName() + "/" + subRegion.getName() );

    string const & fluidName = solver.fluidModelNames()[targetIndex];
    string const & discretizationName = solver.getDiscretization();
    string const & capName = solver.capPresModelNames()[targetIndex];

    MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
    CapillaryPressureBase & cap = subRegion.getConstitutiveModel< CapillaryPressureBase >( capName );
    // reset the solver state to zero out variable updates
    solver.resetStateToBeginningOfStep( domain );

    FluxApproximationBase const & fluxApprox =
      domain.getNumericalMethodManager().getFiniteVolumeManager().getFluxApproximation( discretizationName );
    fluxApprox.forStencils< CellElementStencilTPFA >( mesh, [&]( auto const & stencil )
    {

      auto const & weights = stencil.getWeights();
      auto const & seri = stencil.getElementRegionIndices();
      auto const & sesri = stencil.getElementSubRegionIndices();
      auto const & sei = stencil.getElementIndices();
      localIndex constexpr NUM_ELEMS = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
      localIndex constexpr MAX_STENCIL = CellElementStencilTPFA::MAX_STENCIL_SIZE;

      for( auto iconn = 0; iconn < stencil.size(); ++iconn )
      {
        DEFINE_SUBR_FIELDS( isGrav )
        DEFINE_CAP_FIELDS( isCapillary )
        // -- fluid fetched fields
        DEFINE_FLUID_FIELDS()


        real64 totFlux = 0;
        localIndex k_up_hu{};
        localIndex k_up_pu{};
        real64 fflow[2]{};
        real64 dFflow_dP[2][MAX_STENCIL]{};
        real64 dFflow_dC[2][MAX_STENCIL][NC]{};

        localIndex k_up{};
        real64 phaseFlux[NP]{};
        real64 dPhaseFlux_dP[NP][MAX_STENCIL]{};
        real64 dPhaseFlux_dC[NP][MAX_STENCIL][NC]{};
        for( localIndex ip = 0; ip < NP; ++ip )
        {
          UHelpers::formPPUVelocity< NC, NUM_ELEMS, MAX_STENCIL >( NP,
                                                                   ip,
                                                                   MAX_STENCIL,
                                                                   seri[iconn],
                                                                   sesri[iconn],
                                                                   sei[iconn],
                                                                   weights[iconn],
                                                                   presView.toNestedViewConst(),
                                                                   dPresView.toNestedViewConst(),
                                                                   gravCoefView.toNestedViewConst(),
                                                                   phaseMobView.toNestedViewConst(),
                                                                   dPhaseMob_dPView.toNestedViewConst(),
                                                                   dPhaseMob_dCView.toNestedViewConst(),
                                                                   dPhaseVolFrac_dPView.toNestedViewConst(),
                                                                   dPhaseVolFrac_dCView.toNestedViewConst(),
                                                                   dCompFrac_dCompDensView.toNestedViewConst(),
                                                                   phaseMassDensView.toNestedViewConst(),
                                                                   dPhaseMassDens_dPView.toNestedViewConst(),
                                                                   dPhaseMassDens_dCView.toNestedViewConst(),
                                                                   phaseCapPressureView.toNestedViewConst(),
                                                                   dPhaseCapPressure_dPhaseVolFracView.toNestedViewConst(),
                                                                   capPressureFlag,
                                                                   k_up,
                                                                   phaseFlux[ip],
                                                                   dPhaseFlux_dP[ip],
                                                                   dPhaseFlux_dC[ip] );

          totFlux += phaseFlux[ip];

        }

        //standard loop
        for( localIndex ip = 0; ip < NP; ++ip )
        {
          UHelpers::formFracFlow< NC, NUM_ELEMS, MAX_STENCIL,
                                  T,
                                  HybridUpwind >( NP,
                                                  ip,
                                                  MAX_STENCIL,
                                                  seri[iconn],
                                                  sesri[iconn],
                                                  sei[iconn],
                                                  weights[iconn],
                                                  (isViscous) ? totFlux : 0.0,
                                                  presView.toNestedViewConst(),
                                                  dPresView.toNestedViewConst(),
                                                  gravCoefView.toNestedViewConst(),
                                                  dCompFrac_dCompDensView.toNestedViewConst(),
                                                  phaseMassDensView.toNestedViewConst(),
                                                  dPhaseMassDens_dPView.toNestedViewConst(),
                                                  dPhaseMassDens_dCView.toNestedViewConst(),
                                                  phaseMobView.toNestedViewConst(),
                                                  dPhaseMob_dPView.toNestedViewConst(),
                                                  dPhaseMob_dCView.toNestedViewConst(),
                                                  dPhaseVolFrac_dPView.toNestedViewConst(),
                                                  dPhaseVolFrac_dCView.toNestedViewConst(),
                                                  phaseCapPressureView.toNestedViewConst(),
                                                  dPhaseCapPressure_dPhaseVolFracView.toNestedViewConst(),
                                                  capPressureFlag,
                                                  k_up_hu,
                                                  fflow[0],
                                                  dFflow_dP[0],
                                                  dFflow_dC[0] );

          UHelpers::formFracFlow< NC, NUM_ELEMS, MAX_STENCIL,
                                  T,
                                  PhaseUpwind >( NP,
                                                 ip,
                                                 MAX_STENCIL,
                                                 seri[iconn],
                                                 sesri[iconn],
                                                 sei[iconn],
                                                 weights[iconn],
                                                 (isViscous) ? totFlux : 0.0,
                                                 presView.toNestedViewConst(),
                                                 dPresView.toNestedViewConst(),
                                                 gravCoefView.toNestedViewConst(),
                                                 dCompFrac_dCompDensView.toNestedViewConst(),
                                                 phaseMassDensView.toNestedViewConst(),
                                                 dPhaseMassDens_dPView.toNestedViewConst(),
                                                 dPhaseMassDens_dCView.toNestedViewConst(),
                                                 phaseMobView.toNestedViewConst(),
                                                 dPhaseMob_dPView.toNestedViewConst(),
                                                 dPhaseMob_dCView.toNestedViewConst(),
                                                 dPhaseVolFrac_dPView.toNestedViewConst(),
                                                 dPhaseVolFrac_dCView.toNestedViewConst(),
                                                 phaseCapPressureView.toNestedViewConst(),
                                                 dPhaseCapPressure_dPhaseVolFracView.toNestedViewConst(),
                                                 capPressureFlag,
                                                 k_up_pu,
                                                 fflow[1],
                                                 dFflow_dP[1],
                                                 dFflow_dC[1] );

          EXPECT_EQ( k_up_hu, k_up_pu );
          EXPECT_EQ( fflow[0], fflow[1] );
        }

      }
    } );
  } );
}

class CompositionalMultiphaseFlowUpwindHelperKernelsTest : public ::testing::Test
{
public:

  CompositionalMultiphaseFlowUpwindHelperKernelsTest()
    :
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  { }

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< CompositionalMultiphaseFVM >( "compflow" );

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
  CompositionalMultiphaseFVM * solver;
};

real64 constexpr CompositionalMultiphaseFlowUpwindHelperKernelsTest::time;
real64 constexpr CompositionalMultiphaseFlowUpwindHelperKernelsTest::dt;
real64 constexpr CompositionalMultiphaseFlowUpwindHelperKernelsTest::eps;

TEST_F( CompositionalMultiphaseFlowUpwindHelperKernelsTest, test_standardPPU )
{
  localIndex constexpr NP = 2;
  localIndex constexpr NC = 4;
  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  testCompositionalStandardUpwind< NC, NP >( *solver, domain );
}

TEST_F( CompositionalMultiphaseFlowUpwindHelperKernelsTest, test_HUPUViscous )
{
  localIndex constexpr NP = 2;
  localIndex constexpr NC = 4;
  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  testCompositionalUpwindHUPU< NC, NP, term::Viscous >( *solver, domain );
}

TEST_F( CompositionalMultiphaseFlowUpwindHelperKernelsTest, test_HUPUGravity )
{
  localIndex constexpr NP = 2;
  localIndex constexpr NC = 4;
  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  testCompositionalUpwindHUPU< NC, NP, term::Gravity >( *solver, domain );
}

TEST_F( CompositionalMultiphaseFlowUpwindHelperKernelsTest, test_HUPUCapillary )
{
  localIndex constexpr NP = 2;
  localIndex constexpr NC = 4;
  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  testCompositionalUpwindHUPU< NC, NP, term::Capillary >( *solver, domain );
}

int main( int argc,
          char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
