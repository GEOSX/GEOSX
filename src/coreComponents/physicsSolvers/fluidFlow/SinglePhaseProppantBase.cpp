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

/**
 * @file SinglePhaseProppantBase.cpp
 */

#include "SinglePhaseProppantBase.hpp"

#include "constitutive/fluid/slurryFluidSelector.hpp"
#include "physicsSolvers/fluidFlow/ProppantTransport.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseProppantBaseKernels.hpp"

namespace geosx
{
using constitutive::SlurryFluidBase;

SinglePhaseProppantBase::SinglePhaseProppantBase( const std::string & name,
                                                  Group * const parent ) :
  SinglePhaseBase( name, parent )
{}

SinglePhaseProppantBase::~SinglePhaseProppantBase()
{
  // TODO Auto-generated destructor stub
}

void
SinglePhaseProppantBase::ValidateFluidModels( DomainPartition const & domain ) const
{
  // Validate fluid models in regions
  for( auto & mesh : domain.getMeshBodies()->GetSubGroups() )
  {
    MeshLevel const & meshLevel =
      *Group::group_cast< MeshBody const * >( mesh.second )->getMeshLevel( 0 );
    ValidateModelMapping< SlurryFluidBase >( *meshLevel.getElemManager(),
                                             m_fluidModelNames );
  }
}

SinglePhaseBase::FluidPropViews
SinglePhaseProppantBase::getFluidProperties(
  constitutive::ConstitutiveBase const & fluid ) const
{
  SlurryFluidBase const & slurryFluid = dynamicCast< SlurryFluidBase const & >( fluid );
  return {
    slurryFluid.density(),
    slurryFluid.dDensity_dPressure(),
    slurryFluid.viscosity(),
    slurryFluid.dViscosity_dPressure(),
    slurryFluid
      .getWrapper< array2d< real64 > >( SlurryFluidBase::viewKeyStruct::densityString )
      ->getDefaultValue(),
    slurryFluid
      .getWrapper< array2d< real64 > >( SlurryFluidBase::viewKeyStruct::viscosityString )
      ->getDefaultValue() };
}

arrayView1d< real64 const > const &
SinglePhaseProppantBase::getPoreVolumeMult(
  ElementSubRegionBase const & subRegion ) const
{
  return subRegion.getReference< array1d< real64 > >(
    ProppantTransport::viewKeyStruct::poroMultiplierString );
}

void
SinglePhaseProppantBase::UpdateFluidModel( Group & dataGroup,
                                           localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const & pres =
    dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString );

  arrayView1d< real64 const > const & dPres =
    dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  arrayView1d< real64 const > const & proppantConcentration =
    dataGroup.getReference< array1d< real64 > >(
      ProppantTransport::viewKeyStruct::proppantConcentrationString );

  arrayView1d< real64 const > const & dProppantConcentration =
    dataGroup.getReference< array1d< real64 > >(
      ProppantTransport::viewKeyStruct::deltaProppantConcentrationString );

  arrayView2d< real64 const > const & componentConcentration =
    dataGroup.getReference< array2d< real64 > >(
      ProppantTransport::viewKeyStruct::componentConcentrationString );

  arrayView1d< R1Tensor const > const & cellBasedFlux =
    dataGroup.getReference< array1d< R1Tensor > >(
      ProppantTransport::viewKeyStruct::cellBasedFluxString );

  arrayView1d< integer const > const & isProppantBoundaryElement =
    dataGroup.getReference< array1d< integer > >(
      ProppantTransport::viewKeyStruct::isProppantBoundaryString );

  SlurryFluidBase & fluid =
    GetConstitutiveModel< SlurryFluidBase >( dataGroup,
                                             m_fluidModelNames[targetIndex] );

  constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid ) {
    typename TYPEOFREF( castedFluid )::KernelWrapper fluidWrapper =
      castedFluid.createKernelWrapper();
    SinglePhaseProppantBaseKernels::FluidUpdateKernel::Launch(
      fluidWrapper,
      pres,
      dPres,
      proppantConcentration,
      dProppantConcentration,
      componentConcentration,
      cellBasedFlux,
      isProppantBoundaryElement );
  } );
}

void
SinglePhaseProppantBase::ResetViewsPrivate( ElementRegionManager const & elemManager )
{
  m_density.clear();
  m_density = elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >(
    SlurryFluidBase::viewKeyStruct::densityString,
    targetRegionNames(),
    fluidModelNames() );
  m_density.setName( getName() + "/accessors/" +
                     SlurryFluidBase::viewKeyStruct::densityString );

  m_dDens_dPres.clear();
  m_dDens_dPres = elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >(
    SlurryFluidBase::viewKeyStruct::dDens_dPresString,
    targetRegionNames(),
    fluidModelNames() );
  m_dDens_dPres.setName( getName() + "/accessors/" +
                         SlurryFluidBase::viewKeyStruct::dDens_dPresString );

  m_viscosity.clear();
  m_viscosity = elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >(
    SlurryFluidBase::viewKeyStruct::viscosityString,
    targetRegionNames(),
    fluidModelNames() );
  m_viscosity.setName( getName() + "/accessors/" +
                       SlurryFluidBase::viewKeyStruct::viscosityString );

  m_dVisc_dPres.clear();
  m_dVisc_dPres = elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >(
    SlurryFluidBase::viewKeyStruct::dVisc_dPresString,
    targetRegionNames(),
    fluidModelNames() );
  m_dVisc_dPres.setName( getName() + "/accessors/" +
                         SlurryFluidBase::viewKeyStruct::dVisc_dPresString );

  m_transTMultiplier.clear();
  m_transTMultiplier = elemManager.ConstructArrayViewAccessor< R1Tensor, 1 >(
    ProppantTransport::viewKeyStruct::transTMultiplierString );
  m_transTMultiplier.setName(
    getName() + "/accessors/" +
    ProppantTransport::viewKeyStruct::transTMultiplierString );
}

}  // namespace geosx
