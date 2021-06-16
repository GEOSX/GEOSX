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
 * @file SolidBase.cpp
 */

#include "SolidBase.hpp"

namespace geosx
{
using namespace dataRepository;

namespace constitutive
{

SolidBase::SolidBase( string const & name, Group * const parent ):
  RockBase( name, parent ),
  m_newStress( 0, 0, 6 ),
  m_oldStress( 0, 0, 6 )
{
  registerWrapper( viewKeyStruct::stressString(), &m_newStress ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setApplyDefaultValue( 0 ). // default to zero initial stress
    setDescription( "Current Material Stress" );

  registerWrapper( viewKeyStruct::oldStressString(), &m_oldStress ).
    setApplyDefaultValue( 0 ). // default to zero initial stress
    setDescription( "Previous Material Stress" );
}


SolidBase::~SolidBase()
{}


void SolidBase::postProcessInput()
{}


void SolidBase::allocateConstitutiveData( dataRepository::Group & parent,
                                          localIndex const numConstitutivePointsPerParentIndex )
{
  m_newStress.resize( 0, numConstitutivePointsPerParentIndex, 6 );
  m_oldStress.resize( 0, numConstitutivePointsPerParentIndex, 6 );

  RockBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}


void SolidBase::saveConvergedState() const
{
  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView3d< real64 const, solid::STRESS_USD > newStress = m_newStress;
  arrayView3d< real64, solid::STRESS_USD > oldStress = m_oldStress;

  arrayView2d< real64 const > newPorosity = m_newPorosity;
  arrayView2d< real64 > oldPorosity = m_oldPorosity;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      LvArray::tensorOps::copy< 6 >( oldStress[k][q], newStress[k][q] );
      oldPorosity[k][q] = newPorosity[k][q];
    }
  } );
}


} /* namespace constitutive */
} /* namespace geosx */
