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
 *  @file ElasticOrthotropic.cpp
 */

#include "ElasticOrthotropic.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

ElasticOrthotropic::ElasticOrthotropic( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_defaultE1(),
  m_defaultE2(),
  m_defaultE3(),
  m_defaultNu12(),
  m_defaultNu13(),
  m_defaultNu23(),
  m_defaultG12(),
  m_defaultG13(),
  m_defaultG23(),
  m_c11(),
  m_c12(),
  m_c13(),
  m_c22(),
  m_c23(),
  m_c33(),
  m_c44(),
  m_c55(),
  m_c66()
{
  registerWrapper( viewKeyStruct::defaultE1String(), &m_defaultE1 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Young's Modulus E1" );

  registerWrapper( viewKeyStruct::defaultE2String(), &m_defaultE2 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Young's Modulus E2" );

  registerWrapper( viewKeyStruct::defaultE3String(), &m_defaultE3 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Young's Modulus E3" );

  registerWrapper( viewKeyStruct::defaultNu12String(), &m_defaultNu12 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Poission's Ratio Nu12" );

  registerWrapper( viewKeyStruct::defaultNu13String(), &m_defaultNu13 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Poission's Ratio Nu13" );

  registerWrapper( viewKeyStruct::defaultNu23String(), &m_defaultNu23 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Poission's Ratio Nu23" );

  registerWrapper( viewKeyStruct::defaultG12String(), &m_defaultG12 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Shear Modulus G12" );

  registerWrapper( viewKeyStruct::defaultG13String(), &m_defaultG13 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Shear Modulus G13" );

  registerWrapper( viewKeyStruct::defaultG23String(), &m_defaultG23 ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Shear Modulus G23" );

  registerWrapper< real64 >( viewKeyStruct::defaultC11String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C11 Component of Voigt Stiffness Tensor" );

  registerWrapper< real64 >( viewKeyStruct::defaultC12String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C12 Component of Voigt Stiffness Tensor" );

  registerWrapper< real64 >( viewKeyStruct::defaultC13String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C13 Component of Voigt Stiffness Tensor" );

  registerWrapper< real64 >( viewKeyStruct::defaultC22String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C22 Component of Voigt Stiffness Tensor" );

  registerWrapper< real64 >( viewKeyStruct::defaultC23String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C23 Component of Voigt Stiffness Tensor" );

  registerWrapper< real64 >( viewKeyStruct::defaultC33String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C33 Component of Voigt Stiffness Tensor" );

  registerWrapper< real64 >( viewKeyStruct::defaultC44String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C44 Component of Voigt Stiffness Tensor" );

  registerWrapper< real64 >( viewKeyStruct::defaultC55String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C55 Component of Voigt Stiffness Tensor" );

  registerWrapper< real64 >( viewKeyStruct::defaultC66String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default C66 Component of Voigt Stiffness Tensor" );

  registerWrapper( viewKeyStruct::c11String(), &m_c11 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C11" );

  registerWrapper( viewKeyStruct::c12String(), &m_c12 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C12" );

  registerWrapper( viewKeyStruct::c13String(), &m_c13 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C13" );

  registerWrapper( viewKeyStruct::c22String(), &m_c22 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C22" );

  registerWrapper( viewKeyStruct::c23String(), &m_c23 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C23" );

  registerWrapper( viewKeyStruct::c33String(), &m_c33 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C33" );

  registerWrapper( viewKeyStruct::c44String(), &m_c44 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C44" );

  registerWrapper( viewKeyStruct::c55String(), &m_c55 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C55" );

  registerWrapper( viewKeyStruct::c66String(), &m_c66 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C66" );
}

ElasticOrthotropic::~ElasticOrthotropic()
{}

void ElasticOrthotropic::postProcessInput()
{
  SolidBase::postProcessInput();

  real64 & c11 = getReference< real64 >( viewKeyStruct::defaultC11String() );
  real64 & c12 = getReference< real64 >( viewKeyStruct::defaultC12String() );
  real64 & c13 = getReference< real64 >( viewKeyStruct::defaultC13String() );
  real64 & c22 = getReference< real64 >( viewKeyStruct::defaultC22String() );
  real64 & c23 = getReference< real64 >( viewKeyStruct::defaultC23String() );
  real64 & c33 = getReference< real64 >( viewKeyStruct::defaultC33String() );
  real64 & c44 = getReference< real64 >( viewKeyStruct::defaultC44String() );
  real64 & c55 = getReference< real64 >( viewKeyStruct::defaultC55String() );
  real64 & c66 = getReference< real64 >( viewKeyStruct::defaultC66String() );

  real64 const E1   = m_defaultE1;
  real64 const E2   = m_defaultE2;
  real64 const E3   = m_defaultE3;
  real64 const Nu12 = m_defaultNu12;
  real64 const Nu13 = m_defaultNu13;
  real64 const Nu23 = m_defaultNu23;
  real64 const G12  = m_defaultG12;
  real64 const G13  = m_defaultG13;
  real64 const G23  = m_defaultG23;

  if( E1 > 0.0 && E2 > 0.0 && E3 > 0.0 && G12 > 0.0 && G13 > 0.0 && G23 > 0.0 )
  {
    real64 const Nu21 = Nu12 * E2 / E1;
    real64 const Nu31 = Nu13 * E3 / E1;
    real64 const Nu32 = Nu23 * E3 / E2;

    real64 const delta = 1.0 - Nu12 * Nu21 - Nu13 * Nu31 - Nu23 * Nu32 - 2.0 * Nu12 * Nu23 * Nu31;

    if( delta > 0.0 && Nu23 * Nu32 < 1.0 && Nu31 * Nu13 < 1.0 && Nu12 * Nu21 < 1.0 )
    {
      c11 = ( 1.0 - Nu23 * Nu32 )  * E1 / delta;
      c12 = ( Nu21 + Nu31 * Nu23 ) * E1 / delta;
      c13 = ( Nu31 + Nu21 * Nu32 ) * E1 / delta;

      c22 = ( 1.0 - Nu31 * Nu13 )  * E2 / delta;
      c23 = ( Nu32 + Nu31 * Nu12 ) * E2 / delta;

      c33 = ( 1.0 - Nu12 * Nu21 )  * E3 / delta;

      c44 = G23;
      c55 = G13;
      c66 = G12;
    }
    else
    {
      GEOSX_ERROR( "Invalid specification for default elastic constants." );
    }
  }
  else if( c11 <= 0.0 || c22 <= 0.0 || c33 <= 0.0 || c44 <= 0.0 || c55 <= 0.0 || c66 <= 0.0 )
  {
    GEOSX_ERROR( "Invalid specification for default elastic stiffnesses." );
  }

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c11String() ).
    setApplyDefaultValue( c11 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c12String() ).
    setApplyDefaultValue( c12 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c13String() ).
    setApplyDefaultValue( c13 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c22String() ).
    setApplyDefaultValue( c22 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c23String() ).
    setApplyDefaultValue( c23 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c33String() ).
    setApplyDefaultValue( c33 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c44String() ).
    setApplyDefaultValue( c44 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c55String() ).
    setApplyDefaultValue( c55 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c66String() ).
    setApplyDefaultValue( c66 );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ElasticOrthotropic, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geosx */
