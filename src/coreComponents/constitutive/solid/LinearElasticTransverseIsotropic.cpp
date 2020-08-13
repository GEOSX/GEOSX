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
 *  @file LinearElasticTransverseIsotropic.cpp
 */

#include "LinearElasticTransverseIsotropic.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{
LinearElasticTransverseIsotropic::LinearElasticTransverseIsotropic(
  std::string const & name,
  Group * const parent ) :
  SolidBase( name, parent ),
  m_defaultYoungsModulusTransverse(),
  m_defaultYoungsModulusAxial(),
  m_defaultPoissonTransverse(),
  m_defaultPoissonAxialTransverse(),
  m_defaultShearModulusAxialTransverse(),
  m_c11(),
  m_c13(),
  m_c33(),
  m_c44(),
  m_c66()
{
  registerWrapper( viewKeyStruct::defaultYoungsModulusTransverse,
                   &m_defaultYoungsModulusTransverse )
    ->setApplyDefaultValue( -1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Elastic Bulk Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultYoungsModulusAxial,
                   &m_defaultYoungsModulusAxial )
    ->setApplyDefaultValue( -1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Elastic Shear Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultPoissonRatioTransverse,
                   &m_defaultPoissonTransverse )
    ->setApplyDefaultValue( -1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Elastic Shear Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultPoissonRatioAxialTransverse,
                   &m_defaultPoissonAxialTransverse )
    ->setApplyDefaultValue( -1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Elastic Shear Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultShearModulusAxialTransverse,
                   &m_defaultShearModulusAxialTransverse )
    ->setApplyDefaultValue( -1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Elastic Shear Modulus Parameter" );

  registerWrapper( viewKeyStruct::c11, &m_c11 )
    ->setApplyDefaultValue( -1 )
    ->setDescription( "Elastic Bulk Modulus Field" );

  registerWrapper( viewKeyStruct::c13, &m_c13 )
    ->setApplyDefaultValue( -1 )
    ->setDescription( "Elastic Bulk Modulus Field" );

  registerWrapper( viewKeyStruct::c33, &m_c33 )
    ->setApplyDefaultValue( -1 )
    ->setDescription( "Elastic Bulk Modulus Field" );

  registerWrapper( viewKeyStruct::c44, &m_c44 )
    ->setApplyDefaultValue( -1 )
    ->setDescription( "Elastic Bulk Modulus Field" );

  registerWrapper( viewKeyStruct::c66, &m_c66 )
    ->setApplyDefaultValue( -1 )
    ->setDescription( "Elastic Bulk Modulus Field" );
}

LinearElasticTransverseIsotropic::~LinearElasticTransverseIsotropic()
{}

void
LinearElasticTransverseIsotropic::DeliverClone(
  string const & name,
  Group * const parent,
  std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< LinearElasticTransverseIsotropic >( name, parent );
  }
  SolidBase::DeliverClone( name, parent, clone );
  LinearElasticTransverseIsotropic * const newConstitutiveRelation =
    dynamic_cast< LinearElasticTransverseIsotropic * >( clone.get() );

  newConstitutiveRelation->m_defaultYoungsModulusTransverse =
    m_defaultYoungsModulusTransverse;
  newConstitutiveRelation->m_defaultYoungsModulusAxial =
    m_defaultYoungsModulusAxial;
  newConstitutiveRelation->m_defaultPoissonTransverse =
    m_defaultPoissonTransverse;
  newConstitutiveRelation->m_defaultPoissonAxialTransverse =
    m_defaultPoissonAxialTransverse;
  newConstitutiveRelation->m_defaultShearModulusAxialTransverse =
    m_defaultShearModulusAxialTransverse;

  newConstitutiveRelation->m_c11 = m_c11;
  newConstitutiveRelation->m_c13 = m_c13;
  newConstitutiveRelation->m_c33 = m_c33;
  newConstitutiveRelation->m_c44 = m_c44;
  newConstitutiveRelation->m_c66 = m_c66;
  newConstitutiveRelation->m_stress = m_stress;
}

void
LinearElasticTransverseIsotropic::AllocateConstitutiveData(
  dataRepository::Group * const parent,
  localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::AllocateConstitutiveData( parent,
                                       numConstitutivePointsPerParentIndex );
  this->resize( parent->size() );
}

void
LinearElasticTransverseIsotropic::PostProcessInput()
{
  if( !m_postProcessed )
  {
    real64 const Et = m_defaultYoungsModulusTransverse;
    real64 const Ea = m_defaultYoungsModulusAxial;
    real64 const Nut = m_defaultPoissonTransverse;
    real64 const Nuat = m_defaultPoissonAxialTransverse;
    real64 const Gat = m_defaultShearModulusAxialTransverse;
    real64 const Nuta = Nuat * ( Et / Ea );

    real64 const delta = ( 1 + Nut ) * ( 1 - Nut - 2 * Nuta * Nuat ) / ( Et * Et * Ea );
    real64 const c11Default = ( 1.0 - Nuta * Nuat ) / ( Ea * Et * delta );
    real64 const c13Default = Nuat * ( 1.0 + Nut ) / ( Ea * Et * delta );
    real64 const c33Default = ( 1 - Nut * Nut ) / ( Ea * Ea * delta );
    real64 const c44Default = Gat;
    real64 const c66Default = 0.5 * Ea / ( 1 + Nut );

    this->getWrapper< array1d< real64 > >( viewKeyStruct::c11 )
      ->setApplyDefaultValue( c11Default );

    this->getWrapper< array1d< real64 > >( viewKeyStruct::c13 )
      ->setApplyDefaultValue( c13Default );

    this->getWrapper< array1d< real64 > >( viewKeyStruct::c33 )
      ->setApplyDefaultValue( c33Default );

    this->getWrapper< array1d< real64 > >( viewKeyStruct::c44 )
      ->setApplyDefaultValue( c44Default );

    this->getWrapper< array1d< real64 > >( viewKeyStruct::c66 )
      ->setApplyDefaultValue( c66Default );
  }
  m_postProcessed = true;
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        LinearElasticTransverseIsotropic,
                        std::string const &,
                        Group * const )
}  // namespace constitutive
} /* namespace geosx */
