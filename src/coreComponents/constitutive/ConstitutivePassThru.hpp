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
 * @file solidSelector.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_SOLIDSELECTOR_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_SOLIDSELECTOR_HPP_

#include "NullModel.hpp"
#include "solid/LinearElasticIsotropic.hpp"
#include "solid/LinearElasticAnisotropic.hpp"
#include "solid/LinearElasticTransverseIsotropic.hpp"
#include "solid/PoroElastic.hpp"

namespace geosx
{
namespace constitutive
{



template< typename BASETYPE >
struct ConstitutivePassThru;

//template< template< typename BASE > class DERIVED >
//struct ConstitutivePassThru;




template<>
struct ConstitutivePassThru< SolidBase >
{
  template< typename LAMBDA >
  static
  GEOSX_FORCE_INLINE
  void Execute( ConstitutiveBase * const constitutiveRelation,
                LAMBDA && lambda )
  {
    if( dynamic_cast< LinearElasticIsotropic * >( constitutiveRelation ) )
    {
      lambda( static_cast< LinearElasticIsotropic * >( constitutiveRelation) );
    }
    else if( dynamic_cast< LinearElasticTransverseIsotropic * >( constitutiveRelation ) )
    {
      lambda( static_cast< LinearElasticTransverseIsotropic * >( constitutiveRelation) );
    }
    else if( dynamic_cast< LinearElasticAnisotropic * >( constitutiveRelation ) )
    {
      lambda( static_cast< LinearElasticAnisotropic * >( constitutiveRelation) );
    }
    else
    {
      string name;
      if( constitutiveRelation !=nullptr )
      {
        name = constitutiveRelation->getName();
      }
      GEOSX_ERROR( "ConstitutivePassThru<SolidBase>::Execute( "<<
                   constitutiveRelation<<" ) failed. ( "<<
                   constitutiveRelation<<" ) is named "<<name );
    }
  }
};



template<>
struct ConstitutivePassThru< NullModel >
{
  template< typename LAMBDA >
  static
  GEOSX_FORCE_INLINE
  void Execute( ConstitutiveBase * const constitutiveRelation,
                LAMBDA && lambda )
  {
    if( dynamic_cast< NullModel * >( constitutiveRelation ) )
    {
      lambda( static_cast< NullModel * >( constitutiveRelation ) );
    }
    else
    {
//      lambda( constitutiveRelation );

      string name;
      if( constitutiveRelation !=nullptr )
      {
        name = constitutiveRelation->getName();
      }
      GEOSX_ERROR( "ConstitutivePassThru<NullModel>::Execute( "<<
                   constitutiveRelation<<" ) failed. ( "<<
                   constitutiveRelation<<" ) is named "<<name );

    }
  }
};



template<>
struct ConstitutivePassThru< PoroElasticBase >
{
  template< typename LAMBDA >
  static
  GEOSX_FORCE_INLINE
  void Execute( ConstitutiveBase * const constitutiveRelation,
                LAMBDA && lambda )
  {
    if( dynamic_cast< PoroElastic< LinearElasticIsotropic > * >( constitutiveRelation ) )
    {
      lambda( static_cast< PoroElastic< LinearElasticIsotropic > * >( constitutiveRelation) );
    }
    else if( dynamic_cast< PoroElastic< LinearElasticTransverseIsotropic >  * >( constitutiveRelation ) )
    {
      lambda( static_cast< PoroElastic< LinearElasticTransverseIsotropic >  * >( constitutiveRelation) );
    }
    else if( dynamic_cast< PoroElastic< LinearElasticAnisotropic > * >( constitutiveRelation ) )
    {
      lambda( static_cast< PoroElastic< LinearElasticAnisotropic > * >( constitutiveRelation) );
    }
    else
    {
      string name;
      if( constitutiveRelation !=nullptr )
      {
        name = constitutiveRelation->getName();
      }
      GEOSX_ERROR( "ConstitutivePassThru<SolidBase>::Execute( "<<
                   constitutiveRelation<<" ) failed. ( "<<
                   constitutiveRelation<<" ) is named "<<name );
    }
  }
};


}
}

#endif /* GEOSX_CONSTITUTIVE_SOLID_SOLIDSELECTOR_HPP_ */