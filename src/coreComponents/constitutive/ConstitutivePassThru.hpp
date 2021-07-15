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
 * @file ConstitutivePassThru.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CONSTITUTIVEPASSTHRU_HPP_
#define GEOSX_CONSTITUTIVE_CONSTITUTIVEPASSTHRU_HPP_


#include "ConstitutivePassThruHandler.hpp"
#include "NullModel.hpp"
#include "solid/DamageVolDev.hpp"
#include "solid/DamageSpectral.hpp"
#include "solid/DruckerPrager.hpp"
#include "solid/DruckerPragerExtended.hpp"
#include "solid/ElasticIsotropic.hpp"
#include "solid/ElasticTransverseIsotropic.hpp"
#include "solid/PorousSolid.hpp"
#include "solid/PoroElastic.hpp"


namespace geosx
{
namespace constitutive
{

/**
 * @struct ConstitutivePassThru
 * @brief Struct to facilitate launching of lambda functions with a compile
 *   time knowledge of what constitutive model is used.
 *
 * This struct works by implementing an if-else or switch-case block for a
 * specific constitutive base type, and executing the lambda passing it a
 * casted pointer to the constitutive relation.
 */
template< typename BASETYPE >
struct ConstitutivePassThru;

/**
 * Specialization for models that derive from SolidBase.
 */
template<>
struct ConstitutivePassThru< SolidBase >
{

  // NOTE: The switch order here can be fragile if a model derives from another
  //       model, as the dynamic_cast will also cast to a base version.
  //       Models should be ordered such that children come before parents.
  //       For example, DruckerPrager before ElasticIsotropic, DamageVolDev before
  //       Damage, etc.

  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< DamageSpectral< ElasticIsotropic >,
                                 DamageVolDev< ElasticIsotropic >,
                                 Damage< ElasticIsotropic >,
                                 DruckerPragerExtended,
                                 DruckerPrager,
                                 ElasticIsotropic,
                                 ElasticTransverseIsotropic >::execute( constitutiveRelation,
                                                                        std::forward< LAMBDA >( lambda ) );
  }
};


/**
 * Specialization for the PoroElastic models.
 */
template<>
struct ConstitutivePassThru< PoroElasticBase >
{
  template< typename LAMBDA >
  static void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< PoroElastic< DruckerPragerExtended >,
                                 PoroElastic< DruckerPrager >,
                                 PoroElastic< ElasticIsotropic >,
                                 PoroElastic< ElasticTransverseIsotropic > >::execute( constitutiveRelation,
                                                                                       std::forward< LAMBDA >( lambda ) );
  }
};



/**
 * Specialization for the NullModel.
 */
template<>
struct ConstitutivePassThru< NullModel >
{
  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    if( auto * const ptr = dynamic_cast< NullModel * >( &constitutiveRelation ) )
    {
      lambda( *ptr );
    }
    else
    {
      GEOSX_ERROR( "ConstitutivePassThru< NullModel >::execute failed. The constitutive relation is named "
                   << constitutiveRelation.getName() << " with type "
                   << LvArray::system::demangleType( constitutiveRelation ) );
    }
  }
};


/**
 * Specialization for the Damage models.
 */
template<>
struct ConstitutivePassThru< DamageBase >
{
  template< typename LAMBDA >
  static void execute( ConstitutiveBase & constitutiveRelation,
                       LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< DamageSpectral< ElasticIsotropic >,
                                 DamageVolDev< ElasticIsotropic >,
                                 Damage< ElasticIsotropic > >::execute( constitutiveRelation,
                                                                        std::forward< LAMBDA >( lambda ) );
  }
};


/**
 * Specialization for the PorousSolid models.
 */
template<>
struct ConstitutivePassThru< PorousSolidBase >
{
  template< typename LAMBDA >
  static void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< PorousSolid< DruckerPragerExtended >,
                                 PorousSolid< DruckerPrager >,
                                 PorousSolid< ElasticIsotropic >,
                                 PorousSolid< ElasticTransverseIsotropic >,
                                 PorousSolid< DamageSpectral< ElasticIsotropic > >,
                                 PorousSolid< DamageVolDev< ElasticIsotropic > >,
                                 PorousSolid< Damage< ElasticIsotropic > > >::execute( constitutiveRelation,
                                                                                       std::forward< LAMBDA >( lambda ) );
  }
};


}
}

#endif /* GEOSX_CONSTITUTIVE_CONSTITUTIVEPASSTHRU_HPP_ */
