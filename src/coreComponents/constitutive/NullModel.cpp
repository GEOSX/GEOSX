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

#include "NullModel.hpp"

namespace geosx
{
namespace constitutive
{
NullModel::NullModel( string const & name, Group * const parent ) :
  ConstitutiveBase( name, parent )
{}

NullModel::~NullModel()
{}

void
NullModel::DeliverClone( string const & name,
                         Group * const parent,
                         std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< NullModel >( name, parent );
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        NullModel,
                        std::string const &,
                        dataRepository::Group * const )

}  // namespace constitutive
} /* namespace geosx */
