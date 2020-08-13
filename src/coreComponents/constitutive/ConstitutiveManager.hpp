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
 * @file ConstitutiveManager.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_
#define GEOSX_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "dataRepository/ReferenceWrapper.hpp"
#include "ConstitutiveBase.hpp"

namespace geosx
{
namespace constitutive
{
/**
 * @typedef This is an alias for array of reference wrappers which are used to provide a
 * multidimensional array like interface via the operator[].
 */
template< typename VIEWTYPE >
using ViewAccessor = array1d< VIEWTYPE >;

/**
 * @class ConstitutiveManager
 * @brief Class to manage the allocation and access to constitutive relations
 */
class ConstitutiveManager : public dataRepository::Group
{
public:
  ConstitutiveManager() = delete;

  ConstitutiveManager( string const & name, Group * const parent );

  virtual Group *
  CreateChild( string const & childKey,
               string const & childName ) override final;

  /// This function is used to expand any catalogs in the data structure
  virtual void
  ExpandObjectCatalogs() override;

  ConstitutiveBase *
  HangConstitutiveRelation(
    string const & constitutiveRelationInstanceName,
    dataRepository::Group * const parent,
    localIndex const numConstitutivePointsPerParentIndex ) const;

  ~ConstitutiveManager() override;

  template< typename T = ConstitutiveBase >
  T const *
  GetConstitutiveRelation( string const & constitutiveRelationInstanceName ) const
  {
    return this->GetGroup< T >( constitutiveRelationInstanceName );
  }

  template< typename T = ConstitutiveBase >
  T *
  GetConstitutiveRelation( string const & constitutiveRelationInstanceName )
  {
    return this->GetGroup< T >( constitutiveRelationInstanceName );
  }

  template< typename T = ConstitutiveBase >
  T const *
  GetConstitutiveRelation( localIndex const index ) const
  {
    return this->GetGroup< T >( index );
  }

  template< typename T = ConstitutiveBase >
  T *
  GetConstitutiveRelation( localIndex const index )
  {
    return this->GetGroup< T >( index );
  }

  // template< typename T >
  // ViewAccessor< T >
  // GetConstitutiveData( string const & name,
  //                      dataRepository::Group * const relationGroup );

  template< typename T >
  ViewAccessor< T >
  GetConstitutiveData(
    string const & name,
    dataRepository::Group const * const relationGroup ) const;

  struct groupKeyStruct
  {
    static constexpr auto constitutiveModelsString = "ConstitutiveModels";
  } m_ConstitutiveManagerGroupKeys;
};

template< typename T >
ViewAccessor< T >
ConstitutiveManager::GetConstitutiveData(
  string const & name,
  dataRepository::Group const * const relationGroup ) const
{
  ViewAccessor< T const > rval( relationGroup->numSubGroups() );

  rval.resize( relationGroup->numSubGroups() );
  for( localIndex a = 0; a < this->GetSubGroups().size(); ++a )
  {
    ConstitutiveBase const * const material =
      relationGroup->GetGroup< ConstitutiveBase >( a );
    if( material->hasWrapper( name ) )
    {
      rval[a] = material->getReference< T >( name );
    }
  }
  return rval;
}

// template< typename T >
// ViewAccessor< T >
// ConstitutiveManager::GetConstitutiveData( string const & name,
//                                           dataRepository::Group * const relationGroup )
// {
//   return const_cast< ViewAccessor<T> >(const_cast<ConstitutiveManager const *>(this)->
//                                        GetConstitutiveData<T>( name, relationGroup ) );
// }

}  // namespace constitutive
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_ */
