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
 * @file FieldSpecificationManager.hpp
 */

#ifndef GEOSX_MANAGERS_FIELDSPECIFICATION_FIELDSPECIFICATIONMANAGER_HPP_
#define GEOSX_MANAGERS_FIELDSPECIFICATION_FIELDSPECIFICATIONMANAGER_HPP_

#include "codingUtilities/StringUtilities.hpp"
#include "managers/FieldSpecification/FieldSpecificationBase.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
/**
 * @brief The key for BoundaryConditionManager
 * @return the key
 */
string const boundaryConditionManager( "BoundaryConditionManager" );
}  // namespace keys
}  // namespace dataRepository

/**
 * @class FieldSpecificationManager
 * This class contains the field objects and provides an interface for administering
 * to specify them The class is a singleton.
 */
class FieldSpecificationManager : public dataRepository::Group
{
public:
  /**
   * @brief Singleton getter returns a pointer to the Singleton instance of
   *        FieldSpecificationManager.
   * @return a pointer to the singleton FieldSpecificationManager
   */
  static FieldSpecificationManager &
  get();

  /**
   * @brief Create a new FieldSpecificationBase object as a child of this group.
   * @param childKey the catalog key of the new FieldSpecificationBase derived type to create
   * @param childName the name of the new FieldSpecificationBase object in the repository
   * @return the group child
   */
  virtual Group *
  CreateChild( string const & childKey,
               string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void
  ExpandObjectCatalogs() override;

  /**
   * @brief Function to apply a value directly to a field variable.
   * @tparam LAMBDA The type of the lambda function
   * @param time The time at which the value will be evaluated. For instance if the
   *             field is a time dependent function, this is the evaluation time.
   * @param domain The DomainPartition object.
   * @param fieldPath The path to the object that contains the variable described in fieldName. This
   *                  path need not be the complete path, but rather the check that is performed is
   *                  that the fieldPath specified is contained in the string that specifies the
   *                  path from the actual field specification. In other words, the fieldPath variable
   *                  can be a substring of the path specified, and that will be
   *                  sufficient to proceed with the application of the.
   * @param fieldName The name of the field/variable that the value will be applied to.
   *                  It may not be necessary that this name is in the data repository, as the user
   *                  supplied lambda may apply whatever it condition it would like. However, this
   *                  name is used for comparing against the value given in the field specification.
   *
   * Interface function that applies a boundary condition directly to a field variable by looping
   * through all FieldSpecificationBase objects present in the FieldSpecificationManager. Searches
   * for the string specified in fieldPath as a substring in the objectPath specified in the input
   * file, checks if fieldName is equal to the fieldName specified in the input file, and check if
   * the time parameter falls within the beginTime and endTime of the FieldSpecificationBase object,
   * and calls FieldSpecificationBase::ApplyFieldValue().
   *
   */
  template< typename POLICY = parallelHostPolicy >
  void
  ApplyFieldValue( real64 const time,
                   dataRepository::Group * domain,
                   string const & fieldPath,
                   string const & fieldName ) const
  {
    GEOSX_MARK_FUNCTION;

    ApplyFieldValue< POLICY >( time,
                               domain,
                               fieldPath,
                               fieldName,
                               [&]( FieldSpecificationBase const * const,
                                    SortedArrayView< localIndex const > const & ) {} );
  }

  /**
   * @brief Function to apply a value directly to a field variable and applies a lambda
   *        for any post operations that are needed.
   * @tparam LAMBDA The type of the lambda function
   * @param time The time at which the field will be evaluated. For instance if the
   *             field is a time dependent function, this is the evaluation time.
   * @param domain The DomainPartition object.
   * @param fieldPath The path to the object that contains the variable described in fieldName. This
   *                  path need not be the complete path, but rather the check that is performed is
   *                  that the fieldPath specified is contained in the string that specifies the
   *                  path from the actual field specification. In other words, the fieldPath variable
   *                  can be a substring of the path specified , and that will be
   *                  sufficient to proceed with the application of the value.
   * @param fieldName The name of the field/variable that the value will be applied to.
   *                  It may not be necessary that this name is in the data repository, as the user
   *                  supplied lambda may apply whatever it condition it would like. However, this
   *                  name is used for comparing against the value given in the field specification.
   * @param lambda A lambda function that defines any operations that should be performed
   *               after application of value to the field.
   *
   * Interface function that applies a value directly to a field by looping
   * through all FieldSpecificationBase objects present in the FieldSpecificationManager. Searches
   * for the string specified in fieldPath as a substring in the objectPath specified in the input
   * file, checks if fieldName is equal to the fieldName specified in the input file, and check if
   * the time parameter falls within the beginTime and endTime of the FieldSpecificationBase object,
   * and calls FieldSpecificationBase::ApplyFieldValue(), and calls the lambda function
   * to apply any operations required for completing the application of the value to the field in addition to
   * setting the target field.
   */
  template< typename POLICY = parallelHostPolicy, typename LAMBDA = void >
  void
  ApplyFieldValue( real64 const time,
                   dataRepository::Group * domain,
                   string const & fieldPath,
                   string const & fieldName,
                   LAMBDA && lambda ) const;

  /**
   * @brief Function to apply a value directly to a field variable and applies a lambda
   *        for any post operations that are needed.
   * @tparam PRELAMBDA The type of the lambda function to be called before the ApplyFieldValue
   * @tparam POSTLAMBDA The type of the lambda function to be called after the ApplyFieldValue
   * @param time The time at which the field will be evaluated. For instance if the
   *             field is a time dependent function, this is the evaluation time.
   * @param domain The DomainPartition object.
   * @param fieldPath The path to the object that contains the variable described in fieldName. This
   *                  path need not be the complete path, but rather the check that is performed is
   *                  that the fieldPath specified is contained in the string that specifies the
   *                  path from the actual field specification. In other words, the fieldPath variable
   *                  can be a substring of the path specified , and that will be
   *                  sufficient to proceed with the application of the value.
   * @param fieldName The name of the field/variable that the value will be applied to.
   *                  It may not be necessary that this name is in the data repository, as the user
   *                  supplied lambda may apply whatever it condition it would like. However, this
   *                  name is used for comparing against the value given in the field specification.
   * @param preLambda A lambda function that defines any operations that should be performed
   *                  before application of value to the field.
   * @param postLambda A lambda function that defines any operations that should be performed
   *                   after application of value to the field.
   *
   * Interface function that applies a value directly to a field by looping
   * through all FieldSpecificationBase objects present in the FieldSpecificationManager. Calls the
   * preLambda function to apply any operations required before the application of the value to the
   * field in addition to setting the target field. Searches for the string specified in fieldPath
   * as a substring in the objectPath specified in the input file, checks if fieldName is equal to
   * the fieldName specified in the input file, and check if the time parameter falls within the
   * beginTime and endTime of the FieldSpecificationBase object, and calls
   * FieldSpecificationBase::ApplyFieldValue(), and calls the postLambda function to apply any
   * operations required for completing the application of the value to the field in addition to
   * setting the target field.
   */
  template< typename POLICY = parallelHostPolicy,
            typename PRELAMBDA = void,
            typename POSTLAMBDA = void >
  void
  ApplyFieldValue( real64 const time,
                   dataRepository::Group * domain,
                   string const & fieldPath,
                   string const & fieldName,
                   PRELAMBDA && preLambda,
                   POSTLAMBDA && postLambda ) const;

  /**
   * @brief function to apply initial conditions
   * @param domain the DomainParition object
   */
  void
  ApplyInitialConditions( dataRepository::Group * domain ) const;

  /**
   * @brief This function is the main driver for the field applications
   * @tparam LAMBDA The type of the lambda function
   * @param time The time at which the field will be evaluated. For instance if the
   *             field is a time dependent function, this is the evaluation time.
   * @param domain The DomainPartition object.
   * @param fieldPath The path to the object that contains the variable described in fieldName. This
   *                  path need not be the complete path, but rather the check that is performed is
   *                  that the fieldPath specified is contained in the string that specifies the
   *                  path from the actual field specification. In other words, the fieldPath variable
   *                  can be a substring of the path specified, and that will be
   *                  sufficient to proceed with the application of the value to the field.
   * @param fieldName The name of the field/variable that the value will be applied to.
   *                  It may not be necessary that this name is in the data repository, as the user
   *                  supplied lambda may apply whatever it condition it would like. However, this
   *                  name is used for comparing against the value given in the specification.
   * @param lambda A lambda function that defines the application of the field.
   *
   * This function loops through all available fields, checks to see if they
   * should be applied, and applies them. More specifically, this function simply checks
   * values of fieldPath,fieldName, against each FieldSpecificationBase object contained in the
   * FieldSpecificationManager and decides on whether or not to call the user defined lambda.
   */
  template< typename LAMBDA >
  void
  Apply( real64 const time,
         dataRepository::Group * domain,
         string const & fieldPath,
         string const & fieldName,
         LAMBDA && lambda ) const
  {
    GEOSX_MARK_FUNCTION;
    // loop over all FieldSpecificationBase objects
    for( auto & subGroup : this->GetSubGroups() )
    {
      FieldSpecificationBase const * fs =
        subGroup.second->group_cast< FieldSpecificationBase const * >();
      int const isInitialCondition = fs->initialCondition();

      if( ( isInitialCondition && fieldPath == "" ) ||
          ( !isInitialCondition &&
            fs->GetObjectPath().find( fieldPath ) != string::npos ) )
      {
        string_array const targetPath =
          stringutilities::Tokenize( fs->GetObjectPath(), "/" );
        localIndex const targetPathLength =
          LvArray::integerConversion< localIndex >( targetPath.size() );
        string const targetName = fs->GetFieldName();

        if( ( isInitialCondition && fieldName == "" ) ||
            ( !isInitialCondition && time >= fs->GetStartTime() &&
              time < fs->GetEndTime() && targetName == fieldName ) )
        {
          MeshLevel * const meshLevel =
            domain->group_cast< DomainPartition * >()->getMeshBody( 0 )->getMeshLevel(
              0 );

          dataRepository::Group * targetGroup = meshLevel;

          string processedPath;
          for( localIndex pathLevel = 0; pathLevel < targetPathLength; ++pathLevel )
          {
            dataRepository::Group * const elemRegionSubGroup =
              targetGroup->GetGroup(
                ElementRegionManager::groupKeyStruct::elementRegionsGroup );
            if( elemRegionSubGroup != nullptr )
            {
              targetGroup = elemRegionSubGroup;
            }

            dataRepository::Group * const elemSubRegionSubGroup =
              targetGroup->GetGroup(
                ElementRegionBase::viewKeyStruct::elementSubRegions );
            if( elemSubRegionSubGroup != nullptr )
            {
              targetGroup = elemSubRegionSubGroup;
            }

            if( targetPath[pathLevel] ==
                  ElementRegionManager::groupKeyStruct::elementRegionsGroup ||
                targetPath[pathLevel] ==
                  ElementRegionBase::viewKeyStruct::elementSubRegions )
            {
              continue;
            }

            targetGroup = targetGroup->GetGroup( targetPath[pathLevel] );
            processedPath += "/" + targetPath[pathLevel];

            GEOSX_ERROR_IF(
              targetGroup == nullptr,
              "ApplyBoundaryCondition(): Last entry in objectPath ("
                << processedPath << ") is not found" );
          }
          ApplyOnTargetRecursive( targetGroup, fs, targetName, lambda );
        }
      }
    }
  }

private:
  /**
   * @brief private constructor for the singleton BoundaryConditionManager.
   * @param name The name of the BoundaryConditionManager in the data repository.
   * @param parent The parent of BoundaryConditionManager in the data repository.
   */
  FieldSpecificationManager( string const & name,
                             dataRepository::Group * const parent );
  virtual ~FieldSpecificationManager() override;

  template< typename LAMBDA >
  void
  ApplyOnTargetRecursive( Group * target,
                          FieldSpecificationBase const * fs,
                          string const & targetName,
                          LAMBDA && lambda ) const
  {
    if( ( target->getParent()->getName() ==
            ElementRegionBase::viewKeyStruct::elementSubRegions ||
          target->getName() == "nodeManager" || target->getName() == "FaceManager" ||
          target->getName() ==
            "edgeManager" )  // TODO these 3 strings are harcoded because for the moment, there are
                             // inconsistencies with the name of the Managers...
        && target->getName() != ObjectManagerBase::groupKeyStruct::setsString &&
        target->getName() != ObjectManagerBase::groupKeyStruct::neighborDataString )
    {
      dataRepository::Group const * setGroup =
        target->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
      string_array setNames = fs->GetSetNames();
      for( auto & setName : setNames )
      {
        dataRepository::Wrapper< SortedArray< localIndex > > const * const setWrapper =
          setGroup->getWrapper< SortedArray< localIndex > >( setName );
        if( setWrapper != nullptr )
        {
          SortedArrayView< localIndex const > const & targetSet =
            setWrapper->reference();
          lambda( fs, setName, targetSet, target, targetName );
        }
      }
    }
    else
    {
      target->forSubGroups( [&]( Group & subTarget ) {
        ApplyOnTargetRecursive( &subTarget, fs, targetName, lambda );
      } );
    }
  }
};

template< typename POLICY, typename LAMBDA >
void
FieldSpecificationManager::ApplyFieldValue( real64 const time,
                                            dataRepository::Group * domain,
                                            string const & fieldPath,
                                            string const & fieldName,
                                            LAMBDA && lambda ) const
{
  GEOSX_MARK_FUNCTION;

  Apply( time,
         domain,
         fieldPath,
         fieldName,
         [&]( FieldSpecificationBase const * const fs,
              string const &,
              SortedArrayView< localIndex const > const & targetSet,
              Group * const targetGroup,
              string const & targetField ) {
           fs->ApplyFieldValue< FieldSpecificationEqual, POLICY >( targetSet,
                                                                   time,
                                                                   targetGroup,
                                                                   targetField );
           lambda( fs, targetSet );
         } );
}

template< typename POLICY, typename PRELAMBDA, typename POSTLAMBDA >
void
FieldSpecificationManager::ApplyFieldValue( real64 const time,
                                            dataRepository::Group * domain,
                                            string const & fieldPath,
                                            string const & fieldName,
                                            PRELAMBDA && preLambda,
                                            POSTLAMBDA && postLambda ) const
{
  GEOSX_MARK_FUNCTION;

  Apply( time,
         domain,
         fieldPath,
         fieldName,
         [&]( FieldSpecificationBase const * const fs,
              string const &,
              SortedArrayView< localIndex const > const & targetSet,
              Group * const targetGroup,
              string const & targetField ) {
           preLambda( fs, targetSet );
           fs->ApplyFieldValue< FieldSpecificationEqual, POLICY >( targetSet,
                                                                   time,
                                                                   targetGroup,
                                                                   targetField );
           postLambda( fs, targetSet );
         } );
}

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_FIELDSPECIFICATION_FIELDSPECIFICATIONMANAGER_HPP_ */
