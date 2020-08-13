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
 * @file ContactRelationBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CONTACTRELATIONS_CONTACTRELATIONBASE_HPP_
#define GEOSX_CONSTITUTIVE_CONTACTRELATIONS_CONTACTRELATIONBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"
#include "managers/Functions/FunctionBase.hpp"

namespace geosx
{
namespace constitutive
{
/**
 * @class ContactRelationBase
 *
 * This class will serve as the interface for implementing contact enforcement constitutive relations.
 * This does not include the actual enforcement algorithm, but only the constitutive relations that
 * govern the behavior of the contact. So things like penalty, or friction, or kinematic constraint.
 *
 * Currently this class actually implements something. It should be converted to a base class when
 * it grows up.
 */
class ContactRelationBase : public ConstitutiveBase
{
public:
  /**
   * @brief The standard data repository constructor
   * @param name The name of the relation in the data repository
   * @param parent The name of the parent Group that holds this relation object.
   */
  ContactRelationBase( string const & name, Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~ContactRelationBase() override;

  /**
   * @brief Name that is used to register this a type of "ContactRelationBase" in the object catalog
   * @return See description
   */
  static string
  CatalogName()
  {
    return "Contact";
  }

  virtual string
  GetCatalogName() override
  {
    return CatalogName();
  }

  virtual void
  DeliverClone(
    string const & GEOSX_UNUSED_PARAM( name ),
    Group * const GEOSX_UNUSED_PARAM( parent ),
    std::unique_ptr< ConstitutiveBase > & GEOSX_UNUSED_PARAM( clone ) ) const override
  {}

  virtual Group *
  CreateChild( string const & catalogKey,
               string const & name ) override;

  /**
   * This function is used to inform the schema generator
   * that table functions are allowed as children.
   */
  virtual void
  SetSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
                       xmlWrapper::xmlNode schemaParent,
                       integer documentationType ) override;

  virtual void
  InitializePreSubGroups( Group * const ) override;

  virtual real64
  limitTangentialTractionNorm( real64 const normalTraction ) const;

  virtual real64
  dLimitTangentialTractionNorm_dNormalTraction(
    real64 const normalTraction ) const;

  /// accessor for penalty stiffness
  inline real64
  stiffness() const
  {
    return m_penaltyStiffness;
  }

  /**
   * @brief evaluation of effective aperture given in input model aperture/gap
   * @param[in] aperture the model aperture/gap
   * @return And effective physical aperture that is always > 0
   */
  inline real64
  effectiveAperture( real64 const aperture ) const
  {
    return m_apertureFunction->Evaluate( &aperture );
  }

  /**
   * @brief evaluation of the derivative of the effective physical aperture
   * @param[in] aperture the model aperture/gap
   * @return
   */
  inline real64
  dEffectiveAperture_dAperture( real64 const aperture ) const
  {
    real64 aperPlus = aperture;
    real64 aperMinus = aperture - 1.0e-6;
    real64 slope = ( m_apertureFunction->Evaluate( &aperPlus ) -
                     m_apertureFunction->Evaluate( &aperMinus ) ) /
      1.0e-6;
    return slope;
  }

  inline real64
  apertureTolerance() const
  {
    return m_apertureTolerance;
  }

  /**
   * @struct Structure to hold scoped key names
   */
  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto penaltyStiffnessString = "penaltyStiffness";
    static constexpr auto apertureToleranceString = "apertureTolerance";
  };

private:
  /// the value of penalty to penetration
  real64 m_penaltyStiffness;

  /// pointer to the function that limits the model aperture to a physically admissible value.
  FunctionBase * m_apertureFunction;

  real64 m_apertureTolerance;
};

}  // namespace constitutive
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_CONTACTRELATIONS_CONTACTRELATIONBASE_HPP_ */
