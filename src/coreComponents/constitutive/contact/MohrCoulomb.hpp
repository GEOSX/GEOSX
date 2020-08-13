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
 *  @file MohrCoulomb.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CONTACTRELATIONS_MOHRCOULOMB_HPP_
#define GEOSX_CONSTITUTIVE_CONTACTRELATIONS_MOHRCOULOMB_HPP_
#include "ContactRelationBase.hpp"

namespace geosx
{
namespace constitutive
{
/**
 * @class LinearElasticIsotropic
 *
 * Class to provide a linear elastic isotropic material response.
 */
class MohrCoulomb : public ContactRelationBase
{
public:
  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  MohrCoulomb( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~MohrCoulomb() override;

  virtual void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  virtual real64
  limitTangentialTractionNorm(
    real64 const normalTraction ) const override final;

  virtual real64
  dLimitTangentialTractionNorm_dNormalTraction(
    real64 const normalTraction ) const override final;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "MohrCoulomb";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string
  CatalogName()
  {
    return m_catalogNameString;
  }

  virtual string
  GetCatalogName() override
  {
    return CatalogName();
  }

  ///@}

  bool m_postProcessed = false;

  /**
   * @struct Set of "char const *" and keys for data specified in this class.
   */
  struct viewKeyStruct : public ContactRelationBase::viewKeyStruct
  {
    /// string/key for default cohesion
    static constexpr auto defaultCohesionString = "defaultCohesion";

    /// string/key for default friction angle
    static constexpr auto defaultFrictionAngleString = "defaultFrictionAngle";

    /// string/key for cohesion
    static constexpr auto cohesionString = "cohesion";

    /// string/key for friction angle input (in radians)
    static constexpr auto frictionAngleString = "frictionAngle";

    /// string/key for friction coefficient
    static constexpr auto frictionCoefficientString = "frictionCoefficient";
  };

  /**
   * @brief Accessor for cohesion
   * @return A const reference to arrayView1d<real64> containing the cohesions
   *         (at every element).
   */
  real64 const &
  cohesion()
  {
    return m_cohesion;
  }

  /**
   * @brief Const accessor for cohesion
   * @return A const reference to arrayView1d<real64 const> containing the
   *         cohesions (at every element).
   */
  real64 const &
  cohesion() const
  {
    return m_cohesion;
  }

  /**
   * @brief Accessor for friction angle
   * @return A const reference to arrayView1d<real64> containing the friction
   *         coefficient (at every element).
   */
  real64 const &
  frictionCoefficient()
  {
    return m_frictionCoefficient;
  }

  /**
   * @brief Const accessor for friction angle
   * @return A const reference to arrayView1d<real64 const> containing the
   *         friction coefficient (at every element).
   */
  real64 const &
  frictionCoefficient() const
  {
    return m_frictionCoefficient;
  }

protected:
  virtual void
  PostProcessInput() override;

private:
  /// The cohesion for each upper level dimension (i.e. cell) of *this
  real64 m_cohesion;

  /// The friction angle for each upper level dimension (i.e. cell) of *this
  real64 m_frictionAngle;

  /// The friction coefficient for each upper level dimension (i.e. cell) of *this
  real64 m_frictionCoefficient;
};

}  // namespace constitutive

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_CONTACTRELATIONS_MOHRCOULOMB_HPP_ */
