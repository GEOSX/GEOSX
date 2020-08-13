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
 *  @file LinearElasticTransverseIsotropic.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_LINEARELASTICTRANSVERSEISOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_LINEARELASTICTRANSVERSEISOTROPIC_HPP_
#include "SolidBase.hpp"
#include "constitutive/ExponentialRelation.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{
namespace constitutive
{
/**
 * @class LinearElasticTransverseIsotropicUpdates
 *
 * Class to provide linear elastic transverse isotropic material updates that
 * may be called from a kernel function.
 *
 * @note The "transverse" directions are 1 and 2 (or 0 and 1 in C-index)
 */
class LinearElasticTransverseIsotropicUpdates : public SolidBaseUpdates
{
public:
  /**
   * @brief Constructor
   * @param[in] c11 The 11 component of the Voigt stiffness tensor.
   * @param[in] c13 The 13 component of the Voigt stiffness tensor.
   * @param[in] c33 The 33 component of the Voigt stiffness tensor.
   * @param[in] c44 The 44 component of the Voigt stiffness tensor.
   * @param[in] c66 The 66 component of the Voigt stiffness tensor.
   * @param[in] stress The ArrayView holding the stress data for each quadrature
   *                   point.
   */
  LinearElasticTransverseIsotropicUpdates(
    arrayView1d< real64 const > const & c11,
    arrayView1d< real64 const > const & c13,
    arrayView1d< real64 const > const & c33,
    arrayView1d< real64 const > const & c44,
    arrayView1d< real64 const > const & c66,
    arrayView3d< real64, solid::STRESS_USD > const & stress ) :
    SolidBaseUpdates( stress ),
    m_c11( c11 ),
    m_c13( c13 ),
    m_c33( c33 ),
    m_c44( c44 ),
    m_c66( c66 )
  {}

  /// Deleted default constructor
  LinearElasticTransverseIsotropicUpdates() = delete;

  /// Default copy constructor
  LinearElasticTransverseIsotropicUpdates(
    LinearElasticTransverseIsotropicUpdates const & ) = default;

  /// Default move constructor
  LinearElasticTransverseIsotropicUpdates(
    LinearElasticTransverseIsotropicUpdates && ) = default;

  /// Deleted copy assignment operator
  LinearElasticTransverseIsotropicUpdates &
  operator=(
    LinearElasticTransverseIsotropicUpdates const & ) = delete;

  /// Deleted move assignment operator
  LinearElasticTransverseIsotropicUpdates &
  operator=(
    LinearElasticTransverseIsotropicUpdates && ) = delete;

  GEOSX_HOST_DEVICE
  virtual void
  SmallStrainNoState( localIndex const k,
                      real64 const ( &voigtStrain )[6],
                      real64 ( &stress )[6] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void
  SmallStrain( localIndex const k,
               localIndex const q,
               real64 const ( &voigtStrainInc )[6] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void
  HypoElastic( localIndex const k,
               localIndex const q,
               real64 const ( &Ddt )[6],
               real64 const ( &Rot )[3][3] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void
  HyperElastic( localIndex const k,
                real64 const ( &FmI )[3][3],
                real64 ( &stress )[6] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void
  HyperElastic( localIndex const k,
                localIndex const q,
                real64 const ( &FmI )[3][3] ) const override final;

  GEOSX_HOST_DEVICE inline virtual void
  GetStiffness(
    localIndex const k,
    real64 ( &c )[6][6] ) const override final
  {
    memset( c, 0, sizeof( c ) );
    c[0][0] = m_c11[k];
    c[0][1] = m_c11[k] - 2 * m_c66[k];
    c[0][2] = m_c13[k];
    c[1][0] = c[0][1];
    c[1][1] = m_c11[k];
    c[1][2] = m_c13[k];
    c[2][0] = c[0][2];
    c[2][1] = c[1][2];
    c[2][2] = m_c33[k];
    c[3][3] = m_c44[k];
    c[4][4] = m_c44[k];
    c[5][5] = m_c66[k];
  }

private:
  /// A reference to the ArrayView holding c11 for each element.
  arrayView1d< real64 const > const m_c11;

  /// A reference to the ArrayView holding c13 for each element.
  arrayView1d< real64 const > const m_c13;

  /// A reference to the ArrayView holding c33 for each element.
  arrayView1d< real64 const > const m_c33;

  /// A reference to the ArrayView holding c44 for each element.
  arrayView1d< real64 const > const m_c44;

  /// A reference to the ArrayView holding c66 for each element.
  arrayView1d< real64 const > const m_c66;
};

GEOSX_FORCE_INLINE
GEOSX_HOST_DEVICE
void
LinearElasticTransverseIsotropicUpdates::SmallStrainNoState(
  localIndex const k,
  real64 const ( &voigtStrain )[6],
  real64 ( &stress )[6] ) const
{
  real64 const c12temp = ( m_c11[k] - 2.0 * m_c66[k] );
  stress[0] = m_c11[k] * voigtStrain[0] + c12temp * voigtStrain[1] +
    m_c13[k] * voigtStrain[2];
  stress[1] = c12temp * voigtStrain[0] + m_c11[k] * voigtStrain[1] +
    m_c13[k] * voigtStrain[2];
  stress[2] = m_c13[k] * voigtStrain[0] + m_c13[k] * voigtStrain[1] +
    m_c33[k] * voigtStrain[2];

  stress[3] = m_c44[k] * voigtStrain[3];
  stress[4] = m_c44[k] * voigtStrain[4];
  stress[5] = m_c66[k] * voigtStrain[5];
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
LinearElasticTransverseIsotropicUpdates::SmallStrain(
  localIndex const k,
  localIndex const q,
  real64 const ( &voigtStrainInc )[6] ) const
{
  real64 const temp = m_c11[k] * ( voigtStrainInc[0] + voigtStrainInc[1] ) +
    m_c13[k] * voigtStrainInc[2];
  m_stress( k, q, 0 ) += -2.0 * m_c66[k] * voigtStrainInc[1] + temp;
  m_stress( k, q, 1 ) += -2.0 * m_c66[k] * voigtStrainInc[0] + temp;
  m_stress( k, q, 2 ) = m_stress( k, q, 2 ) +
    m_c13[k] * ( voigtStrainInc[0] + voigtStrainInc[1] ) +
    m_c33[k] * voigtStrainInc[2];
  m_stress( k, q, 3 ) = m_stress( k, q, 3 ) + m_c44[k] * voigtStrainInc[3];
  m_stress( k, q, 4 ) = m_stress( k, q, 4 ) + m_c44[k] * voigtStrainInc[4];
  m_stress( k, q, 5 ) = m_stress( k, q, 5 ) + m_c66[k] * voigtStrainInc[5];
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
LinearElasticTransverseIsotropicUpdates::HypoElastic(
  localIndex const k,
  localIndex const q,
  real64 const ( &Ddt )[6],
  real64 const ( &Rot )[3][3] ) const
{
  SmallStrain( k, q, Ddt );
  real64 temp[6];
  LvArray::tensorOps::AikSymBklAjl< 3 >( temp, Rot, m_stress[k][q] );
  LvArray::tensorOps::copy< 6 >( m_stress[k][q], temp );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
LinearElasticTransverseIsotropicUpdates::HyperElastic(
  localIndex const GEOSX_UNUSED_PARAM( k ),
  real64 const ( &GEOSX_UNUSED_PARAM( FmI ) )[3][3],
  real64 ( & )[6] ) const
{
  GEOSX_ERROR(
    "LinearElasticTransverseIsotropicKernelWrapper::HyperElastic() is not "
    "implemented!" );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
LinearElasticTransverseIsotropicUpdates::HyperElastic(
  localIndex const GEOSX_UNUSED_PARAM( k ),
  localIndex const GEOSX_UNUSED_PARAM( q ),
  real64 const ( &GEOSX_UNUSED_PARAM( FmI ) )[3][3] ) const
{
  GEOSX_ERROR(
    "LinearElasticTransverseIsotropicKernelWrapper::HyperElastic() is not "
    "implemented!" );
}

/**
 * @class LinearElasticTransverseIsotropic
 *
 * Class to provide a linear elastic transverse isotropic material response.
 */
class LinearElasticTransverseIsotropic : public SolidBase
{
public:
  /// @typedef Alias for LinearElasticTransverseIsotropicUpdates
  using KernelWrapper = LinearElasticTransverseIsotropicUpdates;

  /**
   * @brief constructor
   * @param[in]name name of the instance in the catalog
   * @param[in]parent the group which contains this instance
   */
  LinearElasticTransverseIsotropic( string const & name, Group * const parent );

  /**
   * Destructor
   */
  virtual ~LinearElasticTransverseIsotropic() override;

  virtual void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  virtual void
  AllocateConstitutiveData(
    dataRepository::Group * const parent,
    localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString =
    "LinearElasticTransverseIsotropic";

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

  /**
   * @struct Set of "char const *" and keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for transverse youngs modulus
    static constexpr auto defaultYoungsModulusTransverse =
      "defaultYoungsModulusTransverse";

    /// string/key for axial Young's modulus
    static constexpr auto defaultYoungsModulusAxial =
      "defaultYoungsModulusAxial";

    /// string/key for transverse Poisson's Ratio
    static constexpr auto defaultPoissonRatioTransverse =
      "defaultPoissonRatioTransverse";

    /// string/key for axial Poisson's Ratio
    static constexpr auto defaultPoissonRatioAxialTransverse =
      "defaultPoissonRatioAxialTransverse";

    /// string/key for transverse shear modulus
    static constexpr auto defaultShearModulusAxialTransverse =
      "defaultShearModulusAxialTransverse";

    /// string/key for c11 component of Voigt stiffness tensor
    static constexpr auto c11 = "c11";

    /// string/key for c13 component of Voigt stiffness tensor
    static constexpr auto c13 = "c13";

    /// string/key for c33 component of Voigt stiffness tensor
    static constexpr auto c33 = "c33";

    /// string/key for c44 component of Voigt stiffness tensor
    static constexpr auto c44 = "c44";

    /// string/key for c66 component of Voigt stiffness tensor
    static constexpr auto c66 = "c66";
  };

  /**
   * @brief Getter for default transverse Young's modulus
   * @return The value of the default transverse Young's modulus.
   */
  real64
  getDefaultYoungsModulusTransverse()
  {
    return m_defaultYoungsModulusTransverse;
  }

  /**
   * @brief Setter for the default transverse Young's modulus.
   * @param[in] input New value for the default transverse Young's modulus
   */
  void
  setDefaultYoungsModulusTransverse( real64 const input )
  {
    m_defaultYoungsModulusTransverse = input;
  }

  /**
   * @brief Getter for default axial Young's modulus
   * @return The value of the default axial Young's modulus.
   */
  real64
  getDefaultYoungsModulusAxial()
  {
    return m_defaultYoungsModulusAxial;
  }

  /**
   * @brief Setter for the default axial Young's modulus.
   * @param[in] input New value for the default axial Young's modulus
   */
  void
  setDefaultYoungsModulusAxial( real64 const input )
  {
    m_defaultYoungsModulusAxial = input;
  }

  /**
   * @brief Getter for default transverse Poisson's ratio
   * @return The value of the default transverse Poisson's ratio.
   */
  real64
  getDefaultPoissonsRatioTransverse()
  {
    return m_defaultPoissonTransverse;
  }

  /**
   * @brief Setter for the default transverse Poisson's ratio.
   * @param[in] input New value for the default transverse Poisson's ratio
   */
  void
  setDefaultPoissonsRatioTransverse( real64 const input )
  {
    m_defaultPoissonTransverse = input;
  }

  /**
   * @brief Getter for default axial Poisson's ratio
   * @return The value of the default axial/transverse Poisson's modulus.
   */
  real64
  getDefaultPoissonsRatioAxialTransverse()
  {
    return m_defaultPoissonAxialTransverse;
  }

  /**
   * @brief Setter for the default axial Poisson's modulus.
   * @param[in] input New value for the default axial/transverse Poisson's
   *             modulus
   */
  void
  setDefaultPoissonsRatioAxialTransverse( real64 const input )
  {
    m_defaultPoissonAxialTransverse = input;
  }

  /**
   * @brief Getter for default axial/transverse Shear modulus
   * @return The value of the default axial/transverse Shear modulus.
   */
  real64
  getDefaultShearModulusAxialTransverse()
  {
    return m_defaultShearModulusAxialTransverse;
  }

  /**
   * @brief Setter for the default axial/transverse Shear modulus.
   * @param[in] input New value for the default axial/transverse Shear modulus
   */
  void
  setDefaultShearModulusAxialTransverse( real64 const input )
  {
    m_defaultShearModulusAxialTransverse = input;
  }

  /**
   * @brief Const-Getter for 11 component of Voigt stiffness tensor.
   * @return reference to immutable 11 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > const &
  getC11() const
  {
    return m_c11;
  }

  /**
   * @brief Getter for 11 component of Voigt stiffness tensor.
   * @return reference to mutable 11 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > const &
  getC11()
  {
    return m_c11;
  }

  /**
   * @brief Const-Getter for 13 component of Voigt stiffness tensor.
   * @return reference to immutable 13 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > const &
  getC13() const
  {
    return m_c13;
  }

  /**
   * @brief Getter for 13 component of Voigt stiffness tensor.
   * @return reference to mutable 13 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > const &
  getC13()
  {
    return m_c13;
  }

  /**
   * @brief Const-Getter for 33 component of Voigt stiffness tensor.
   * @return reference to immutable 33 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > const &
  getC33() const
  {
    return m_c33;
  }

  /**
   * @brief Getter for 33 component of Voigt stiffness tensor.
   * @return reference to mutable 33 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > const &
  getC33()
  {
    return m_c33;
  }

  /**
   * @brief Const-Getter for 44 component of Voigt stiffness tensor.
   * @return reference to immutable 44 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > const &
  getC44() const
  {
    return m_c44;
  }

  /**
   * @brief Getter for 44 component of Voigt stiffness tensor.
   * @return reference to mutable 44 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > const &
  getC44()
  {
    return m_c44;
  }

  /**
   * @brief Const-Getter for 66 component of Voigt stiffness tensor.
   * @return reference to immutable 66 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > const &
  getC66() const
  {
    return m_c66;
  }

  /**
   * @brief Getter for 66 component of Voigt stiffness tensor.
   * @return reference to mutable 66 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > const &
  getC66()
  {
    return m_c66;
  }

  /**
   * @brief Create a instantiation of the
   *        LinearElasticTransverseIsotropicUpdates class that refers to the
   *        data in this.
   * @return An instantiation of LinearElasticTransverseIsotropicUpdates.
   */
  LinearElasticTransverseIsotropicUpdates
  createKernelUpdates()
  {
    return LinearElasticTransverseIsotropicUpdates( m_c11,
                                                    m_c13,
                                                    m_c33,
                                                    m_c44,
                                                    m_c66,
                                                    m_stress );
  }

  /**
   * @brief Construct an update kernel for a derived type.
   * @tparam UPDATE_KERNEL The type of update kernel from the derived type.
   * @tparam PARAMS The parameter pack to hold the constructor parameters for
   *   the derived update kernel.
   * @param constructorParams The constructor parameter for the derived type.
   * @return An @p UPDATE_KERNEL object.
   */
  template< typename UPDATE_KERNEL, typename... PARAMS >
  UPDATE_KERNEL
  createDerivedKernelUpdates( PARAMS &&... constructorParams )
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_c11,
                          m_c13,
                          m_c33,
                          m_c44,
                          m_c66,
                          m_stress );
  }

protected:
  virtual void
  PostProcessInput() override;

private:
  /// The default value of the transverse Young's modulus for any new
  /// allocations.
  real64 m_defaultYoungsModulusTransverse;

  /// The default value of the axial Young's modulus for any new
  /// allocations.
  real64 m_defaultYoungsModulusAxial;

  /// The default value of the transverse Poisson's ratio for any new
  /// allocations.
  real64 m_defaultPoissonTransverse;

  /// The default value of the axial/transverse Poisson's ratio for any new
  /// allocations.
  real64 m_defaultPoissonAxialTransverse;

  /// The default value of the axial/transverse Shear modulus for any new
  /// allocations.
  real64 m_defaultShearModulusAxialTransverse;

  /// The 11 component of the Voigt stiffness tensor.
  array1d< real64 > m_c11;

  /// The 13 component of the Voigt stiffness tensor.
  array1d< real64 > m_c13;

  /// The 33 component of the Voigt stiffness tensor.
  array1d< real64 > m_c33;

  /// The 44 component of the Voigt stiffness tensor.
  array1d< real64 > m_c44;

  /// The 66 component of the Voigt stiffness tensor.
  array1d< real64 > m_c66;
};
}  // namespace constitutive

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_LINEARELASTICTRANSVERSEISOTROPIC_HPP_ */
