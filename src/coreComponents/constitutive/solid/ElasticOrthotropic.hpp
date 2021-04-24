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
 *  @file ElasticOrthotropic.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_ELASTICORTHOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_ELASTICORTHOTROPIC_HPP_

#include "SolidBase.hpp"
#include "constitutive/ExponentialRelation.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "SolidModelDiscretizationOpsOrthotropic.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class ElasticOrthotropicUpdates
 *
 * Class to provide elastic orthotropic material updates that
 * may be called from a kernel function.
 */
class ElasticOrthotropicUpdates : public SolidBaseUpdates
{
public:
  /**
   * @brief Constructor
   * @param[in] c11 The 11 component of the Voigt stiffness tensor.
   * @param[in] c12 The 12 component of the Voigt stiffness tensor.
   * @param[in] c13 The 13 component of the Voigt stiffness tensor.
   * @param[in] c22 The 22 component of the Voigt stiffness tensor.
   * @param[in] c23 The 23 component of the Voigt stiffness tensor.
   * @param[in] c33 The 33 component of the Voigt stiffness tensor.
   * @param[in] c44 The 44 component of the Voigt stiffness tensor.
   * @param[in] c55 The 55 component of the Voigt stiffness tensor.
   * @param[in] c66 The 66 component of the Voigt stiffness tensor.
   * @param[in] stress The ArrayView holding the stress data for each quadrature
   *                   point.
   */
  ElasticOrthotropicUpdates( arrayView1d< real64 const > const & c11,
                             arrayView1d< real64 const > const & c12,
                             arrayView1d< real64 const > const & c13,
                             arrayView1d< real64 const > const & c22,
                             arrayView1d< real64 const > const & c23,
                             arrayView1d< real64 const > const & c33,
                             arrayView1d< real64 const > const & c44,
                             arrayView1d< real64 const > const & c55,
                             arrayView1d< real64 const > const & c66,
                             arrayView3d< real64, solid::STRESS_USD > const & newStress,
                             arrayView3d< real64, solid::STRESS_USD > const & oldStress ):
    SolidBaseUpdates( newStress, oldStress ),
    m_c11( c11 ),
    m_c12( c12 ),
    m_c13( c13 ),
    m_c22( c22 ),
    m_c23( c23 ),
    m_c33( c33 ),
    m_c44( c44 ),
    m_c55( c55 ),
    m_c66( c66 )
  {}

  /// Deleted default constructor
  ElasticOrthotropicUpdates() = delete;

  /// Default copy constructor
  ElasticOrthotropicUpdates( ElasticOrthotropicUpdates const & ) = default;

  /// Default move constructor
  ElasticOrthotropicUpdates( ElasticOrthotropicUpdates && ) = default;

  /// Deleted copy assignment operator
  ElasticOrthotropicUpdates & operator=( ElasticOrthotropicUpdates const & ) = delete;

  /// Deleted move assignment operator
  ElasticOrthotropicUpdates & operator=( ElasticOrthotropicUpdates && ) =  delete;

  // Use orthotropic form of inner product compression

  using DiscretizationOps = SolidModelDiscretizationOpsOrthotropic;

  // total strain interfaces

  GEOSX_HOST_DEVICE
  virtual void smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                    localIndex const q,
                                                    real64 const ( &totalStrain )[6],
                                                    real64 ( &stress )[6] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void smallStrainNoStateUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const ( &totalStrain )[6],
                                         real64 ( &stress )[6],
                                         real64 ( &stiffness )[6][6] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void smallStrainNoStateUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const ( &totalStrain )[6],
                                         real64 ( &stress )[6],
                                         DiscretizationOps & stiffness ) const final;

  // incremental strain interfaces

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate_StressOnly( localIndex const k,
                                             localIndex const q,
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( &stress )[6] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  real64 ( &stiffness )[6][6] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  DiscretizationOps & stiffness ) const final;

  // hypo interface

  GEOSX_HOST_DEVICE
  virtual void hypoUpdate( localIndex const k,
                           localIndex const q,
                           real64 const ( &Ddt )[6],
                           real64 const ( &Rot )[3][3],
                           real64 ( &stress )[6],
                           real64 ( &stiffness )[6][6] ) const override final;

  // miscellaneous getters

  GEOSX_HOST_DEVICE
  virtual void getElasticStiffness( localIndex const k, real64 ( &stiffness )[6][6] ) const override final;


private:
  /// A reference to the ArrayView holding c11 for each element.
  arrayView1d< real64 const > const m_c11;

  /// A reference to the ArrayView holding c12 for each element.
  arrayView1d< real64 const > const m_c12;

  /// A reference to the ArrayView holding c13 for each element.
  arrayView1d< real64 const > const m_c13;

  /// A reference to the ArrayView holding c22 for each element.
  arrayView1d< real64 const > const m_c22;

  /// A reference to the ArrayView holding c23 for each element.
  arrayView1d< real64 const > const m_c23;

  /// A reference to the ArrayView holding c33 for each element.
  arrayView1d< real64 const > const m_c33;

  /// A reference to the ArrayView holding c44 for each element.
  arrayView1d< real64 const > const m_c44;

  /// A reference to the ArrayView holding c55 for each element.
  arrayView1d< real64 const > const m_c55;

  /// A reference to the ArrayView holding c66 for each element.
  arrayView1d< real64 const > const m_c66;
};


GEOSX_FORCE_INLINE
GEOSX_HOST_DEVICE
void ElasticOrthotropicUpdates::getElasticStiffness( localIndex const k,
                                                     real64 ( & stiffness )[6][6] ) const
{
  LvArray::tensorOps::fill< 6, 6 >( stiffness, 0 );

  stiffness[0][0] = m_c11[k];
  stiffness[0][1] = m_c12[k];
  stiffness[0][2] = m_c13[k];

  stiffness[1][0] = stiffness[0][1];
  stiffness[1][1] = m_c22[k];
  stiffness[1][2] = m_c23[k];

  stiffness[2][0] = stiffness[0][2];
  stiffness[2][1] = stiffness[1][2];
  stiffness[2][2] = m_c33[k];

  stiffness[3][3] = m_c44[k];
  stiffness[4][4] = m_c55[k];
  stiffness[5][5] = m_c66[k];
}

GEOSX_FORCE_INLINE
GEOSX_HOST_DEVICE
void ElasticOrthotropicUpdates::smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                                     localIndex const q,
                                                                     real64 const ( &totalStrain )[6],
                                                                     real64 ( & stress )[6] ) const
{
  GEOSX_UNUSED_VAR( q );
  stress[0] = m_c11[k] * totalStrain[0] + m_c12[k] * totalStrain[1] + m_c13[k]*totalStrain[2];
  stress[1] = m_c12[k] * totalStrain[0] + m_c22[k] * totalStrain[1] + m_c23[k]*totalStrain[2];
  stress[2] = m_c13[k] * totalStrain[0] + m_c23[k] * totalStrain[1] + m_c33[k]*totalStrain[2];

  stress[3] = m_c44[k] * totalStrain[3];
  stress[4] = m_c55[k] * totalStrain[4];
  stress[5] = m_c66[k] * totalStrain[5];
}


GEOSX_FORCE_INLINE
GEOSX_HOST_DEVICE
void ElasticOrthotropicUpdates::smallStrainNoStateUpdate( localIndex const k,
                                                          localIndex const q,
                                                          real64 const ( &totalStrain )[6],
                                                          real64 ( & stress )[6],
                                                          real64 ( & stiffness )[6][6] ) const
{
  smallStrainNoStateUpdate_StressOnly( k, q, totalStrain, stress );
  getElasticStiffness( k, stiffness );
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ElasticOrthotropicUpdates::smallStrainNoStateUpdate( localIndex const k,
                                                          localIndex const q,
                                                          real64 const ( &totalStrain )[6],
                                                          real64 ( & stress )[6],
                                                          DiscretizationOps & stiffness ) const
{
  smallStrainNoStateUpdate_StressOnly( k, q, totalStrain, stress );
  stiffness.m_c11 = m_c11[k];
  stiffness.m_c12 = m_c12[k];
  stiffness.m_c13 = m_c13[k];
  stiffness.m_c22 = m_c22[k];
  stiffness.m_c23 = m_c23[k];
  stiffness.m_c33 = m_c33[k];
  stiffness.m_c44 = m_c44[k];
  stiffness.m_c55 = m_c55[k];
  stiffness.m_c66 = m_c66[k];
}


GEOSX_FORCE_INLINE
GEOSX_HOST_DEVICE
void ElasticOrthotropicUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                              localIndex const q,
                                                              real64 const ( &strainIncrement )[6],
                                                              real64 ( & stress )[6] ) const
{
  smallStrainNoStateUpdate_StressOnly( k, q, strainIncrement, stress ); // stress =  incrementalStress
  LvArray::tensorOps::add< 6 >( stress, m_oldStress[k][q] );            // stress += m_oldStress
  saveStress( k, q, stress );                                           // m_newStress = stress
}


GEOSX_FORCE_INLINE
GEOSX_HOST_DEVICE
void ElasticOrthotropicUpdates::smallStrainUpdate( localIndex const k,
                                                   localIndex const q,
                                                   real64 const ( &strainIncrement )[6],
                                                   real64 ( & stress )[6],
                                                   real64 ( & stiffness )[6][6] ) const
{
  smallStrainUpdate_StressOnly( k, q, strainIncrement, stress );
  getElasticStiffness( k, stiffness );
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ElasticOrthotropicUpdates::smallStrainUpdate( localIndex const k,
                                                   localIndex const q,
                                                   real64 const ( &strainIncrement )[6],
                                                   real64 ( & stress )[6],
                                                   DiscretizationOps & stiffness ) const
{
  smallStrainUpdate_StressOnly( k, q, strainIncrement, stress );
  stiffness.m_c11 = m_c11[k];
  stiffness.m_c12 = m_c12[k];
  stiffness.m_c13 = m_c13[k];
  stiffness.m_c22 = m_c22[k];
  stiffness.m_c23 = m_c23[k];
  stiffness.m_c33 = m_c33[k];
  stiffness.m_c44 = m_c44[k];
  stiffness.m_c55 = m_c55[k];
  stiffness.m_c66 = m_c66[k];
}


GEOSX_FORCE_INLINE
GEOSX_HOST_DEVICE
void ElasticOrthotropicUpdates::hypoUpdate( localIndex const k,
                                            localIndex const q,
                                            real64 const ( &Ddt )[6],
                                            real64 const ( &Rot )[3][3],
                                            real64 ( & stress )[6],
                                            real64 ( & stiffness )[6][6] ) const
{
  GEOSX_UNUSED_VAR( k );
  GEOSX_UNUSED_VAR( q );
  GEOSX_UNUSED_VAR( Ddt );
  GEOSX_UNUSED_VAR( Rot );
  GEOSX_UNUSED_VAR( stress );
  GEOSX_UNUSED_VAR( stiffness );
  GEOSX_ERROR( "hypoUpdate() disabled for anisotropic models when using Hughes-Winget integration" );
}


/**
 * @class ElasticOrthotropic
 *
 * Class to provide an elastic orthotropic material response.
 */
class ElasticOrthotropic : public SolidBase
{
public:

  /// @typedef Alias for ElasticOrthotropicUpdates
  using KernelWrapper = ElasticOrthotropicUpdates;

  /**
   * @brief constructor
   * @param[in]name name of the instance in the catalog
   * @param[in]parent the group which contains this instance
   */
  ElasticOrthotropic( string const & name, Group * const parent );

  /**
   * Destructor
   */
  virtual ~ElasticOrthotropic() override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "ElasticOrthotropic";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static string catalogName() { return m_catalogNameString; }

  virtual string getCatalogName() const override { return catalogName(); }
  ///@}

  /**
   * Keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for Young's modulus E1
    static constexpr char const * defaultE1String() { return "defaultE1"; };

    /// string/key for Young's modulus E2
    static constexpr char const * defaultE2String() { return "defaultE2"; };

    /// string/key for Young's modulus E3
    static constexpr char const * defaultE3String() { return "defaultE3"; };

    /// string/key for Poisson's ratio Nu12
    static constexpr char const * defaultNu12String() { return "defaultNu12"; };

    /// string/key for Poisson's ratio Nu13
    static constexpr char const * defaultNu13String() { return "defaultNu13"; };

    /// string/key for Poisson's ratio Nu23
    static constexpr char const * defaultNu23String() { return "defaultNu23"; };

    /// string/key for shear modulus G12
    static constexpr char const * defaultG12String() { return "defaultG12"; };

    /// string/key for shear modulus G13
    static constexpr char const * defaultG13String() { return "defaultG13"; };

    /// string/key for shear modulus G23
    static constexpr char const * defaultG23String() { return "defaultG23"; };

    /// string/key for default c11 component of Voigt stiffness tensor
    static constexpr char const * defaultC11String() { return "defaultC11"; };

    /// string/key for default c12 component of Voigt stiffness tensor
    static constexpr char const * defaultC12String() { return "defaultC12"; };

    /// string/key for default c13 component of Voigt stiffness tensor
    static constexpr char const * defaultC13String() { return "defaultC13"; };

    /// string/key for default c22 component of Voigt stiffness tensor
    static constexpr char const * defaultC22String() { return "defaultC22"; };

    /// string/key for default c23 component of Voigt stiffness tensor
    static constexpr char const * defaultC23String() { return "defaultC23"; };

    /// string/key for default c33 component of Voigt stiffness tensor
    static constexpr char const * defaultC33String() { return "defaultC33"; };

    /// string/key for default c44 component of Voigt stiffness tensor
    static constexpr char const * defaultC44String() { return "defaultC44"; };

    /// string/key for default c55 component of Voigt stiffness tensor
    static constexpr char const * defaultC55String() { return "defaultC55"; };

    /// string/key for default c66 component of Voigt stiffness tensor
    static constexpr char const * defaultC66String() { return "defaultC66"; };

    /// string/key for c11 component of Voigt stiffness tensor
    static constexpr char const * c11String() { return "c11"; };

    /// string/key for c12 component of Voigt stiffness tensor
    static constexpr char const * c12String() { return "c12"; };

    /// string/key for c13 component of Voigt stiffness tensor
    static constexpr char const * c13String() { return "c13"; };

    /// string/key for c22 component of Voigt stiffness tensor
    static constexpr char const * c22String() { return "c22"; };

    /// string/key for c23 component of Voigt stiffness tensor
    static constexpr char const * c23String() { return "c23"; };

    /// string/key for c33 component of Voigt stiffness tensor
    static constexpr char const * c33String() { return "c33"; };

    /// string/key for c44 component of Voigt stiffness tensor
    static constexpr char const * c44String() { return "c44"; };

    /// string/key for c55 component of Voigt stiffness tensor
    static constexpr char const * c55String() { return "c55"; };

    /// string/key for c66 component of Voigt stiffness tensor
    static constexpr char const * c66String() { return "c66"; };
  };

  /**
   * @brief Getter for default Young's modulus E1.
   * @return The value of the default Young's modulus E1.
   */

  real64 getDefaultE1() const
  {
    return m_defaultE1;
  }

  /**
   * @brief Setter for the default Young's modulus E1.
   * @param[in] input New value for the default Young's modulus E1.
   */
  void setDefaultE1( real64 const input )
  {
    m_defaultE1 = input;
  }

  /**
   * @brief Getter for default Young's modulus E2.
   * @return The value of the default Young's modulus E2.
   */

  real64 getDefaultE2() const
  {
    return m_defaultE2;
  }

  /**
   * @brief Setter for the default Young's modulus E2.
   * @param[in] input New value for the default Young's modulus E2.
   */
  void setDefaultE2( real64 const input )
  {
    m_defaultE2 = input;
  }

  /**
   * @brief Getter for default Young's modulus E3.
   * @return The value of the default Young's modulus E3.
   */

  real64 getDefaultE3() const
  {
    return m_defaultE3;
  }

  /**
   * @brief Setter for the default Young's modulus E3.
   * @param[in] input New value for the default Young's modulus E3.
   */
  void setDefaultE3( real64 const input )
  {
    m_defaultE3 = input;
  }

  /**
   * @brief Getter for default Poisson's ratio Nu12.
   * @return The value of the default Poisson's ratio Nu12.
   */

  real64 getDefaultNu12() const
  {
    return m_defaultNu12;
  }

  /**
   * @brief Setter for the default Poisson's ratio Nu12.
   * @param[in] input New value for the default Poisson's ratio Nu12.
   */
  void setDefaultNu12( real64 const input )
  {
    m_defaultNu12 = input;
  }

  /**
   * @brief Getter for default Poisson's ratio Nu13.
   * @return The value of the default Poisson's ratio Nu13.
   */

  real64 getDefaultNu13() const
  {
    return m_defaultNu13;
  }

  /**
   * @brief Setter for the default Poisson's ratio Nu13.
   * @param[in] input New value for the default Poisson's ratio Nu13.
   */
  void setDefaultNu13( real64 const input )
  {
    m_defaultNu13 = input;
  }

  /**
   * @brief Getter for default Poisson's ratio Nu23.
   * @return The value of the default Poisson's ratio Nu23.
   */

  real64 getDefaultNu23() const
  {
    return m_defaultNu23;
  }

  /**
   * @brief Setter for the default Poisson's ratio Nu23.
   * @param[in] input New value for the default Poisson's ratio Nu23.
   */
  void setDefaultNu23( real64 const input )
  {
    m_defaultNu23 = input;
  }

  /**
   * @brief Getter for default shear modulus G12.
   * @return The value of the default shear modulus G12.
   */

  real64 getDefaultG12() const
  {
    return m_defaultG12;
  }

  /**
   * @brief Setter for the default shear modulus G12.
   * @param[in] input New value for the default shear modulus G12.
   */
  void setDefaultG12( real64 const input )
  {
    m_defaultG12 = input;
  }

  /**
   * @brief Getter for default shear modulus G13.
   * @return The value of the default shear modulus G13.
   */

  real64 getDefaultG13() const
  {
    return m_defaultG13;
  }

  /**
   * @brief Setter for the default shear modulus G13.
   * @param[in] input New value for the default shear modulus G13.
   */
  void setDefaultG13( real64 const input )
  {
    m_defaultG13 = input;
  }

  /**
   * @brief Getter for default shear modulus G23.
   * @return The value of the default shear modulus G23.
   */

  real64 getDefaultG23() const
  {
    return m_defaultG23;
  }

  /**
   * @brief Setter for the default shear modulus G23.
   * @param[in] input New value for the default shear modulus G23.
   */
  void setDefaultG23( real64 const input )
  {
    m_defaultG23 = input;
  }

  /**
   * @brief Const-Getter for 11 component of Voigt stiffness tensor.
   * @return reference to immutable 11 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > getC11() const { return m_c11; }

  /**
   * @brief Getter for 11 component of Voigt stiffness tensor.
   * @return reference to mutable 11 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > getC11() { return m_c11; }

  /**
   * @brief Const-Getter for 13 component of Voigt stiffness tensor.
   * @return reference to immutable 13 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > getC13() const { return m_c13; }

  /**
   * @brief Getter for 13 component of Voigt stiffness tensor.
   * @return reference to mutable 13 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > getC13() { return m_c13; }

  /**
   * @brief Const-Getter for 33 component of Voigt stiffness tensor.
   * @return reference to immutable 33 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > getC33() const { return m_c33; }

  /**
   * @brief Getter for 33 component of Voigt stiffness tensor.
   * @return reference to mutable 33 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > getC33() { return m_c33; }

  /**
   * @brief Const-Getter for 44 component of Voigt stiffness tensor.
   * @return reference to immutable 44 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > getC44() const { return m_c44; }

  /**
   * @brief Getter for 44 component of Voigt stiffness tensor.
   * @return reference to mutable 44 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > getC44() { return m_c44; }

  /**
   * @brief Const-Getter for 66 component of Voigt stiffness tensor.
   * @return reference to immutable 66 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 const > getC66() const { return m_c66; }

  /**
   * @brief Getter for 66 component of Voigt stiffness tensor.
   * @return reference to mutable 66 component of Voigt stiffness tensor.
   */
  arrayView1d< real64 > getC66() { return m_c66; }

  /**
   * @brief Create a instantiation of the
   *        ElasticOrthotropicUpdates class that refers to the
   *        data in this.
   * @return An instantiation of ElasticOrthotropicUpdates.
   */
  ElasticOrthotropicUpdates createKernelUpdates()
  {
    return ElasticOrthotropicUpdates( m_c11,
                                      m_c12,
                                      m_c13,
                                      m_c22,
                                      m_c23,
                                      m_c33,
                                      m_c44,
                                      m_c55,
                                      m_c66,
                                      m_newStress,
                                      m_oldStress );
  }

  /**
   * @brief Construct an update kernel for a derived type.
   * @tparam UPDATE_KERNEL The type of update kernel from the derived type.
   * @tparam PARAMS The parameter pack to hold the constructor parameters for
   *   the derived update kernel.
   * @param constructorParams The constructor parameter for the derived type.
   * @return An @p UPDATE_KERNEL object.
   */
  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams )
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_c11,
                          m_c12,
                          m_c13,
                          m_c22,
                          m_c23,
                          m_c33,
                          m_c44,
                          m_c55,
                          m_c66,
                          m_newStress,
                          m_oldStress );
  }

protected:
  virtual void postProcessInput() override;

  /// The default value of the Young's modulus E1 for any new
  /// allocations.
  real64 m_defaultE1;

  /// The default value of the Young's modulus E2 for any new
  /// allocations.
  real64 m_defaultE2;

  /// The default value of the Young's modulus E3 for any new
  /// allocations.
  real64 m_defaultE3;

  /// The default value of the Poisson's ratio Nu12 for any new
  /// allocations.
  real64 m_defaultNu12;

  /// The default value of the Poisson's ratio Nu13 for any new
  /// allocations.
  real64 m_defaultNu13;

  /// The default value of the Poisson's ratio Nu23 for any new
  /// allocations.
  real64 m_defaultNu23;

  /// The default value of the Shear modulus G12 for any new
  /// allocations.
  real64 m_defaultG12;

  /// The default value of the Shear modulus G13 for any new
  /// allocations.
  real64 m_defaultG13;

  /// The default value of the Shear modulus G23 for any new
  /// allocations.
  real64 m_defaultG23;

  /// The 11 component of the Voigt stiffness tensor.
  array1d< real64 > m_c11;

  /// The 12 component of the Voigt stiffness tensor.
  array1d< real64 > m_c12;

  /// The 13 component of the Voigt stiffness tensor.
  array1d< real64 > m_c13;

  /// The 22 component of the Voigt stiffness tensor.
  array1d< real64 > m_c22;

  /// The 23 component of the Voigt stiffness tensor.
  array1d< real64 > m_c23;

  /// The 33 component of the Voigt stiffness tensor.
  array1d< real64 > m_c33;

  /// The 44 component of the Voigt stiffness tensor.
  array1d< real64 > m_c44;

  /// The 55 component of the Voigt stiffness tensor.
  array1d< real64 > m_c55;

  /// The 66 component of the Voigt stiffness tensor.
  array1d< real64 > m_c66;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_ELASTICORTHOTROPIC_HPP_ */
