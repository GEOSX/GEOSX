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
 *  @file LinearElasticIsotropic.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_LINEARELASTICISOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_LINEARELASTICISOTROPIC_HPP_
#include "SolidBase.hpp"
#include "constitutive/ExponentialRelation.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "SolidModelDiscretizationOpsIsotropic.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class LinearElasticIsotropicUpdates
 *
 * Class to provide linear elastic isotropic material updates that may be
 * called from a kernel function.
 */
class LinearElasticIsotropicUpdates : public SolidBaseUpdates
{
public:
  using DiscretizationOps = SolidModelDiscretizationOpsIsotropic;

  /**
   * @brief Constructor
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each
   *                        element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each
   *                         element.
   * @param[in] stress The ArrayView holding the stress data for each quadrature
   *                   point.
   */
  LinearElasticIsotropicUpdates( arrayView1d< real64 const > const & bulkModulus,
                                 arrayView1d< real64 const > const & shearModulus,
                                 arrayView3d< real64, solid::STRESS_USD > const & stress ):
    SolidBaseUpdates( stress ),
    m_bulkModulus( bulkModulus ),
    m_shearModulus( shearModulus )
  {}

  /// Deleted default constructor
  LinearElasticIsotropicUpdates() = delete;

  /// Default copy constructor
  LinearElasticIsotropicUpdates( LinearElasticIsotropicUpdates const & ) = default;

  /// Default move constructor
  LinearElasticIsotropicUpdates( LinearElasticIsotropicUpdates && ) = default;

  /// Deleted copy assignment operator
  LinearElasticIsotropicUpdates & operator=( LinearElasticIsotropicUpdates const & ) = delete;

  /// Deleted move assignment operator
  LinearElasticIsotropicUpdates & operator=( LinearElasticIsotropicUpdates && ) =  delete;


  /**
   * @copydoc SolidBase::GetStiffness
   */
  GEOSX_HOST_DEVICE inline
  virtual void GetStiffness( localIndex const k,
                             localIndex const q,
                             real64 (& c)[6][6] ) const override
  {
    GEOSX_UNUSED_VAR( q );
    real64 const G = m_shearModulus[k];
    real64 const Lame = m_bulkModulus[k] - 2.0/3.0 * G;

    LvArray::tensorOps::fill< 6, 6 >( c, 0 );

    c[0][0] = Lame + 2 * G;
    c[0][1] = Lame;
    c[0][2] = Lame;

    c[1][0] = Lame;
    c[1][1] = Lame + 2 * G;
    c[1][2] = Lame;

    c[2][0] = Lame;
    c[2][1] = Lame;
    c[2][2] = Lame + 2 * G;

    c[3][3] = G;

    c[4][4] = G;

    c[5][5] = G;
  }

  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  void setDiscretizationOps( localIndex const k,
                             localIndex const q,
                             DiscretizationOps & discOps ) const
  {
    GEOSX_UNUSED_VAR( q )
    real64 const G = m_shearModulus[k];
    real64 const Lame = m_bulkModulus[k] - 2.0/3.0 * G;
    discOps.m_lambda = Lame;
    discOps.m_shearModulus = G;
  }

  GEOSX_HOST_DEVICE
  virtual void SmallStrainNoState( localIndex const k,
                                   real64 const ( &voigtStrain )[ 6 ],
                                   real64 ( &stress )[ 6 ] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void SmallStrain( localIndex const k,
                            localIndex const q,
                            real64 const ( &voigtStrainInc )[ 6 ] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void HypoElastic( localIndex const k,
                            localIndex const q,
                            real64 const ( &Ddt )[ 6 ],
                            real64 const ( &Rot )[ 3 ][ 3 ] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void HyperElastic( localIndex const k,
                             real64 const (&FmI)[3][3],
                             real64 ( &stress )[ 6 ] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void HyperElastic( localIndex const k,
                             localIndex const q,
                             real64 const (&FmI)[3][3] ) const override final;

  GEOSX_HOST_DEVICE
  virtual real64 calculateStrainEnergyDensity( localIndex const k,
                                               localIndex const q ) const override;

  GEOSX_HOST_DEVICE
  virtual real64 getBulkModulus( localIndex const k ) const
  {
    return m_bulkModulus[k];
  }

  GEOSX_HOST_DEVICE
  virtual real64 getShearModulus( localIndex const k ) const
  {
    return m_shearModulus[k];
  }

private:
  /// A reference to the ArrayView holding the bulk modulus for each element.
  arrayView1d< real64 const > const m_bulkModulus;

  /// A reference to the ArrayView holding the shear modulus for each element.
  arrayView1d< real64 const > const m_shearModulus;

};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void LinearElasticIsotropicUpdates::SmallStrainNoState( localIndex const k,
                                                        real64 const ( &voigtStrain )[ 6 ],
                                                        real64 ( & stress )[ 6 ] ) const
{
  real64 const lambda = m_bulkModulus[k] - 2.0/3.0 * m_shearModulus[k];
  real64 const diag = lambda * ( voigtStrain[0] + voigtStrain[1] + voigtStrain[2] );
  real64 const TwoG = 2.0 * m_shearModulus[k];

  stress[0] = diag + TwoG * voigtStrain[0];
  stress[1] = diag + TwoG * voigtStrain[1];
  stress[2] = diag + TwoG * voigtStrain[2];
  stress[3] = m_shearModulus[k] * voigtStrain[3];
  stress[4] = m_shearModulus[k] * voigtStrain[4];
  stress[5] = m_shearModulus[k] * voigtStrain[5];

}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void LinearElasticIsotropicUpdates::SmallStrain( localIndex const k,
                                                 localIndex const q,
                                                 real64 const ( &voigtStrainInc )[ 6 ] ) const
{
  real64 const lambda = m_bulkModulus[k] - 2.0/3.0 * m_shearModulus[k];
  real64 const volStrain = ( voigtStrainInc[0] + voigtStrainInc[1] + voigtStrainInc[2] );
  real64 const TwoG = 2.0 * m_shearModulus[k];

  m_stress( k, q, 0 ) =  m_stress( k, q, 0 ) + TwoG * voigtStrainInc[0] + lambda * volStrain;
  m_stress( k, q, 1 ) =  m_stress( k, q, 1 ) + TwoG * voigtStrainInc[1] + lambda * volStrain;
  m_stress( k, q, 2 ) =  m_stress( k, q, 2 ) + TwoG * voigtStrainInc[2] + lambda * volStrain;
  m_stress( k, q, 3 ) =  m_stress( k, q, 3 ) + m_shearModulus[k] * voigtStrainInc[3];
  m_stress( k, q, 4 ) =  m_stress( k, q, 4 ) + m_shearModulus[k] * voigtStrainInc[4];
  m_stress( k, q, 5 ) =  m_stress( k, q, 5 ) + m_shearModulus[k] * voigtStrainInc[5];

}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void LinearElasticIsotropicUpdates::HypoElastic( localIndex const k,
                                                 localIndex const q,
                                                 real64 const ( &Ddt )[ 6 ],
                                                 real64 const ( &Rot )[ 3 ][ 3 ] ) const
{
  real64 const lambda = m_bulkModulus[ k ] - 2.0 / 3.0 * m_shearModulus[ k ];
  real64 const volStrain = ( Ddt[ 0 ] + Ddt[ 1 ] + Ddt[ 2 ] );
  real64 const TwoG = 2.0 * m_shearModulus[ k ];

  m_stress( k, q, 0 ) =  m_stress( k, q, 0 ) + TwoG * Ddt[ 0 ] + lambda * volStrain;
  m_stress( k, q, 1 ) =  m_stress( k, q, 1 ) + TwoG * Ddt[ 1 ] + lambda * volStrain;
  m_stress( k, q, 2 ) =  m_stress( k, q, 2 ) + TwoG * Ddt[ 2 ] + lambda * volStrain;
  m_stress( k, q, 3 ) =  m_stress( k, q, 3 ) + TwoG * Ddt[ 3 ];
  m_stress( k, q, 4 ) =  m_stress( k, q, 4 ) + TwoG * Ddt[ 4 ];
  m_stress( k, q, 5 ) =  m_stress( k, q, 5 ) + TwoG * Ddt[ 5 ];

  real64 temp[ 6 ] = { 0 };
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( temp, Rot, m_stress[ k ][ q ] );
  LvArray::tensorOps::copy< 6 >( m_stress[ k ][ q ], temp );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void LinearElasticIsotropicUpdates::HyperElastic( localIndex const k,
                                                  real64 const (&FmI)[3][3],
                                                  real64 ( & stress )[ 6 ] ) const
{
  real64 const C1 = 0.5 * m_shearModulus[k];
  real64 const D1 = 0.5 * m_bulkModulus[k];
  real64 const detFm1 = FmI[0][0] + FmI[1][1] + FmI[2][2]
                        - FmI[1][2]*FmI[2][1] + FmI[1][1]*FmI[2][2]
                        + FmI[0][2]*(-FmI[2][0] - FmI[1][1]*FmI[2][0] + FmI[1][0]*FmI[2][1])
                        + FmI[0][1]*(-FmI[1][0] + FmI[1][2]*FmI[2][0] - FmI[1][0]*FmI[2][2])
                        + FmI[0][0]*( FmI[1][1] - FmI[1][2]*FmI[2][1] + FmI[2][2] + FmI[1][1]*FmI[2][2]);


  real64 const p = -2 * D1 * ( detFm1 + 1.0 ) * detFm1;
  real64 devB[6] = { 1/3 * (2 * FmI[0][0] * (2 + FmI[0][0]) - FmI[1][1] * (2 + FmI[1][1]) - FmI[2][2] * (2 + FmI[2][2]) +
                            2 * FmI[0][1]*FmI[0][1] + 2 * FmI[0][2] * FmI[0][2] - FmI[1][0] * FmI[1][0] - FmI[1][2] * FmI[1][2] -
                            FmI[2][0] * FmI[2][0] - FmI[2][1] * FmI[2][1]),
                     1/3 * (-FmI[0][0] * (2 + FmI[0][0]) + 2 * FmI[1][1] * ( 2 + FmI[1][1]) - FmI[2][2] * (2 + FmI[2][2]) -
                            FmI[0][1]*FmI[0][1] - FmI[0][2]*FmI[0][2] + 2 * FmI[1][0]*FmI[1][0] + 2 * FmI[1][2]*FmI[1][2] - FmI[2][0]*FmI[2][0] - FmI[2][1]*
                            FmI[2][1]),
                     1/3 *(-FmI[0][0] * (2 + FmI[0][0]) - FmI[1][1] * (2 + FmI[1][1]) + 2 * FmI[2][2] * (2 + FmI[2][2]) -
                           FmI[0][1]*FmI[0][1] - FmI[0][2]*FmI[0][2] - FmI[1][0]*FmI[1][0] - FmI[1][2]*FmI[1][2] + 2 * FmI[2][0]*FmI[2][0] + 2 * FmI[2][1]*
                           FmI[2][1]),
                     FmI[1][2] + FmI[1][0] * FmI[2][0] + FmI[2][1] + FmI[1][1]*FmI[2][1] + FmI[1][2]*FmI[2][2],
                     FmI[0][2] + FmI[2][0] + FmI[0][0] * FmI[2][0] + FmI[0][1]*FmI[2][1] + FmI[0][2]*FmI[2][2],
                     FmI[0][1] + FmI[1][0] + FmI[0][0] * FmI[1][0] + FmI[0][1]*FmI[1][1] + FmI[0][2]*FmI[1][2]
  };

  real64 const C = 2 * C1 / pow( detFm1 + 1, 2.0/3.0 );
  stress[0] = -p + C * devB[0];
  stress[1] = -p + C * devB[1];
  stress[2] = -p + C * devB[2];
  stress[3] = C * devB[3];
  stress[4] = C * devB[4];
  stress[5] = C * devB[5];
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void LinearElasticIsotropicUpdates::HyperElastic( localIndex const k,
                                                  localIndex const q,
                                                  real64 const (&FmI)[3][3] ) const
{
  real64 stress[ 6 ];
  HyperElastic( k, FmI, stress );
  LvArray::tensorOps::copy< 6 >( m_stress[ k ][ q ], stress );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 LinearElasticIsotropicUpdates::calculateStrainEnergyDensity( localIndex const k,
                                                                    localIndex const q ) const
{
  real64 const invE = ( 3.0 * m_bulkModulus[k] + m_shearModulus[k] ) / ( 9.0 * m_bulkModulus[k] * m_shearModulus[k] );
  real64 const nu = ( 1.5 * m_bulkModulus[k] - m_shearModulus[k] ) / ( 3.0 * m_bulkModulus[k] + m_shearModulus[k] );

  auto const & stress = m_stress[k][q];
  real64 const newStrainEnergyDensity = ( stress[0]*stress[0] + stress[1]*stress[1] + stress[2]*stress[2] -
                                          2 * ( nu       * ( stress[1]*stress[2] + stress[0]*stress[1] + stress[0]*stress[2] ) -
                                                (1 + nu) * ( stress[3]*stress[3] + stress[4]*stress[4] + stress[5]*stress[5] )
                                                )
                                          ) * invE * 0.5;
  // Make sure strain energy is always non-negative
  GEOSX_ASSERT_MSG( newStrainEnergyDensity >= 0.0,
                    "negative strain energy density" );

  //make adjustments for Volumetric Split
  /////////////////////////////////
  //real64 traceOfStress = this->m_stress(k,q,0) + this->m_stress(k,q,1) + this->m_stress(k,q,2);
  //real64 compressionIndicator = 0;
  // if (traceOfStress < 0.0)
  // {
  //   compressionIndicator = 1;
  // }
  // real64 const activeStrainEnergyDensity = newStrainEnergyDensity -
  // compressionIndicator*(traceOfStress/3.0)*(traceOfStress/3.0)/(2*m_bulkModulus[k]);
  // return activeStrainEnergyDensity;
  /////////////////////////////////
  //end of adjustments

  //regular return
  return newStrainEnergyDensity;
}


/**
 * @class LinearElasticIsotropic
 *
 * Class to provide a linear elastic isotropic material response.
 */
class LinearElasticIsotropic : public SolidBase
{
public:

  /// @typedef Alias for LinearElasticIsotropicUpdates
  using KernelWrapper = LinearElasticIsotropicUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  LinearElasticIsotropic( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~LinearElasticIsotropic() override;
  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "LinearElasticIsotropic";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string CatalogName() { return m_catalogNameString; }

  virtual string getCatalogName() const override { return CatalogName(); }

  ///@}

  /**
   * @struct Set of "char const *" and keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for default bulk modulus
    static constexpr auto defaultBulkModulusString  = "defaultBulkModulus";

    /// string/key for default poisson ratio
    static constexpr auto defaultPoissonRatioString =  "defaultPoissonRatio";

    /// string/key for default shear modulus
    static constexpr auto defaultShearModulusString = "defaultShearModulus";

    /// string/key for default youngs modulus
    static constexpr auto defaultYoungsModulusString =  "defaultYoungsModulus";

    /// string/key for bulk modulus
    static constexpr auto bulkModulusString  = "BulkModulus";
    /// string/key for shear modulus
    static constexpr auto shearModulusString = "ShearModulus";
  };

  /**
   * @brief Setter for the default bulk modulus.
   * @param[in] bulkModulus The value that m_defaultBulkModulus will be set to.
   */
  void setDefaultBulkModulus( real64 const bulkModulus );

  /**
   * @brief Setter for the default shear modulus.
   * @param[in] bulkModulus The value that m_defaultShearModulus will be set to.
   */
  void setDefaultShearModulus( real64 const shearModulus );

  /**
   * @brief Accessor for bulk modulus
   * @return A const reference to arrayView1d<real64> containing the bulk
   *         modulus (at every element).
   */
  arrayView1d< real64 > bulkModulus() { return m_bulkModulus; }

  /**
   * @brief Const accessor for bulk modulus
   * @return A const reference to arrayView1d<real64 const> containing the bulk
   *         modulus (at every element).
   */
  arrayView1d< real64 const > bulkModulus() const { return m_bulkModulus; }

  /**
   * @brief Accessor for shear modulus
   * @return A const reference to arrayView1d<real64> containing the shear
   *         modulus (at every element).
   */
  arrayView1d< real64 > shearModulus() { return m_shearModulus; }

  /**
   * @brief Const accessor for shear modulus
   * @return A const reference to arrayView1d<real64 const> containing the
   *         shear modulus (at every element).
   */
  arrayView1d< real64 const > shearModulus() const { return m_shearModulus; }

  /**
   * @brief Create a instantiation of the LinearElasticIsotropicUpdate class
   *        that refers to the data in this.
   * @return An instantiation of LinearElasticIsotropicUpdate.
   */
  LinearElasticIsotropicUpdates createKernelUpdates( bool const includeState = true ) const
  {
    if( includeState )
    {
      return LinearElasticIsotropicUpdates( m_bulkModulus, m_shearModulus, m_stress );
    }
    else
    {
      return LinearElasticIsotropicUpdates( m_bulkModulus,
                                            m_shearModulus,
                                            arrayView3d< real64, solid::STRESS_USD >() );
    }
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
                          m_bulkModulus,
                          m_shearModulus,
                          m_stress );
  }


protected:
  virtual void PostProcessInput() override;

private:

  /// The default value of the bulk modulus for any new allocations.
  real64 m_defaultBulkModulus;

  /// The default value of the shear modulus for any new allocations.
  real64 m_defaultShearModulus;

  /// The bulk modulus for each upper level dimension (i.e. cell) of *this
  array1d< real64 > m_bulkModulus;

  /// The shear modulus for each upper level dimension (i.e. cell) of *this
  array1d< real64 > m_shearModulus;

};

}

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_LINEARELASTICISOTROPIC_HPP_ */
