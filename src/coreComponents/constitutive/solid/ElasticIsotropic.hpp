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
 *  @file ElasticIsotropic.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_ELASTICISOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_ELASTICISOTROPIC_HPP_

#include "SolidBase.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsIsotropic.hpp"
#include "constitutive/ExponentialRelation.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class ElasticIsotropicUpdates
 *
 * Class to provide elastic isotropic material updates that may be
 * called from a kernel function.
 */
class ElasticIsotropicUpdates : public SolidBaseUpdates
{
public:
  /**
   * @brief Constructor
   * @param[in] bulkModulus  The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] newStress    The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress    The ArrayView holding the old stress data for each quadrature point.
   */
  ElasticIsotropicUpdates( arrayView1d< real64 const > const & bulkModulus,
                           arrayView1d< real64 const > const & shearModulus,
                           arrayView3d< real64, solid::STRESS_USD > const & newStress,
                           arrayView3d< real64, solid::STRESS_USD > const & oldStress ):
    SolidBaseUpdates( newStress, oldStress ),
    m_bulkModulus( bulkModulus ),
    m_shearModulus( shearModulus )
  {}

  /// Deleted default constructor
  ElasticIsotropicUpdates() = delete;

  /// Default copy constructor
  ElasticIsotropicUpdates( ElasticIsotropicUpdates const & ) = default;

  /// Default move constructor
  ElasticIsotropicUpdates( ElasticIsotropicUpdates && ) = default;

  /// Deleted copy assignment operator
  ElasticIsotropicUpdates & operator=( ElasticIsotropicUpdates const & ) = delete;

  /// Deleted move assignment operator
  ElasticIsotropicUpdates & operator=( ElasticIsotropicUpdates && ) =  delete;

  /// Use the "isotropic" form of inner product compression
  using DiscretizationOps = SolidModelDiscretizationOpsIsotropic;

  /// Use base version of saveConvergedState
  using SolidBaseUpdates::saveConvergedState;

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

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate_StressOnly( localIndex const k,
                                             localIndex const q,
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( &stress )[6] ) const override;

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  real64 ( &stiffness )[6][6] ) const override;

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  DiscretizationOps & stiffness ) const;

  GEOSX_HOST_DEVICE
  virtual void getElasticStiffness( localIndex const k,
                                    real64 ( &stiffness )[6][6] ) const override;

  GEOSX_HOST_DEVICE
  virtual void getElasticStrain( localIndex const k,
                                 localIndex const q,
                                 real64 ( &elasticStrain )[6] ) const override final;

  GEOSX_HOST_DEVICE
  virtual real64 getBulkModulus( localIndex const k ) const override final
  {
    return m_bulkModulus[k];
  }


  // TODO: confirm hyper stress/strain measures before activatiing

  /*
     GEOSX_HOST_DEVICE
     virtual void hyperUpdate( localIndex const k,
                            localIndex const q,
                            real64 const ( & FminusI )[3][3],
                            real64 ( & stress )[6] ) const override final;

     GEOSX_HOST_DEVICE
     virtual void hyperUpdate( localIndex const k,
                            localIndex const q,
                            real64 const ( & FminusI )[3][3],
                            real64 ( & stress )[6],
                            real64 ( & stiffness )[6][6] ) const override final;
   */

protected:

  /// A reference to the ArrayView holding the bulk modulus for each element.
  arrayView1d< real64 const > const m_bulkModulus;

  /// A reference to the ArrayView holding the shear modulus for each element.
  arrayView1d< real64 const > const m_shearModulus;
};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ElasticIsotropicUpdates::getElasticStiffness( localIndex const k,
                                                   real64 ( & stiffness )[6][6] ) const
{
  real64 const G = m_shearModulus[k];
  real64 const lambda = conversions::BulkModAndShearMod::toFirstLame( m_bulkModulus[k], G );

  LvArray::tensorOps::fill< 6, 6 >( stiffness, 0 );

  stiffness[0][0] = lambda + 2*G;
  stiffness[0][1] = lambda;
  stiffness[0][2] = lambda;

  stiffness[1][0] = lambda;
  stiffness[1][1] = lambda + 2*G;
  stiffness[1][2] = lambda;

  stiffness[2][0] = lambda;
  stiffness[2][1] = lambda;
  stiffness[2][2] = lambda + 2*G;

  stiffness[3][3] = G;
  stiffness[4][4] = G;
  stiffness[5][5] = G;
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ElasticIsotropicUpdates::getElasticStrain( localIndex const k,
                                                localIndex const q,
                                                real64 ( & elasticStrain)[6] ) const
{
  real64 const E = conversions::BulkModAndShearMod::toYoungsMod( m_bulkModulus[k], m_shearModulus[k] );
  real64 const nu = conversions::BulkModAndShearMod::toPoissonRatio( m_bulkModulus[k], m_shearModulus[k] );

  elasticStrain[0] = (    m_newStress[k][q][0] - nu*m_newStress[k][q][1] - nu*m_newStress[k][q][2])/E;
  elasticStrain[1] = (-nu*m_newStress[k][q][0] +    m_newStress[k][q][1] - nu*m_newStress[k][q][2])/E;
  elasticStrain[2] = (-nu*m_newStress[k][q][0] - nu*m_newStress[k][q][1] +    m_newStress[k][q][2])/E;

  elasticStrain[3] = m_newStress[k][q][3] / m_shearModulus[k];
  elasticStrain[4] = m_newStress[k][q][4] / m_shearModulus[k];
  elasticStrain[5] = m_newStress[k][q][5] / m_shearModulus[k];
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ElasticIsotropicUpdates::smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                                   localIndex const q,
                                                                   real64 const ( &totalStrain )[6],
                                                                   real64 ( & stress )[6] ) const
{
  GEOSX_UNUSED_VAR( q );

  real64 const twoG   = 2 * m_shearModulus[k];
  real64 const lambda = conversions::BulkModAndShearMod::toFirstLame( m_bulkModulus[k], m_shearModulus[k] );
  real64 const vol    = lambda * ( totalStrain[0] + totalStrain[1] + totalStrain[2] );

  stress[0] = vol + twoG * totalStrain[0];
  stress[1] = vol + twoG * totalStrain[1];
  stress[2] = vol + twoG * totalStrain[2];

  stress[3] = m_shearModulus[k] * totalStrain[3];
  stress[4] = m_shearModulus[k] * totalStrain[4];
  stress[5] = m_shearModulus[k] * totalStrain[5];
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ElasticIsotropicUpdates::smallStrainNoStateUpdate( localIndex const k,
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
void ElasticIsotropicUpdates::smallStrainNoStateUpdate( localIndex const k,
                                                        localIndex const q,
                                                        real64 const ( &totalStrain )[6],
                                                        real64 ( & stress )[6],
                                                        DiscretizationOps & stiffness ) const
{
  smallStrainNoStateUpdate_StressOnly( k, q, totalStrain, stress );
  stiffness.m_bulkModulus = m_bulkModulus[k];
  stiffness.m_shearModulus = m_shearModulus[k];
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ElasticIsotropicUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                            localIndex const q,
                                                            real64 const ( &strainIncrement )[6],
                                                            real64 ( & stress )[6] ) const
{
  smallStrainNoStateUpdate_StressOnly( k, q, strainIncrement, stress ); // stress  = incrementalStress
  LvArray::tensorOps::add< 6 >( stress, m_oldStress[k][q] );            // stress += m_oldStress
  saveStress( k, q, stress );                                           // m_newStress = stress
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ElasticIsotropicUpdates::smallStrainUpdate( localIndex const k,
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
void ElasticIsotropicUpdates::smallStrainUpdate( localIndex const k,
                                                 localIndex const q,
                                                 real64 const ( &strainIncrement )[6],
                                                 real64 ( & stress )[6],
                                                 DiscretizationOps & stiffness ) const
{
  smallStrainUpdate_StressOnly( k, q, strainIncrement, stress );
  stiffness.m_bulkModulus = m_bulkModulus[k];
  stiffness.m_shearModulus = m_shearModulus[k];
}


// TODO: need to confirm stress / strain measures before activating hyper inferface
/*
   GEOSX_HOST_DEVICE
   GEOSX_FORCE_INLINE
   void ElasticIsotropicUpdates::hyperUpdate( localIndex const k,
                                           localIndex const q,
                                           real64 const (&FmI)[3][3],
                                           real64 ( & stress )[ 6 ] ) const
   {
   GEOSX_UNUSED_VAR(q);

   real64 const C1 = 0.5 * m_shearModulus[k];
   real64 const D1 = 0.5 * m_bulkModulus[k];
   real64 const detFm1 = FmI[0][0] + FmI[1][1] + FmI[2][2]
                        - FmI[1][2]*FmI[2][1] + FmI[1][1]*FmI[2][2]
 + FmI[0][2]*(-FmI[2][0] - FmI[1][1]*FmI[2][0] + FmI[1][0]*FmI[2][1])
 + FmI[0][1]*(-FmI[1][0] + FmI[1][2]*FmI[2][0] - FmI[1][0]*FmI[2][2])
 + FmI[0][0]*( FmI[1][1] - FmI[1][2]*FmI[2][1] + FmI[2][2] + FmI[1][1]*FmI[2][2]);


   real64 const p = -2 * D1 * ( detFm1 + 1.0 ) * detFm1;
   real64 devB[6] = { 1/3 * (2 * FmI[0][0] * (2 + FmI[0][0]) - FmI[1][1] * (2 + FmI[1][1]) - FmI[2][2] * (2 + FmI[2][2])
 ++
                            2 * FmI[0][1]*FmI[0][1] + 2 * FmI[0][2] * FmI[0][2] - FmI[1][0] * FmI[1][0] - FmI[1][2] *
 +FmI[1][2] -
                            FmI[2][0] * FmI[2][0] - FmI[2][1] * FmI[2][1]),
                     1/3 * (-FmI[0][0] * (2 + FmI[0][0]) + 2 * FmI[1][1] * ( 2 + FmI[1][1]) - FmI[2][2] * (2 +
 +FmI[2][2]) -
                            FmI[0][1]*FmI[0][1] - FmI[0][2]*FmI[0][2] + 2 * FmI[1][0]*FmI[1][0] + 2 *
 +FmI[1][2]*FmI[1][2] - FmI[2][0]*FmI[2][0] - FmI[2][1]*
                            FmI[2][1]),
                     1/3 *(-FmI[0][0] * (2 + FmI[0][0]) - FmI[1][1] * (2 + FmI[1][1]) + 2 * FmI[2][2] * (2 + FmI[2][2])
 +-
                           FmI[0][1]*FmI[0][1] - FmI[0][2]*FmI[0][2] - FmI[1][0]*FmI[1][0] - FmI[1][2]*FmI[1][2] + 2 *
 +FmI[2][0]*FmI[2][0] + 2 * FmI[2][1]*
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
 */


/**
 * @class ElasticIsotropic
 *
 * Class to provide an elastic isotropic material response.
 */
class ElasticIsotropic : public SolidBase
{
public:

  /// Alias for ElasticIsotropicUpdates
  using KernelWrapper = ElasticIsotropicUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  ElasticIsotropic( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~ElasticIsotropic() override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "ElasticIsotropic";

  /**
   * @brief Static catalog string
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string catalogName() { return m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  /// Keys for data specified in this class.
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for default bulk modulus
    static constexpr char const * defaultBulkModulusString() { return "defaultBulkModulus"; }

    /// string/key for default poisson ratio
    static constexpr char const * defaultPoissonRatioString() { return "defaultPoissonRatio"; }

    /// string/key for default shear modulus
    static constexpr char const * defaultShearModulusString() { return "defaultShearModulus"; }

    /// string/key for default Young's modulus
    static constexpr char const * defaultYoungsModulusString() { return "defaultYoungsModulus"; }

    /// string/key for bulk modulus
    static constexpr char const * bulkModulusString() { return "bulkModulus"; }

    /// string/key for shear modulus
    static constexpr char const * shearModulusString() { return "shearModulus"; }
  };

  /**
   * @brief Accessor for bulk modulus
   * @return A const reference to arrayView1d<real64> containing the bulk
   *         modulus (at every element).
   */
  arrayView1d< real64 > const bulkModulus() { return m_bulkModulus; }

  /**
   * @brief Const accessor for bulk modulus
   * @return A const reference to arrayView1d<real64 const> containing the bulk
   *         modulus (at every element).
   */
  arrayView1d< real64 const > const bulkModulus() const { return m_bulkModulus; }

  /**
   * @brief Accessor for shear modulus
   * @return A const reference to arrayView1d<real64> containing the shear
   *         modulus (at every element).
   */
  arrayView1d< real64 > const shearModulus() { return m_shearModulus; }

  /**
   * @brief Const accessor for shear modulus
   * @return A const reference to arrayView1d<real64 const> containing the
   *         shear modulus (at every element).
   */
  arrayView1d< real64 const > const shearModulus() const { return m_shearModulus; }

  /**
   * @brief Create a instantiation of the ElasticIsotropicUpdate class
   *        that refers to the data in this.
   * @param includeState Flag whether to pass state arrays that may not be needed for "no-state" updates
   * @return An instantiation of ElasticIsotropicUpdate.
   */
  ElasticIsotropicUpdates createKernelUpdates( bool const includeState = true ) const
  {
    if( includeState )
    {
      return ElasticIsotropicUpdates( m_bulkModulus,
                                      m_shearModulus,
                                      m_newStress,
                                      m_oldStress );
    }
    else // for "no state" updates, pass empty views to avoid transfer of stress data to device
    {
      return ElasticIsotropicUpdates( m_bulkModulus,
                                      m_shearModulus,
                                      arrayView3d< real64, solid::STRESS_USD >(),
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
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams ) const
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_bulkModulus,
                          m_shearModulus,
                          m_newStress,
                          m_oldStress );
  }


protected:

  /// Post-process XML data
  virtual void postProcessInput() override;

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

#endif /* GEOSX_CONSTITUTIVE_SOLID_ELASTICISOTROPIC_HPP_ */
