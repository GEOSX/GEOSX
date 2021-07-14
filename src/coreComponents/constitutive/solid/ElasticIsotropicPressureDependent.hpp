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
 *  @file ElasticIsotropicPressureDependent.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_ELASTICISOTROPICPRESSUREDEPENDENT_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_ELASTICISOTROPICPRESSUREDEPENDENT_HPP_

#include "SolidBase.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsIsotropic.hpp"
#include "constitutive/ExponentialRelation.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class ElasticIsotropicPressureDependentUpdates
 *
 * Class to provide elastic isotropic material updates that may be
 * called from a kernel function.
 */
class ElasticIsotropicPressureDependentUpdates : public SolidBaseUpdates
{
public:
  /**
   * @brief Constructor
   * @param[in] bulkModulus  The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] newStress    The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress    The ArrayView holding the old stress data for each quadrature point.
   */
  ElasticIsotropicPressureDependentUpdates( real64 const & refPressure,
                                            real64 const & refStrainVol,
                                            arrayView1d< real64 const > const & recompressionIndex,
                                            arrayView1d< real64 const > const & shearModulus,
                                            arrayView3d< real64, solid::STRESS_USD > const & newStress,
                                            arrayView3d< real64, solid::STRESS_USD > const & oldStress ):
    SolidBaseUpdates( newStress, oldStress ),
    m_refPressure( refPressure ),
    m_refStrainVol( refStrainVol ),
    m_recompressionIndex( recompressionIndex ),
    m_shearModulus( shearModulus )
  {}

  /// Deleted default constructor
  ElasticIsotropicPressureDependentUpdates() = delete;

  /// Default copy constructor
  ElasticIsotropicPressureDependentUpdates( ElasticIsotropicPressureDependentUpdates const & ) = default;

  /// Default move constructor
  ElasticIsotropicPressureDependentUpdates( ElasticIsotropicPressureDependentUpdates && ) = default;

  /// Deleted copy assignment operator
  ElasticIsotropicPressureDependentUpdates & operator=( ElasticIsotropicPressureDependentUpdates const & ) = delete;

  /// Deleted move assignment operator
  ElasticIsotropicPressureDependentUpdates & operator=( ElasticIsotropicPressureDependentUpdates && ) =  delete;

  /// Use the "isotropic" form of inner product compression
  using DiscretizationOps = SolidModelDiscretizationOpsIsotropic;

  /// Use base version of saveConvergedState
  using SolidBaseUpdates::saveConvergedState;

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
                                    localIndex const q,
                                    real64 ( &stiffness )[6][6] ) const override;

  GEOSX_HOST_DEVICE
  virtual void getElasticStrain( localIndex const k,
                                 localIndex const q,
                                 real64 ( &elasticStrain )[6] ) const override final;


protected:

  /// A reference to the ArrayView holding the reference pressure for each element.
  real64 const m_refPressure;

  /// A reference to the ArrayView holding the reference volumetric strain for each element.
  real64 const m_refStrainVol;

  /// A reference to the ArrayView holding the recompression index for each element.
  arrayView1d< real64 const > const m_recompressionIndex;

  /// A reference to the ArrayView holding the shear modulus for each element.
  arrayView1d< real64 const > const m_shearModulus;
};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ElasticIsotropicPressureDependentUpdates::getElasticStiffness( localIndex const k,
                                                                    localIndex const q,
                                                                    real64 ( & stiffness )[6][6] ) const
{
  real64 const mu     = m_shearModulus[k];
  real64 const Cr     = m_recompressionIndex[k];

  real64 bulkModulus;

  real64 deviator[6];
  real64 stress[6];
  real64 P;
  real64 Q;

  for( localIndex i=0; i<6; ++i )
  {
    stress[i] = m_newStress[k][q][i];
  }

  twoInvariant::stressDecomposition( stress,
                                     P,
                                     Q,
                                     deviator );

  bulkModulus = -P/Cr;


  LvArray::tensorOps::fill< 6, 6 >( stiffness, 0 );



  real64 const lambda = bulkModulus - 2./3. * mu;

  stiffness[0][0] = lambda + 2*mu;
  stiffness[0][1] = lambda;
  stiffness[0][2] = lambda;

  stiffness[1][0] = lambda;
  stiffness[1][1] = lambda + 2*mu;
  stiffness[1][2] = lambda;

  stiffness[2][0] = lambda;
  stiffness[2][1] = lambda;
  stiffness[2][2] = lambda + 2*mu;

  stiffness[3][3] = mu;
  stiffness[4][4] = mu;
  stiffness[5][5] = mu;
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ElasticIsotropicPressureDependentUpdates::getElasticStrain( localIndex const k,
                                                                 localIndex const q,
                                                                 real64 ( & elasticStrain)[6] ) const
{
  real64 const mu     = m_shearModulus[k];
  real64 const p0     = m_refPressure;
  real64 const eps_v0 = m_refStrainVol;
  real64 const Cr     = m_recompressionIndex[k];
  real64 deviator[6];
  real64 stress[6];
  real64 P;
  real64 Q;
  real64 elasticStrainVol;
  real64 elasticStrainDev;

  for( localIndex i=0; i<6; ++i )
  {
    stress[i] = m_newStress[k][q][i];
  }

  twoInvariant::stressDecomposition( stress,
                                     P,
                                     Q,
                                     deviator );

  elasticStrainVol = std::log( P/p0 ) * Cr * (-1.0) + eps_v0;


  elasticStrainDev = Q/3./mu;

  twoInvariant::strainRecomposition( elasticStrainVol,
                                     elasticStrainDev,
                                     deviator,
                                     elasticStrain );

}



GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ElasticIsotropicPressureDependentUpdates::smallStrainUpdate( localIndex const k,
                                                                  localIndex const q,
                                                                  real64 const ( &strainIncrement )[6],
                                                                  real64 ( & stress )[6],
                                                                  real64 ( & stiffness )[6][6] ) const
{
  //smallStrainUpdate_StressOnly( k, q, strainIncrement, stress );
  // Rename variables for easier implementation

  real64 const mu     = m_shearModulus[k];
  real64 const p0     = m_refPressure;
  real64 const eps_v0 = m_refStrainVol;
  real64 const Cr     = m_recompressionIndex[k];

  // two-invariant decomposition of old stress in P-Q space (mean & deviatoric stress)

  real64 oldP;
  real64 oldQ;
  real64 P;
  real64 Q;
  real64 oldDeviator[6];
  real64 deviator[6];
  real64 oldStrainElastic[6];
  real64 strainElasticTotal[6];
  real64 eps_s_elastic;
  real64 eps_v_elastic;
  real64 oldElasticStrainVol;
  real64 oldElasticStrainDev;


  for( localIndex i=0; i<6; ++i )
  {
    stress[i] = m_oldStress[k][q][i];
  }

  twoInvariant::stressDecomposition( stress,
                                     oldP,
                                     oldQ,
                                     oldDeviator );

  // Recover elastic strains from the previous step, based on stress from the previous step
  // [Note: in order to minimize data transfer, we are not storing and passing elastic strains]

  oldElasticStrainVol = std::log( oldP/p0 ) * Cr * (-1.0) + eps_v0;



  oldElasticStrainDev = oldQ/3./mu;

  // Now recover the old strain tensor from the strain invariants.
  // Note that we need the deviatoric direction (n-hat) from the previous step.

  twoInvariant::strainRecomposition( oldElasticStrainVol,
                                     oldElasticStrainDev,
                                     oldDeviator,
                                     oldStrainElastic );

  // Total elastic strain

  for( localIndex i=0; i<6; ++i )
  {
    strainElasticTotal[i] = oldStrainElastic[i] + strainIncrement[i];
  }
  // two-invariant decomposition of trial elastic strain

  twoInvariant::strainDecomposition( strainElasticTotal,
                                     eps_v_elastic,
                                     eps_s_elastic,
                                     deviator );

  // Calculate trial mean and deviatoric stress

  P = p0 * std::exp( -1./Cr* (eps_v_elastic-eps_v0));


  Q = 3. * mu * eps_s_elastic;

  twoInvariant::stressRecomposition( P,
                                     Q,
                                     deviator,
                                     stress );
  saveStress( k, q, stress );
  getElasticStiffness( k, q, stiffness );
}

//TODO: implement the discretizationOps version of smallStrainUpdate
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ElasticIsotropicPressureDependentUpdates::smallStrainUpdate( localIndex const k,
                                                                  localIndex const q,
                                                                  real64 const ( &strainIncrement )[6],
                                                                  real64 ( & stress )[6],
                                                                  DiscretizationOps & stiffness ) const
{
  //smallStrainUpdate_StressOnly( k, q, strainIncrement, stress );
  // Rename variables for easier implementation

  real64 const mu     = m_shearModulus[k];
  real64 const p0     = m_refPressure;
  real64 const eps_v0 = m_refStrainVol;
  real64 const Cr     = m_recompressionIndex[k];

  // two-invariant decomposition of old stress in P-Q space (mean & deviatoric stress)

  real64 oldP;
  real64 oldQ;
  real64 P;
  real64 Q;
  real64 oldDeviator[6];
  real64 deviator[6];
  real64 oldStrainElastic[6];
  real64 strainElasticTotal[6];
  real64 eps_s_elastic;
  real64 eps_v_elastic;
  real64 oldElasticStrainVol;
  real64 oldElasticStrainDev;

  real64 bulkModulus = -p0/Cr;

  for( localIndex i=0; i<6; ++i )
  {
    stress[i] = m_oldStress[k][q][i];
  }

  twoInvariant::stressDecomposition( stress,
                                     oldP,
                                     oldQ,
                                     oldDeviator );

  // Recover elastic strains from the previous step, based on stress from the previous step
  // [Note: in order to minimize data transfer, we are not storing and passing elastic strains]

  oldElasticStrainVol = std::log( oldP/p0 ) * Cr * (-1.0) + eps_v0;


  oldElasticStrainDev = oldQ/3./mu;

  // Now recover the old strain tensor from the strain invariants.
  // Note that we need the deviatoric direction (n-hat) from the previous step.

  twoInvariant::strainRecomposition( oldElasticStrainVol,
                                     oldElasticStrainDev,
                                     oldDeviator,
                                     oldStrainElastic );

  // Total elastic strain

  for( localIndex i=0; i<6; ++i )
  {
    strainElasticTotal[i] = oldStrainElastic[i] + strainIncrement[i];
  }
  // two-invariant decomposition of elastic strain

  twoInvariant::strainDecomposition( strainElasticTotal,
                                     eps_v_elastic,
                                     eps_s_elastic,
                                     deviator );

  // Calculate mean and deviatoric stress

  P = p0 * std::exp( -1./Cr* (eps_v_elastic-eps_v0));


  Q = 3. * mu * eps_s_elastic;

  twoInvariant::stressRecomposition( P,
                                     Q,
                                     deviator,
                                     stress );

  bulkModulus = -P/Cr;


  saveStress( k, q, stress );
  stiffness.m_bulkModulus = bulkModulus;
  stiffness.m_shearModulus = m_shearModulus[k];
}



/**
 * @class ElasticIsotropicPressureDependent
 *
 * Class to provide an elastic isotropic material response.
 */
class ElasticIsotropicPressureDependent : public SolidBase
{
public:

  /// Alias for ElasticIsotropicPressureDependentUpdates
  using KernelWrapper = ElasticIsotropicPressureDependentUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  ElasticIsotropicPressureDependent( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~ElasticIsotropicPressureDependent() override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "ElasticIsotropicPressureDependent";

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
//    /// string/key for default poisson ratio
//    static constexpr char const * defaultPoissonRatioString() { return "defaultPoissonRatio"; }

    /// string/key for default shear modulus
    static constexpr char const * defaultShearModulusString() { return "defaultShearModulus"; }

//    /// string/key for default Young's modulus
//    static constexpr char const * defaultYoungsModulusString() { return "defaultYoungsModulus"; }

    /// string/key for default reference pressure
    static constexpr char const * defaultRefPressureString() { return "defaultRefPressure"; }

    /// string/key for default reference volumetric strain
    static constexpr char const * defaultRefStrainVolString() { return "defaultRefStrainVol"; }

    /// string/key for default recompression index
    static constexpr char const * defaultRecompressionIndexString() { return "defaultRecompressionIndex"; }

    /// string/key for reference pressure
    static constexpr char const * refPressureString() { return "refPressure"; }

    /// string/key for reference volumetric strain
    static constexpr char const * refStrainVolString() { return "refStrainVol"; }

    /// string/key for recompression index
    static constexpr char const * recompressionIndexString() { return "recompressionIndex"; }

    /// string/key for shear modulus
    static constexpr char const * shearModulusString() { return "shearModulus"; }
  };

  /**
   * @brief Accessor for recompresion index
   * @return A const reference to arrayView1d<real64> containing the bulk
   *         modulus (at every element).
   */
  arrayView1d< real64 > const recompressionIndex() { return m_recompressionIndex; }

  /**
   * @brief Const accessor for recompression index
   * @return A const reference to arrayView1d<real64 const> containing the bulk
   *         modulus (at every element).
   */
  arrayView1d< real64 const > const recompressionIndex() const { return m_recompressionIndex; }

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
   * @brief Create a instantiation of the ElasticIsotropicPressureDependentUpdate class
   *        that refers to the data in this.
   * @param includeState Flag whether to pass state arrays that may not be needed for "no-state" updates
   * @return An instantiation of ElasticIsotropicPressureDependentUpdate.
   */
  ElasticIsotropicPressureDependentUpdates createKernelUpdates( bool const includeState = true ) const
  {
    if( includeState )
    {
      return ElasticIsotropicPressureDependentUpdates( m_refPressure,
                                                       m_refStrainVol,
                                                       m_recompressionIndex,
                                                       m_shearModulus,
                                                       m_newStress,
                                                       m_oldStress );
    }
    else // for "no state" updates, pass empty views to avoid transfer of stress data to device
    {
      return ElasticIsotropicPressureDependentUpdates( m_refPressure,
                                                       m_refStrainVol,
                                                       m_recompressionIndex,
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
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams )
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_refPressure,
                          m_refStrainVol,
                          m_recompressionIndex,
                          m_shearModulus,
                          m_newStress,
                          m_oldStress );
  }


protected:

  /// Post-process XML data
  virtual void postProcessInput() override;

  /// The default value of the bulk modulus for any new allocations.
  real64 m_defaultRefPressure;

  /// The default value of the shear modulus for any new allocations.
  real64 m_defaultRefStrainVol;

  /// The default value of the bulk modulus for any new allocations.
  real64 m_defaultRecompressionIndex;

  /// The default value of the shear modulus for any new allocations.
  real64 m_defaultShearModulus;

  /// The bulk modulus for each upper level dimension (i.e. cell) of *this
  real64 m_refPressure;

  /// The shear modulus for each upper level dimension (i.e. cell) of *this
  real64 m_refStrainVol;

  /// The bulk modulus for each upper level dimension (i.e. cell) of *this
  array1d< real64 > m_recompressionIndex;

  /// The shear modulus for each upper level dimension (i.e. cell) of *this
  array1d< real64 > m_shearModulus;

};

}

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_ELASTICISOTROPICPRESSUREDEPENDENT_HPP_ */
