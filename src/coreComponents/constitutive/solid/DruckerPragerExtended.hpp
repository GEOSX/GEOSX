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
 *  @file DruckerPragerExtended.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGEREXTENDED_HPP
#define GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGEREXTENDED_HPP

#include "ElasticIsotropic.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class DruckerPragerExtendedUpdates
 *
 * Class to provide material updates that may be called from a kernel function.
 */
class DruckerPragerExtendedUpdates : public ElasticIsotropicUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] stress The ArrayView holding the stress data for each quadrature point.
   */
  DruckerPragerExtendedUpdates( arrayView1d< real64 const > const & initialFriction,
                                arrayView1d< real64 const > const & residualFriction,
                                arrayView1d< real64 const > const & dilationRatio,
                                arrayView1d< real64 const > const & pressureIntercept,
                                arrayView1d< real64 const > const & hardening,
                                arrayView2d< real64 > const & newState,
                                arrayView2d< real64 > const & oldState,
                                arrayView2d< real64 > const & newPorosity,
                                arrayView2d< real64 > const & oldPorosity,
                                arrayView2d< real64 > const & dPorosity_dPressure,
                                arrayView1d< real64 > const & referencePorosity,
                                real64 const & grainBulkModulus,
                                arrayView2d< real64 > const & grainDensity,
                                arrayView3d< real64, solid::STRESS_USD > const & newStress,
                                arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                                arrayView1d< real64 const > const & bulkModulus,
                                arrayView1d< real64 const > const & shearModulus ):
    ElasticIsotropicUpdates( newPorosity,
                             oldPorosity,
                             dPorosity_dPressure,
                             referencePorosity,
                             grainBulkModulus,
                             grainDensity,
                             newStress,
                             oldStress,
                             bulkModulus,
                             shearModulus ),
    m_initialFriction( initialFriction ),
    m_residualFriction( residualFriction ),
    m_dilationRatio( dilationRatio ),
    m_pressureIntercept( pressureIntercept ),
    m_hardening( hardening ),
    m_newState( newState ),
    m_oldState( oldState )
  {}

  /// Default copy constructor
  DruckerPragerExtendedUpdates( DruckerPragerExtendedUpdates const & ) = default;

  /// Default move constructor
  DruckerPragerExtendedUpdates( DruckerPragerExtendedUpdates && ) = default;

  /// Deleted default constructor
  DruckerPragerExtendedUpdates() = delete;

  /// Deleted copy assignment operator
  DruckerPragerExtendedUpdates & operator=( DruckerPragerExtendedUpdates const & ) = delete;

  /// Deleted move assignment operator
  DruckerPragerExtendedUpdates & operator=( DruckerPragerExtendedUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotroipic; // TODO: typo in anistropic (fix in DiscOps PR)

  // Bring in base implementations to prevent hiding warnings
  using ElasticIsotropicUpdates::smallStrainUpdate;
  using ElasticIsotropicUpdates::smallStrainUpdate_porosity;

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

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate_porosity( localIndex const k,
                                           localIndex const q,
                                           real64 const & pressure,
                                           real64 const & deltaPressure,
                                           real64 const ( &strainIncrement )[6],
                                           real64 & porosity,
                                           real64 & dPorosity_dPressure,
                                           real64 & dPorosity_dVolStrainIncrement,
                                           real64 ( &stress )[6],
                                           real64 & dTotalStress_dPressure,
                                           DiscretizationOps & stiffness ) const final;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const override final
  {
    ElasticIsotropicUpdates::saveConvergedState( k, q );
    m_oldState[k][q] = m_newState[k][q];
  }

private:
  /// A reference to the ArrayView holding the initial friction param for each element.
  arrayView1d< real64 const > const m_initialFriction;

  /// A reference to the ArrayView holding the residual friction param for each element.
  arrayView1d< real64 const > const m_residualFriction;

  /// A reference to the ArrayView holding the dilation ratio for each element.
  arrayView1d< real64 const > const m_dilationRatio;

  /// A reference to the ArrayView holding the pressure intercept for each element.
  arrayView1d< real64 const > const m_pressureIntercept;

  /// A reference to the ArrayView holding the hardening parameter for each element.
  arrayView1d< real64 const > const m_hardening;

  /// A reference to the ArrayView holding the new state variable for each integration point
  arrayView2d< real64 > const m_newState;

  /// A reference to the ArrayView holding the old state variable for each integration point
  arrayView2d< real64 > const m_oldState;

  /// Hyperbolic model for friction hardening
  GEOSX_HOST_DEVICE
  void hyperbolicModel( real64 const y1,
                        real64 const y2,
                        real64 const m,
                        real64 const x,
                        real64 & y,
                        real64 & dy_dx ) const
  {
    if( x<0 )
    {
      y = y1;
      dy_dx = 0;
    }
    else
    {
      y = y1 + (y2-y1)*x/(m+x);
      dy_dx = (y2-y1)*m/((m+x)*(m+x));
    }
  };

};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerExtendedUpdates::smallStrainUpdate( localIndex const k,
                                                      localIndex const q,
                                                      real64 const ( &strainIncrement )[6],
                                                      real64 ( & stress )[6],
                                                      real64 ( & stiffness )[6][6] ) const
{
  // elastic predictor (assume strainIncrement is all elastic)

  ElasticIsotropicUpdates::smallStrainUpdate( k, q, strainIncrement, stress, stiffness );

  // decompose into mean (P) and von mises (Q) stress invariants

  real64 trialP;
  real64 trialQ;
  real64 deviator[6];

  twoInvariant::stressDecomposition( stress,
                                     trialP,
                                     trialQ,
                                     deviator );

  // check yield function F <= 0, using old state

  real64 friction, dfriction_dstate;

  hyperbolicModel( m_initialFriction[k],
                   m_residualFriction[k],
                   m_hardening[k],
                   m_oldState[k][q],
                   friction,
                   dfriction_dstate );

  real64 yield = trialQ + friction * (trialP - m_pressureIntercept[k]);

  if( yield < 1e-9 ) // elasticity
  {
    return;
  }

  // else, plasticity (trial stress point lies outside yield surface)
  // the return mapping can in general be written as a newton iteration.

  real64 solution[3], residual[3], delta[3];
  real64 jacobian[3][3] = {{}}, jacobianInv[3][3] = {{}};

  solution[0] = trialP; // initial guess for newP
  solution[1] = trialQ; // initial guess for newQ
  solution[2] = 1e-5;   // initial guess for plastic multiplier

  real64 norm, normZero = 1e30;

  // begin newton loop

  for( localIndex iter=0; iter<20; ++iter )
  {
    m_newState[k][q] = m_oldState[k][q] + solution[2];

    hyperbolicModel( m_initialFriction[k],
                     m_residualFriction[k],
                     m_hardening[k],
                     m_newState[k][q],
                     friction,
                     dfriction_dstate );

    // assemble residual system
    // resid1 = P - trialP + dlambda*bulkMod*dG/dP = 0
    // resid2 = Q - trialQ + dlambda*3*shearMod*dG/dQ = 0
    // resid3 = F = 0

    residual[0] = solution[0] - trialP + solution[2] * m_bulkModulus[k] * m_dilationRatio[k] * friction;
    residual[1] = solution[1] - trialQ + solution[2] * 3 * m_shearModulus[k];
    residual[2] = solution[1] + friction * (solution[0] - m_pressureIntercept[k]);

    // check for convergence

    norm = LvArray::tensorOps::l2Norm< 3 >( residual );

    if( iter==0 )
    {
      normZero = norm;
    }

    if( norm < 1e-8*(normZero+1))
    {
      break;
    }

    // solve Newton system

    jacobian[0][0] = 1;
    jacobian[0][2] = m_bulkModulus[k] * m_dilationRatio[k] * friction + solution[2] * m_bulkModulus[k] * m_dilationRatio[k] * dfriction_dstate;
    jacobian[1][1] = 1;
    jacobian[1][2] = 3 * m_shearModulus[k];
    jacobian[2][0] = friction;
    jacobian[2][1] = 1;
    jacobian[2][2] = dfriction_dstate * (solution[0]-m_pressureIntercept[k]);

    LvArray::tensorOps::invert< 3 >( jacobianInv, jacobian );
    LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( delta, jacobianInv, residual );

    for( localIndex i=0; i<3; ++i )
    {
      solution[i] -= delta[i];
    }
  }

  // re-construct stress = P*eye + sqrt(2/3)*Q*nhat

  twoInvariant::stressRecomposition( solution[0],
                                     solution[1],
                                     deviator,
                                     stress );

  // construct consistent tangent stiffness
  // note: if trialQ = 0, we will get a divide by zero error below,
  // but this is an unphysical (zero-strength) state anyway

  LvArray::tensorOps::fill< 6, 6 >( stiffness, 0 );

  real64 c1 = 2 * m_shearModulus[k] * solution[1] / trialQ;
  real64 c2 = jacobianInv[0][0] * m_bulkModulus[k] - c1 / 3;
  real64 c3 = sqrt( 2./3 ) * 3 * m_shearModulus[k] * jacobianInv[0][1];
  real64 c4 = sqrt( 2./3 ) * m_bulkModulus[k] * jacobianInv[1][0];
  real64 c5 = 2 * jacobianInv[1][1] * m_shearModulus[k] - c1;

  real64 identity[6];

  for( localIndex i=0; i<3; ++i )
  {
    stiffness[i][i] = c1;
    stiffness[i+3][i+3] = 0.5 * c1;
    identity[i] = 1.0;
    identity[i+3] = 0.0;
  }

  for( localIndex i=0; i<6; ++i )
  {
    for( localIndex j=0; j<6; ++j )
    {
      stiffness[i][j] +=   c2 * identity[i] * identity[j]
                         + c3 * identity[i] * deviator[j]
                         + c4 * deviator[i] * identity[j]
                         + c5 * deviator[i] * deviator[j];
    }
  }

  // save new stress and return
  saveStress( k, q, stress );
  return;
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerExtendedUpdates::smallStrainUpdate( localIndex const k,
                                                      localIndex const q,
                                                      real64 const ( &strainIncrement )[6],
                                                      real64 ( & stress )[6],
                                                      DiscretizationOps & stiffness ) const
{
  smallStrainUpdate( k, q, strainIncrement, stress, stiffness.m_c );
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerExtendedUpdates::smallStrainUpdate_porosity( localIndex const k,
                                                               localIndex const q,
                                                               real64 const & pressure,
                                                               real64 const & deltaPressure,
                                                               real64 const ( &strainIncrement )[6],
                                                               real64 & porosity,
                                                               real64 & dPorosity_dPressure,
                                                               real64 & dPorosity_dVolStrainIncrement,
                                                               real64 ( & stress )[6],
                                                               real64 & dTotalStress_dPressure,
                                                               DiscretizationOps & stiffness ) const
{
  GEOSX_UNUSED_VAR( pressure );
  GEOSX_UNUSED_VAR( deltaPressure );
  GEOSX_UNUSED_VAR( porosity );
  GEOSX_UNUSED_VAR( dPorosity_dPressure );
  GEOSX_UNUSED_VAR( dPorosity_dVolStrainIncrement );
  GEOSX_UNUSED_VAR( dTotalStress_dPressure );
  smallStrainUpdate( k, q, strainIncrement, stress, stiffness );

  //TODO fill in this.
}


/**
 * @class DruckerPragerExtended
 *
 * Extended Drucker-Prager material model.
 */
class DruckerPragerExtended : public ElasticIsotropic
{
public:

  /// @typedef Alias for DruckerPragerExtendedUpdates
  using KernelWrapper = DruckerPragerExtendedUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  DruckerPragerExtended( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~DruckerPragerExtended() override;


  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "ExtendedDruckerPrager";

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
    /// string/key for default initial friction angle
    static constexpr char const * defaultInitialFrictionAngleString() { return "defaultInitialFrictionAngle"; }

    /// string/key for default initial friction angle
    static constexpr char const * defaultResidualFrictionAngleString() { return "defaultResidualFrictionAngle"; }

    /// string/key for default dilation angle
    static constexpr char const * defaultDilationRatioString() { return "defaultDilationRatio"; }

    /// string/key for default hardening rate
    static constexpr char const * defaultHardeningString() { return "defaultHardening"; }

    /// string/key for default cohesion
    static constexpr char const * defaultCohesionString() { return "defaultCohesion"; }

    /// string/key for initial friction angle
    static constexpr char const * initialFrictionString() { return "initialFriction"; }

    /// string/key for final friction angle
    static constexpr char const * residualFrictionString() { return "residualFriction"; }

    /// string/key for dilation angle
    static constexpr char const * dilationRatioString() { return "dilationRatio"; }

    /// string/key for pressure intercept
    static constexpr char const * pressureInterceptString() { return "pressureIntercept"; }

    /// string/key for cohesion
    static constexpr char const * hardeningString() { return "hardening"; }

    /// string/key for state variable
    static constexpr char const * newStateString() { return "stateVariable"; }

    /// string/key for state variable
    static constexpr char const * oldStateString() { return "oldStateVariable"; }
  };

  /**
   * @brief Create a instantiation of the DruckerPragerExtendedUpdate class that refers to the data in this.
   * @return An instantiation of DruckerPragerExtendedUpdate.
   */
  DruckerPragerExtendedUpdates createKernelUpdates() const
  {
    return DruckerPragerExtendedUpdates( m_initialFriction,
                                         m_residualFriction,
                                         m_dilationRatio,
                                         m_pressureIntercept,
                                         m_hardening,
                                         m_newState,
                                         m_oldState,
                                         m_newPorosity,
                                         m_oldPorosity,
                                         m_dPorosity_dPressure,
                                         m_referencePorosity,
                                         m_grainBulkModulus,
                                         m_grainDensity,
                                         m_newStress,
                                         m_oldStress,
                                         m_bulkModulus,
                                         m_shearModulus );
  }

  /**
   * @brief Construct an update kernel for a derived type.
   * @tparam UPDATE_KERNEL The type of update kernel from the derived type.
   * @tparam PARAMS The parameter pack to hold the constructor parameters for the derived update kernel.
   * @param constructorParams The constructor parameter for the derived type.
   * @return An @p UPDATE_KERNEL object.
   */
  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams )
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_initialFriction,
                          m_residualFriction,
                          m_dilationRatio,
                          m_pressureIntercept,
                          m_hardening,
                          m_newState,
                          m_oldState,
                          m_newPorosity,
                          m_oldPorosity,
                          m_dPorosity_dPressure,
                          m_referencePorosity,
                          m_grainBulkModulus,
                          m_grainDensity,
                          m_newStress,
                          m_oldStress,
                          m_bulkModulus,
                          m_shearModulus );
  }


protected:
  virtual void postProcessInput() override;

  /// Material parameter: The default value of the initial yield surface slope
  real64 m_defaultInitialFrictionAngle;

  /// Material parameter: The default value of the final yield surface slope
  real64 m_defaultResidualFrictionAngle;

  /// Material parameter: The default value of the plastic potential slope ratio
  real64 m_defaultDilationRatio;

  /// Material parameter: The default value of the initial cohesion
  real64 m_defaultCohesion;

  /// Material parameter: The default value of the hardening rate
  real64 m_defaultHardening;

  /// Material parameter: The initial yield surface slope param for each element
  array1d< real64 > m_initialFriction;

  /// Material parameter: The final yield surface slope param for each element
  array1d< real64 > m_residualFriction;

  /// Material parameter: The plastic potential slope param (a ratio w.r.t. current yield surface)
  array1d< real64 > m_dilationRatio;

  /// Material parameter: The pressure intercept (location of cone vertex) for each element
  array1d< real64 > m_pressureIntercept;

  /// Material parameter: The hyperbolic hardening parameter for each element
  array1d< real64 > m_hardening;

  /// State variable: The current equivalent plastic shear strain for each quadrature point
  array2d< real64 > m_newState;

  /// State variable: The previous equivalent plastic shear strain for each quadrature point
  array2d< real64 > m_oldState;
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGEREXTENDED_HPP_ */
