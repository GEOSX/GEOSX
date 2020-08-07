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
 * @file SolidBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_SOLIDBASE_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_SOLIDBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

namespace geosx
{
namespace constitutive
{


/**
 * @class SolidBaseUpdates
 * This class serves as a base class for all solid constitutive kernel wrapper
 * classes. The responsibility of this base is to:
 *
 * 1) Contain views to the constitutive state and parameter data for solid
 *    constitutive relations.
 * 2) Specify an interface for state update functions.
 *
 * In general, the ArrayView data in the wrapper is specified to be of type
 * "arrayView<T> const" or "arrayView<T const> const". The "const-ness"
 * of the data indicates the distinction from a parameter and a state variable,
 * with the parameters being "T const" and the state variables being "T".
 *
 * @note If an allocation occurs on the underlying Array after a KernelWrapper
 * is created, then the ArrayView members of that KernelWrapper are silently
 * invalid.
 */
class SolidBaseUpdates
{
protected:
  /**
   * @brief constructor
   * @param[in] stress The stress data from the constitutive model class.
   */
  SolidBaseUpdates( arrayView3d< real64, solid::STRESS_USD > const & stress ):
    m_stress( stress )
  {}


  /// Deleted default constructor
  SolidBaseUpdates() = delete;


  /**
   * @brief Copy Constructor
   * @param source Object to copy
   */
  SolidBaseUpdates( SolidBaseUpdates const & source ) = default;

  /**
   * @brief Move Constructor
   * @param source Object to move resources from
   */
  SolidBaseUpdates( SolidBaseUpdates && source ) = default;

  /// Deleted copy assignment operator
  SolidBaseUpdates & operator=( SolidBaseUpdates const & ) = delete;

  /// Deleted move assignment operator
  SolidBaseUpdates & operator=( SolidBaseUpdates && ) =  delete;

public:

  /// A reference the material stress at quadrature points.
  arrayView3d< real64, solid::STRESS_USD > const m_stress;

private:
  /**
   * @brief New interface proposal for small strain update
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] strainIncrement Strain increment in voight notation (linearized strain)
   * @param[out] stress New stress value (Cauchy stress)
   * @param[out] stiffness New tangent stiffness value
   */
  GEOSX_HOST_DEVICE
  virtual void SmallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( & strainIncrement )[6],
                                  real64 ( & stress )[6],
                                  real64 ( & stiffness )[6][6] )
  {
    GEOSX_UNUSED_VAR(k);
    GEOSX_UNUSED_VAR(q);
    GEOSX_UNUSED_VAR(strainIncrement);
    GEOSX_UNUSED_VAR(stress);
    GEOSX_UNUSED_VAR(stiffness);
    GEOSX_ERROR("SolidBase::SmallStrainUpdate() not implemented");
  }
  
  /**
   * @brief Save history variables in preparation for next timestep
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   */
  GEOSX_HOST_DEVICE
  virtual void SaveConvergedState( localIndex const k, localIndex const q )
  {
    GEOSX_UNUSED_VAR(k);
    GEOSX_UNUSED_VAR(q);
  }
  
  //////////////// LEGACY INTERFACE BELOW -- TO REVISIT ////////////////////////
  
  /**
   * accessor to return the stiffness at a given element
   * @param k the element number
   * @param c the stiffness array
   */
  GEOSX_HOST_DEVICE
  virtual void GetStiffness( localIndex const k, real64 ( &c )[6][6] ) const
  {
    GEOSX_UNUSED_VAR(k);
    GEOSX_UNUSED_VAR(c);
    GEOSX_ERROR("SolidBase::GetStiffness() not implemented");
  }

  /**
   * @brief Calculate stress using input generated under small strain
   *        assumptions.
   * @param[in] k The element index.
   * @param[in] voigtStrain The total strain tensor in Voigt notation.
   * @param[out] stress Pointer to the stress data in Voigt notation.
   */
  GEOSX_HOST_DEVICE
  virtual void SmallStrainNoState( localIndex const k,
                                   real64 const ( &voigtStrain )[ 6 ],
                                   real64 ( &stress )[ 6 ] ) const
  {
    GEOSX_UNUSED_VAR(k);
    GEOSX_UNUSED_VAR(voigtStrain);
    GEOSX_UNUSED_VAR(stress);
    GEOSX_ERROR("SolidBase::SmallStrainNoState not implemented");
  }

  /**
   * @brief Update the constitutive state using input generated under small
   *        strain assumptions.
   * @param[in] k The element index.
   * @param[in] q The quadrature point index.
   * @param[in] voigtStrainIncrement The increment in strain expressed in Voigt
   *                                 notation.
   */
  GEOSX_HOST_DEVICE
  virtual void SmallStrain( localIndex const k,
                            localIndex const q,
                            real64 const ( &voigtStrainInc )[ 6 ] ) const
  {
    GEOSX_UNUSED_VAR(k);
    GEOSX_UNUSED_VAR(q);
    GEOSX_UNUSED_VAR(voigtStrainInc);
    GEOSX_ERROR("SolidBase::SmallStrain not implemented");
  }

  /**
   * @brief Hypoelastic update to the constitutive state using input generated
   *        under finite strain assumptions.
   * @param[in] k The element index.
   * @param[in] q The quadrature point index.
   * @param[in] Ddt The incremental deformation tensor
   *                (rate of deformation tensor * dt)
   * @param[in] Rot The incremental rotation tensor
   */
  GEOSX_HOST_DEVICE
  virtual void HypoElastic( localIndex const k,
                            localIndex const q,
                            real64 const ( &Ddt )[ 6 ],
                            real64 const ( &Rot )[ 3 ][ 3 ] ) const
  {
    GEOSX_UNUSED_VAR(k);
    GEOSX_UNUSED_VAR(q);
    GEOSX_UNUSED_VAR(Ddt);
    GEOSX_UNUSED_VAR(Rot);
    GEOSX_ERROR("SolidBase::HypoElastic not implemented");
  }

  /**
   * @brief Hyper-elastic stress update
   * @param[in] k The element index.
   * @param[in] FmI The deformation gradient minus Identity
   * @param[out] stress Pointer to the stress data in Voigt notation.
   */
  GEOSX_HOST_DEVICE
  virtual void HyperElastic( localIndex const k,
                             real64 const (&FmI)[3][3],
                             real64 ( &stress )[ 6 ] ) const
  {
    GEOSX_UNUSED_VAR(k);
    GEOSX_UNUSED_VAR(FmI);
    GEOSX_UNUSED_VAR(stress);
    GEOSX_ERROR("SolidBase::HyperElastic() not implemented");
  }

  /**
   * @brief Hyper-elastic state update
   * @param[in] k The element index.
   * @param[in] q The quadrature point index.
   * @param[in] FmI The deformation gradient minus Identity
   */
  GEOSX_HOST_DEVICE
  virtual void HyperElastic( localIndex const k,
                             localIndex const q,
                             real64 const (&FmI)[3][3] ) const
  {
    GEOSX_UNUSED_VAR(k);
    GEOSX_UNUSED_VAR(q);
    GEOSX_UNUSED_VAR(FmI);
    GEOSX_ERROR("SolidBase::HyperElastic() not implemented");
  }

};


/**
 * @class SolidBase
 * This class serves as the base class for solid constitutive models.
 */
class SolidBase : public constitutive::ConstitutiveBase
{
public:
  /**
   * @brief Constructor
   * @param name Name of the SolidBase object in the repository.
   * @param parent The parent group of the SolidBase object.
   */
  SolidBase( string const & name,
             Group * const parent );

  /**
   * Destructor
   */
  virtual ~SolidBase() override;

  virtual void DeliverClone( string const & name,
                             Group * const parent,
                             std::unique_ptr< ConstitutiveBase > & clone ) const override;

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;
  
  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto defaultDensityString  = "defaultDensity";
    static constexpr auto densityString  = "density";
    static constexpr auto stressString = "stress";
  };

  /**
   * @name Accessors
   */
  ///@{

  /**
   * Getter for default density
   * @return The default density
   */
  real64 getDefaultDensity() const
  {
    return m_defaultDensity;
  }

  /**
   * Setter for default density
   * @param density The value to set m_defaultDensity equal to.
   */
  void setDefaultDensity( real64 const density )
  {
    m_defaultDensity = density;
  }

  /// Non-const/Mutable accessor for density.
  arrayView2d< real64 > const & getDensity()
  {
    return m_density;
  }

  /// Const/non-mutable accessor for density
  arrayView2d< real64 const > const & getDensity() const
  {
    return m_density;
  }

  /// Non-const/mutable accessor for stress
  arrayView3d< real64, solid::STRESS_USD > const & getStress()
  {
    return m_stress;
  }

  /// Const/non-mutable accessor for stress
  arrayView3d< real64 const, solid::STRESS_USD > const & getStress() const
  {
    return m_stress;
  }

  ///@}

protected:
  /// The default density for new allocations.
  real64 m_defaultDensity = 0;

  /// The material density at a quadrature point.
  array2d< real64 > m_density;

  /// The material stress at a quadrature point.
  array3d< real64, solid::STRESS_PERMUTATION > m_stress;
  
  /// band-aid fix...going to have to remove this after we clean up
  /// initialization for constitutive models.
  bool m_postProcessed = false;
};

} // namespace constitutive
} // namespace geosx

#endif /* GEOSX_CONSTITUTIVE_SOLID_SOLIDBASE_HPP_ */
