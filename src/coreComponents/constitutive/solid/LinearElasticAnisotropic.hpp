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
 *  @file LinearElasticAnisotropic.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_LINEARELASTICANISOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_LINEARELASTICANISOTROPIC_HPP_
#include "SolidBase.hpp"
#include "constitutive/ExponentialRelation.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{


namespace constitutive
{


/**
 * @class LinearElasticAnisotropicUpdates
 *
 * Class to provide fully anisotropic material updates that may be called from
 * a kernel function.
 *
 */
class LinearElasticAnisotropicUpdates : public SolidBaseUpdates
{
public:

  /**
   * @brief Constructor
   * @param [in] C The Voigt stiffness tensor
   * @param[in] stress The ArrayView holding the stress data for each quadrature
   *                   point.
   */
  LinearElasticAnisotropicUpdates( arrayView3d< real64 const, solid::STIFFNESS_USD > const & C,
                                   arrayView3d< real64, solid::STRESS_USD > const & stress ):
    SolidBaseUpdates( stress ),
    m_stiffnessView( C )
  {}

  /// Deleted default constructor
  LinearElasticAnisotropicUpdates() = delete;

  /// Default copy constructor
  LinearElasticAnisotropicUpdates( LinearElasticAnisotropicUpdates const & ) = default;

  /// Default move constructor
  LinearElasticAnisotropicUpdates( LinearElasticAnisotropicUpdates && ) = default;

  /// Deleted copy assignment operator
  LinearElasticAnisotropicUpdates & operator=( LinearElasticAnisotropicUpdates const & ) = delete;

  /// Deleted move assignment operator
  LinearElasticAnisotropicUpdates & operator=( LinearElasticAnisotropicUpdates && ) =  delete;
                            
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

  /**
   * accessor to return the stiffness at a given element
   * @param k the element number
   * @param c the stiffness array
   */
  GEOSX_HOST_DEVICE inline
  virtual void GetStiffness( localIndex const k, real64 (& c)[6][6] ) const override final
  {
    LvArray::tensorOps::copy< 6, 6 >( c, m_stiffnessView[ k ] );
  }

  /// A reference to the ArrayView holding the Voigt Stiffness tensor in each
  /// element.
  arrayView3d< real64 const, solid::STIFFNESS_USD > const m_stiffnessView;
};

GEOSX_FORCE_INLINE
GEOSX_HOST_DEVICE
void
LinearElasticAnisotropicUpdates::
  SmallStrainNoState( localIndex const k,
                      real64 const ( &voigtStrain )[ 6 ],
                      real64 ( & stress )[ 6 ] ) const
{
  for( localIndex i=0; i<6; ++i )
  {
    for( localIndex j=0; j<6; ++j )
    {
      stress[i] = stress[i] + m_stiffnessView( k, i, j ) * voigtStrain[j];
    }
  }
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
LinearElasticAnisotropicUpdates::
  SmallStrain( localIndex const k,
               localIndex const q,
               real64 const ( &voigtStrainInc )[ 6 ] ) const
{
  LvArray::tensorOps::plusAijBj< 6, 6 >( m_stress[ k ][ q ], m_stiffnessView[ k ], voigtStrainInc );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
LinearElasticAnisotropicUpdates::
  HypoElastic( localIndex const k,
               localIndex const q,
               real64 const ( &Ddt )[ 6 ],
               real64 const ( &Rot )[ 3 ][ 3 ] ) const
{
  for( localIndex j=0; j<3; ++j )
  {
    m_stress( k, q, 0 ) = m_stress( k, q, 0 ) + m_stiffnessView( k, 0, j ) * Ddt[ j ];
    m_stress( k, q, 1 ) = m_stress( k, q, 1 ) + m_stiffnessView( k, 1, j ) * Ddt[ j ];
    m_stress( k, q, 2 ) = m_stress( k, q, 2 ) + m_stiffnessView( k, 2, j ) * Ddt[ j ];
    m_stress( k, q, 3 ) = m_stress( k, q, 3 ) + m_stiffnessView( k, 3, j ) * Ddt[ j ];
    m_stress( k, q, 4 ) = m_stress( k, q, 4 ) + m_stiffnessView( k, 4, j ) * Ddt[ j ];
    m_stress( k, q, 5 ) = m_stress( k, q, 5 ) + m_stiffnessView( k, 5, j ) * Ddt[ j ];
  }

  for( localIndex j=3; j<6; ++j )
  {
    m_stress( k, q, 0 ) = m_stress( k, q, 0 ) + m_stiffnessView( k, 0, j ) * 2 * Ddt[ j ];
    m_stress( k, q, 1 ) = m_stress( k, q, 1 ) + m_stiffnessView( k, 1, j ) * 2 * Ddt[ j ];
    m_stress( k, q, 2 ) = m_stress( k, q, 2 ) + m_stiffnessView( k, 2, j ) * 2 * Ddt[ j ];
    m_stress( k, q, 3 ) = m_stress( k, q, 3 ) + m_stiffnessView( k, 3, j ) * 2 * Ddt[ j ];
    m_stress( k, q, 4 ) = m_stress( k, q, 4 ) + m_stiffnessView( k, 4, j ) * 2 * Ddt[ j ];
    m_stress( k, q, 5 ) = m_stress( k, q, 5 ) + m_stiffnessView( k, 5, j ) * 2 * Ddt[ j ];
  }

  real64 temp[ 6 ];
  LvArray::tensorOps::AikSymBklAjl< 3 >( temp, Rot, m_stress[ k ][ q ] );
  LvArray::tensorOps::copy< 6 >( m_stress[ k ][ q ], temp );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
LinearElasticAnisotropicUpdates::
  HyperElastic( localIndex const GEOSX_UNUSED_PARAM( k ),
                real64 const (&GEOSX_UNUSED_PARAM( FmI ))[3][3],
                real64 ( & )[ 6 ] ) const
{
  GEOSX_ERROR( "LinearElasticAnisotropicKernelWrapper::HyperElastic() is not implemented!" );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
LinearElasticAnisotropicUpdates::
  HyperElastic( localIndex const GEOSX_UNUSED_PARAM( k ),
                localIndex const GEOSX_UNUSED_PARAM( q ),
                real64 const (&GEOSX_UNUSED_PARAM( FmI ))[3][3] ) const
{
  GEOSX_ERROR( "LinearElasticAnisotropicKernelWrapper::HyperElastic() is not implemented!" );
}

/**
 * @class LinearElasticAnisotropic
 *
 * Class to provide a linear elastic anisotropic material response.
 */
class LinearElasticAnisotropic : public SolidBase
{
public:
  /// @typedef Alias for LinearElasticAnisotropicUpdates
  using KernelWrapper = LinearElasticAnisotropicUpdates;

  /**
   * @brief constructor
   * @param[in]name name of the instance in the catalog
   * @param[in]parent the group which contains this instance
   */
  LinearElasticAnisotropic( string const & name, Group * const parent );

  /**
   * Destructor
   */
  virtual ~LinearElasticAnisotropic() override;

  virtual void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "LinearElasticAnisotropic";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string CatalogName() { return m_catalogNameString; }

  virtual string GetCatalogName() override { return CatalogName(); }
  ///@}

  /**
   * @struct Set of "char const *" and keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for default stiffness tensor
    static constexpr auto defaultStiffnessString  = "defaultStiffness";

    /// string/key for stiffness tensor
    static constexpr auto stiffnessString= "stiffness";
  };

  /**
   * @brief Setter for the default stiffness
   * @param c Input c-array holding the components of the new default stiffness
   *          tensor.
   */
  void setDefaultStiffness( real64 const c[6][6] );


  /**
   * @brief Create a instantiation of the LinearElasticAnisotropicUpdate class
   *        that refers to the data in this.
   * @return An instantiation of LinearElasticAnisotropicUpdate.
   */
  LinearElasticAnisotropicUpdates createKernelUpdates() const
  {
    return LinearElasticAnisotropicUpdates( m_stiffness.toViewConst(),
                                            m_stress.toView() );
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
                          m_stiffness.toViewConst(),
                          m_stress.toView() );
  }

  /**
   * @brief Const Getter for stiffness tensor
   * @return ArrayView to the stiffness tensor
   */
  arrayView3d< real64 const, solid::STIFFNESS_USD > const & getStiffness() const
  {
    return m_stiffness;
  }

  /**
   * @brief Non-const Getter for stiffness tensor
   * @return ArrayView to the stiffness tensor
   */
  arrayView3d< real64, solid::STIFFNESS_USD > const & getStiffness()
  {
    return m_stiffness;
  }

protected:
  virtual void PostProcessInput() override;

private:
  /// default value for stiffness tensor
  array2d< real64 > m_defaultStiffness;

  /// stiffness tensor in Voigt notation.
  array3d< real64, solid::STIFFNESS_PERMUTATION > m_stiffness;

};


}

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_LINEARELASTICANISOTROPIC_HPP_ */
