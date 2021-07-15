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
 * @file CoupledSolid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_COUPLEDSOLID_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_COUPLEDSOLID_HPP_

#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/solid/SolidBase.hpp"

namespace geosx
{
namespace constitutive
{

template< typename SOLID_TYPE,
          typename PORO_TYPE,
          typename PERM_TYPE >
class CoupledSolidUpdates
{
public:
  /**
   * @brief Constructor
   */
  CoupledSolidUpdates( SOLID_TYPE const & solidModel,
                       PORO_TYPE const & porosityModel,
                       PERM_TYPE const & permModel ):
    m_solidUpdate( solidModel.createKernelUpdates() ),
    m_porosityUpdate( porosityModel.createKernelUpdates() ),
    m_permUpdate( permModel.createKernelWrapper() )
  {}

  GEOSX_HOST_DEVICE
  real64 getOldPorosity( localIndex const k,
                         localIndex const q ) const
  {
    return m_porosityUpdate.getOldPorosity( k, q );
  }

  GEOSX_HOST_DEVICE
  real64 getPorosity( localIndex const k,
                      localIndex const q ) const
  {
    return m_porosityUpdate.getPorosity( k, q );
  }

protected:
  typename SOLID_TYPE::KernelWrapper const m_solidUpdate;

  typename PORO_TYPE::KernelWrapper const m_porosityUpdate;

  typename PERM_TYPE::KernelWrapper const m_permUpdate;
};



/**
 * @brief Class to represent a coupled solid model
 */
template< typename SOLID_TYPE,
          typename PORO_TYPE,
          typename PERM_TYPE >
class CoupledSolid : public ConstitutiveBase
{
public:

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  CoupledSolid( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~CoupledSolid() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "Coupled" ) + SOLID_TYPE::m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

//  /// Post-process XML input
//  virtual void postProcessInput() override;

  struct viewKeyStruct
  {
    static constexpr char const * solidModelNameString() { return "solidModelName"; }
    static constexpr char const * porosityModelNameString() { return "porosityModelName"; }
    static constexpr char const * permeabilityModelNameString() { return "permeabilityModelName"; }
  };

  /**
   * @brief Create a instantiation of the PorousSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of PorousSolidUpdates.
   */
  CoupledSolidUpdates< SOLID_TYPE, PORO_TYPE, PERM_TYPE > createKernelUpdates() const
  {

    return CoupledSolidUpdates< SOLID_TYPE, PORO_TYPE, PERM_TYPE >( getSolidModel(),
                                                                    getPorosityModel(),
                                                                    getPermModel() );
  }

  arrayView2d< real64 const > const  dPorosity_dPressure() const
  { return getPorosityModel().dPorosity_dPressure(); }

  virtual void saveConvergedState() const override final
  {
    // getSolidModel().saveConvergedState();

    getPorosityModel().saveConvergedState();
  }

protected:

  SOLID_TYPE const & getSolidModel() const
  { return this->getParent().template getGroup< SOLID_TYPE >( m_solidModelName ); }

  PORO_TYPE const & getPorosityModel() const
  { return this->getParent().template getGroup< PORO_TYPE >( m_porosityModelName ); }

  PERM_TYPE const & getPermModel() const
  { return this->getParent().template getGroup< PERM_TYPE >( m_permeabilityModelName ); }

  // the name of the solid model
  string m_solidModelName;

  // the name of the porosity model
  string m_porosityModelName;

  // the name of the porosity model
  string m_permeabilityModelName;
};


template< typename SOLID_TYPE,
          typename PORO_TYPE,
          typename PERM_TYPE >
CoupledSolid< SOLID_TYPE, PORO_TYPE, PERM_TYPE >::CoupledSolid( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_solidModelName(),
  m_porosityModelName()
{
  registerWrapper( viewKeyStruct::solidModelNameString(), &m_solidModelName ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the solid model." );

  registerWrapper( viewKeyStruct::porosityModelNameString(), &m_porosityModelName ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the porosity model." );

  registerWrapper( viewKeyStruct::permeabilityModelNameString(), &m_permeabilityModelName ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the permeability model." );
}

template< typename SOLID_TYPE,
          typename PORO_TYPE,
          typename PERM_TYPE >
CoupledSolid< SOLID_TYPE, PORO_TYPE, PERM_TYPE >::~CoupledSolid()
{}


}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_COUPLEDSOLID_HPP_ */
