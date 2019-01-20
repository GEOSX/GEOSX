/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
  * @file ThreePhaseBakerRelativePermeability.cpp
  */

#include "ThreePhaseBakerRelativePermeability.hpp"

#include <cmath>

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{


ThreePhaseBakerRelativePermeability::ThreePhaseBakerRelativePermeability( std::string const & name,
                                                                          ManagedGroup * const parent )
  : RelativePermeabilityBase( name, parent )
{
  RegisterViewWrapper( viewKeyStruct::phaseMinVolumeFractionString, &m_phaseMinVolumeFraction, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Minimum volume fraction value for each phase");

  
  RegisterViewWrapper( viewKeyStruct::waterOilRelPermExponentString,   &m_waterOilRelPermExponent,   false )->
    setApplyDefaultValue(1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Rel perm power law exponent for the pair (water phase, oil phase) at residual gas saturation");

  RegisterViewWrapper( viewKeyStruct::waterOilRelPermMaxValueString,   &m_waterOilRelPermMaxValue,   false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Maximum rel perm value for the pair (water phase, oil phase) at residual gas saturation");


  RegisterViewWrapper( viewKeyStruct::gasOilRelPermExponentString,   &m_gasOilRelPermExponent,   false )->
    setApplyDefaultValue(1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Rel perm power law exponent for the pair (gas phase, oil phase) at residual water saturation");
  
  RegisterViewWrapper( viewKeyStruct::gasOilRelPermMaxValueString,   &m_gasOilRelPermMaxValue,   false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Maximum rel perm value for the pair (gas phase, oil phase) at residual water saturation");

}

ThreePhaseBakerRelativePermeability::~ThreePhaseBakerRelativePermeability()
{

}

std::unique_ptr<ConstitutiveBase>
ThreePhaseBakerRelativePermeability::DeliverClone(string const & name, ManagedGroup * const parent) const
{
  std::unique_ptr< ThreePhaseBakerRelativePermeability > clone = std::make_unique<ThreePhaseBakerRelativePermeability>( name, parent );

  clone->m_phaseNames = this->m_phaseNames;
  clone->m_phaseTypes = this->m_phaseTypes;
  clone->m_phaseOrder = this->m_phaseOrder;

  clone->m_phaseMinVolumeFraction = this->m_phaseMinVolumeFraction;
  
  clone->m_waterOilRelPermExponent = this->m_waterOilRelPermExponent;
  clone->m_waterOilRelPermMaxValue = this->m_waterOilRelPermMaxValue;
  
  clone->m_gasOilRelPermExponent   = this->m_gasOilRelPermExponent;
  clone->m_gasOilRelPermMaxValue   = this->m_gasOilRelPermMaxValue;

  clone->m_volFracScale = this->m_volFracScale;

  std::unique_ptr<ConstitutiveBase> rval = std::move( clone );
  return rval;
}


void ThreePhaseBakerRelativePermeability::PostProcessInput()
{
  RelativePermeabilityBase::PostProcessInput();

  localIndex const NP = numFluidPhases();

#define COREY_CHECK_INPUT_LENGTH( data, expected, attr ) \
  if (integer_conversion<localIndex>((data).size()) != integer_conversion<localIndex>(expected)) \
  { \
    GEOS_ERROR( "ThreePhaseBakerRelativePermeability: invalid number of entries in " \
                << (attr) << " attribute (" \
                << (data).size() << " given, " \
                << (expected) << " expected)"); \
  }

  COREY_CHECK_INPUT_LENGTH( m_phaseMinVolumeFraction,  NP,   viewKeyStruct::phaseMinVolumeFractionString );
    
  COREY_CHECK_INPUT_LENGTH( m_waterOilRelPermExponent, NP-1, viewKeyStruct::waterOilRelPermExponentString );
  COREY_CHECK_INPUT_LENGTH( m_waterOilRelPermMaxValue, NP-1, viewKeyStruct::waterOilRelPermMaxValueString );
  
  COREY_CHECK_INPUT_LENGTH( m_gasOilRelPermExponent,   NP-1, viewKeyStruct::gasOilRelPermExponentString );
  COREY_CHECK_INPUT_LENGTH( m_gasOilRelPermMaxValue,   NP-1, viewKeyStruct::gasOilRelPermMaxValueString );

#undef COREY_CHECK_INPUT_LENGTH

  m_volFracScale = 1.0;
  for (localIndex ip = 0; ip < NP; ++ip)
  {
    GEOS_ERROR_IF( m_phaseMinVolumeFraction[ip] < 0.0 || m_phaseMinVolumeFraction[ip] > 1.0,
                   "ThreePhaseBakerRelativePermeability: invalid phase min volume fraction value: " << m_phaseMinVolumeFraction[ip] );
    m_volFracScale -= m_phaseMinVolumeFraction[ip];
  }
  GEOS_ERROR_IF( m_volFracScale < 0.0, "ThreePhaseBakerRelativePermeability: sum of min volume fractions exceeds 1.0" );

  
  for (localIndex ip = 0; ip < NP-1; ++ip)
  {
    GEOS_ERROR_IF( m_waterOilRelPermExponent[ip] < 0.0,
	  	   "ThreePhaseBakerRelativePermeability: invalid water-oil exponent value: " << m_waterOilRelPermExponent[ip] );
    GEOS_ERROR_IF( m_waterOilRelPermMaxValue[ip] < 0.0 || m_waterOilRelPermMaxValue[ip] > 1.0,
  		   "ThreePhaseBakerRelativePermeability: invalid maximum value: " << m_waterOilRelPermMaxValue[ip] );

    GEOS_ERROR_IF( m_gasOilRelPermExponent[ip] < 0.0,
	  	   "ThreePhaseBakerRelativePermeability: invalid gas-oil exponent value: " << m_gasOilRelPermExponent[ip] );
    GEOS_ERROR_IF( m_gasOilRelPermMaxValue[ip] < 0.0 || m_gasOilRelPermMaxValue[ip] > 1.0,
  		   "ThreePhaseBakerRelativePermeability: invalid maximum value: " << m_gasOilRelPermMaxValue[ip] );
  }

}


void ThreePhaseBakerRelativePermeability::BatchUpdate( arrayView2d<real64 const> const & phaseVolumeFraction )
{

  arrayView1d<real64 const> const & phaseMinVolumeFraction = m_phaseMinVolumeFraction;
  
  arrayView1d<real64 const> const & waterOilRelPermExponent = m_waterOilRelPermExponent;
  arrayView1d<real64 const> const & waterOilRelPermMaxValue = m_waterOilRelPermMaxValue;
  
  arrayView1d<real64 const> const & gasOilRelPermExponent   = m_gasOilRelPermExponent;
  arrayView1d<real64 const> const & gasOilRelPermMaxValue   = m_gasOilRelPermMaxValue;

  RelativePermeabilityBase::BatchUpdateKernel<ThreePhaseBakerRelativePermeability>( phaseVolumeFraction,
										    m_phaseOrder,
                                                                                    phaseMinVolumeFraction,
                                                                                    waterOilRelPermExponent,
										    waterOilRelPermMaxValue,
										    gasOilRelPermExponent,
										    gasOilRelPermMaxValue,
                                                                                    m_volFracScale);
}


void ThreePhaseBakerRelativePermeability::PointUpdate( arraySlice1d<real64 const> const & phaseVolFraction,
                                                       localIndex const k,
                                                       localIndex const q )
{
  arraySlice1d<real64> const relPerm           = m_phaseRelPerm[k][q];
  arraySlice2d<real64> const dRelPerm_dVolFrac = m_dPhaseRelPerm_dPhaseVolFrac[k][q];

  localIndex const NP = numFluidPhases();

  Compute( NP,
           phaseVolFraction,
           relPerm,
           dRelPerm_dVolFrac,
	   m_phaseOrder,
           m_phaseMinVolumeFraction,
           m_waterOilRelPermExponent,
	   m_waterOilRelPermMaxValue,
	   m_gasOilRelPermExponent,
	   m_gasOilRelPermMaxValue,
           m_volFracScale );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThreePhaseBakerRelativePermeability, std::string const &, ManagedGroup * const )
} // namespace constitutive

} // namespace geosx
