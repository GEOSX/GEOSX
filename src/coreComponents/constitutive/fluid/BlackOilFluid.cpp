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
 * @file BlackOilFluid.cpp
 */

#include "BlackOilFluid.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/Path.hpp"

// PVTPackage includes
#include "MultiphaseSystem/BlackOilMultiphaseSystem.hpp"
#include "MultiphaseSystem/DeadOilMultiphaseSystem.hpp"

using namespace PVTPackage;

namespace geosx
{
using namespace dataRepository;

namespace constitutive
{
namespace
{
BlackOilFluid::FluidType
getBlackOilFluidType( string const & name )
{
  static std::map< string, BlackOilFluid::FluidType > const fluidTypes = {
    { "LiveOil", BlackOilFluid::FluidType::LiveOil },
    { "DeadOil", BlackOilFluid::FluidType::DeadOil },
  };
  auto const it = fluidTypes.find( name );
  GEOSX_ERROR_IF( it == fluidTypes.end(),
                  "Black-oil fluid type not supported by PVTPackage: " << name );
  return it->second;
}

}  // namespace

BlackOilFluid::BlackOilFluid( std::string const & name, Group * const parent ) :
  MultiFluidPVTPackageWrapper( name, parent )
{
  getWrapperBase( viewKeyStruct::componentMolarWeightString )
    ->setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::phaseNamesString )->setInputFlag( InputFlags::REQUIRED );

  registerWrapper( viewKeyStruct::surfaceDensitiesString, &m_surfaceDensities )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "List of surface densities for each phase" );

  registerWrapper( viewKeyStruct::tableFilesString, &m_tableFiles )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setRestartFlags( RestartFlags::NO_WRITE )
    ->setDescription( "List of filenames with input PVT tables" );

  registerWrapper( viewKeyStruct::fluidTypeString, &m_fluidTypeString )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Type of black-oil fluid (LiveOil/DeadOil)" );
}

BlackOilFluid::~BlackOilFluid()
{}

void
BlackOilFluid::DeliverClone( string const & name,
                             Group * const parent,
                             std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< BlackOilFluid >( name, parent );
  }

  MultiFluidPVTPackageWrapper::DeliverClone( name, parent, clone );
  BlackOilFluid & fluid = dynamicCast< BlackOilFluid & >( *clone );

  fluid.m_surfaceDensities = m_surfaceDensities;
  fluid.m_tableFiles = m_tableFiles;
  fluid.m_fluidTypeString = m_fluidTypeString;
  fluid.m_fluidType = m_fluidType;

  fluid.createFluid();
}

void
BlackOilFluid::PostProcessInput()
{
  // TODO maybe use different names?
  m_componentNames = m_phaseNames;

  MultiFluidPVTPackageWrapper::PostProcessInput();

  localIndex const NP = numFluidPhases();

#define BOFLUID_CHECK_INPUT_LENGTH( data, expected, attr )                     \
  if( LvArray::integerConversion< localIndex >( ( data ).size() ) !=           \
      LvArray::integerConversion< localIndex >( expected ) )                   \
  {                                                                            \
    GEOSX_ERROR( "BlackOilFluid: invalid number of entries in "                \
                 << ( attr ) << " attribute (" << ( data ).size() << "given, " \
                 << ( expected ) << " expected)" );                            \
  }

  BOFLUID_CHECK_INPUT_LENGTH( m_surfaceDensities,
                              NP,
                              viewKeyStruct::surfaceDensitiesString )
  BOFLUID_CHECK_INPUT_LENGTH( m_tableFiles, NP, viewKeyStruct::surfaceDensitiesString )

#undef BOFLUID_CHECK_INPUT_LENGTH

  m_fluidType = getBlackOilFluidType( m_fluidTypeString );
}

void
BlackOilFluid::createFluid()
{
  std::vector< PVTPackage::PHASE_TYPE > phases( m_phaseTypes.begin(),
                                                m_phaseTypes.end() );
  std::vector< std::string > tableFiles( m_tableFiles.begin(), m_tableFiles.end() );
  std::vector< double > densities( m_surfaceDensities.begin(),
                                   m_surfaceDensities.end() );
  std::vector< double > molarWeights( m_componentMolarWeight.begin(),
                                      m_componentMolarWeight.end() );

  switch( m_fluidType )
  {
    case FluidType::LiveOil:
    {
      m_fluid = std::make_unique< BlackOilMultiphaseSystem >( phases,
                                                              tableFiles,
                                                              densities,
                                                              molarWeights );
      break;
    }
    case FluidType::DeadOil:
    {
      m_fluid = std::make_unique< DeadOilMultiphaseSystem >( phases,
                                                             tableFiles,
                                                             densities,
                                                             molarWeights );
      break;
    }
    default:
    {
      GEOSX_ERROR( "Unknown fluid type" );
    }
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        BlackOilFluid,
                        std::string const &,
                        Group * const )
}  // namespace constitutive

}  // namespace geosx
