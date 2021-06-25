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

/*
 * @file WellControls.cpp
 */

#include "WellControls.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"

namespace geosx
{

using namespace dataRepository;

WellControls::WellControls( string const & name, Group * const parent )
  : Group( name, parent ),
  m_type( Type::PRODUCER ),
  m_refElevation( 0.0 ),
  m_refGravCoef( 0.0 ),
  m_currentControl( Control::BHP ),
  m_targetBHP( 0.0 ),
  m_targetTotalRate( 0.0 ),
  m_targetPhaseRate( 0.0 ),
  m_targetPhaseName( "" ),
  m_useSurfaceConditions( 0 ),
  m_surfacePres( 0.0 ),
  m_surfaceTemp( 0.0 ),
  m_targetTotalRateTableName( "" ),
  m_targetPhaseRateTableName( "" ),
  m_targetBHPTableName( "" ),
  m_targetTotalRateTable( nullptr ),
  m_targetPhaseRateTable( nullptr ),
  m_targetBHPTable( nullptr )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::typeString(), &m_type ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Well type. Valid options:\n* " + EnumStrings< Type >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::controlString(), &m_currentControl ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Well control. Valid options:\n* " + EnumStrings< Control >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::targetBHPString(), &m_targetBHP ).
    setDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Target bottom-hole pressure" );

  registerWrapper( viewKeyStruct::targetTotalRateString(), &m_targetTotalRate ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Target total volumetric rate" );

  registerWrapper( viewKeyStruct::targetPhaseRateString(), &m_targetPhaseRate ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Target phase volumetric rate" );

  registerWrapper( viewKeyStruct::targetPhaseNameString(), &m_targetPhaseName ).
    setDefaultValue( "" ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the target phase" );

  registerWrapper( viewKeyStruct::refElevString(), &m_refElevation ).
    setDefaultValue( -1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference elevation where BHP control is enforced" );

  registerWrapper( viewKeyStruct::injectionStreamString(), &m_injectionStream ).
    setDefaultValue( -1 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Global component densities for the injection stream" );

  registerWrapper( viewKeyStruct::useSurfaceConditionsString(), &m_useSurfaceConditions ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to specify whether rates are checked at surface or reservoir conditions.\n"
                    "Equal to 1 for surface conditions, and to 0 for reservoir conditions" );

  registerWrapper( viewKeyStruct::surfacePressureString(), &m_surfacePres ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Surface pressure used to compute volumetric rates when surface conditions are used" );

  registerWrapper( viewKeyStruct::surfaceTemperatureString(), &m_surfaceTemp ).
    setDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Surface temperature used to compute volumetric rates when surface conditions are used" );

  registerWrapper( viewKeyStruct::targetBHPTableNameString(), &m_targetBHPTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the BHP table when the rate is a time dependent function" );

  registerWrapper( viewKeyStruct::targetTotalRateTableNameString(), &m_targetTotalRateTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the rate table when the rate is a time dependent function" );

  registerWrapper( viewKeyStruct::targetPhaseRateTableNameString(), &m_targetPhaseRateTableName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the phase rate table when the rate is a time dependent function" );
}


WellControls::~WellControls()
{}

void WellControls::switchToBHPControl( real64 const & val )
{
  m_currentControl = Control::BHP;
  m_targetBHP = val;
}

void WellControls::switchToTotalRateControl( real64 const & val )
{
  m_currentControl = Control::TOTALVOLRATE;
  m_targetTotalRate = val;
}

void WellControls::switchToPhaseRateControl( real64 const & val )
{
  m_currentControl = Control::PHASEVOLRATE;
  m_targetPhaseRate = val;
}

void WellControls::postProcessInput()
{
  // 3.a) check target BHP
  if( m_targetBHP < 0 )
  {
    GEOSX_ERROR( "Target bottom-hole pressure for well "<< getName() << " is negative" );
  }

  // 3.b) check target rates
  if( m_targetTotalRate < 0 )
  {
    GEOSX_ERROR( "Target rate for well "<< getName() << " is negative" );
  }
  if( m_targetPhaseRate < 0 )
  {
    GEOSX_ERROR( "Target oil rate for well "<< getName() << " is negative" );
  }

  // 4) check injection stream
  if( !m_injectionStream.empty())
  {
    real64 sum = 0.0;
    for( localIndex ic = 0; ic < m_injectionStream.size(); ++ic )
    {
      GEOSX_ERROR_IF( m_injectionStream[ic] < 0.0 || m_injectionStream[ic] > 1.0,
                      "Invalid injection stream for well " << getName() );
      sum += m_injectionStream[ic];
    }
    GEOSX_ERROR_IF( std::abs( 1.0 - sum ) > std::numeric_limits< real64 >::epsilon(),
                    "Invalid injection stream for well " << getName() );
  }

  // 5) check the flag for surface / reservoir conditions
  GEOSX_ERROR_IF( m_useSurfaceConditions != 0 && m_useSurfaceConditions != 1,
                  "The flag to select surface/reservoir conditions must be equal to 0 or 1" );

  // 6) check the flag for surface / reservoir conditions
  GEOSX_ERROR_IF( m_useSurfaceConditions == 1 && m_surfacePres <= 0,
                  "When useSurfaceConditions == 1, the surface pressure must be defined" );

  // 7) check that at least one rate constraint has been defined
  GEOSX_ERROR_IF( ((m_targetPhaseRate <= 0.0 && m_targetPhaseRateTableName == "") &&
                   (m_targetTotalRate <= 0.0 && m_targetTotalRateTableName == "")),
                  "You need to specify a phase rate constraint or a total rate constraint for injectors" );

  // 8) check whether redundant information has been provided
  GEOSX_ERROR_IF( ((m_targetPhaseRate > 0.0 && m_targetPhaseRateTableName != "")),
                  "You are provided redundant information for well phase rate" );

  GEOSX_ERROR_IF( ((m_targetTotalRate > 0.0 && m_targetTotalRateTableName != "")),
                  "You are provided redundant information for well total rate" );

  GEOSX_ERROR_IF( ((m_targetBHP > 0.0 && m_targetBHPTableName != "")),
                  "You are provided redundant information for well BHP" );

  GEOSX_ERROR_IF( ((m_targetBHP <= 0.0 && m_targetBHPTableName == "")),
                  "You have to provide well BHP by specifying either " << viewKeyStruct::targetBHPString() << " or " << viewKeyStruct::targetBHPTableNameString() );

  //  9) Create time-dependent BHP table
  FunctionManager & functionManager = FunctionManager::getInstance();
  if( m_targetBHPTableName != "" )
  {
    m_targetBHPTable = &(functionManager.getGroup< TableFunction >( m_targetBHPTableName ));
  }
  else
  {
    array1d< array1d< real64 > > timeCoord;
    timeCoord.resize( 1 );
    timeCoord[0].emplace_back( 0 );
    array1d< real64 > constantBHPValue;
    constantBHPValue.emplace_back( m_targetBHP );

    m_targetBHPTable = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", getName()+"constantBHPTable" ));
    m_targetBHPTable->setTableCoordinates( timeCoord );
    m_targetBHPTable->setTableValues( constantBHPValue );
    m_targetBHPTable->reInitializeFunction();
    m_targetBHPTable->setInterpolationMethod( TableFunction::InterpolationType::Lower );
  }
  //  10) Create time-dependent total rate table
  if( m_targetTotalRateTableName != "" )
  {
    m_targetTotalRateTable = &(functionManager.getGroup< TableFunction >( m_targetTotalRateTableName ));
    GEOSX_THROW_IF( m_targetTotalRateTable->getInterpolationMethod() != TableFunction::InterpolationType::Lower,
                    "The interpolation method for the time-dependent rate table " << m_targetTotalRateTable->getName() << " should be TableFunction::InterpolationType::Lower",
                    InputError );
  }
  else
  {
    array1d< array1d< real64 > > timeCoord;
    timeCoord.resize( 1 );
    timeCoord[0].emplace_back( 0 );
    array1d< real64 > constantTotalRateValue;
    constantTotalRateValue.emplace_back( m_targetTotalRate );

    m_targetTotalRateTable = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", getName()+"constantTotalRateTable" ));
    m_targetTotalRateTable->setTableCoordinates( timeCoord );
    m_targetTotalRateTable->setTableValues( constantTotalRateValue );
    m_targetTotalRateTable->reInitializeFunction();
    m_targetTotalRateTable->setInterpolationMethod( TableFunction::InterpolationType::Lower );

  }

  //  11) Create time-dependent phase rate table
  if( m_targetPhaseRateTableName != "" )
  {
    m_targetPhaseRateTable = &(functionManager.getGroup< TableFunction >( m_targetPhaseRateTableName ));
    GEOSX_THROW_IF( m_targetPhaseRateTable->getInterpolationMethod() != TableFunction::InterpolationType::Lower,
                    "The interpolation method for the time-dependent rate table " << m_targetPhaseRateTable->getName() << " should be TableFunction::InterpolationType::Lower",
                    InputError );
  }
  else
  {
    array1d< array1d< real64 > > timeCoord;
    timeCoord.resize( 1 );
    timeCoord[0].emplace_back( 0 );
    array1d< real64 > constantPhaseRateValue;
    constantPhaseRateValue.emplace_back( m_targetPhaseRate );

    m_targetPhaseRateTable = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", getName()+"constantPhaseRateTable" ));
    m_targetPhaseRateTable->setTableCoordinates( timeCoord );
    m_targetPhaseRateTable->setTableValues( constantPhaseRateValue );
    m_targetPhaseRateTable->reInitializeFunction();
    m_targetPhaseRateTable->setInterpolationMethod( TableFunction::InterpolationType::Lower );
  }

  // For producer, the rate is timed by -1 to reflect mass extraction from the grid cell
  if( getType() == Type::PRODUCER )
  {
    array1d< real64 > & PhaseRatetableValues = m_targetPhaseRateTable->getValues();
    for( localIndex i = 0; i < PhaseRatetableValues.size(); ++i )
    {
      PhaseRatetableValues( i ) *= -1;
    }

    array1d< real64 > & TotalRatetableValues = m_targetTotalRateTable->getValues();
    for( localIndex i = 0; i < TotalRatetableValues.size(); ++i )
    {
      TotalRatetableValues( i ) *= -1;
    }
  }

}

void WellControls::initializePostInitialConditionsPreSubGroups()
{
//   // for a producer, the solvers compute negative rates, so we adjust the input here
//   if( getType() == Type::PRODUCER
//       && (m_targetPhaseRate > 0.0 || m_targetTotalRate > 0.0) )
//   {
//     m_targetTotalRate *= -1;
//     m_targetPhaseRate *= -1;
//   }
}

} //namespace geosx
