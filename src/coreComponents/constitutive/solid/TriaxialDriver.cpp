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
 * @file TriaxialDriver.cpp
 */

#include "TriaxialDriver.hpp"
#include "fileIO/Outputs/OutputBase.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;


TriaxialDriver::TriaxialDriver( const string & name,
                                Group * const parent ):
  TaskBase( name, parent )
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::solidMaterialNameString(), &m_solidMaterialName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Solid material to test" );

  registerWrapper( viewKeyStruct::modeString(), &m_mode ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Test mode [stressControl, strainControl, mixedControl]" );

  registerWrapper( viewKeyStruct::axialFunctionString(), &m_axialFunctionName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Function controlling axial stress or strain (depending on test mode)" );

  registerWrapper( viewKeyStruct::radialFunctionString(), &m_radialFunctionName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Function controlling radial stress or strain (depending on test mode)" );

  registerWrapper( viewKeyStruct::initialStressString(), &m_initialStress ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Initial stress (scalar used to set an isotropic stress state)" );

  registerWrapper( viewKeyStruct::numStepsString(), &m_numSteps ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Number of load steps to take" );

  registerWrapper( viewKeyStruct::outputString(), &m_outputFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "none" ).
    setDescription( "Output file" );

  registerWrapper( viewKeyStruct::baselineString(), &m_baselineFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "none" ).
    setDescription( "Baseline file" );
}


TriaxialDriver::~TriaxialDriver()
{}


void TriaxialDriver::postProcessInput()
{

  GEOSX_THROW_IF( m_mode != "stressControl" && m_mode != "strainControl" && m_mode != "mixedControl",
                  "Test mode \'" << m_mode << "\' not recognized.",
                  InputError );

  // initialize table functions

  FunctionManager & functionManager = FunctionManager::getInstance();

  TableFunction & axialFunction = functionManager.getGroup< TableFunction >( m_axialFunctionName );
  TableFunction & radialFunction = functionManager.getGroup< TableFunction >( m_radialFunctionName );

  axialFunction.initializeFunction();
  radialFunction.initializeFunction();

  // determine time increment

  ArrayOfArraysView< real64 > coordinates = axialFunction.getCoordinates();
  real64 const minTime = coordinates[0][0];
  real64 const maxTime = coordinates[0][coordinates.sizeOfArray( 0 )-1];
  real64 const dt = (maxTime-minTime) / m_numSteps;

  // resize data arrays

  integer const length = m_numSteps+1;
  m_table.resize( length, m_numColumns );

  // set time column

  for( integer n=0; n<length; ++n )
  {
    m_table( n, TIME ) = minTime + n*dt;
  }

  // initial stress is always isotropic

  m_table( 0, SIG0 ) = m_initialStress;
  m_table( 0, SIG1 ) = m_initialStress;
  m_table( 0, SIG2 ) = m_initialStress;

  // preset certain columns depending on testing mode:
  //   mixedControl .... specified axial strain and radial stress
  //   strainControl ... specified axial and radial strain
  //   stressControl ... specified axial and radial stress

  for( integer n=0; n<length; ++n )
  {
    real64 axi = axialFunction.evaluate( &m_table( n, TIME ) );
    real64 rad = radialFunction.evaluate( &m_table( n, TIME ) );

    if( m_mode == "mixedControl" )
    {
      m_table( n, EPS0 ) = axi;
      m_table( n, SIG1 ) = rad;
      m_table( n, SIG2 ) = rad;
    }
    else if( m_mode == "strainControl" )
    {
      m_table( n, EPS0 ) = axi;
      m_table( n, EPS1 ) = rad;
      m_table( n, EPS2 ) = rad;
    }
    else if( m_mode == "stressControl" )
    {
      m_table( n, SIG0 ) = axi;
      m_table( n, SIG1 ) = rad;
      m_table( n, SIG2 ) = rad;
    }
  }

  // double check the initial stress value is consistent with any function values that
  // may overwrite it.

  GEOSX_THROW_IF( !isEqual( m_initialStress, m_table( 0, SIG0 ), 1e-6 ),
                  "Initial stress values indicated by initialStress and axialFunction(time=0) appear inconsistent",
                  InputError );

  GEOSX_THROW_IF( !isEqual( m_initialStress, m_table( 0, SIG1 ), 1e-6 ),
                  "Initial stress values indicated by initialStress and radialFunction(time=0) appear inconsistent",
                  InputError );
}


template< typename SOLID_TYPE >
void TriaxialDriver::runStrainControlTest( SOLID_TYPE & solid, arrayView2d< real64 > & table )
{
  typename SOLID_TYPE::KernelWrapper updates = solid.createKernelUpdates();

  forAll< parallelDevicePolicy<> >( 1, [=]  GEOSX_HOST_DEVICE ( integer const ei )
  {
    real64 stress[6] = {};
    real64 strainIncrement[6] = {};
    real64 stiffness[6][6] = {{}};

    for( integer n=1; n<=m_numSteps; ++n )
    {
      strainIncrement[0] = table( n, EPS0 )-table( n-1, EPS0 );
      strainIncrement[1] = table( n, EPS1 )-table( n-1, EPS1 );
      strainIncrement[2] = table( n, EPS2 )-table( n-1, EPS2 );

      updates.smallStrainUpdate( ei, 0, strainIncrement, stress, stiffness );
      updates.saveConvergedState ( ei, 0 );

      table( n, SIG0 ) = stress[0];
      table( n, SIG1 ) = stress[1];
      table( n, SIG2 ) = stress[2];

      table( n, ITER ) = 0;
    }
  } );
}


template< typename SOLID_TYPE >
void TriaxialDriver::runMixedControlTest( SOLID_TYPE & solid, arrayView2d< real64 > & table )
{
  typename SOLID_TYPE::KernelWrapper updates = solid.createKernelUpdates();

  forAll< parallelDevicePolicy<> >( 1, [=]  GEOSX_HOST_DEVICE ( integer const ei )
  {
    real64 stress[6] = {};
    real64 strainIncrement[6] = {};
    real64 deltaStrainIncrement = 0;
    real64 stiffness[6][6] = {{}};

    real64 scale = 0;
    for( integer n=1; n<=m_numSteps; ++n )
    {
      scale += fabs( table( n, SIG0 )) + fabs( table( n, SIG1 )) + fabs( table( n, SIG2 ));
    }
    scale = 3*m_numSteps / scale;

    for( integer n=1; n<=m_numSteps; ++n )
    {
      strainIncrement[0] = table( n, EPS0 )-table( n-1, EPS0 );
      strainIncrement[1] = 0;
      strainIncrement[2] = 0;

      real64 norm, normZero = 1e30;
      integer k = 0;
      integer cuts = 0;

      for(; k<m_maxIter; ++k )
      {
        updates.smallStrainUpdate( ei, 0, strainIncrement, stress, stiffness );

        norm = scale*fabs( stress[1]-table( n, SIG1 ) );

        if( k == 0 )
        {
          normZero = norm;
        }

        if( norm < m_newtonTol ) // success
        {
          break;
        }
        else if( k > 0 && norm > normZero && cuts < m_maxCuts ) // backtrack by half delta
        {
          cuts++;
          deltaStrainIncrement *= 0.5;
          strainIncrement[1] += deltaStrainIncrement;
          strainIncrement[2]  = strainIncrement[1];
        }
        else // newton update
        {
          deltaStrainIncrement  = (stress[1]-table( n, SIG1 )) / (stiffness[1][1]+stiffness[1][2]);
          strainIncrement[1]   -= deltaStrainIncrement;
          strainIncrement[2]    = strainIncrement[1];
        }
      }

      updates.saveConvergedState ( ei, 0 );

      table( n, SIG0 ) = stress[0];
      table( n, EPS1 ) = table( n-1, EPS1 )+strainIncrement[1];
      table( n, EPS2 ) = table( n, EPS1 );

      table( n, ITER ) = k;
      table( n, NORM ) = norm;

      if( norm > m_newtonTol )
      {
        break;
      }
    }
  } );
}


template< typename SOLID_TYPE >
void TriaxialDriver::runStressControlTest( SOLID_TYPE & solid, arrayView2d< real64 > & table )
{
  typename SOLID_TYPE::KernelWrapper updates = solid.createKernelUpdates();

  forAll< parallelDevicePolicy<> >( 1, [=]  GEOSX_HOST_DEVICE ( integer const ei )
  {
    real64 stress[6] = {};
    real64 strainIncrement[6] = {};
    real64 deltaStrainIncrement[6] = {};
    real64 stiffness[6][6] = {{}};

    real64 resid[2] = {};
    real64 jacobian[2][2] = {{}};

    real64 scale = 0;
    for( integer n=1; n<=m_numSteps; ++n )
    {
      scale += fabs( table( n, SIG0 )) + fabs( table( n, SIG1 )) + fabs( table( n, SIG2 ));
    }
    scale = 3*m_numSteps / scale;

    for( integer n=1; n<=m_numSteps; ++n )
    {
      //   std::cout<<"time step="<<n<<std::endl;
      strainIncrement[0] = 0;
      strainIncrement[1] = 0;
      strainIncrement[2] = 0;

      real64 norm, normZero = 1e30, det;
      integer k = 0;
      integer cuts = 0;

      for(; k<m_maxIter; ++k )
      {
        updates.smallStrainUpdate( ei, 0, strainIncrement, stress, stiffness );

        resid[0] = scale * (stress[0]-table( n, SIG0 ));
        resid[1] = scale * (stress[1]-table( n, SIG1 ));

        norm = sqrt( resid[0]*resid[0] + resid[1]*resid[1] );
        //  std::cout<<"k= "<<k<<std::endl;
        // std::cout<<"norm ="<<norm<<std::endl;

        if( k == 0 )
        {
          normZero = norm;
        }

        if( norm < m_newtonTol ) // success
        {
          break;
        }
        else if( k > 0 && norm > normZero && cuts < m_maxCuts ) // backtrack by half delta
        {
          cuts++;
          deltaStrainIncrement[0] *= 0.5;
          deltaStrainIncrement[1] *= 0.5;
          strainIncrement[0] += deltaStrainIncrement[0];
          strainIncrement[1] += deltaStrainIncrement[1];
          strainIncrement[2]  = strainIncrement[1];
          //  std::cout<<"k="<<k<<" , cuts="<<cuts<<std::endl;
        }
        else // newton update
        {
          cuts = 0;
          jacobian[0][0] = scale*stiffness[0][0];
          jacobian[1][0] = scale*stiffness[1][0];
          jacobian[0][1] = scale*(stiffness[0][1] + stiffness[0][2]);
          jacobian[1][1] = scale*(stiffness[1][1] + stiffness[1][2]);

          det = jacobian[0][0]*jacobian[1][1]-jacobian[0][1]*jacobian[1][0];

          deltaStrainIncrement[0] = (jacobian[1][1]*resid[0]-jacobian[0][1]*resid[1] ) / det;
          deltaStrainIncrement[1] = (jacobian[0][0]*resid[1]-jacobian[1][0]*resid[0] ) / det;

          strainIncrement[0] -= deltaStrainIncrement[0];
          strainIncrement[1] -= deltaStrainIncrement[1];
          strainIncrement[2]  = strainIncrement[1];
        }
      }

      updates.saveConvergedState ( ei, 0 );

      table( n, EPS0 ) = table( n-1, EPS0 )+strainIncrement[0];
      table( n, EPS1 ) = table( n-1, EPS1 )+strainIncrement[1];
      table( n, EPS2 ) = table( n, EPS1 );

      table( n, ITER ) = k;
      table( n, NORM ) = norm;

      if( norm > m_newtonTol )
      {
        break;
      }
    }
  } );
}


bool TriaxialDriver::execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                              real64 const GEOSX_UNUSED_PARAM( dt ),
                              integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                              integer const GEOSX_UNUSED_PARAM( eventCounter ),
                              real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                              DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  // this code only makes sense in serial

  GEOSX_THROW_IF( MpiWrapper::commRank() > 0, "Triaxial Driver should only be run in serial", std::runtime_error );

  // get the solid out of the constitutive manager.
  // for the moment it is of type SolidBase.

  ConstitutiveManager & constitutiveManager = this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );

  SolidBase & baseSolid = constitutiveManager.getGroup< SolidBase >( m_solidMaterialName );

  // depending on logLevel, print some useful info

  if( getLogLevel() > 0 )
  {
    GEOSX_LOG_RANK_0( "Launching Triaxial Driver" );
    GEOSX_LOG_RANK_0( "  Material .......... " << m_solidMaterialName );
    GEOSX_LOG_RANK_0( "  Type .............. " << baseSolid.getCatalogName() );
    GEOSX_LOG_RANK_0( "  Mode .............. " << m_mode );
    GEOSX_LOG_RANK_0( "  Axial Control ..... " << m_axialFunctionName );
    GEOSX_LOG_RANK_0( "  Radial Control .... " << m_radialFunctionName );
    GEOSX_LOG_RANK_0( "  Initial Stress .... " << m_initialStress );
    GEOSX_LOG_RANK_0( "  Steps ............. " << m_numSteps );
    GEOSX_LOG_RANK_0( "  Output ............ " << m_outputFile );
    GEOSX_LOG_RANK_0( "  Baseline .......... " << m_baselineFile );
  }

  // create a dummy discretization with one quadrature point for
  // storing constitutive data

  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  dataRepository::Group discretization( "discretization", &rootGroup );

  discretization.resize( 1 );   // one element
  baseSolid.allocateConstitutiveData( discretization, 1 );   // one quadrature point

  // set the initial stress state using the data table

  arrayView3d< real64, solid::STRESS_USD > stressArray = baseSolid.getStress();

  stressArray( 0, 0, 0 ) = m_table( 0, SIG0 );
  stressArray( 0, 0, 1 ) = m_table( 0, SIG1 );
  stressArray( 0, 0, 2 ) = m_table( 0, SIG2 );

  baseSolid.saveConvergedState();

  // pass the solid through the ConstitutivePassThru to downcast from the
  // base type to a known model type.  the lambda here then executes the
  // appropriate test driver. note that these calls will move data to device if available.

  ConstitutivePassThru< SolidBase >::execute( baseSolid, [&]( auto & selectedSolid )
  {
    using SOLID_TYPE = TYPEOFREF( selectedSolid );

    if( m_mode == "mixedControl" )
    {
      runMixedControlTest< SOLID_TYPE >( selectedSolid, m_table );
    }
    else if( m_mode == "strainControl" )
    {
      runStrainControlTest< SOLID_TYPE >( selectedSolid, m_table );
    }
    else if( m_mode == "stressControl" )
    {
      runStressControlTest< SOLID_TYPE >( selectedSolid, m_table );
    }
  } );

  // move table back to host for output
  m_table.move( LvArray::MemorySpace::host );

  validateResults();

  if( m_outputFile != "none" )
  {
    outputResults();
  }

  if( m_baselineFile != "none" )
  {
    compareWithBaseline();
  }

  return false;
}


void TriaxialDriver::validateResults()
{
  for( integer n=0; n<m_numSteps; ++n )
  {
    if( m_table( n, NORM ) > m_newtonTol )
    {
      GEOSX_LOG_RANK_0( "WARNING: Material driver failed to converge at loadstep " << n << "." );
      GEOSX_LOG_RANK_0( "         This usually indicates the material has completely failed and/or the loading state is inadmissible." );
      GEOSX_LOG_RANK_0( "         In rare cases, it may indicate a problem in the material model implementation." );

      for( integer col=EPS0; col<ITER; ++col )
      {
        m_table( n, col ) = 0;
      }
    }
  }
}


void TriaxialDriver::outputResults()
{
  // TODO: improve file path output to grab command line -o directory
  //       for the moment, we just use the specified m_outputFile directly

  FILE * fp = fopen( m_outputFile.c_str(), "w" );

  /*
     string const outputDir = OutputBase::getOutputDirectory();
     string const outputPath = joinPath( outputDir, m_outputFile );
     FILE * fp = fopen( outputPath.c_str(), "w" );
   */

  fprintf( fp, "# column 1 = time\n" );
  fprintf( fp, "# column 2 = axial_strain\n" );
  fprintf( fp, "# column 3 = radial_strain_1\n" );
  fprintf( fp, "# column 4 = radial_strain_2\n" );
  fprintf( fp, "# column 5 = axial_stress\n" );
  fprintf( fp, "# column 6 = radial_stress_1\n" );
  fprintf( fp, "# column 7 = radial_stress_2\n" );
  fprintf( fp, "# column 8 = newton_iter\n" );
  fprintf( fp, "# column 9 = residual_norm\n" );

  for( integer n=0; n<=m_numSteps; ++n )
  {
    for( integer col=0; col<m_numColumns; ++col )
    {
      fprintf( fp, "%.4e ", m_table( n, col ) );
    }
    fprintf( fp, "\n" );
  }
  fclose( fp );
}


void TriaxialDriver::compareWithBaseline()
{
  // open baseline file

  std::ifstream file( m_baselineFile.c_str() );
  GEOSX_THROW_IF( !file.is_open(), "Can't seem to open the baseline file " << m_baselineFile, InputError );

  // discard file header

  string line;
  for( integer row=0; row < m_numColumns; ++row )
  {
    getline( file, line );
  }

  // read data block.  we assume the file size is consistent with m_table,
  // but check for a premature end-of-file. we then compare results value by value.
  // we ignore the newton iteration and residual columns, as those may be platform
  // specific.

  real64 value;
  real64 error;

  for( integer row=0; row < m_table.size( 0 ); ++row )
  {
    for( integer col=0; col < m_table.size( 1 ); ++col )
    {
      GEOSX_THROW_IF( file.eof(), "Baseline file appears shorter than internal results", std::runtime_error );
      file >> value;

      if( col < ITER ) // only compare "real" data columns
      {
        error = fabs( m_table[row][col]-value ) / ( fabs( value )+1 );
        GEOSX_THROW_IF( error > m_baselineTol, "Results do not match baseline at data row " << row+1
                                                                                            << " (row " << row+10 << " with header)"
                                                                                            << " and column " << col+1, std::runtime_error );
      }
    }
  }

  // check we actually reached the end of the baseline file

  file >> value;
  GEOSX_THROW_IF( !file.eof(), "Baseline file appears longer than internal results", std::runtime_error );

  // success

  if( getLogLevel() > 0 )
  {
    GEOSX_LOG_RANK_0( "  Comparison ........ Internal results consistent with baseline." );
  }

  file.close();
}


REGISTER_CATALOG_ENTRY( TaskBase,
                        TriaxialDriver,
                        string const &, dataRepository::Group * const )

} /* namespace geosx */
