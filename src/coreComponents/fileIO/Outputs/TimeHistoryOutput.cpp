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

#include "TimeHistoryOutput.hpp"

namespace geosx
{
TimeHistoryOutput::TimeHistoryOutput( string const & name,
                                      Group * const parent ):
  OutputBase( name, parent ),
  m_collectorPaths( ),
  m_format( ),
  m_filename( ),
  m_recordCount( 0 ),
  m_io( )
{
  registerWrapper( viewKeys::timeHistoryOutputTargetString(), &m_collectorPaths ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "A list of collectors from which to collect and output time history information." );

  registerWrapper( viewKeys::timeHistoryOutputFilenameString(), &m_filename ).
    setApplyDefaultValue( "TimeHistory" ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The filename to which to write time history output." );

  registerWrapper( viewKeys::timeHistoryOutputFormatString(), &m_format ).
    setApplyDefaultValue( "hdf" ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The output file format for time history output." );

  registerWrapper( viewKeys::timeHistoryRestartString(), &m_recordCount ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "The current history record to be written, on restart from an earlier time allows use to remove invalid future history." );

}

void TimeHistoryOutput::initCollectorParallel( DomainPartition & domain, HistoryCollection & collector )
{
  bool const freshInit = ( m_recordCount == 0 );

  string const outputDirectory = getOutputDirectory();
  string const outputFile = joinPath( outputDirectory, m_filename );

  localIndex ioCount = 0;
  for( localIndex ii = 0; ii < collector.getCollectionCount( ); ++ii )
  {
    HistoryMetadata metadata = collector.getMetadata( domain, ii );
    m_io.emplace_back( std::make_unique< HDFHistIO >( outputFile, metadata, m_recordCount ) );
    collector.registerBufferCall( ii, [this, ii, ioCount, &domain, &collector]()
    {
      collector.updateSetsIndices( domain );
      HistoryMetadata md = collector.getMetadata( domain, ii );
      m_io[ioCount]->updateCollectingCount( md.getDims( )[0] );
      return m_io[ioCount]->getBufferHead( );
    } );
    m_io.back()->init( !freshInit );
    ++ioCount;
  }

  localIndex metaCollectorCount = collector.getNumMetaCollectors( );
  for( localIndex metaIdx = 0; metaIdx < metaCollectorCount; ++metaIdx )
  {
    HistoryCollection & metaCollector = collector.getMetaCollector( metaIdx );
    for( localIndex ii = 0; ii < metaCollector.getCollectionCount( ); ++ii )
    {
      HistoryMetadata metaMetadata = metaCollector.getMetadata( domain, ii );
      metaMetadata.setName( collector.getTargetName() + " " + metaMetadata.getName( ) );
      m_io.emplace_back( std::make_unique< HDFHistIO >( outputFile, metaMetadata, m_recordCount ) );
      metaCollector.registerBufferCall( ii, [this, ii, ioCount, &domain, &metaCollector] ()
      {
        metaCollector.updateSetsIndices( domain );
        HistoryMetadata md = metaCollector.getMetadata( domain, ii );
        m_io[ ioCount ]->updateCollectingCount( md.getDims()[0] );
        return m_io[ ioCount ]->getBufferHead( );
      } );
      m_io.back()->init( !freshInit );
      ++ioCount;
    }

  }

  // do the time output last so its at the end of the m_io list, since writes are parallel we need
  //  the rest of the collectors to share position in the list across the world comm

  // rank == 0 does time output for the collector
  if( MpiWrapper::commRank( ) == 0 )
  {
    HistoryMetadata timeMetadata = collector.getTimeMetadata( );
    m_io.emplace_back( std::make_unique< HDFHistIO >( outputFile, timeMetadata, m_recordCount, 1, 2, MPI_COMM_SELF ) );
    collector.registerTimeBufferCall( [this]() { return m_io.back()->getBufferHead( ); } );
    m_io.back()->init( !freshInit );
  }

  MpiWrapper::barrier( MPI_COMM_GEOSX );
}

void TimeHistoryOutput::initializePostSubGroups()
{
  {
    // check whether to truncate or append to the file up front so we don't have to bother during later accesses
    string const outputDirectory = getOutputDirectory();
    if( MpiWrapper::commRank( MPI_COMM_GEOSX ) == 0 )
    {
      makeDirsForPath( outputDirectory );
    }
    MpiWrapper::barrier( MPI_COMM_GEOSX );
    string const outputFile = joinPath( outputDirectory, m_filename );
    HDFFile( outputFile, (m_recordCount == 0), true, MPI_COMM_GEOSX );
  }

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  for( auto collectorPath : m_collectorPaths )
  {
    HistoryCollection & collector = this->getGroupByPath< HistoryCollection >( collectorPath );
    collector.initializePostSubGroups();
    initCollectorParallel( domain, collector );
  }
}

bool TimeHistoryOutput::execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                 real64 const GEOSX_UNUSED_PARAM( dt ),
                                 integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                 integer const GEOSX_UNUSED_PARAM( eventCounter ),
                                 real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                                 DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_MARK_FUNCTION;
  localIndex newBuffered = m_io.front()->getBufferedCount( );
  for( auto & th_io : m_io )
  {
    th_io->write( );
  }
  m_recordCount += newBuffered;
  return false;
}

void TimeHistoryOutput::cleanup( real64 const time_n,
                                 integer const cycleNumber,
                                 integer const eventCounter,
                                 real64 const eventProgress,
                                 DomainPartition & domain )
{
  execute( time_n, 0.0, cycleNumber, eventCounter, eventProgress, domain );
  // remove any unused trailing space reserved to write additional histories
  for( auto & th_io : m_io )
  {
    th_io->compressInFile();
  }
}

REGISTER_CATALOG_ENTRY( OutputBase, TimeHistoryOutput, string const &, Group * const )
}
