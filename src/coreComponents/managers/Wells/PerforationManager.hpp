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

/*
 * @file PerforationManager.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATIONMANAGER_HPP
#define GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATIONMANAGER_HPP

#include "dataRepository/ManagedGroup.hpp"
#include "managers/ObjectManagerBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
static constexpr auto perforations = "Perforations";
}
}

class DomainPartition;
class Perforation;

class PerforationManager : public ObjectManagerBase
{
public:

  explicit PerforationManager( string const & name, dataRepository::ManagedGroup * const parent );
  ~PerforationManager() override;

  PerforationManager() = delete;
  PerforationManager( PerforationManager const &) = delete;
  PerforationManager( PerforationManager && ) = delete;

  virtual void FillDocumentationNode() override;
  virtual void CreateChild( string const & childKey, string const & childName ) override;

  virtual void FinalInitializationPreSubGroups(ManagedGroup * const problemManager) override;

  virtual const string getCatalogName() const override;

  localIndex numConnectionsGlobal() const { return integer_conversion<localIndex>(m_allPerfList.size()); }
  localIndex numConnectionsLocal()  const { return integer_conversion<localIndex>(size());         }

  Perforation const * getPerforation( localIndex iperf ) const;
  Perforation *       getPerforation( localIndex iperf );

  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {

    static constexpr auto connectionElementRegionString    = "connectionElementRegion";
    static constexpr auto connectionElementSubregionString = "connectionElementSubregion";
    static constexpr auto connectionElementIndexString     = "connectionElementIndex";
    static constexpr auto connectionPerforationIndexString = "connectionPerforationIndex";

    static constexpr auto gravityDepthString               = "gravityDepth";

    dataRepository::ViewKey connectionElementRegion    = { connectionElementRegionString    };
    dataRepository::ViewKey connectionElementSubregion = { connectionElementSubregionString };
    dataRepository::ViewKey connectionElementIndex     = { connectionElementIndexString     };
    dataRepository::ViewKey connectionPerforationIndex = { connectionPerforationIndexString };

    dataRepository::ViewKey gravityDepth               = { gravityDepthString     };

  } viewKeysPerfManager;

struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    static constexpr auto perforationString = "Perforation";

    dataRepository::GroupKey perforation = { perforationString };

  } groupKeysPerfManager;

private:

  void ConnectToCells( DomainPartition const * domain );
  void PrecomputeData( DomainPartition const * domain );

  array1d<localIndex> m_connectionElementRegion;
  array1d<localIndex> m_connectionElementSubregion;
  array1d<localIndex> m_connectionElementIndex;
  array1d<localIndex> m_connectionPerforationIndex;

  string_array m_allPerfList;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATIONMANAGER_HPP
