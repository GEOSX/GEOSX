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
 * @file TimeHistoryCollector.hpp
 */

#ifndef GEOSX_SetIndexCollection_HPP_
#define GEOSX_SetIndexCollection_HPP_

#include "TimeHistoryCollection.hpp"

namespace geosx
{

/**
 * @class SetIndexCollection
 *
 * A task class for serializing set collection-index information (indices into HistoryOutput related to a Set,
 *   not the simulation indices of the Set).
 */
class SetIndexCollection : public HistoryCollection
{
public:
/**
 * @brief Constructor
 * @param set_index_offset An offset constituting all 'prior' rank-sets + local-set index counts to correctly provide history output set
 * indices.
 * @copydoc geosx::dataRepository::Group::Group(std::string const & name, Group * const parent)
 */
  SetIndexCollection( string const & object_path, string const & set_name, globalIndex set_index_offset );

/// @copydoc HistoryCollection::GetMetadata
  virtual HistoryMetadata GetMetadata( Group * problem_group ) override;

/// @copydoc HistoryCollection::Collect
  virtual void Collect( Group * domain_group,
                        real64 const GEOSX_UNUSED_PARAM( time_n ),
                        real64 const GEOSX_UNUSED_PARAM( dt ),
                        buffer_unit_type * & buffer ) override;

protected:
  string m_object_path;
  string m_set_name;
  globalIndex m_set_index_offset;
};
}
#endif
