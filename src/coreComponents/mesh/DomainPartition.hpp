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
 * @file DomainPartition.hpp
 */

#ifndef GEOSX_MESH_DOMAINPARTITION_HPP_
#define GEOSX_MESH_DOMAINPARTITION_HPP_

#include "common/MpiWrapper.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "dataRepository/Group.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"

namespace geosx
{

class SiloFile;
namespace dataRepository
{
namespace keys
{
/// @return PartitionManager string key
string const partitionManager( "partitionManager" );
}
}

class ObjectManagerBase;
class PartitionBase;

/**
 * @brief Partition of the decomposed physical domain. It also manages the connexion information to its neighbors.
 *
 * Two types of neighors are being managed. One MPI cartesian communicator that DomainPartition shall build.
 * Or through `metis`, but it is not the responsibility of DomainPartition to build the decomposition in that case.
 */
class DomainPartition : public dataRepository::Group
{
public:
  /**
   * @brief Constructor.
   * @param[in] name Name of this object manager
   * @param[in] parent Parent Group
   */
  DomainPartition( string const & name,
                   Group * const parent );

  /**
   * @brief Destructor.
   */
  ~DomainPartition() override;

  /**
   * @name Default, copy and assignment constructors are deleted.
   */
  ///@{
  /// @cond DO_NOT_DOCUMENT
  DomainPartition() = delete;

  DomainPartition( DomainPartition const & ) = delete;

  DomainPartition( DomainPartition && ) = delete;

  DomainPartition & operator=( DomainPartition const & ) = delete;

  DomainPartition & operator=( DomainPartition && ) = delete;
  /// @endcond
  ///@}

  void initializationOrder( string_array & order ) override final;

  /**
   * @brief Build all the sets of the DomainPartition.
   *
   * A domain contain sets of nodes or elements (that can be used to defined boundary conditions, etc.).
   * This member functions build those sets.
   */
  void generateSets();

  /**
   * @name MPI functionality
   */
  ///@{
  /**
   * @brief Constructs the communications between this DomainPartition and its neighbors.
   * @param use_nonblocking If true complete the communications of each phase in the order they are received.
   */
  void setupCommunications( bool use_nonblocking );

  /**
   * @brief Recursively builds neighbors if an MPI cartesian topology is used (i.e. not metis).
   * @param idim Dimension index in the cartesian.
   * @param cartcomm Communicator with cartesian structure.
   * @param ncoords Cartesian coordinates of a process (assumed to be of length 3).
   *
   * This recursive function builds the neighbors recursively by increasing
   * the dimension index of the current DomainPartition until all the dimensions (3) are done.
   * The relevant values for initiating this functions are therefore @p ibim = 0
   * and a non-initialized vector @p ncoords of length 3.
   *
   * This functions should have been implemented `private`
   * and an additional functions to initiate the recursion could have been implemented.
   */
  void addNeighbors( const unsigned int idim,
                     MPI_Comm & cartcomm,
                     int * ncoords );
  ///@}


  /**
   * @brief struct to serve as a container for group strings and keys
   * @struct groupKeysStruct
   */
  struct groupKeysStruct
  {
    /// @return String key to the Group holding the MeshBodies
    static constexpr char const * meshBodiesString() { return "MeshBodies"; }
    /// @return String key to the Group holding the ConstitutiveManager
    static constexpr char const * constitutiveManagerString() { return "Constitutive"; }

    /// View key to the Group holding the MeshBodies
    dataRepository::GroupKey meshBodies = { meshBodiesString() };
    /// View key to the Group holding the ConstitutiveManager
    dataRepository::GroupKey constitutiveManager = { constitutiveManagerString() };
    /// View key to the Group holding the CommunicationManager
    dataRepository::GroupKey communicationManager = { "communicationManager" };
  }
  /// groupKey struct for the DomainPartition class
  groupKeys;

  /**
   * @brief Get the constitutive manager, const version.
   * @return Pointer to a const instance of a ConstitutiveManager.
   */
  constitutive::ConstitutiveManager const & getConstitutiveManager() const
  { return this->getGroup< constitutive::ConstitutiveManager >( groupKeys.constitutiveManager ); }

  /**
   * @brief Get the constitutive manager.
   * @return Pointer to an instance of a ConstitutiveManager.
   */
  constitutive::ConstitutiveManager & getConstitutiveManager()
  { return this->getGroup< constitutive::ConstitutiveManager >( groupKeys.constitutiveManager ); }

  /**
   * @brief @return Return a reference to const NumericalMethodsManager from ProblemManager
   */
  NumericalMethodsManager const & getNumericalMethodManager() const
  { return this->getParent().getGroup< NumericalMethodsManager >( "NumericalMethods" ); }

  /**
   * @brief @return Return a reference to NumericalMethodsManager from ProblemManager
   */
  NumericalMethodsManager & getNumericalMethodManager()
  { return this->getParent().getGroup< NumericalMethodsManager >( "NumericalMethods" ); }

  /**
   * @brief Get the mesh bodies, const version.
   * @return Reference to a const instance of a Group that contains MeshBody instances.
   */
  Group const & getMeshBodies() const
  { return this->getGroup( groupKeys.meshBodies ); }

  /**
   * @brief Get the mesh bodies.
   * @return Reference to a instance of a Group that contains MeshBody instances.
   */
  Group & getMeshBodies()
  { return this->getGroup( groupKeys.meshBodies ); }

  /**
   * @brief Get a MeshBody by name, const version.
   * @tparam KEY_TYPE The type of the key used to look up the MeshBody.
   * @param key The key to the MeshBody.
   * @return Reference to a const MeshBody instance matching @p key.
   */
  template< typename KEY_TYPE >
  MeshBody const & getMeshBody( KEY_TYPE const & key ) const
  { return getMeshBodies().getGroup< MeshBody >( key ); }

  /**
   * @brief Get a MeshBody by name.
   * @tparam KEY_TYPE The type of the key used to look up the MeshBody.
   * @param key The key to the MeshBody.
   * @return Reference to a const MeshBody instance matching @p key.
   */
  template< typename KEY_TYPE >
  MeshBody & getMeshBody( KEY_TYPE const & key )
  { return getMeshBodies().getGroup< MeshBody >( key ); }

  /**
   * @brief Get the metis neighbors indices.  @see DomainPartition#m_metisNeighborList
   * @return Container of global indices.
   */
  std::set< int > & getMetisNeighborList()
  { return m_metisNeighborList; }

  /**
   * @brief Get the metis neighbors indices, const version. @see DomainPartition#m_metisNeighborList
   * @return Container of global indices.
   */
  std::set< int > const & getMetisNeighborList() const
  { return m_metisNeighborList; }

  /**
   * @brief Get the neighbor communicators. @see DomainPartition#m_neighbors.
   * @return Container of communicators.
   */
  std::vector< NeighborCommunicator > & getNeighbors()
  { return m_neighbors; }

  /**
   * @brief Get the neighbor communicators, const version. @see DomainPartition#m_neighbors.
   * @return Container of communicators.
   */
  std::vector< NeighborCommunicator > const & getNeighbors() const
  { return m_neighbors; };

  /**
   * @brief Get a list of neighbor ranks.
   * @return Container of neighbor ranks.
   */
  std::vector< int > getNeighborRanks() const
  {
    std::vector< int > ranks;
    ranks.reserve( m_neighbors.size() );
    for( NeighborCommunicator const & neighbor : m_neighbors )
    {
      ranks.push_back( neighbor.neighborRank() );
    }
    return ranks;
  }

private:

  /**
   * @brief Contains the global indices of the metis neighbors in case `metis` is used. Empty otherwise.
   */
  std::set< int > m_metisNeighborList;
  /**
   * @brief Contains all the communicators from this DomainPartition to its neighbors.
   */
  std::vector< NeighborCommunicator > m_neighbors;
};

} /* namespace geosx */

#endif /* GEOSX_MESH_DOMAINPARTITION_HPP_ */
