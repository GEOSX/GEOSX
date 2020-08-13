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
 * @file DomainPartition.hpp
 */

#ifndef GEOSX_MANAGERS_DOMAINPARTITION_HPP_
#define GEOSX_MANAGERS_DOMAINPARTITION_HPP_

#include "dataRepository/Group.hpp"
#include "mesh/MeshBody.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"

namespace geosx
{
class SiloFile;
namespace dataRepository
{
namespace keys
{
/// @return PartitionManager string key
string const partitionManager( "partitionManager" );
}  // namespace keys
}  // namespace dataRepository

class ObjectManagerBase;
class PartitionBase;

/**
 * @brief Partition of the decomposed physical domain. It also manages the connexion information to its neighbors.
 *
 * Two types of neighors are being managed. One MPI cartesian communicator that DomainPartition shall build.
 * Or through `metis`, but it is not the responsibility of DomainParition to build the decomposition in that case.
 */
class DomainPartition : public dataRepository::Group
{
public:
  /**
   * @brief Constructor.
   * @param[in] name Name of this object manager
   * @param[in] parent Parent Group
   */
  DomainPartition( std::string const & name, Group * const parent );

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

  DomainPartition &
  operator=( DomainPartition const & ) = delete;

  DomainPartition &
  operator=( DomainPartition && ) = delete;
  /// @endcond
  ///@}

  /**
   * @copydoc dataRepository::Group::RegisterDataOnMeshRecursive( Group * const )
   */
  virtual void
  RegisterDataOnMeshRecursive( Group * const meshBodies ) override final;

  void
  InitializationOrder( string_array & order ) override final;

  /**
   * @brief Build all the sets of the DomainPartition.
   *
   * A domain contain sets of nodes or elements (that can be used to defined boundary conditions, etc.).
   * This member functions build those sets.
   */
  void
  GenerateSets();

  /**
   * @name MPI functionality
   */
  ///@{
  /**
   * @brief Constructs the communications between this DomainPartition and its neighbors.
   * @param use_nonblocking If true complete the communications of each phase in the order they are received.
   */
  void
  SetupCommunications( bool use_nonblocking );

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
  void
  AddNeighbors( const unsigned int idim, MPI_Comm & cartcomm, int * ncoords );
  ///@}

  /**
   * @brief struct to serve as a container for group strings and keys
   * @struct groupKeysStruct
   */
  struct groupKeysStruct
  {
    /// String key to the Group holding the MeshBodies
    static constexpr auto meshBodiesString = "MeshBodies";
    /// String key to the Group holding the ConstitutiveManager
    static constexpr auto constitutiveManagerString = "Constitutive";

    /// View key to the Group holding the MeshBodies
    dataRepository::GroupKey meshBodies = { meshBodiesString };
    /// View key to the Group holding the ConstitutiveManager
    dataRepository::GroupKey constitutiveManager = { constitutiveManagerString };
    /// View key to the Group holding the CommunicationManager
    dataRepository::GroupKey communicationManager = { "communicationManager" };
  }
  /// groupKey struct for the DomainPartition class
  groupKeys;

  /**
   * @brief Get the constitutive manager, const version.
   * @return Pointer to a const instance of a ConstitutiveManager.
   */
  constitutive::ConstitutiveManager const *
  getConstitutiveManager() const
  {
    return this->GetGroup< constitutive::ConstitutiveManager >(
      groupKeys.constitutiveManager );
  }

  /**
   * @brief Get the constitutive manager.
   * @return Pointer to an instance of a ConstitutiveManager.
   */
  constitutive::ConstitutiveManager *
  getConstitutiveManager()
  {
    return this->GetGroup< constitutive::ConstitutiveManager >(
      groupKeys.constitutiveManager );
  }

  /**
   * @brief @return Return a reference to const NumericalMethodsManager from ProblemManager
   */
  NumericalMethodsManager const &
  getNumericalMethodManager() const
  {
    return *(
      this->getParent()->GetGroup< NumericalMethodsManager >( "NumericalMethods" ) );
  }

  /**
   * @brief @return Return a reference to NumericalMethodsManager from ProblemManager
   */
  NumericalMethodsManager &
  getNumericalMethodManager()
  {
    return *(
      this->getParent()->GetGroup< NumericalMethodsManager >( "NumericalMethods" ) );
  }

  /**
   * @brief Get the mesh bodies, const version.
   * @return Pointer to a const instance of a Group that contains MeshBody instances.
   */
  Group const *
  getMeshBodies() const
  {
    return this->GetGroup( groupKeys.meshBodies );
  }

  /**
   * @brief Get the mesh bodies.
   * @return Pointer to a instance of a Group that contains MeshBody instances.
   */
  Group *
  getMeshBodies()
  {
    return this->GetGroup( groupKeys.meshBodies );
  }

  /**
   * @brief Get a MeshBody by name, const version.
   * @param meshName The name of the MeshBody.
   * @return Pointer to a const MeshBody instance matching @p meshName.
   */
  MeshBody const *
  getMeshBody( string const & meshName ) const
  {
    return this->GetGroup( groupKeys.meshBodies )->GetGroup< MeshBody >( meshName );
  }

  /**
   * @brief Get a MeshBody by name.
   * @param meshName The name of the MeshBody.
   * @return Pointer to a const MeshBody instance matching @p meshName.
   */
  MeshBody *
  getMeshBody( string const & meshName )
  {
    return this->GetGroup( groupKeys.meshBodies )->GetGroup< MeshBody >( meshName );
  }

  /**
   * @brief Get a MeshBody by index, const version.
   * @param index The index of the MeshBody.
   * @return Pointer to a const MeshBody instance at @p index position.
   */
  MeshBody const *
  getMeshBody( localIndex const index ) const
  {
    return this->GetGroup( groupKeys.meshBodies )->GetGroup< MeshBody >( index );
  }

  /**
   * @brief Get MeshBody by index.
   * @param index The index of the MeshBody.
   * @return Pointer to a MeshBody instance at @p index position.
   */
  MeshBody *
  getMeshBody( localIndex const index )
  {
    return this->GetGroup( groupKeys.meshBodies )->GetGroup< MeshBody >( index );
  }

  /**
   * @brief Get the metis neighbors indices.  @see DomainPartition#m_metisNeighborList
   * @return Container of global indices.
   */
  std::set< int > &
  getMetisNeighborList()
  {
    return m_metisNeighborList;
  }

  /**
   * @brief Get the metis neighbors indices, const version. @see DomainPartition#m_metisNeighborList
   * @return Container of global indices.
   */
  std::set< int > const &
  getMetisNeighborList() const
  {
    return m_metisNeighborList;
  }

  /**
   * @brief Get the neighbor communicators. @see DomainPartition#m_neighbors.
   * @return Container of communicators.
   */
  std::vector< NeighborCommunicator > &
  getNeighbors()
  {
    return m_neighbors;
  }

  /**
   * @brief Get the neighbor communicators, const version. @see DomainPartition#m_neighbors.
   * @return Container of communicators.
   */
  std::vector< NeighborCommunicator > const &
  getNeighbors() const
  {
    return m_neighbors;
  };

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

#endif /* GEOSX_MANAGERS_DOMAINPARTITION_HPP_ */
