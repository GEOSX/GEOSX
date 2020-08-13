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
 * @file DofManager.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_DOFMANAGER_HPP_
#define GEOSX_LINEARALGEBRA_DOFMANAGER_HPP_

#include "LvArray/src/SparsityPattern.hpp"
#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/NodeManager.hpp"

namespace geosx
{
class DomainPartition;
class MeshLevel;
class ObjectManagerBase;
class FluxApproximationBase;

/**
 * @class DofManager
 * @brief The DoFManager is responsible for allocating global dofs, constructing
 * sparsity patterns, and generally simplifying the interaction between
 * PhysicsSolvers and linear algebra operations.
 */
class DofManager
{
public:
  /**
   * @brief Enumeration of geometric objects for support location. Note that this
   * enum is nearly identical to Connectivity, but we keep both for code readability
   * in function calls.
   */
  enum class Location
  {
    Elem,  //!< location is element (like pressure in finite volumes)
    Face,  //!< location is face (like flux in mixed finite elements)
    Edge,  //!< location is edge (like flux between fracture elements)
    Node   //!< location is node (like displacements in finite elements)
  };

  /**
   * @brief Enumeration of geometric objects for connectivity type. Note that this
   * enum is nearly identical to Location, but we keep both for code readability
   * in function calls.
   */
  enum class Connector
  {
    Elem,    //!< connectivity is element (like in finite elements)
    Face,    //!< connectivity is face (like in finite volumes TPFA)
    Edge,    //!< connectivity is edge (like fracture element connectors)
    Node,    //!< connectivity is node (like in finite volumes MPFA)
    None,    //!< there is no connectivity (self connected field, like a lumped mass matrix)
    Stencil  //!< connectivity is through a (set of) user-provided stencil(s)
  };

  /**
   * Limit on max number of fields
   */
  static localIndex constexpr MAX_FIELDS = 10;

  /**
   * Field description
   */
  struct FieldDescription
  {
    string name;                    //!< field name
    Location location;              //!< support location
    std::vector< string > regions;  //!< list of support region names
    localIndex numComponents;       //!< number of vector components
    string key;                     //!< string key for index array
    string docstring;               //!< documentation string
    localIndex numLocalDof;         //!< number of local rows
    globalIndex numGlobalDof;       //!< number of global rows
    globalIndex blockOffset;        //!< offset of this field's block in a block-wise ordered system
    globalIndex rankOffset;         //!< field's first DoF on current processor (within its block, ignoring other fields)
    globalIndex globalOffset;       //!< global offset of field's DOFs on current processor for multi-field problems
  };

  /**
   * Coupling description
   */
  struct CouplingDescription
  {
    Connector connector =
      Connector::None;                   //!< geometric object defining dof connections
    std::vector< std::string > regions;  //!< list of region names
    FluxApproximationBase const * stencils =
      nullptr;  //!< pointer to flux stencils for stencil based connections
  };

  /**
   * @brief Constructor.
   *
   * @param [in] name a unique name for this DoF manager
   */
  explicit DofManager( string name );

  /**
   * @brief Deleted copy constructor.
   */
  DofManager( DofManager const & ) = delete;

  /**
   * @brief Move constructor.
   */
  DofManager( DofManager && ) = default;

  /**
   * @brief Deleted copy assignment.
   * @return
   */
  DofManager &
  operator=( DofManager const & ) = delete;

  /**
   * @brief Deleted move assignment.
   * @return
   */
  DofManager &
  operator=( DofManager && ) = delete;

  /**
   * @brief Destructor.
   */
  ~DofManager() = default;

  /**
   * Remove all fields
   */
  void
  clear();

  /**
   * @brief Assign a mesh.
   *
   * @param [in] domain DomainPartition the input domain.
   * @param [in] meshLevelIndex Optional localIndex the mesh level.
   * @param [in] meshBodyIndex Optional localIndex the body level.
   */
  void
  setMesh( DomainPartition & domain,
           localIndex const meshLevelIndex = 0,
           localIndex const meshBodyIndex = 0 );

  /**
   * @brief Just an interface to allow only three parameters.
   *
   * @param [in] fieldName string the name of the field.
   * @param [in] location Location where it is defined.
   */
  void
  addField( string const & fieldName, Location const location );

  /**
   * @brief Just another interface to allow four parameters (no regions, default is everywhere).
   *
   * @param [in] fieldName string the name of the field.
   * @param [in] location Location where it is defined.
   * @param [in] components localIndex number of components (for vector fields).
   */
  void
  addField( string const & fieldName,
            Location const location,
            localIndex const components );

  /**
   * @brief Just another interface to allow four parameters (no components, default is 1).
   *
   * @param [in] fieldName string the name of the field.
   * @param [in] location Location where it is defined.
   * @param [in] regions names of regions where this field is defined.
   */
  void
  addField( string const & fieldName,
            Location const location,
            arrayView1d< string const > const & regions );

  /**
   * @brief The user can add a field with a support location, connectivity type, string key, number of scalar
   * components, and a list of element regions over which the field is active. If the region list is empty,
   * it is assumed the field exists on all regions.
   *
   * The connectivity type is used to infer the sparsity pattern that connects degrees of freedom.
   * If LC denotes a boolean connectivity graph between support locations L and connectors C, the desired sparsity
   * pattern will be computed as LC*CL. For example, for a TPFA discretization we have dofs located at cell centers,
   * and connected through adjacent faces. In this example, LC is the cell-to-face connectivity, and LC*CL is the
   * desired TPFA sparsity pattern. More generally,
   *
   * - Example 1 = ("displacement",NODE,ELEM,3) for a Q1 finite-element interpolation for elasticity
   * - Example 2 = ("pressure",ELEM,FACE,1) for a scalar TPFA-type approximation
   * - Example 3 = ("pressure",ELEM,NODE,1) for a scalar MPFA-type approximation
   * - Example 4 = ("mass",ELEM,NONE,1) for a diagonal-only sparsity pattern (no connectivitys)
   *
   * When the number of components is greater than one, we always assume they are tightly coupled to one another
   * and form a dense block. The sparsity pattern LC*CL is then interpreted as the super-node pattern, containing
   * dense sub-blocks.
   *
   * @param [in] fieldName string the name of the field.
   * @param [in] location Location where it is defined.
   * @param [in] components localIndex number of components (for vector fields).
   * @param [in] regions names of regions where this field is defined.
   */
  void
  addField( string const & fieldName,
            Location const location,
            localIndex const components,
            arrayView1d< string const > const & regions );

  /**
   * @brief Just an interface to allow only three parameters.
   *
   * @param [in] rowFieldName string the name of the row field.
   * @param [in] colFieldName string the name of the col field.
   * @param [in] connectivity Connectivity through what they are connected.
   */
  void
  addCoupling( string const & rowFieldName,
               string const & colFieldName,
               Connector const connectivity );

  /**
   * @brief Just another interface to allow four parameters (no symmetry, default is true).
   *
   * @param [in] rowFieldName string the name of the row field.
   * @param [in] colFieldName string the name of the col field.
   * @param [in] connectivity Connectivity through what they are connected.
   * @param [in] regions names of regions where this coupling is defined.
   */
  void
  addCoupling( string const & rowFieldName,
               string const & colFieldName,
               Connector const connectivity,
               arrayView1d< string const > const & regions );

  /**
   * @brief Just another interface to allow four parameters (no regions, default is everywhere).
   *
   * @param [in] rowFieldName string the name of the row field.
   * @param [in] colFieldName string the name of the col field.
   * @param [in] connectivity Connectivity through what they are connected.
   * @param [in] symmetric bool is it symmetric, i.e., both row-col and col-row?
   */
  void
  addCoupling( string const & rowFieldName,
               string const & colFieldName,
               Connector const connectivity,
               bool const symmetric );

  /**
   * @brief Add coupling between two fields.
   * The connectivity argument defines how the two fields couple. If the first field has support location A,
   * the second field has support location B, and the connecting object is C, the sparsity pattern will be
   * defined as (AC)(CB). The final argument indicates if the coupling is symmetric, in the sense that there
   * is a two-way coupling between the fields. Without this argument, a nonzero block will be added to the
   * system matrix for block AB, but block BA will remain zero (one-way coupling).
   *
   * - Example 1 = ("node_field","elem_field", ELEM, true) couples all dofs sharing a common element (two-way coupling)
   * - Example 2 = ("node_field_1","node_field_2", NODE, true) couples all dofs sharing a common node (two-way coupling)
   * - Example 3 = ("node_field_1","face_field", NODE, false) couples nodal dofs to adjacent faces (one-way coupling)
   *
   * @param [in] rowFieldName string the name of the row field.
   * @param [in] colFieldName string the name of the col field.
   * @param [in] connectivity Connectivity through what they are connected.
   * @param [in] regions names of regions where this coupling is defined.
   * @param [in] symmetric bool is it symmetric, i.e., both row-col and col-row?
   */
  void
  addCoupling( string const & rowFieldName,
               string const & colFieldName,
               Connector const connectivity,
               arrayView1d< string const > const & regions,
               bool const symmetric );

  /**
   * @brief Special interface for self-connectivity through a stencil.
   * @param [in] fieldName name of the field (this is only for diagonal blocks)
   * @param [in] stencils a pointer to FluxApproximation storing the stencils
   *
   * The field must be defined on element support. The set of regions is taken
   * automatically from the field definition.
   */
  void
  addCoupling( string const & fieldName,
               FluxApproximationBase const & stencils );

  /**
   * @brief Finish populating fields and apply appropriate dof renumbering
   *
   * This function must be called after all field and coupling information has been added.
   * It adjusts DoF index arrays to account for presence of other fields (in a global monolithic fashion).
   *
   * @note After DofManager has been closed, new fields and coupling cannot be added, until
   *       @ref clear or @ref setMesh is called.
   *
   * @note After reorderByRank() is called, the meaning of FieldDescription::globalOffset changes from
   *       "global offset of field's block in a global field-wise ordered (block) system" to
   *       "global offset of field's block on current processor in a rank-wise ordered system".
   *       This meaning is consistent with its use throughout. For example, this is the row/col
   *       global offset used to insert the field's sparsity block into a global coupled system.
   */
  void
  reorderByRank();

  /**
   * @brief Check if string key is already being used
   *
   * @param name field key to check
   * @return flag true if exists
   */
  bool
  fieldExists( string const & name ) const;

  /**
   * @brief Return the key used to record the field in the DofManager.
   *
   * @param [in] fieldName string the name of the field.
   * @return string indicating name of the field.
   */
  string const &
  getKey( string const & fieldName ) const;

  /**
   * @brief Return global number of dofs across all processors. If field argument is empty, return
   * monolithic size.
   *
   * @param [in] fieldName Optional string the name of the field.
   * @return     number of global dofs
   */
  globalIndex
  numGlobalDofs( string const & fieldName = "" ) const;

  /**
   * @brief Return local number of dofs on this processor. If field argument is empty, return
   * monolithic size.
   *
   * @param [in] fieldName Optional string the name of the field.
   * @return     number of local dofs
   */
  localIndex
  numLocalDofs( string const & fieldName = "" ) const;

  /**
   * @brief Return an array of local number of dofs on this processor
   * sorted by field registration order.
   *
   * @return     array of number of local dofs
   */
  array1d< localIndex >
  numLocalDofsPerField() const;

  /**
   * @brief Computes an array of size equal to sum of all field local number of dofs containing
   * unique integer labels associated to components stored in the field descriptions.
   *
   * @return array1d of localIndex labels
   */
  array1d< localIndex >
  getLocalDofComponentLabels() const;

  /**
   * @brief Return the sum of local dofs across all previous processors w.r.t. to the calling one for
   * the specified field.
   *
   * @param [in] fieldName Optional string the name of the field.
   * @return     the rank offset
   */
  globalIndex
  rankOffset( string const & fieldName = "" ) const;

  /**
   * @brief Get the number of components in a field. If field argument is empty, return
   * the total number of components across fields.
   *
   * @param [in] fieldName Optional string the name of the field.
   * @return     the number of dof components
   */
  localIndex
  numComponents( string const & fieldName = "" ) const;

  /**
   * @brief Return an array of number of components per field, sorted by field
   * registration order.
   *
   * @return     array of number of components
   */
  array1d< localIndex >
  numComponentsPerField() const;

  /**
   * @brief Get the local number of support points on this processor.
   * @param [in] fieldName the name of the field
   * @return number of local support points
   */
  localIndex
  numLocalSupport( string const & fieldName ) const;

  /**
   * @brief Get the local number of support points across all processors.
   * @param [in] fieldName name of the field
   * @return number of global support points
   */
  globalIndex
  numGlobalSupport( string const & fieldName ) const;

  /**
   * @brief Get the support location type of the field.
   * @param [in] fieldName name of the field
   * @return support location type
   */
  Location
  getLocation( string const & fieldName ) const;

  /**
   * @brief Get global offset of field's block on current processor in the system matrix.
   * @param [in] fieldName name of the field.
   * @return global offset of the field
   */
  globalIndex
  globalOffset( string const & fieldName ) const;

  /**
   * @brief Populate sparsity pattern of the entire system matrix.
   *
   * @param [out] matrix the target parallel matrix
   * @param [in]  closePattern whether to reorderByRank the matrix upon pattern assembly
   */
  template< typename MATRIX >
  void
  setSparsityPattern( MATRIX & matrix, bool closePattern = true ) const;

  /**
   * @brief Populate sparsity pattern for one block of the system matrix.
   *
   * @param [out] matrix the target parallel matrix
   * @param [in]  rowFieldName the name of the row field.
   * @param [in]  colFieldName the name of the col field.
   * @param [in]  closePattern whether to reorderByRank the matrix upon pattern assembly
   */
  template< typename MATRIX >
  void
  setSparsityPattern( MATRIX & matrix,
                      string const & rowFieldName,
                      string const & colFieldName,
                      bool closePattern = true ) const;

  /**
   * @brief Populate sparsity pattern of the entire system matrix.
   * @param [out] pattern the target sparsity pattern
   */
  void
  setSparsityPattern( SparsityPattern< globalIndex > & pattern ) const;

  /**
   * @brief Populate sparsity pattern for one block of the system matrix.
   * @param [out] pattern the target sparsity pattern
   * @param [in]  rowFieldName name of the row field
   * @param [in]  colFieldName name of the col field
   */
  void
  setSparsityPattern( SparsityPattern< globalIndex > & pattern,
                      string const & rowFieldName,
                      string const & colFieldName ) const;

  /**
   * @brief Copy values from LA vectors to simulation data arrays.
   *
   * @tparam FIELD_OP operation to perform (see FieldSpecificationOps.hpp)
   * @param vector source LA vector
   * @param srcFieldName name of the source field (as defined in DofManager)
   * @param dstFieldName name of the destination field (view wrapper key on the manager)
   * @param scalingFactor a factor to scale vector values by
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  template< typename VECTOR >
  void
  copyVectorToField( VECTOR const & vector,
                     string const & srcFieldName,
                     string const & dstFieldName,
                     real64 const scalingFactor,
                     localIndex const loCompIndex = 0,
                     localIndex const hiCompIndex = -1 ) const;

  /**
   * @brief Copy values from LA vectors to simulation data arrays.
   *
   * @param localVector source local vector
   * @param srcFieldName name of the source field (as defined in DofManager)
   * @param dstFieldName name of the destination field (view wrapper key on the manager)
   * @param scalingFactor a factor to scale vector values by
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  void
  copyVectorToField( arrayView1d< real64 const > const & localVector,
                     string const & srcFieldName,
                     string const & dstFieldName,
                     real64 const scalingFactor,
                     localIndex const loCompIndex = 0,
                     localIndex const hiCompIndex = -1 ) const;

  /**
   * @brief Add values from LA vectors to simulation data arrays.
   *
   * @tparam FIELD_OP operation to perform (see FieldSpecificationOps.hpp)
   * @param vector source LA vector
   * @param srcFieldName name of the source field (as defined in DofManager)
   * @param dstFieldName name of the destination field (view wrapper key on the manager)
   * @param scalingFactor a factor to scale vector values by
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  template< typename VECTOR >
  void
  addVectorToField( VECTOR const & vector,
                    string const & srcFieldName,
                    string const & dstFieldName,
                    real64 const scalingFactor,
                    localIndex const loCompIndex = 0,
                    localIndex const hiCompIndex = -1 ) const;

  /**
   * @brief Add values from LA vectors to simulation data arrays.
   *
   * @param localVector source local vector
   * @param srcFieldName name of the source field (as defined in DofManager)
   * @param dstFieldName name of the destination field (view wrapper key on the manager)
   * @param scalingFactor a factor to scale vector values by
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  void
  addVectorToField( arrayView1d< real64 const > const & localVector,
                    string const & srcFieldName,
                    string const & dstFieldName,
                    real64 const scalingFactor,
                    localIndex const loCompIndex = 0,
                    localIndex const hiCompIndex = -1 ) const;

  /**
   * @brief Copy values from nodes to DOFs.
   *
   * @tparam FIELD_OP operation to perform (see FieldSpecificationOps.hpp)
   * @param vector target LA vector
   * @param srcFieldName name of the source field (view wrapper key on the manager)
   * @param dstFieldName name of the destination field (as defined in DofManager)
   * @param scalingFactor a factor to scale vector values by
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  template< typename VECTOR >
  void
  copyFieldToVector( VECTOR & vector,
                     string const & srcFieldName,
                     string const & dstFieldName,
                     real64 const scalingFactor,
                     localIndex const loCompIndex = 0,
                     localIndex const hiCompIndex = -1 ) const;

  /**
   * @brief Copy values from simulation data arrays to vectors.
   *
   * @tparam FIELD_OP operation to perform (see FieldSpecificationOps.hpp)
   * @param localVector target LA vector
   * @param srcFieldName name of the source field (view wrapper key on the manager)
   * @param dstFieldName name of the destination field (as defined in DofManager)
   * @param scalingFactor a factor to scale vector values by
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  void
  copyFieldToVector( arrayView1d< real64 > const & localVector,
                     string const & srcFieldName,
                     string const & dstFieldName,
                     real64 const scalingFactor,
                     localIndex const loCompIndex = 0,
                     localIndex const hiCompIndex = -1 ) const;

  /**
   * @brief Add values from a simulation data array to a DOF vector.
   *
   * @tparam FIELD_OP operation to perform (see FieldSpecificationOps.hpp)
   * @param vector target LA vector
   * @param srcFieldName name of the source field (view wrapper key on the manager)
   * @param dstFieldName name of the destination field (as defined in DofManager)
   * @param scalingFactor a factor to scale vector values by
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  template< typename VECTOR >
  void
  addFieldToVector( VECTOR & vector,
                    string const & srcFieldName,
                    string const & dstFieldName,
                    real64 const scalingFactor,
                    localIndex const loCompIndex = 0,
                    localIndex const hiCompIndex = -1 ) const;

  /**
   * @brief Add values from a simulation data array to a DOF vector.
   *
   * @tparam FIELD_OP operation to perform (see FieldSpecificationOps.hpp)
   * @param localVector target vector
   * @param srcFieldName name of the source field (view wrapper key on the manager)
   * @param dstFieldName name of the destination field (as defined in DofManager)
   * @param scalingFactor a factor to scale vector values by
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  void
  addFieldToVector( arrayView1d< real64 > const & localVector,
                    string const & srcFieldName,
                    string const & dstFieldName,
                    real64 const scalingFactor,
                    localIndex const loCompIndex = 0,
                    localIndex const hiCompIndex = -1 ) const;

  /**
   * @brief Describes a selection of components from a DoF field.
   *
   * A half-open range [@p loComp, @p hiComp) is selected.
   */
  struct SubComponent
  {
    string fieldName;   //!< Name of the DOF field in DofManager
    localIndex loComp;  //!< Low component index (included in selection)
    localIndex hiComp;  //!< High component index (excluded from selection)
  };

  /**
   * @brief Create a dof selection by filtering out excluded components
   * @param excluded a list of dof components to exclude
   * @return a vector of remaining dof components
   *
   * @note Removed components must not have repeats, and each entry must either have
   *       loComp = 0 or hiComp = numComponents(fieldName) (or both). In other words,
   *       filtered out components must not leave "holes" in DOFs.
   */
  std::vector< SubComponent >
  filterDofs(
    std::vector< SubComponent > const & excluded ) const;

  /**
   * @brief Create a matrix that restricts vectors and matrices to a subset of DOFs
   * @tparam MATRIX type of matrix used for restrictor
   * @param selection a list of fields to select; each entry is a struct containing
   *                  the name of the field and low and high selected component indices
   * @param comm the MPI communicator to use in the operator
   * @param transpose if @p true, the transpose (prolongation) operator will be created
   * @param restrictor resulting operator
   *
   * @note Can only be called after reorderByRank(), since global DOF indexing is required
   *       for the restrictor to make sense.
   */
  template< typename MATRIX >
  void
  makeRestrictor( std::vector< SubComponent > const & selection,
                  MPI_Comm const & comm,
                  bool transpose,
                  MATRIX & restrictor ) const;

  /**
   * @brief Print the summary of declared fields and coupling.
   *
   * @param os output stream
   */
  void
  printFieldInfo( std::ostream & os = std::cout ) const;

private:
  /**
   * @brief Initialize data structure for connectivity and sparsity pattern
   */
  void
  initializeDataStructure();

  /**
   * @brief Get field index from string key
   */
  localIndex
  getFieldIndex( string const & name ) const;

  /**
   * @brief Create index array for the field
   */
  void
  createIndexArray( FieldDescription & field );

  /**
   * @brief Remove an index array for the field
   */
  void
  removeIndexArray( FieldDescription const & field );

  /**
   * @brief Populate the sparsity pattern for a coupling block between given fields.
   * @param pattern the sparsity to be filled
   * @param rowFieldIndex index of row field (must be non-negative)
   * @param colFieldIndex index of col field (must be non-negative)
   *
   * This private function is used as a building block by higher-level SetSparsityPattern()
   */
  template< typename MATRIX >
  void
  setSparsityPatternOneBlock( MATRIX & pattern,
                              localIndex const rowFieldIndex,
                              localIndex const colFieldIndex ) const;

  template< typename MATRIX >
  void
  setSparsityPatternFromStencil( MATRIX & pattern,
                                 localIndex const fieldIndex ) const;

  /**
   * @brief Calculate or estimate the number of nonzero entries in each local row
   * @param rowLengths array of row lengths (values are be incremented, not overwritten)
   * @param rowFieldIndex index of row field (must be non-negative)
   * @param colFieldIndex index of col field (must be non-negative)
   */
  void
  countRowLengthsOneBlock( arrayView1d< localIndex > const & rowLengths,
                           localIndex const rowFieldIndex,
                           localIndex const colFieldIndex ) const;

  void
  countRowLengthsFromStencil( arrayView1d< localIndex > const & rowLengths,
                              localIndex const fieldIndex ) const;

  /**
   * @brief Populate the sparsity pattern for a coupling block between given fields.
   * @param pattern the sparsity to be filled
   * @param rowFieldIndex index of row field (must be non-negative)
   * @param colFieldIndex index of col field (must be non-negative)
   *
   * This private function is used as a building block by higher-level SetSparsityPattern()
   */
  void
  setSparsityPatternOneBlock( SparsityPattern< globalIndex > & pattern,
                              localIndex const rowFieldIndex,
                              localIndex const colFieldIndex ) const;

  void
  setSparsityPatternFromStencil( SparsityPattern< globalIndex > & pattern,
                                 localIndex const fieldIndex ) const;

  template< int DIMS_PER_DOF >
  void
  setFiniteElementSparsityPattern( SparsityPattern< globalIndex > & pattern,
                                   localIndex const fieldIndex ) const;

  /**
   * @brief Generic implementation for @ref copyVectorToField and @ref addVectorToField
   * @tparam FIELD_OP operation to perform (see FieldSpecificationOps.hpp)
   * @tparam POLICY execution policy for the kernel
   * @param vector source LA vector
   * @param srcFieldName name of the source field (as defined in DofManager)
   * @param scalingFactor a factor to scale vector values by
   * @param manager mesh object manager that contains the target field (subregion for elements)
   * @param dstFieldName name of the destination field (view wrapper key on the manager)
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  template< typename FIELD_OP, typename POLICY, typename LOCAL_VECTOR >
  void
  vectorToField( LOCAL_VECTOR const localVector,
                 string const & srcFieldName,
                 string const & dstFieldName,
                 real64 const scalingFactor,
                 localIndex const loCompIndex,
                 localIndex const hiCompIndex ) const;

  /**
   * @brief Generic implementation for @ref copyFieldToVector and @ref addFieldToVector
   * @tparam FIELD_OP operation to perform (see FieldSpecificationOps.hpp)
   * @tparam POLICY execution policy for the kernel
   * @param manager mesh object manager that contains the target field (subregion for elements)
   * @param srcFieldName name of the source field (view wrapper key on the manager)
   * @param scalingFactor a factor to scale vector values by
   * @param vector ponter to target vector local data (host or device)
   * @param dstFieldName name of the destination field (as defined in DofManager)
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  template< typename FIELD_OP, typename POLICY, typename LOCAL_VECTOR >
  void
  fieldToVector( LOCAL_VECTOR localVector,
                 string const & srcFieldName,
                 string const & dstFieldName,
                 real64 const scalingFactor,
                 localIndex const loCompIndex,
                 localIndex const hiCompIndex ) const;

  /// Name of the manager (unique, for unique identification of index array keys)
  string m_name;

  /// Pointer to domain manager
  DomainPartition * m_domain = nullptr;

  /// Pointer to corresponding MeshLevel
  MeshLevel * m_mesh = nullptr;

  /// Array of field descriptions
  std::vector< FieldDescription > m_fields;

  /// Table of connector types within and between fields
  std::vector< std::vector< CouplingDescription > > m_coupling;

  /// Flag indicating that DOFs have been reordered rank-wise.
  bool m_reordered;
};

} /* namespace geosx */

#endif /*GEOSX_LINEARALGEBRA_DOFMANAGER_HPP_*/
