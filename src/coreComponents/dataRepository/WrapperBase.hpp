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

/** @file */

#ifndef GEOSX_DATAREPOSITORY_WRAPPERBASE_HPP_
#define GEOSX_DATAREPOSITORY_WRAPPERBASE_HPP_

#include "common/DataTypes.hpp"
#include "InputFlags.hpp"
#include "xmlWrapper.hpp"
#include "RestartFlags.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include <string>
#include <memory>
#include <set>

namespace conduit
{
class Node;
}

namespace geosx
{
namespace dataRepository
{
class Group;

/**
 * @class WrapperBase
 * @brief Base class for all wrappers containing common operations
 */
class WrapperBase
{
public:
  /**
   * @name Constructor, desctructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param[in] name name of the object
   * @param[in] parent pointer to Group that holds this WrapperBase
   */
  explicit WrapperBase( string const & name, Group * const parent );

  /// @cond DO_NOT_DOCUMENT
  WrapperBase() = delete;
  WrapperBase( WrapperBase const & ) = delete;
  WrapperBase( WrapperBase && ) = delete;
  WrapperBase &
  operator=( WrapperBase const & ) = delete;
  WrapperBase &
  operator=( WrapperBase && ) = delete;
  /// @endcond

  /**
   * @brief Default destructor.
   */
  virtual ~WrapperBase();

  ///@}

  /**
   * @name Methods that delegate to the wrapped type
   *
   * These functions will call the corresponding method on the wrapped
   * object, if such method is declared in wrapped type.
   */
  ///@{

  /**
   * @brief Calls T::size()
   * @return result of T::size()
   */
  virtual localIndex
  size() const = 0;

  /**
   * @brief @return a const void pointer to the data.
   */
  virtual void const *
  voidPointer() const = 0;

  /**
   * @brief @return the number of bytes in an object of unit size.
   *        Ie elementByteSize() * size() gives the size of memory pointed to by voidPointer().
   */
  virtual localIndex
  elementByteSize() const = 0;

  /**
   * @brief Calls T::resize( num_dims, dims )
   * @param[in] num_dims number of dimensions in T
   * @param[in] dims pointer to the new dims
   */
  virtual void
  resize( int num_dims, localIndex const * const dims ) = 0;

  /**
   * @brief Calls T::reserve( newCapacity ) if it exists, otherwise a no-op.
   * @param[in] newCapacity the new capacity of the T.
   */
  virtual void
  reserve( localIndex const newCapacity ) = 0;

  /**
   * @brief @return T::capacity() if it exists, other wise calls size().
   */
  virtual localIndex
  capacity() const = 0;

  /**
   * @brief Calls T::resize(newsize) if it exists.
   * @param[in] newsize parameter to pass to T::resize(newsize)
   */
  virtual void
  resize( localIndex newsize ) = 0;

  /**
   * @brief Calls resize(newsize) where newsize is taken from the parent Group.
   */
  void
  resize();

  /**
   * @brief Calls T::copy(sourceIndex, destIndex)
   * @param[in] sourceIndex index to copy from
   * @param[in] destIndex index to copy to
   */
  virtual void
  copy( localIndex const sourceIndex, localIndex const destIndex ) = 0;

  /**
   * @brief Calls T::move(space, touch)
   * @param[in] space A CHAI execution space to move the data into
   * @param[in] touch whether to register a touch in target space
   */
  virtual void
  move( LvArray::MemorySpace const space, bool const touch ) const = 0;

  ///@}

  /**
   * @brief Return true iff this wrapper has a valid default value.
   * @return True iff this wrapper has a valid default value.
   */
  virtual bool
  hasDefaultValue() const = 0;

  /**
   * @brief Return a string representing the default value.
   * @return A string representing the default value.
   */
  virtual std::string
  getDefaultValueString() const = 0;

  /**
   * @brief Initialize the wrapper from the input xml node.
   * @param targetNode the xml node to initialize from.
   * @return True iff the wrapper initialized itself from the file.
   */
  virtual bool
  processInputFile( xmlWrapper::xmlNode const & targetNode ) = 0;

  /**
   * @brief Push the data in the wrapper into a Conduit blueprint field.
   * @param fields The Conduit Node containg the blueprint fields.
   * @param name The name of the field.
   * @param topology The topology associated with the field.
   * @param componentNames The name of the components, if not specified they are auto generated.
   * @note This wrapper must hold an LvArray::Array.
   */
  virtual void
  addBlueprintField(
    conduit::Node & fields,
    std::string const & name,
    std::string const & topology,
    std::vector< std::string > const & componentNames = {} ) const = 0;

  /**
   * @brief Push the data in the wrapper into a Conduit Blueprint mcarray.
   * @param node The Conduit Node to put the data into.
   * @param componentNames The names of the components, if not specified they are auto generated.
   * @note This wrapper must hold an LvArray::Array.
   */
  virtual void
  populateMCArray(
    conduit::Node & node,
    std::vector< std::string > const & componentNames = {} ) const = 0;

  /**
   * @brief Create a new Wrapper with values averaged over the second dimension.
   * @param name The name to give the new wrapper.
   * @param group The group to hang the new Wrapper from.
   * @return The newly created wrapper.
   * @note This Wrapper must hold an LvArray::Array of dimension 2 or greater.
   * @note The new Wrapper is not registered with @p group.
   */
  virtual std::unique_ptr< WrapperBase >
  averageOverSecondDim( std::string const & name,
                        Group & group ) const = 0;

  /**
   * @name Restart output methods
   */
  ///@{

  /**
   * @brief Register the wrapper's data for writing with Conduit.
   */
  virtual void
  registerToWrite() const = 0;

  /**
   * @brief Write the wrapped data into Conduit.
   */
  virtual void
  finishWriting() const = 0;

  /**
   * @brief Read the wrapped data from Conduit.
   * @return True iff the Wrapper read in data.
   */
  virtual bool
  loadFromConduit() = 0;

  ///@}

  /**
   * @name Methods for buffer packing/unpacking
   *
   * This group of functions is used to pack/unpack wrapped object to/from binary buffers
   */
  ///@{

  /**
   * @brief Check whether wrapped type is can be packed into a buffer on host or device.
   * @param[in] on_device    determine whether the wrapper is packable on host vs device
   * @return @p true if @p T is packable, @p false otherwise
   */
  virtual bool
  isPackable( bool on_device = false ) const = 0;

  /**
   * @brief Pack the entire wrapped object into a buffer.
   * @param[in,out] buffer the binary buffer pointer, advanced upon completion
   * @param[in] on_device    whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @return               the number of @p buffer_unit_type units packed
   */
  virtual localIndex
  Pack( buffer_unit_type *& buffer,
        bool on_device = false ) const = 0;

  /**
   * @brief For indexable types, pack selected indices of wrapped object into a buffer.
   * @param[in,out] buffer the binary buffer pointer, advanced upon completion
   * @param[in] packList   the list of indices to pack
   * @param[in] on_device    whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @return               the number of @p buffer_unit_type units packed
   */
  virtual localIndex
  PackByIndex( buffer_unit_type *& buffer,
               arrayView1d< localIndex const > const & packList,
               bool on_device = false ) const = 0;

  /**
   * @brief Get the buffer size needed to pack the entire wrapped object.
   * @param[in] on_device    whether to use device-based packing functions
   *                         this matters as the size on device differs from the size on host
   *                         as we pack less metadata on device
   * @return the number of @p buffer_unit_type units needed to pack
   */
  virtual localIndex
  PackSize( bool on_device = false ) const = 0;

  /**
   * @brief Get the buffer size needed to pack the selected indices wrapped object.
   * @param[in] packList the list of indices to pack
   * @param[in] on_device    whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @return             the number of @p buffer_unit_type units needed to pack
   */
  virtual localIndex
  PackByIndexSize( arrayView1d< localIndex const > const & packList,
                   bool on_device = false ) const = 0;

  /**
   * @brief Unpack the entire wrapped object from a buffer.
   * @param[in,out] buffer the binary buffer pointer, advanced upon completion
   * @param[in] on_device    whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @return               the number of @p buffer_unit_type units unpacked
   */
  virtual localIndex
  Unpack( buffer_unit_type const *& buffer,
          bool on_device = false ) = 0;

  /**
   * @brief For indexable types, unpack selected indices of wrapped object from a buffer.
   * @param[in,out] buffer    the binary buffer pointer, advanced upon completion
   * @param[in] unpackIndices the list of indices to pack
   * @param[in] on_device    whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @return                  the number of @p buffer_unit_type units unpacked
   */
  virtual localIndex
  UnpackByIndex( buffer_unit_type const *& buffer,
                 arrayView1d< localIndex const > const & unpackIndices,
                 bool on_device = false ) = 0;

  ///@}

  /**
   * @name Wrapper attribute getters and setters
   */
  ///@{

  /**
   * @brief Check whether this wrapper is resized when its parent is resized.
   * @return @p true if wrapper is resized with parent group, @p false otherwise
   */
  int
  sizedFromParent() const
  {
    return m_sizedFromParent;
  }

  /**
   * @brief Set whether this wrapper is resized when its parent is resized.
   * @param val an int that is converted into a bool
   * @return a pointer to this wrapper
   */
  WrapperBase *
  setSizedFromParent( int val )
  {
    m_sizedFromParent = val;
    return this;
  }

  /**
   * @brief Get the RestartFlags of the wrapper.
   * @return this wrapper's restart flags
   */
  RestartFlags
  getRestartFlags() const
  {
    return m_restart_flags;
  }

  /**
   * @brief Set the RestartFlags of the wrapper.
   * @param flags the new RestartFlags value
   * @return a pointer to this wrapper
   */
  WrapperBase *
  setRestartFlags( RestartFlags flags )
  {
    m_restart_flags = flags;
    return this;
  }

  /**
   * @brief Get PlotLevel for this wrapper.
   * @return this wrapper's plot level
   */
  PlotLevel
  getPlotLevel() const
  {
    return m_plotLevel;
  }

  /**
   * @brief Set the PlotLevel of the wrapper.
   * @param flag the new PlotLevel value
   * @return a pointer to this wrapper
   */
  WrapperBase *
  setPlotLevel( PlotLevel const flag )
  {
    m_plotLevel = flag;
    return this;
  }

  /**
   * @brief Get name of the wrapper.
   * @return name of the wrapper
   */
  string const &
  getName() const
  {
    return m_name;
  }

  /**
   * @brief Set the InputFlag of the wrapper.
   * @param input the new InputFlags value
   * @return a pointer to this wrapper
   */
  WrapperBase *
  setInputFlag( InputFlags const input )
  {
    if( input == InputFlags::OPTIONAL || input == InputFlags::REQUIRED )
    {
      this->setSizedFromParent( 0 );
    }
    m_inputFlag = input;
    return this;
  }

  /**
   * @brief Get the InputFlag of the wrapper.
   * @return this wrapper's input flags
   */
  InputFlags
  getInputFlag() const
  {
    return m_inputFlag;
  }

  /**
   * @brief Set the description string of the wrapper.
   * @param description the description
   * @return a pointer to this wrapper
   */
  WrapperBase *
  setDescription( string const & description )
  {
    m_description = description;
    return this;
  }

  /**
   * @brief Get the description string of the wrapper.
   * @return this wrapper's description string
   */
  string const &
  getDescription() const
  {
    return m_description;
  }

  /**
   * @brief @return a table formatted string containing the input options.
   * @param outputHeader If true outputs the table header, otherwise just
   *                     outputs a row.
   */
  string
  dumpInputOptions( bool const outputHeader ) const;

  /**
   * @brief Get the list of names of groups that registered this wrapper.
   * @return vector of object names
   */
  std::set< string > const &
  getRegisteringObjects() const
  {
    return m_registeringObjects;
  }

  /**
   * @brief Add a new name to the list of groups that register this wrapper.
   * @param objectName name of the registering object
   * @return pointer to this wrapper
   */
  WrapperBase *
  setRegisteringObjects( string const & objectName )
  {
    m_registeringObjects.insert( objectName );
    return this;
  }

  ///@}

  /**
   * @name Miscellaneous
   */
  ///@{

  /**
   * @brief Copy attributes from another wrapper
   * @param[in] source the source wrapper, must wrap the same type @p T
   */
  virtual void
  CopyWrapperAttributes( WrapperBase const & source );

  /**
   * @brief Creates a clone of @p *this WrapperBase
   * @param[in] name name of the clone
   * @param[in] parent parent Group that will hold this clone
   * @return
   *
   * The overridden function will create a copy of the derived Wrapper<T> the using the provided
   * values of name and parent to differentiate itself from the source.
   */
  virtual std::unique_ptr< WrapperBase >
  clone( string const & name,
         Group * const parent ) = 0;

  /**
   * @brief Get the typeid of T.
   * @return type_info of the wrapped type "typeid(T)"
   */
  virtual std::type_info const &
  get_typeid() const = 0;

  ///@}

#if defined( USE_TOTALVIEW_OUTPUT )
  /**
   * @brief Virtual function to return the the typename for a Wrapper derived type that is
   *                represented by a WrapperBase *.
   * @return A string that contains the typename for use in totalview.
   */
  virtual string
  totalviewTypeName() const = 0;

  /**
   * @brief Function to execute the TV_tff_add_row() calls that represent each data member that
   *        will be displayed.
   * @return 0
   */
  virtual int
  setTotalviewDisplay() const;
//  static int TV_ttf_display_type( const WrapperBase * wrapper);
#endif

protected:
  /// @cond DO_NOT_DOCUMENT

  conduit::Node &
  getConduitNode()
  {
    return m_conduitNode;
  }

  /// @endcond

protected:
  /// Name of the object that is being wrapped
  string m_name;

  /// Pointer to Group that holds this WrapperBase
  Group * m_parent;

  /// Integer to indicate whether or not this wrapped object should be resized when m_parent is resized
  int m_sizedFromParent;

  /// Flag to determine the restart behavior for this wrapped object
  RestartFlags m_restart_flags;

  /// Flag to store the plotLevel
  PlotLevel m_plotLevel;

  /// Flag to store if this wrapped object should be read from input
  InputFlags m_inputFlag;

  /// A string description of the wrapped object
  string m_description;

  /// A vector of the names of the objects that created this Wrapper.
  std::set< string > m_registeringObjects;

  /// A reference to the corresponding conduit::Node.
  conduit::Node & m_conduitNode;
};

}  // namespace dataRepository
}  // namespace geosx

#endif /* GEOSX_DATAREPOSITORY_WRAPPERBASE_HPP_ */
