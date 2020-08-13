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
 * @file BufferAllocator.hpp
 */

#ifndef GEOSX_BUFFER_ALLOCATOR_HPP
#define GEOSX_BUFFER_ALLOCATOR_HPP

#include "common/GeosxConfig.hpp"

#ifdef GEOSX_USE_CHAI
  #include <umpire/ResourceManager.hpp>
  #include <umpire/TypedAllocator.hpp>

namespace geosx
{
/**
 * @brief Set the current desired behaviour of the BufferAllocator
 * @param p Whether or not BufferAllocators should be instantiated
 *          with a preference for using pinned memory.
 */
void
setPreferPinned( bool p );

/**
 * @brief Get the current desired behaviour of the BufferAllocator
 * @return Whether or not BufferAllocators should be instantiated
 *         with a preference for using pinned memory.
 */
bool
getPreferPinned();

/**
 * @brief Wrapper class for umpire allocator, only used to determine which umpire allocator to use based on
 * availability.
 * @note It would be better in some ways to have this be a singleton simply because the underlying umpire allocators
 *       are also singletons. It currently isn't mostly because we use the default constructor in the buffer_type
 *       instantiations in the neighbor comm (which can be refactored to use a singleton instance-getter fairly easily).
 *       In general though it would mean any trivial usage of the buffer_type this is used in more cumbersome.
 *       At the moment our usage looks like:
 *         buffer_type(size)
 *       but if this was a singleton it would necessitate that we create buffers as: (outside of neighbor comm)
 *         buffer_type(size,BufferAllocator<buffer_type>::getInstance())
 *       which is messy. Changing the buffer_type can fix the issue but would introduce additional refactoring so for
 *       the moment this implementation suffices.
 */
template< typename T >
class BufferAllocator
{
public:
  // The type used to instantiate the class, and the underlying umpire allocator.
  using value_type = T;

private:
  // An umpire allocator allocating the type for which this class is instantiated.
  umpire::TypedAllocator< value_type > m_alloc;
  bool m_prefer_pinned_l;

public:
  /**
   * @brief Default behavior is to allocate host memory, if there is a pinned memory allocator
   *        provided by umpire for the target platform, and getPreferPinned returns true,
   *        use that instead.
   */
  BufferAllocator() :
    m_alloc( umpire::TypedAllocator< T >(
      umpire::ResourceManager::getInstance().getAllocator(
        umpire::resource::Host ) ) ),
    m_prefer_pinned_l( getPreferPinned() )
  {
    auto & rm = umpire::ResourceManager::getInstance();
    if( rm.isAllocator( "PINNED" ) && m_prefer_pinned_l )
      m_alloc =
        umpire::TypedAllocator< T >( rm.getAllocator( umpire::resource::Pinned ) );
  }
  /**
   * @brief Allocate a buffer.
   * @param sz The number of elements of type value_type to allocate a buffer for.
   * @return A pointer to the allocated buffer.
   */
  value_type *
  allocate( size_t sz )
  {
    return m_alloc.allocate( sz );
  }
  /**
   * @brief Deallocate a buffer.
   * @param buffer A pointer to the buffer to deallocate
   * @param sz The size of the buffer to deallocate.
   */
  void
  deallocate( value_type * buffer, size_t sz )
  {
    if( buffer != nullptr ) m_alloc.deallocate( buffer, sz );
  }
  /**
   * @brief Inequality operator.
   * @param The other BufferAllocator to test against this buffer allocator for inequality.
   * @return Always false. Since the actual umpire allocator is a singleton, so any properly-typed
   *         buffer can be deallocated from any properly-typed BufferAllocator.
   */
  bool
  operator!=( const BufferAllocator & )
  {
    return false;
  }
  /**
   * @brief Equality operator.
   * @param other The other BufferAllocator to test against this buffer allocator for equality.
   * return Always true. Since the actual umpire allocator is a singleton, so any properly-typed
   *        buffer can be deallocated from any properly-typed BufferAllocator.
   */
  bool
  operator==( const BufferAllocator & other )
  {
    return !operator!=( other );
  }
};

}  // namespace geosx
#endif

#endif
