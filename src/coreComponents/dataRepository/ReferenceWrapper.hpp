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
 * @file ReferenceWrapper.hpp
 * This file contains the class definition of ReferenceWrapper.
 */
#ifndef GEOSX_DATAREPOSITORY_REFERENCEWRAPPER_HPP_
#define GEOSX_DATAREPOSITORY_REFERENCEWRAPPER_HPP_

// System includes
#include <type_traits>

namespace geosx
{
/**
 * @class ReferenceWrapper
 * @tparam Type that is wrapped
 *
 * This class manages a pointer to the templated type, but provides a reference-like interface,
 * thus negating the requirement to dereference the object. The primary use for this is for nested
 * arrays that hold pointers at the last level, but allows for reference-like usage. For instance,
 * consider a collection of object that you would like to refer to thought an array of pointers.
 *
 * <tt>array1d< ReferenceWrapper< array1d< double > > > arr;</tt>
 *
 * where the <tt>array::operator[]</tt> exists. The ReferenceWrapper allows
 *
 * <tt>arr[index1][index2]</tt>
 *
 * note: this is really only useful for heavy array that hold their own data. For light arrays that
 * hold pointers to their data, then this is unnecessary as a copy of the array does not trigger a
 * deep copy.
 */
template< typename T >
class ReferenceWrapper
{
public:
  /**
   * @brief Default constructor sets m_ref to nullptr.
   */
  ReferenceWrapper() :
    m_ref( nullptr )
  {}

  /**
   * @brief Constructor that sets m_ref to address of input.
   * @param[in] source object to wrap
   */
  ReferenceWrapper( T & source ) noexcept :
    m_ref( &source )
  {}

  /**
   * @brief Default destructor.
   */
  ~ReferenceWrapper() = default;

  /**
   * @brief Copy constructor copies the source m_ref to the new m_ref.
   * @param[in] source object to copy
   */
  ReferenceWrapper( ReferenceWrapper const & source ) :
    m_ref( source.m_ref )
  {}

  /**
   * @brief Move constructor copies the source m_ref to the new m_ref.
   * @param[in,out] source object to move from
   */
  ReferenceWrapper( ReferenceWrapper && source ) :
    m_ref( source.m_ref )
  {
    source.m_ref = nullptr;
  }

  /**
   * @brief Copy assignment operator
   * @param[in] source object to copy
   * @return
   */
  ReferenceWrapper &
  operator=( ReferenceWrapper const & source )
  {
    m_ref = source.m_ref;
    return *this;
  }

  /**
   * @tparam T_RHS type of the rhs
   * @brief Assignment operator.
   * @tparam U dummy template parameter to enable SFINAE stuff
   * @param[in] rhs value to be copied
   * @return *this
   *
   * Calls m_ref->operator=() to allow for any type on the rhs
   * if m_ref->operator=() has a valid overload for T_RHS.
   */
  template< typename T_RHS, typename U = T >
  inline
    typename std::enable_if< !std::is_const< U >::value, ReferenceWrapper & >::type
    operator=( T_RHS const & rhs )
  {
    *m_ref = rhs;
    return *this;
  }

  /**
   * @brief Move assignment operator.
   * @param[in,out] source the rhs value to be moved
   * @return reference to this object
   *
   * Sets the value that m_ref refers to to the value of the rhs.
   */
  inline ReferenceWrapper &
  operator=( T && source )
  {
    *m_ref = std::move( source );
    return *this;
  }

  /**
   * @brief User defined conversion to <tt>T &</tt>.
   */
  inline
  operator T &()
  {
    return *m_ref;
  }

  /**
   * @brief User defined conversion to <tt>T const &</tt>
   */
  inline operator T const &() const
  {
    return *m_ref;
  }

  /**
   * @brief Set the address that m_ref points to.
   * @param[in] source reference to object that wrapper will refer to
   */
  inline void
  set( T & source )
  {
    m_ref = &source;
  }

  /**
   * @brief Set the address that m_ref points to.
   * @param[in] source pointer to object that wrapper will refer to
   */
  inline void
  set( T * source )
  {
    m_ref = source;
  }

  /**
   * @brief Accessor for m_ref.
   * @return reference to wrapped value
   */
  inline T &
  get()
  {
    return *m_ref;
  }

  /**
   * @brief Const accessor for m_ref.
   * @return const reference to wrapped value
   */
  inline T const &
  get() const
  {
    return *m_ref;
  }

  /**
   * @brief Check if reference is initialized.
   * @return @p true if the object has been initialized with a value, @p false otherwise
   */
  inline bool
  isValid() const
  {
    return m_ref;
  }

  /**
   * @brief Const accessor for m_ref.
   * @return const reference to wrapped value
   */
  inline T const *
  getPtr() const
  {
    return m_ref;
  }

  /*
   * Unfortunately, Doxygen does not understand decltype in function return types.
   * It does not generate documentation for the following two functions, but still
   * emits a warning. So we just disable docs for these until issue fixed in Doxygen.
   */
  /// @cond DO_NOT_DOCUMENT

  /**
   * @brief A pass thru square bracket operator which calls underlying <tt>T::operator[]</tt>.
   * @tparam U dummy type to allow for SFINAE evaluation of availability of <tt>T::operator[]</tt>
   * @param[in] i index to pass into the <tt>T::operator[]</tt>
   * @return the return type of <tt>T::operator[]</tt>
   */
  template< typename INDEX_TYPE, typename U = T >
  inline decltype( std::declval< U >()[1] )
  operator[]( INDEX_TYPE const i )
  {
    return ( *m_ref )[i];
  }

  /**
   * @brief A const pass thru square bracket operator which calls underlying <tt>T::operator[]</tt>.
   * @tparam U dummy type to allow for SFINAE evaluation of availability of <tt>T::operator[]</tt>
   * @param[in] i index to pass into the <tt>T::operator[]</tt> const
   * @return the return type of <tt>T::operator[]</tt> const
   */
  template< typename INDEX_TYPE, typename U = T >
  inline decltype( std::declval< U const >()[1] )
  operator[]( INDEX_TYPE const i ) const
  {
    return ( *m_ref )[i];
  }

  /// @endcond

  /**
   * @brief A pass thru parenthesis  operator which calls underlying <tt>T::operator()</tt>.
   * @tparam variadic types to pass through to <tt>T::operator()</tt>
   * @param args variadic params to pass through to <tt>T::operator()</tt>
   * @return the return type of <tt>T::operator()</tt>
   */
  template< typename... ARGS >
  inline typename std::result_of< T &( ARGS &&... ) >::type
  operator()( ARGS &&... args )
  {
    return m_ref->operator()( std::forward< ARGS >( args )... );
  }

  /**
   * @brief A pass thru parenthesis operator which calls underlying <tt>T::operator()</tt>.
   * @tparam variadic types to pass through to <tt>T::operator()</tt>
   * @param args variadic params to pass through to <tt>T::operator()</tt>
   * @return the return type of <tt>T::operator()</tt> const
   */
  template< typename... ARGS >
  inline typename std::result_of< T const &( ARGS &&... ) >::type
  operator()(
    ARGS &&... args ) const
  {
    return m_ref->operator()( std::forward< ARGS >( args )... );
  }

private:
  /// pointer to the address of the object that we would like to wrap
  T * m_ref;
};

} /* namespace geosx */

#endif /* GEOSX_DATAREPOSITORY_REFERENCEWRAPPER_HPP_ */
