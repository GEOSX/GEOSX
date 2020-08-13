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
 * @file MappedVector.hpp
 */

#ifndef GEOSX_DATAREPOSITORY_MAPPEDVECTOR_HPP_
#define GEOSX_DATAREPOSITORY_MAPPEDVECTOR_HPP_

// Source includes
#include "KeyIndexT.hpp"
#include "common/GeosxMacros.hpp"
#include "common/Logger.hpp"
#include "LvArray/src/limits.hpp"

// System includes
#include <vector>

namespace geosx
{
/**
 * @class MappedVector
 *
 * This class defines a stl-like container that stores values in an stl vector,
 * and has a map lookup table to access the values by a key. It combines the
 * random access performance of a vector when the index is known, the flexibility
 * of a mapped key lookup O(n) if only the key is known.
 *
 * In addition, a keyIndex can be used for lookup, which will give similar
 * performance to an index lookup after the first use of a keyIndex.
 */
template< typename T,
          typename T_PTR = T *,
          typename KEY_TYPE = std::string,
          typename INDEX_TYPE = int >
class MappedVector
{
public:
  static_assert( std::is_same< T_PTR, T * >::value ||
                   std::is_same< T_PTR, std::unique_ptr< T > >::value,
                 "invalid second template argument for "
                 "MappedVector<T,T_PTR,KEY_TYPE,INDEX_TYPE>. Allowable types "
                 "are T * and std::unique_ptr<T>." );

  /// The type used for the key of the map
  using key_type = KEY_TYPE;

  /// pointer to the value type
  using mapped_type = T_PTR;

  /// the type of the lookup map
  using LookupMapType = std::unordered_map< KEY_TYPE, INDEX_TYPE >;

  /// the type of the values held in the vector
  using value_type = typename std::pair< KEY_TYPE, T_PTR >;

  /// the type of the values with const keys held in the vector
  using const_key_value_type = typename std::pair< KEY_TYPE const, T * >;

  /// a const type of the values held in the vector
  using const_value_type = typename std::pair< KEY_TYPE const, T const * >;

  /// the type of the vector container
  using valueContainer = std::vector< value_type >;

  /// a const type of the vector container
  using constKeyValueContainer = std::vector< const_key_value_type >;

  /// a const type of the vector container
  using constValueContainer = std::vector< const_value_type >;

  /// the pointer type of the value container
  using pointer = typename valueContainer::pointer;

  /// the pointer to const type of the value container
  using const_pointer = typename valueContainer::const_pointer;

  /// reference type of the value container
  using reference = typename valueContainer::reference;

  /// reference to const type of the value container
  using const_reference = typename valueContainer::const_reference;

  /// the size_type of the value container
  using size_type = typename valueContainer::size_type;

  /// the iterator type of the value container
  using iterator = typename constKeyValueContainer::iterator;

  /// the iterator to const type of the value container
  using const_iterator = typename constValueContainer::const_iterator;

  /// the reverse iterator type of the value container
  using reverse_iterator = typename constKeyValueContainer::reverse_iterator;

  /// the reverse iterator to const type of the value container
  using const_reverse_iterator =
    typename constValueContainer::const_reverse_iterator;

  /// alias for the KeyIndex itself
  using KeyIndex = KeyIndexT< KEY_TYPE const, INDEX_TYPE >;

  /// deleted default constructor
  MappedVector() = default;

  /// default destructor
  ~MappedVector()
  {
    clear();
  }

  /**
   * @brief Default copy constructor.
   */
  MappedVector( MappedVector const & ) = default;

  /**
   * @brief Default copy assignment operator.
   * @return
   */
  MappedVector &
  operator=( MappedVector const & ) = default;

  /**
   * @brief Default move operator.
   */
  MappedVector( MappedVector && ) = default;

  /**
   * @brief Default move assignment operator.
   * @return
   */
  MappedVector &
  operator=( MappedVector && ) = default;

  /**
   * @name element access functions
   */
  ///@{

  /**
   * @param index
   * @return pointer to const T
   */
  inline T const *
  operator[]( INDEX_TYPE index ) const
  {
    return ( index > KeyIndex::invalid_index &&
             index < static_cast< INDEX_TYPE >( m_values.size() ) )
      ? const_cast< T const * >( &( *( m_values[index].second ) ) )
      : nullptr;
  }

  /**
   *
   * @param index
   * @return pointer to T
   */
  inline T *
  operator[]( INDEX_TYPE index )
  {
    return const_cast< T * >(
      const_cast< MappedVector< T, T_PTR, KEY_TYPE, INDEX_TYPE > const * >( this )->
      operator[]( index ) );
  }

  /**
   *
   * @param keyName
   * @return pointer to const T
   */
  inline T const *
  operator[]( KEY_TYPE const & keyName ) const
  {
    typename LookupMapType::const_iterator iter = m_keyLookup.find( keyName );
    return ( iter != m_keyLookup.end() ? this->operator[]( iter->second ) : nullptr );
  }

  /**
   *
   * @param keyName
   * @return pointer to T
   */
  inline T *
  operator[]( KEY_TYPE const & keyName )
  {
    return const_cast< T * >(
      const_cast< MappedVector< T, T_PTR, KEY_TYPE, INDEX_TYPE > const * >( this )->
      operator[]( keyName ) );
  }

  /**
   *
   * @param keyIndex
   * @return pointer to const T
   */
  inline T const *
  operator[]( KeyIndex const & keyIndex ) const
  {
    INDEX_TYPE index = keyIndex.Index();

    if( index == KeyIndex::invalid_index )
    {
      index = getIndex( keyIndex.Key() );
      keyIndex.setIndex( index );
    }
#ifdef MAPPED_VECTOR_RANGE_CHECKING
    else if( m_values[index].first != keyIndex.Key() )
    {
      index = getIndex( keyIndex.Key() );
      keyIndex.setIndex( index );
    }
#endif

    return this->operator[]( index );
  }

  /**
   *
   * @param keyIndex
   * @return pointer to T
   */
  inline T *
  operator[]( KeyIndex const & keyIndex )
  {
    return const_cast< T * >(
      const_cast< MappedVector< T, T_PTR, KEY_TYPE, INDEX_TYPE > const * >( this )->
      operator[]( keyIndex ) );
  }

  ///@}

  /**
   * @name iterator functions
   */
  ///@{

  /**
   * @return  a read/write iterator that points to the first
   *  element in the in m_objects.
   */
  iterator
  begin()
  {
    return m_constKeyValues.begin();
  }

  /**
   * @return  a read-only iterator that points to the first
   *  element in the in m_objects.
   */
  const_iterator
  begin() const
  {
    return m_constValues.begin();
  }

  /**
   * @return  a read-only iterator that points to the first
   *  element in the in m_objects.
   */
  const_iterator
  cbegin()
  {
    return m_constValues.begin();
  }

  /**
   * @return  a read/write iterator that points to the last
   *  element in the in m_objects.
   */
  iterator
  end()
  {
    return m_constKeyValues.end();
  }

  /**
   * @return  a read-only iterator that points to the last
   *  element in the in m_objects.
   */
  const_iterator
  end() const
  {
    return m_constValues.end();
  }

  /**
   * @return  a read-only iterator that points to the last
   *  element in the in m_objects.
   */
  const_iterator
  cend()
  {
    return m_constValues.end();
  }
  ///@}

  /**
   *
   * @param key value of the key to use in the lookup
   * @return index associated with key
   */
  inline INDEX_TYPE
  getIndex( KEY_TYPE const & key ) const
  {
    typename LookupMapType::const_iterator iter = m_keyLookup.find( key );
    return ( iter != m_keyLookup.end() ? iter->second : KeyIndex::invalid_index );
  }

  /**
   * @name modifier functions
   */
  ///@{

  /**
   * @brief insert new entry into MappedVector
   * @param keyName key name to assocaite with the new object
   * @param source pointer to object
   * @param takeOwnership whether or not to take ownership of the object
   * @param overwrite if the key already exists, overwrite indicates whether or not overwrite the existing entry
   * @return pointer to the object that is held in the MappedVector
   */
  T *
  insert( KEY_TYPE const & keyName,
          T_PTR source,
          bool takeOwnership,
          bool overwrite = false );

  /**
   * @brief  Remove element at given index
   * @param  index  index of element to remove.
   *
   * Completely remove element at given index and corresponding key lookup.
   * If pointed-to object is owned, it is deleted.
   */
  void
  erase( INDEX_TYPE index )
  {
    // delete the pointed-to value, if owned
    deleteValue( index );

    // delete lookup entry
    m_keyLookup.erase( m_values[index].first );

    // delete and shift vector entries
    m_values.erase( m_values.begin() + index );
    m_ownsValues.erase( m_ownsValues.begin() + index );

    // rebuild parts of const key vectors after deleted entry
    m_constKeyValues.resize( index );
    m_constValues.resize( index );
    for( typename valueContainer::size_type i = index; i < m_values.size(); ++i )
    {
      m_constKeyValues.emplace_back( m_values[i].first, rawPtr( index ) );
      m_constValues.emplace_back( m_values[i].first, rawPtr( index ) );
    }

    // adjust lookup map indices
    for( typename valueContainer::size_type i = index; i < m_values.size(); ++i )
    {
      m_keyLookup[m_values[i].first] = i;
    }
  }

  /**
   *  @brief  Remove element at given key
   *  @param  key  key of element to remove.
   *
   *  This function will set the element at the given key to nullptr.
   */
  void
  erase( KEY_TYPE const & key )
  {
    typename LookupMapType::const_iterator iter = m_keyLookup.find( key );
    if( iter != m_keyLookup.end() )
    {
      erase( iter->second );
    }
  }

  /**
   *  @brief  Remove element at given key
   *  @param  keyIndex  key of element to remove.
   *
   *  This function will set the element at the given key to nullptr.
   */
  void
  erase( KeyIndex & keyIndex )
  {
    INDEX_TYPE index = keyIndex.Index();

    if( ( index == KeyIndex::invalid_index ) ||
        ( m_values[index].first != keyIndex.Key() ) )
    {
      index = getIndex( keyIndex.Key() );
      keyIndex.setIndex( index );
    }
    erase( index );
  }

  /**
   * @brief function to clear the MappedVector
   */
  void
  clear()
  {
    for( typename valueContainer::size_type a = 0; a < m_values.size(); ++a )
    {
      deleteValue( LvArray::integerConversion< INDEX_TYPE >( a ) );
    }
    m_constKeyValues.clear();
    m_constValues.clear();
    m_values.clear();
    m_ownsValues.clear();
    m_keyLookup.clear();
  }

  ///@}

  /**
   * @brief function to return the number of entries stored
   * @return number of entries in MappedVector
   */
  inline INDEX_TYPE
  size() const
  {
    return LvArray::integerConversion< INDEX_TYPE >( m_values.size() );
  }

  /**
   * @brief access for value container
   * @return reference to valueContainer
   */
  inline valueContainer const &
  values()
  {
    return this->m_values;
  }

  /**
   * @brief access for value container
   * @return reference to const valueContainer
   */
  inline constValueContainer const &
  values() const
  {
    return this->m_constValues;
  }

  /**
   * @brief access for key lookup
   * @return reference lookup map
   */
  inline LookupMapType const &
  keys() const
  {
    return m_keyLookup;
  }

private:
  T *
  rawPtr( INDEX_TYPE index )
  {
    return &( *( m_values[index].second ) );
  }

  template< typename U = T_PTR >
  typename std::enable_if< std::is_same< U, T * >::value, void >::type
  deleteValue(
    INDEX_TYPE index )
  {
    if( m_ownsValues[index] )
    {
      delete m_values[index].second;
    }
  }

  template< typename U = T_PTR >
  typename std::enable_if< !std::is_same< U, T * >::value, void >::type
  deleteValue(
    INDEX_TYPE GEOSX_UNUSED_PARAM( index ) )
  {}

  /// random access container that holds the values
  valueContainer m_values;

  /// clone of random access container that holds const keys
  constKeyValueContainer m_constKeyValues;

  /// clone of random access container that holds const keys and pointer to const values
  constValueContainer m_constValues;

  /// map lookup to go from key to index into the m_values container
  LookupMapType m_keyLookup;

  /// flag to indicate whether or not the values in m_values are owned by the container.
  std::vector< int > m_ownsValues;
};

template< typename T, typename T_PTR, typename KEY_TYPE, typename INDEX_TYPE >
T *
MappedVector< T, T_PTR, KEY_TYPE, INDEX_TYPE >::insert( KEY_TYPE const & keyName,
                                                        T_PTR source,
                                                        bool takeOwnership,
                                                        bool overwrite )
{
  INDEX_TYPE index = KeyIndex::invalid_index;
  typename LookupMapType::iterator iterKeyLookup = m_keyLookup.find( keyName );

  // if the key was not found, make DataObject<T> and insert
  if( iterKeyLookup == m_keyLookup.end() )
  {
    value_type newEntry = std::make_pair( keyName, std::move( source ) );
    m_values.push_back( std::move( newEntry ) );
    //TODO this needs to be a safe conversion
    index = static_cast< INDEX_TYPE >( m_values.size() ) - 1;
    m_ownsValues.resize( index + 1 );
    if( takeOwnership )
    {
      m_ownsValues[index] = true;
    }

    m_keyLookup.insert( std::make_pair( keyName, index ) );
    m_constKeyValues.emplace_back( keyName, rawPtr( index ) );
    m_constValues.emplace_back( keyName, rawPtr( index ) );
  }
  // if key was found
  else
  {
    index = iterKeyLookup->second;

    if( takeOwnership )
    {
      m_ownsValues[index] = true;
    }

    // if value is empty, then move source into value slot
    if( m_values[index].second == nullptr )
    {
      m_values[index].second = std::move( source );
      m_constKeyValues[index].second = rawPtr( index );
      m_constValues[index].second = rawPtr( index );
    }
    else
    {
      if( overwrite )
      {
        deleteValue( index );
        m_values[index].second = std::move( source );
        m_constKeyValues[index].second = rawPtr( index );
        m_constValues[index].second = rawPtr( index );
      }
      else if( source->get_typeid() != m_values[index].second->get_typeid() )
      {
        GEOSX_ERROR( "MappedVector::insert(): Tried to insert existing key ("
                     << keyName << ") with a different type without overwrite flag\n"
                     << " " << source->get_typeid().name()
                     << " != " << m_values[index].second->get_typeid().name() );
      }
      else
      {
        delete source;
      }
    }
  }

  return &( *( m_values[index].second ) );
}
}  // namespace geosx

#endif /* GEOSX_DATAREPOSITORY_MAPPEDVECTOR_HPP_ */
