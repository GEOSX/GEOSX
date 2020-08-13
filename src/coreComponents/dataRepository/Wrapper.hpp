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
 * @file Wrapper.hpp
 */

#ifndef GEOSX_DATAREPOSITORY_WRAPPER_HPP_
#define GEOSX_DATAREPOSITORY_WRAPPER_HPP_

// Source inclues
#include "wrapperHelpers.hpp"
#include "KeyNames.hpp"
#include "LvArray/src/limits.hpp"
#include "common/DataTypes.hpp"
#include "codingUtilities/SFINAE_Macros.hpp"
#include "LvArray/src/Macros.hpp"
#include "BufferOps.hpp"
#include "BufferOpsDevice.hpp"
#include "RestartFlags.hpp"
#include "codingUtilities/traits.hpp"
#include "common/GeosxConfig.hpp"
#include "DefaultValue.hpp"
#include "LvArray/src/system.hpp"
#include "WrapperBase.hpp"

// System includes
#include <type_traits>
#include <cstdlib>
#include <type_traits>

namespace geosx
{
namespace dataRepository
{
//template< typename U >
//static void totalViewType( char * const dataType );

/**
 * Templated class to serve as a wrapper to arbitrary objects.
 * @tparam T is any type that is to be wrapped by Wrapper
 */
template< typename T >
class Wrapper final : public WrapperBase
{
public:
  /**
   * @brief Alias for the wrapped type @p T
   */
  using TYPE = T;

  /**
   * @name Constructors, destructor, copy/move assignment
   */
  ///@{

  /**
   * @brief Constructor that creates a new instance of wrapped type
   * @param name name of the object
   * @param parent parent group which owns the Wrapper
   */
  explicit Wrapper( std::string const & name, Group * const parent ) :
    WrapperBase( name, parent ),
    m_ownsData( true ),
    m_data( new T() ),
    m_default()
  {
    if( traits::is_tensorT< T > || std::is_arithmetic< T >::value ||
        traits::is_string< T > )
    {
      this->setSizedFromParent( 0 );
    }

    setName();
  }

  /**
   * @brief Constructor that takes ownership of an existing object.
   * @param name name of the object
   * @param parent parent group that owns the Wrapper
   * @param object object that is being wrapped by the Wrapper
   */
  explicit Wrapper( std::string const & name,
                    Group * const parent,
                    std::unique_ptr< T > object ) :
    WrapperBase( name, parent ),
    m_ownsData( true ),
    m_data( object.release() ),
    m_default()
  {
    if( traits::is_tensorT< T > || std::is_arithmetic< T >::value ||
        traits::is_string< T > )
    {
      this->setSizedFromParent( 0 );
    }

    setName();
  }

  /**
   * @brief Constructor that does not take ownership of an existing object.
   * @param name name of the object
   * @param parent parent group that owns the Wrapper
   * @param object object that is being wrapped by the Wrapper
   */
  explicit Wrapper( std::string const & name, Group * const parent, T * object ) :
    WrapperBase( name, parent ),
    m_ownsData( false ),
    m_data( object ),
    m_default()
  {
    if( traits::is_tensorT< T > || std::is_arithmetic< T >::value ||
        traits::is_string< T > )
    {
      this->setSizedFromParent( 0 );
    }

    setName();
  }

  /**
   * @brief Default destructor
   *
   * Deletes wrapped object if the wrapper is owning
   */
  virtual ~Wrapper() noexcept override
  {
    if( m_ownsData )
    {
      delete m_data;
    }
    //tvTemplateInstantiation();
  }

  /**
   * @brief Copy Assignment Operator
   *
   * @param source rhs
   * @return @p *this
   */
  Wrapper &
  operator=( Wrapper const & source )
  {
    m_data = source.m_data;
    return *this;
  }

  /**
   * @brief Move Assignment Operator
   * @param source
   * @return *this
   */
  Wrapper &
  operator=( Wrapper && source )
  {
    m_data = std::move( source.m_data );
    return *this;
  }

  ///@}

  /**
   * @name Miscellaneous
   */
  ///@{

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual std::unique_ptr< WrapperBase >
  clone( string const & name,
         Group * const parent ) override
  {
    std::unique_ptr< WrapperBase > clonedWrapper =
      std::make_unique< Wrapper< T > >( name, parent, this->m_data );
    clonedWrapper->CopyWrapperAttributes( *this );

    return clonedWrapper;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual void
  CopyWrapperAttributes( WrapperBase const & source ) override
  {
    WrapperBase::CopyWrapperAttributes( source );
    Wrapper< T > const & castedSource = *cast( &source );
    m_ownsData = castedSource.m_ownsData;
    m_default = castedSource.m_default;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual const std::type_info &
  get_typeid() const noexcept override
  {
    return typeid( T );
  }

  ///@}

  /**
   * @name Type-casting static functions
   */
  ///@{

  /**
   * @brief Static function to cast a Wrapper base to a derived Wrapper<T>
   *
   * @param base
   * @return casted Wrapper<T>
   */
  static Wrapper< T > *
  cast( WrapperBase * const base )
  {
    return dynamicCast< Wrapper< T > * >( base );
  }

  /**
   * @brief Static function to cast a Wrapper base to a derived Wrapper<T>
   *
   * @param base
   * @return casted reference to const Wrapper<T>
   */
  static Wrapper< T > const *
  cast( WrapperBase const * const base )
  {
    return dynamicCast< Wrapper< T > const * >( base );
  }

  /**
   * @brief Static function to cast a Wrapper base to a derived Wrapper<T>
   *
   * @param base
   * @return casted Wrapper<T>
   */
  static Wrapper< T > &
  cast( WrapperBase & base )
  {
    return dynamicCast< Wrapper< T > & >( base );
  }

  /**
   * Static function to cast a Wrapper base to a derived Wrapper<T>
   * @param base
   * @return casted reference to const Wrapper<T>
   */
  static Wrapper< T > const &
  cast( WrapperBase const & base )
  {
    return dynamicCast< Wrapper< T > const & >( base );
  }

  ///@}

  /**
   * @name Packing/unpacking methods.
   */
  ///@{

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual bool
  isPackable( bool onDevice ) const override
  {
    if( onDevice )
    {
      return bufferOps::can_memcpy< T >;
    }
    else
    {
      return bufferOps::is_packable< T >;
    }
  }

  /**
   * @brief function to pack T
   * @param buffer the buffer in which to pack T
   * @param[in] onDevice    whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @return number of packed bytes.
   */
  virtual localIndex
  Pack( buffer_unit_type *& buffer, bool onDevice ) const override
  {
    localIndex packedSize = 0;
    packedSize += bufferOps::Pack< true >( buffer, this->getName() );
    if( onDevice )
    {
      packedSize += wrapperHelpers::PackDevice< true >( buffer, reference() );
    }
    else
    {
      packedSize += bufferOps::Pack< true >( buffer, *m_data );
    }
    return packedSize;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual localIndex
  PackByIndex( buffer_unit_type *& buffer,
               arrayView1d< localIndex const > const & packList,
               bool onDevice ) const override
  {
    localIndex packedSize = 0;
    if( sizedFromParent() == 1 )
    {
      packedSize += bufferOps::Pack< true >( buffer, this->getName() );
      if( onDevice )
      {
        packedSize +=
          wrapperHelpers::PackByIndexDevice< true >( buffer, reference(), packList );
      }
      else
      {
        packedSize +=
          wrapperHelpers::PackByIndex< true >( buffer, *m_data, packList );
      }
    }
    return packedSize;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual localIndex
  PackSize( bool onDevice ) const override
  {
    buffer_unit_type * buffer = nullptr;
    localIndex packedSize = 0;
    packedSize += bufferOps::Pack< false >( buffer, this->getName() );
    if( onDevice )
    {
      packedSize += wrapperHelpers::PackDevice< false >( buffer, reference() );
    }
    else
    {
      packedSize += bufferOps::Pack< false >( buffer, *m_data );
    }
    return packedSize;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual localIndex
  PackByIndexSize( arrayView1d< localIndex const > const & packList,
                   bool onDevice ) const override
  {
    localIndex packedSize = 0;
    buffer_unit_type * buffer = nullptr;
    if( sizedFromParent() == 1 )
    {
      packedSize += bufferOps::Pack< false >( buffer, this->getName() );
      if( onDevice )
      {
        packedSize +=
          wrapperHelpers::PackByIndexDevice< false >( buffer, reference(), packList );
      }
      else
      {
        packedSize +=
          wrapperHelpers::PackByIndex< false >( buffer, *m_data, packList );
      }
    }
    return packedSize;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual localIndex
  Unpack( buffer_unit_type const *& buffer, bool onDevice ) override
  {
    localIndex unpackedSize = 0;
    string name;
    unpackedSize += bufferOps::Unpack( buffer, name );
    GEOSX_ERROR_IF( name != this->getName(),
                    "buffer unpack leads to wrapper names that don't match" );
    if( onDevice )
    {
      unpackedSize += wrapperHelpers::UnpackDevice( buffer, referenceAsView() );
    }
    else
    {
      unpackedSize += bufferOps::Unpack( buffer, *m_data );
    }
    return unpackedSize;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual localIndex
  UnpackByIndex( buffer_unit_type const *& buffer,
                 arrayView1d< localIndex const > const & unpackIndices,
                 bool onDevice ) override
  {
    localIndex unpackedSize = 0;
    if( sizedFromParent() == 1 )
    {
      string name;
      unpackedSize += bufferOps::Unpack( buffer, name );
      GEOSX_ERROR_IF( name != this->getName(),
                      "buffer unpack leads to wrapper names that don't match" );
      if( onDevice )
      {
        unpackedSize += wrapperHelpers::UnpackByIndexDevice( buffer,
                                                             referenceAsView(),
                                                             unpackIndices );
      }
      else
      {
        unpackedSize +=
          wrapperHelpers::UnpackByIndex( buffer, *m_data, unpackIndices );
      }
    }
    return unpackedSize;
  }

  ///@}

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  void const *
  voidPointer() const override
  {
    return wrapperHelpers::dataPtr( reference() );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual localIndex
  elementByteSize() const override
  {
    return wrapperHelpers::byteSizeOfElement< T >();
  }

  /**
   * @name Methods that delegate to the wrapped type
   *
   * These functions will call the corresponding method on the wrapped
   * object, if such method is declared in wrapped type.
   */
  ///@{

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual localIndex
  size() const override
  {
    return wrapperHelpers::size( *m_data );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual void
  resize( int ndims, localIndex const * const dims ) override
  {
    wrapperHelpers::move( *m_data, LvArray::MemorySpace::CPU, true );
    wrapperHelpers::resizeDimensions( *m_data, ndims, dims );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual void
  reserve( localIndex const newCapacity ) override
  {
    wrapperHelpers::move( *m_data, LvArray::MemorySpace::CPU, true );
    wrapperHelpers::reserve( reference(), newCapacity );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual localIndex
  capacity() const override
  {
    // We don't use reference() here because that would return an ArrayView which has no capacity method.
    return wrapperHelpers::capacity( *m_data );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual void
  resize( localIndex const newSize ) override
  {
    wrapperHelpers::move( *m_data, LvArray::MemorySpace::CPU, true );
    wrapperHelpers::resizeDefault( reference(), newSize, m_default );
  }

  /// @cond DO_NOT_DOCUMENT
  struct copy_wrapper
  {
    template< typename U, int NDIM, typename PERMUTATION >
    static void
    copy( Array< U, NDIM, PERMUTATION > const & array,
          localIndex const sourceIndex,
          localIndex const destIndex )
    {
      LvArray::forValuesInSliceWithIndices(
        array[sourceIndex],
        [destIndex, &array]( U const & sourceVal, auto const... indices ) {
          array( destIndex, indices... ) = sourceVal;
        } );
    }

    template< typename U >
    static void
    copy( U const &, localIndex const, localIndex const )
    {}
  };
  /// @endcond

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual void
  copy( localIndex const sourceIndex,
        localIndex const destIndex ) override
  {
    if( this->sizedFromParent() )
    {
      copy_wrapper::copy( reference(), sourceIndex, destIndex );
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual void
  move( LvArray::MemorySpace const space,
        bool const touch ) const override
  {
    return wrapperHelpers::move( reference(), space, touch );
  }

  ///@}

  /**
   * @name Methods for wrapped value data access
   */
  ///@{

  /**
   * @brief Accessor for m_data
   * @return reference to T
   */
  T &
  reference()
  {
    return *m_data;
  }

  /**
   * @brief const Accessor for m_data
   * @return reference to T, or in the case of an Array, a reference to an
   *         ArrayView<T const> const.
   */
  traits::ViewTypeConst< T >
  reference() const
  {
    return referenceAsView();
  }

  /**
   * @brief Provide access to wrapped object converted to a view, if possible.
   * @tparam U dummy template parameter to enable SFINAE (do not change)
   * @return the view for wrapped object
   *
   * This is used mainly for @p LvArray classes (arrays, etc.) that can convert
   * themselves into views. For other types, a regular reference is returned.
   */
  template< class U = T >
  std::enable_if_t< traits::HasMemberFunction_toView< U >, traits::ViewType< U > >
  referenceAsView()
  {
    return m_data->toView();
  }

  /**
   * @copydoc referenceAsView()
   */
  template< class U = T >
  std::enable_if_t< !traits::HasMemberFunction_toView< U >, traits::ViewType< U > >
  referenceAsView()
  {
    return *m_data;
  }

  /**
   * @copydoc referenceAsView()
   */
  template< class U = T >
  std::enable_if_t< traits::HasMemberFunction_toView< U >, traits::ViewTypeConst< U > >
  referenceAsView() const
  {
    return m_data->toViewConst();
  }

  /**
   * @copydoc referenceAsView()
   */
  template< class U = T >
  std::enable_if_t< !traits::HasMemberFunction_toView< U >, traits::ViewTypeConst< U > >
  referenceAsView() const
  {
    return *m_data;
  }

  ///@}

  /**
   * @name Methods to manipulate default value for wrapped object
   */
  ///@{

  /**
   * @copydoc WrapperBase::hasDefaultValue()
   */
  virtual bool
  hasDefaultValue() const final override
  {
    return m_default.has_default_value;
  }

  /**
   * @brief Accessor for m_default.
   * @return reference to const m_default member
   */
  template< typename U = T >
  DefaultValue< T > const &
  getDefaultValueStruct() const
  {
    return m_default;
  }

  /**
   * @brief Accessor for default value.
   * @return reference to const T
   */
  template< typename U = T >
  typename std::enable_if< DefaultValue< U >::has_default_value,
                           typename DefaultValue< U >::value_type const & >::type
  getDefaultValue() const
  {
    return m_default.value;
  }

  /**
   * @brief Setter for default value.
   * @param defaultVal the new default value
   * @return pointer to Wrapper<T>
   */
  template< typename U = T >
  typename std::enable_if< DefaultValue< U >::has_default_value, Wrapper< T > * >::type
  setDefaultValue( typename DefaultValue< U >::value_type const & defaultVal )
  {
    m_default.value = defaultVal;
    return this;
  }

  /**
   * @brief Set and apply for default value.
   * @param defaultVal the new default value
   * @return pointer to Wrapper<T>
   */
  template< typename U = T >
  typename std::enable_if< !traits::is_array< U > && DefaultValue< U >::has_default_value,
                           Wrapper< T > * >::type
  setApplyDefaultValue( typename DefaultValue< U >::value_type const & defaultVal )
  {
    m_default.value = defaultVal;
    *m_data = m_default.value;
    return this;
  }

  /**
   * @brief Set and apply for default value.
   * @param defaultVal the new default value
   * @return pointer to Wrapper<T>
   */
  template< typename U = T >
  typename std::enable_if< traits::is_array< U > && DefaultValue< U >::has_default_value,
                           Wrapper< T > * >::type
  setApplyDefaultValue( typename DefaultValue< U >::value_type const & defaultVal )
  {
    m_default.value = defaultVal;
    m_data->template setValues< serialPolicy >( m_default.value );
    return this;
  }

  /**
   * @copydoc WrapperBase::getDefaultValueString()
   */
  virtual std::string
  getDefaultValueString() const override
  {
    // Find the dimensionality of the wrapper value
    std::string wrapper_type = rtTypes::typeNames( std::type_index( get_typeid() ) );
    integer value_dim = 0;
    if( wrapper_type.find( "array3d" ) != std::string::npos )
    {
      value_dim = 3;
    }
    else if( wrapper_type.find( "array2d" ) != std::string::npos )
    {
      value_dim = 2;
    }
    else if( wrapper_type.find( "array" ) != std::string::npos )
    {
      value_dim = 1;
    }

    // Compose the default string
    std::stringstream ss;

    for( integer ii = 0; ii < value_dim; ++ii )
    {
      ss << "{";
    }

    ss << m_default;

    for( integer ii = 0; ii < value_dim; ++ii )
    {
      ss << "}";
    }

    std::string default_string = ss.str();

    // Tensor types will be space-delimited using the << operator
    // Replace these with commas
    if( wrapper_type.find( "Tensor" ) != std::string::npos )
    {
      std::replace( default_string.begin(), default_string.end(), ' ', ',' );
    }

    return default_string;
  }

  virtual bool
  processInputFile( xmlWrapper::xmlNode const & targetNode ) override
  {
    InputFlags const inputFlag = getInputFlag();
    if( inputFlag >= InputFlags::OPTIONAL )
    {
      if( inputFlag == InputFlags::REQUIRED || !hasDefaultValue() )
      {
        bool const readSuccess =
          xmlWrapper::ReadAttributeAsType( reference(),
                                           getName(),
                                           targetNode,
                                           inputFlag == InputFlags::REQUIRED );
        GEOSX_ERROR_IF(
          !readSuccess,
          "Input variable " + getName() + " is required in " + targetNode.path() +
            ". Available options are: \n" + dumpInputOptions( true ) +
            "\nFor more details, please refer to documentation at: \n" +
            "http://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/"
            "userGuide/Index.html \n" );
      }
      else
      {
        xmlWrapper::ReadAttributeAsType( reference(),
                                         getName(),
                                         targetNode,
                                         getDefaultValueStruct() );
      }

      return true;
    }

    return false;
  }

  /// @endcond DO_NOT_DOCUMENT

  /**
   * @brief Associate the path to this wrapper with the object being held.
   * @note This calls T::setName( std::string ) if it exists and is a no-op if not.
   *       LvArray objects implement this method.
   */
  void
  setName()
  {
    wrapperHelpers::setName( reference(), m_conduitNode.path() );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  void
  addBlueprintField( conduit::Node & fields,
                     std::string const & name,
                     std::string const & topology,
                     std::vector< std::string > const & componentNames = {} ) const override
  {
    wrapperHelpers::addBlueprintField( reference(),
                                       fields,
                                       name,
                                       topology,
                                       componentNames );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  void
  populateMCArray( conduit::Node & node,
                   std::vector< std::string > const & componentNames = {} ) const override
  {
    wrapperHelpers::populateMCArray( reference(), node, componentNames );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  std::unique_ptr< WrapperBase >
  averageOverSecondDim( std::string const & name,
                        Group & group ) const override
  {
    auto ptr = wrapperHelpers::averageOverSecondDim( reference() );
    using U = typename decltype( ptr )::element_type;

    GEOSX_ERROR_IF( ptr == nullptr,
                    "Failed to average over the second dimension of." );

    return std::make_unique< Wrapper< U > >( name, &group, std::move( ptr ) );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  void
  registerToWrite() const override
  {
    m_conduitNode.reset();

    if( getRestartFlags() == RestartFlags::NO_WRITE )
    {
      return;
    }

    move( LvArray::MemorySpace::CPU, false );

    m_conduitNode["__sizedFromParent__"].set( sizedFromParent() );

    wrapperHelpers::pushDataToConduitNode( *m_data, m_conduitNode );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  void
  finishWriting() const override
  {
    m_conduitNode.reset();
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  bool
  loadFromConduit() override
  {
    if( getRestartFlags() != RestartFlags::WRITE_AND_READ )
    {
      m_conduitNode.reset();
      return false;
    }

    setSizedFromParent( m_conduitNode["__sizedFromParent__"].value() );

    wrapperHelpers::pullDataFromConduitNode( *m_data, m_conduitNode );

    m_conduitNode.reset();

    return true;
  }

  /**
   *  @name Overridden setters
   *  Setters that replace (hide) corresponding non-virtual functions in WrapperBase
   *  with the goal of returning Wrapper<T> pointers rather than WrapperBase pointers.
   */
  ///@{

  /*
   * @brief Set whether this wrapper is resized when its parent is resized.
   * @param val an int that is converted into a bool
   * @return a pointer to this wrapper
   */

  /**
   * @copydoc WrapperBase::setSizedFromParent(int)
   */
  Wrapper< T > *
  setSizedFromParent( int val )
  {
    WrapperBase::setSizedFromParent( val );
    return this;
  }

  /**
   * @copydoc WrapperBase::setRestartFlags(RestartFlags)
   */
  Wrapper< T > *
  setRestartFlags( RestartFlags flags )
  {
    WrapperBase::setRestartFlags( flags );
    return this;
  }

  /**
   * @copydoc WrapperBase::setPlotLevel(PlotLevel const)
   */
  Wrapper< T > *
  setPlotLevel( PlotLevel const flag )
  {
    WrapperBase::setPlotLevel( flag );
    return this;
  }

  /**
   * @copydoc WrapperBase::setInputFlag(InputFlags const)
   */
  Wrapper< T > *
  setInputFlag( InputFlags const input )
  {
    WrapperBase::setInputFlag( input );
    return this;
  }

  /**
   * @copydoc WrapperBase::setDescription(string const &)
   */
  Wrapper< T > *
  setDescription( string const & description )
  {
    WrapperBase::setDescription( description );
    return this;
  }

  /**
   * @copydoc WrapperBase::setRegisteringObjects(string const &)
   */
  Wrapper< T > *
  setRegisteringObjects( string const & objectName )
  {
    WrapperBase::setRegisteringObjects( objectName );
    return this;
  }

  ///@}

#if defined( USE_TOTALVIEW_OUTPUT )
  virtual string
  totalviewTypeName() const override
  {
    return LvArray::system::demangle( typeid( Wrapper< T > ).name() );
  }

  virtual int
  setTotalviewDisplay() const override
  {
    //std::cout<<"executing Wrapper::setTotalviewDisplay()"<<std::endl;
    WrapperBase::setTotalviewDisplay();
    TV_ttf_add_row( "m_ownsData", "bool", &m_ownsData );
    TV_ttf_add_row( "m_data", LvArray::system::demangle< T >().c_str(), m_data );
    TV_ttf_add_row( "m_default",
                    LvArray::system::demangle< DefaultValue< T > >().c_str(),
                    &m_default );
    return 0;
  }
//  void tvTemplateInstantiation();
#endif

private:
  /// flag to indicate whether or not this wrapper is responsible for allocation/deallocation of the object at the
  /// address of m_data
  bool m_ownsData;

  /// the object being wrapped by this wrapper
  T * m_data;

  /// the default value of the object being wrapped
  DefaultValue< T > m_default;

  Wrapper() = delete;
};

}  // namespace dataRepository
} /* namespace geosx */

//template< typename T >
//int TV_ttf_display_type( geosx::dataRepository::Wrapper<T> const * wrapper)
//{
//  std::cout<<"Executing "<<wrapper->totalviewTypeName()<<"::TV_ttf_display_type()"<<std::endl;
//  return TV_ttf_format_raw;
//}
//
//template int TV_ttf_display_type( geosx::dataRepository::Wrapper<int> const * wrapper );
//
//template< typename T >
//void geosx::dataRepository::Wrapper<T>::tvTemplateInstantiation()
//{
//  TV_ttf_display_type<T>(this);
//}

#endif /* GEOSX_DATAREPOSITORY_WRAPPER_HPP_ */
