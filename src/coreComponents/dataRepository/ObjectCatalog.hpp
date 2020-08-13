#ifndef GEOSX_DATAREPOSITORY_OBJECTCATALOG_HPP_
#define GEOSX_DATAREPOSITORY_OBJECTCATALOG_HPP_
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
 * @file ObjectCatalog.hpp
 * The ObjectCatalog acts as a statically initialized factory. It functions in
 * a similar manner to classic virtual factory method, except that it is no
 * maintained list of derived objects that is required to create new objects.
 * Instead, the ``ObjectCatalog`` creates a "catalog" of derived objects using
 * a ``std::unordered_map``.
 * This ``std::unordered_map`` is then statically initialized through the declaration
 * of a
 */

#include "common/Logger.hpp"
#include "LvArray/src/system.hpp"

#include <unordered_map>
#include <string>
#include <iostream>
#include <memory>

#ifndef OBJECTCATALOGVERBOSE
  /**
 * @brief Enables verbose logging of object catalog
 */
  #define OBJECTCATALOGVERBOSE 0
#endif

#ifndef BASEHOLDSCATALOG
  /**
 * @brief Enables storing catalogs in the base class
 */
  #define BASEHOLDSCATALOG 1
#endif

namespace geosx
{
namespace dataRepository
{
/**
 * @brief This class provides the base class/interface for the catalog value objects.
 * @tparam BASETYPE base class of the objects that the factory produces
 * @tparam ARGS variadic template pack to hold the parameters needed for the
 *              constructor of the @p BASETYPE
 */
//START_SPHINX_0
template< typename BASETYPE, typename... ARGS >
class CatalogInterface
{
public:
  /// This is the type that will be used for the catalog. The catalog is actually instantiated in the @p BASETYPE.
  //START_SPHINX_1
  typedef std::unordered_map< std::string,
                              std::unique_ptr< CatalogInterface< BASETYPE, ARGS... > > >
    CatalogType;
  //STOP_SPHINX

  /**
   * @brief Default constructor.
   */
  CatalogInterface()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOSX_LOG( "Calling constructor for CatalogInterface< "
               << LvArray::system::demangle( typeid( BASETYPE ).name() )
               << " , ... >" );
#endif
  }

  /**
   * @brief Default destructor.
   */
  virtual ~CatalogInterface()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOSX_LOG( "Calling destructor for CatalogInterface< "
               << LvArray::system::demangle( typeid( BASETYPE ).name() )
               << " , ... >" );
#endif
  }

  /**
   * @brief Copy constructor.
   */
  explicit CatalogInterface( CatalogInterface const & ) = default;

  /**
   * @brief Move constructor.
   */
  CatalogInterface( CatalogInterface && ) = default;

  /**
   * @brief Copy assignment operator.
   * @return reference to this object
   */
  CatalogInterface &
  operator=( CatalogInterface const & ) = default;

  /**
   * @brief Move assignment operator.
   * @return reference to this object
   */
  CatalogInterface &
  operator=( CatalogInterface && ) = default;

  /**
   * @brief Get the catalog from that is stored in the target base class.
   * @return returns the catalog for this
   */
  static CatalogType &
  GetCatalog()
  {
#if BASEHOLDSCATALOG == 1
    return BASETYPE::GetCatalog();
#else
    static CatalogType catalog;
    return catalog;
#endif
  }

  /**
   * @brief Create a new object that derives from BASETYPE.
   * @param args arguments to the constructor of the target type
   * @return a unique_ptr<BASETYPE> to the newly allocated class.
   */
  virtual std::unique_ptr< BASETYPE >
  Allocate( ARGS... args ) const = 0;

  /**
   * @brief Check if catalog contains a given key
   * @param objectTypeName name of the type tp look up
   * @return @p true if type has been registered with this catalog, @p false otherwise
   */
  static bool
  hasKeyName( std::string const & objectTypeName )
  {
    return GetCatalog().count( objectTypeName );
  }

  /**
   * @brief Static method to create a new object that derives from BASETYPE
   * @param[in] objectTypeName the key to the catalog entry that is able to create
   * the correct type.
   * @param args these are the arguments to the constructor of the target type
   * @return passes a unique_ptr<BASETYPE> to the newly allocated class.
   */
  //START_SPHINX_2
  static std::unique_ptr< BASETYPE >
  Factory( std::string const & objectTypeName,
           ARGS... args )
  {
    return GetCatalog().at( objectTypeName ).get()->Allocate( args... );
  }
  //STOP_SPHINX

  /**
   * @brief Downcast base type reference to derived type
   * @tparam TYPE type to cast to
   * @param object base type reference to object
   * @return reference to the same object, cast to derived type
   *
   * If @p OBJECTCATALOGVERBOSE is enabled, will check that runtime name of the object
   * is the same as catalog name of the derived type. Therefore may fail for objects
   * that have been assigned a different name (e.g. through XML "name" attribute).
   */
  template< typename TYPE >
  static TYPE &
  catalog_cast( BASETYPE & object )
  {
    std::string castedName = TYPE::CatalogName();
    std::string objectName = object.getName();

    if( castedName != objectName )
    {
#if OBJECTCATALOGVERBOSE > 1
      GEOSX_LOG( "Invalid Cast of " << objectName << " to " << castedName );
#endif
    }

    return static_cast< TYPE & >( object );
  }
};

/**
 * @brief Class to hold allocation capability for specific target derived types.
 * @tparam TYPE this is the derived type
 * @tparam BASETYPE this is the base class that TYPE derives from
 * @tparam ARGS constructor arguments
 */
//START_SPHINX_3
template< typename BASETYPE, typename TYPE, typename... ARGS >
class CatalogEntry final : public CatalogInterface< BASETYPE, ARGS... >
{
public:
  /**
   * @brief Default constructor.
   */
  CatalogEntry() :
    CatalogInterface< BASETYPE, ARGS... >()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOSX_LOG( "Calling constructor for CatalogEntry< "
               << LvArray::system::demangle( typeid( TYPE ).name() ) << " , "
               << LvArray::system::demangle( typeid( BASETYPE ).name() )
               << " , ... >" );
#endif
  }

  /**
   * @brief Default destructor.
   */
  ~CatalogEntry() override
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOSX_LOG( "Calling destructor for CatalogEntry< "
               << LvArray::system::demangle( typeid( TYPE ).name() ) << " , "
               << LvArray::system::demangle( typeid( BASETYPE ).name() )
               << " , ... >" );
#endif
  }

  /**
   * @brief Copy constructor.
   * @param source object to copy
   */
  CatalogEntry( CatalogEntry const & source ) :
    CatalogInterface< BASETYPE, ARGS... >( source )
  {}

  /**
   * @brief Move constructor.
   * @param source object to move from
   */
  CatalogEntry( CatalogEntry && source ) :
    CatalogInterface< BASETYPE, ARGS... >( std::move( source ) )
  {}

  /**
   * @brief Copy assignment operator.
   * @param source object to copy
   * @return reference to this object
   */
  CatalogEntry &
  operator=( CatalogEntry const & source )
  {
    CatalogInterface< BASETYPE, ARGS... >::operator=( source );
  }

  /**
   * @brief Move assignment operator.
   * @param source object to move from
   * @return reference to this object
   */
  CatalogEntry &
  operator=( CatalogEntry && source )
  {
    CatalogInterface< BASETYPE, ARGS... >::operator=( std::move( source ) );
  }

  /**
   * @brief Create a new object that derives from BASETYPE.
   * @param args these are the arguments to the constructor of the target type
   * @return a unique_ptr<BASETYPE> to the newly allocated class.
   */
  //START_SPHINX_4
  virtual std::unique_ptr< BASETYPE >
  Allocate( ARGS... args ) const override
  {
#if OBJECTCATALOGVERBOSE > 0
    GEOSX_LOG( "Creating type "
               << LvArray::system::demangle( typeid( TYPE ).name() )
               << " from catalog of "
               << LvArray::system::demangle( typeid( BASETYPE ).name() ) );
#endif
#if( __cplusplus >= 201402L )
    return std::make_unique< TYPE >( args... );
#else
    return std::unique_ptr< BASETYPE >( new TYPE( args... ) );
#endif
  }
  //STOP_SPHINX
};

/**
 * @brief A class to generate the catalog entry.
 *
 * Once created, instances of this class may be destroyed without consequence.
 */
template< typename BASETYPE, typename TYPE, typename... ARGS >
class CatalogEntryConstructor
{
public:
  /**
   * @brief Constructor creates a catalog entry using the key defined by
   * TYPE::CatalogName(), and value of CatalogEntry<TYPE,BASETYPE,ARGS...>.
   */
  CatalogEntryConstructor()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOSX_LOG( "Calling constructor for CatalogEntryConstructor< "
               << LvArray::system::demangle( typeid( TYPE ).name() ) << " , "
               << LvArray::system::demangle( typeid( BASETYPE ).name() )
               << " , ... >" );
#endif

    std::string name = TYPE::CatalogName();
#if( __cplusplus >= 201402L )
    std::unique_ptr< CatalogEntry< BASETYPE, TYPE, ARGS... > > temp =
      std::make_unique< CatalogEntry< BASETYPE, TYPE, ARGS... > >();
#else
    std::unique_ptr< CatalogEntry< BASETYPE, TYPE, ARGS... > > temp =
      std::unique_ptr< CatalogEntry< BASETYPE, TYPE, ARGS... > >(
        new CatalogEntry< BASETYPE, TYPE, ARGS... >() );
#endif
    ( CatalogInterface< BASETYPE, ARGS... >::GetCatalog() )
      .insert( std::move( std::make_pair( name, std::move( temp ) ) ) );

#if OBJECTCATALOGVERBOSE > 0
    GEOSX_LOG( "Registered "
               << LvArray::system::demangle( typeid( BASETYPE ).name() )
               << " catalog component of derived type "
               << LvArray::system::demangle( typeid( TYPE ).name() ) << " where "
               << LvArray::system::demangle( typeid( TYPE ).name() )
               << "::CatalogName() = " << TYPE::CatalogName() );
#endif
  }

  /**
   * @brief Default destructor.
   */
  ~CatalogEntryConstructor()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOSX_LOG( "Calling destructor for CatalogEntryConstructor< "
               << LvArray::system::demangle( typeid( TYPE ).name() ) << " , "
               << LvArray::system::demangle( typeid( BASETYPE ).name() )
               << " , ... >" );
#endif
  }

  /**
   * @brief Deleted copy constructor.
   */
  CatalogEntryConstructor( CatalogEntryConstructor const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  CatalogEntryConstructor( CatalogEntryConstructor && ) = delete;

  /**
   * @brief Deleted copy assignment operator.
   * @return
   */
  CatalogEntryConstructor &
  operator=( CatalogEntryConstructor const & ) = delete;

  /**
   * @brief Deleted move assignment operator.
   * @return
   */
  CatalogEntryConstructor &
  operator=( CatalogEntryConstructor && ) = delete;
};

/**
 * @brief Specialization of @p CatalogInterface for types with no-argument constructors/
 * @tparam BASETYPE base class that contains the catalog
 */
template< typename BASETYPE >
class CatalogInterface< BASETYPE >
{
public:
  /// This is the type that will be used for the catalog. The catalog is actually instantiated in the @p BASETYPE.
  typedef std::unordered_map< std::string, std::unique_ptr< CatalogInterface< BASETYPE > > >
    CatalogType;

  /**
   * @brief Default constructor.
   */
  CatalogInterface()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOSX_LOG( "Calling constructor for CatalogInterface< "
               << LvArray::system::demangle( typeid( BASETYPE ).name() )
               << " , ... >" );
#endif
  }

  /**
   * @brief Default destructor.
   */
  virtual ~CatalogInterface()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOSX_LOG( "Calling destructor for CatalogInterface< "
               << LvArray::system::demangle( typeid( BASETYPE ).name() )
               << " , ... >" );
#endif
  }

  /**
   * @brief Copy constructor.
   */
  explicit CatalogInterface( CatalogInterface const & ) = default;

  /**
   * @brief Move constructor.
   */
  CatalogInterface( CatalogInterface && ) = default;

  /**
   * @brief Copy assignment operator.
   * @return reference to this object
   */
  CatalogInterface &
  operator=( CatalogInterface const & ) = default;

  /**
   * @brief Move assignment operator.
   * @return reference to this object
   */
  CatalogInterface &
  operator=( CatalogInterface && ) = default;

  /**
   * @brief Get the catalog from that is stored in the target base class.
   * @return returns the catalog for this
   */
  static CatalogType &
  GetCatalog()
  {
#if BASEHOLDSCATALOG == 1
    return BASETYPE::GetCatalog();
#else
    static CatalogType catalog;
    return catalog;
#endif
  }

  /**
   * @brief Create a new object that derives from BASETYPE.
   * @return passes a unique_ptr<BASETYPE> to the newly allocated class.
   */
  virtual std::unique_ptr< BASETYPE >
  Allocate() const = 0;

  /**
   * @brief Create a new object that derives from BASETYPE.
   * @param objectTypeName The key to the catalog entry that is able to create the correct type.
   * @return passes a unique_ptr<BASETYPE> to the newly allocated class.
   */
  static std::unique_ptr< BASETYPE >
  Factory( std::string const & objectTypeName )
  {
    CatalogInterface< BASETYPE > const * const entry =
      GetCatalog().at( objectTypeName ).get();
    return entry->Allocate();
  }

  /**
   * @brief Downcast base type reference to derived type
   * @tparam TYPE type to cast to
   * @param object base type reference to object
   * @return reference to the same object, cast to derived type
   *
   * If @p OBJECTCATALOGVERBOSE is enabled, will check that runtime name of the object
   * is the same as catalog name of the derived type. Therefore may fail for objects
   * that have been assigned a different name (e.g. through XML "name" attribute).
   */
  template< typename TYPE >
  static TYPE &
  catalog_cast( BASETYPE & object )
  {
    std::string castedName = TYPE::CatalogName();
    std::string objectName = object.getName();

    if( castedName != objectName )
    {
#if OBJECTCATALOGVERBOSE > 1
      GEOSX_LOG( "Invalid Cast of " << objectName << " to " << castedName );
#endif
    }

    return static_cast< TYPE & >( object );
  }
};

/**
 * @brief Specialization of @p CatalogEntry for types with no-argument constructors.
 * @tparam BASETYPE
 * @tparam TYPE
 */
template< typename BASETYPE, typename TYPE >
class CatalogEntry< BASETYPE, TYPE > final : public CatalogInterface< BASETYPE >
{
public:
  /**
   * @brief Default constructor.
   */
  CatalogEntry() :
    CatalogInterface< BASETYPE >()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOSX_LOG( "Calling constructor for CatalogEntry< "
               << LvArray::system::demangle( typeid( TYPE ).name() ) << " , "
               << LvArray::system::demangle( typeid( BASETYPE ).name() )
               << " , ... >" );
#endif
  }

  /**
   * @brief Default destructor.
   */
  ~CatalogEntry() override
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOSX_LOG( "Calling destructor for CatalogEntry< "
               << LvArray::system::demangle( typeid( TYPE ).name() ) << " , "
               << LvArray::system::demangle( typeid( BASETYPE ).name() )
               << " , ... >" );
#endif
  }

  /**
   * @brief Copy constructor.
   * @param source object to copy
   */
  CatalogEntry( CatalogEntry const & source ) :
    CatalogInterface< BASETYPE >( source )
  {}

  /**
   * @brief Move constructor.
   * @param source object to move from
   */
  CatalogEntry( CatalogEntry && source ) :
    CatalogInterface< BASETYPE >( std::move( source ) )
  {}

  /**
   * @brief Copy assignment operator.
   * @param source object to copy
   * @return reference to this object
   */
  CatalogEntry &
  operator=( CatalogEntry const & source )
  {
    CatalogInterface< BASETYPE >::operator=( source );
  }

  /**
   * @brief Move assignment operator.
   * @param source object to move from
   * @return reference to this object
   */
  CatalogEntry &
  operator=( CatalogEntry && source )
  {
    CatalogInterface< BASETYPE >::operator=( std::move( source ) );
  }

  /**
   * @brief Create a new instance of @p TYPE.
   * @return a unique_ptr<BASETYPE> that owns the new instance
   */
  virtual std::unique_ptr< BASETYPE >
  Allocate() const override
  {
#if OBJECTCATALOGVERBOSE > 0
    GEOSX_LOG( "Creating type "
               << LvArray::system::demangle( typeid( TYPE ).name() )
               << " from catalog of "
               << LvArray::system::demangle( typeid( BASETYPE ).name() ) );
#endif
#if( __cplusplus >= 201402L )
    return std::make_unique< TYPE >();
#else
    return std::unique_ptr< BASETYPE >( new TYPE() );
#endif
  }
};

/**
 * @brief A specialization of @p CatalogEntryConstructor for types with no-argument constructors.
 */
template< typename BASETYPE, typename TYPE >
class CatalogEntryConstructor< BASETYPE, TYPE >
{
public:
  /**
   * @brief Default constructor.
   */
  CatalogEntryConstructor()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOSX_LOG( "Calling constructor for CatalogEntryConstructor< "
               << LvArray::system::demangle( typeid( TYPE ).name() ) << " , "
               << LvArray::system::demangle( typeid( BASETYPE ).name() )
               << " , ... >" );
#endif

    std::string name = TYPE::CatalogName();
#if( __cplusplus >= 201402L )
    std::unique_ptr< CatalogEntry< BASETYPE, TYPE > > temp =
      std::make_unique< CatalogEntry< BASETYPE, TYPE > >();
#else
    std::unique_ptr< CatalogEntry< BASETYPE, TYPE > > temp =
      std::unique_ptr< CatalogEntry< BASETYPE, TYPE > >(
        new CatalogEntry< BASETYPE, TYPE >() );
#endif
    ( CatalogInterface< BASETYPE >::GetCatalog() )
      .insert( std::move( std::make_pair( name, std::move( temp ) ) ) );

#if OBJECTCATALOGVERBOSE > 0
    GEOSX_LOG( "Registered "
               << LvArray::system::demangle( typeid( BASETYPE ).name() )
               << " catalog component of derived type "
               << LvArray::system::demangle( typeid( TYPE ).name() ) << " where "
               << LvArray::system::demangle( typeid( TYPE ).name() )
               << "::CatalogName() = " << TYPE::CatalogName() );
#endif
  }

  /**
   * @brief Default destuctor.
   */
  ~CatalogEntryConstructor()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOSX_LOG( "Calling destructor for CatalogEntryConstructor< "
               << LvArray::system::demangle( typeid( TYPE ).name() ) << " , "
               << LvArray::system::demangle( typeid( BASETYPE ).name() )
               << " , ... >" );
#endif
  }

  /**
   * @brief Deleted copy constructor.
   */
  CatalogEntryConstructor( CatalogEntryConstructor const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  CatalogEntryConstructor( CatalogEntryConstructor && ) = delete;

  /**
   * @brief Deleted copy assignment operator.
   * @return
   */
  CatalogEntryConstructor &
  operator=( CatalogEntryConstructor const & ) = delete;

  /**
   * @brief Deleted move assignment operator.
   * @return
   */
  CatalogEntryConstructor &
  operator=( CatalogEntryConstructor && ) = delete;
};

}  // namespace dataRepository
}  // namespace geosx

/**
 * @brief Object catalog class registration macro.
 *
 * Macro that takes in the base class of the catalog, the derived class, and the
 * argument types for the constructor of
 * the derived class/base class, and create an object of type
 * CatalogEntryConstructor<ClassName,BaseType,__VA_ARGS__> in
 * an anonymous namespace. This should be called from the source file for the
 * derived class, which will result in the
 * generation of a CatalogEntry<BaseType,ClassName,...> prior to main().
 */
#define REGISTER_CATALOG_ENTRY( BaseType, DerivedType, ... )                           \
  namespace                                                                            \
  {                                                                                    \
  geosx::dataRepository::CatalogEntryConstructor< BaseType, DerivedType, __VA_ARGS__ > \
    catEntry_##DerivedType;                                                            \
  }

/**
 * @brief Same as REGISTER_CATALOG_ENTRY, but for classes with no-argument constructors.
 */
#define REGISTER_CATALOG_ENTRY0( BaseType, DerivedType )                  \
  namespace                                                               \
  {                                                                       \
  geosx::dataRepository::CatalogEntryConstructor< BaseType, DerivedType > \
    catEntry_##DerivedType;                                               \
  }

#endif /* GEOSX_DATAREPOSITORY_OBJECTCATALOG_HPP_ */
