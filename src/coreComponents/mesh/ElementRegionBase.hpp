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

#ifndef GEOSX_MESH_ELEMENTREGIONBASE_HPP
#define GEOSX_MESH_ELEMENTREGIONBASE_HPP

#include "CellElementSubRegion.hpp"
#include "FaceElementSubRegion.hpp"
#include "WellElementSubRegion.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "EmbeddedSurfaceSubRegion.hpp"

namespace geosx
{

class FaceManager;

/**
 * @class ElementRegionBase
 * @brief The ElementRegionBase is the base class to manage the data stored at the element level.
 *
 * The ElementRegion base is the base class for classes such as CellElementRegion, FaceElementRegion,
 * WellElementRegion, and EmbeddedSurfaceRegion.
 */
class ElementRegionBase : public ObjectManagerBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Deleted default constructor.
   */
  ElementRegionBase() = delete;

  /**
   * @brief Main constructor.
   * @param name the name of the element region
   * @param parent the pointer to the parent group
   */
  ElementRegionBase( string const & name, Group * const parent );


  /**
   * @brief Copy constructor.
   * @param init the element region to be copied
   */
  ElementRegionBase( const ElementRegionBase & init );

  /**
   * @brief Default destructor.
   */
  virtual ~ElementRegionBase() override;

  ///@}

  /**
   * @name Generation of the mesh region
   */
  ///@{

  /**
   * @brief Generate mesh.
   * @param cellBlocks cell blocks where the mesh is generated
   */
  virtual void generateMesh( Group & cellBlocks )
  {
    GEOSX_UNUSED_VAR( cellBlocks );
    GEOSX_ERROR( "ElementRegionBase::GenerateMesh() should be overriden if called." );
  }

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Get a collection of the subregions.
   * @return a collection of the subregions
   */
  subGroupMap const & getSubRegions() const
  {
    return getGroup( viewKeyStruct::elementSubRegions() ).getSubGroups();
  }

  /**
   * @brief Get a reference to a subregion.
   * @tparam SUBREGIONTYPE the type that will be used to attempt casting the subregion
   * @tparam KEY_TYPE The type of the key used to lookup the subregion.
   * @param key The key to the subregion.
   * @return A reference to the subregion
   */
  template< typename SUBREGIONTYPE=ElementSubRegionBase, typename KEY_TYPE=void >
  SUBREGIONTYPE const & getSubRegion( KEY_TYPE const & key ) const
  {
    return this->getGroup( viewKeyStruct::elementSubRegions() ).getGroup< SUBREGIONTYPE >( key );
  }

  /**
   * @copydoc getSubRegion( KEY_TYPE const & key ) const
   */
  template< typename SUBREGIONTYPE=ElementSubRegionBase, typename KEY_TYPE=void >
  SUBREGIONTYPE & getSubRegion( KEY_TYPE const & key )
  {
    return this->getGroup( viewKeyStruct::elementSubRegions() ).getGroup< SUBREGIONTYPE >( key );
  }

  /**
   * @brief Get the number of subregions in the region.
   * @return the number of subregions  in the region
   */
  localIndex numSubRegions() const
  {
    return this->getSubRegions().size();
  }

  /**
   * @brief Get the number of elements  in the region
   *        for specific subregion types provided as template arguments.
   * @tparam SUBREGIONTYPE  the first type that will be used in the attempted casting of the subregion
   * @tparam SUBREGIONTYPES a variadic list of types that will be used in the attempted casting of the subregion
   * @return the number of elements contained in the element region
   * @note This function requires that the subRegion types specified
   *       in the variadic template argument can be casted to ElementSubRegionBase
   */
  template< typename SUBREGIONTYPE = ElementSubRegionBase, typename ... SUBREGIONTYPES >
  localIndex getNumberOfElements() const
  {
    localIndex numElem = 0;
    this->forElementSubRegions< SUBREGIONTYPE, SUBREGIONTYPES... >( [&]( Group const & group ) -> void
    {
      numElem += group.size();
    } );
    return numElem;
  }

  /**
   * @copydoc getMaterialList() const
   */
  string_array & getMaterialList() {return m_materialList;}

  /**
   * @brief Get the material list in the element region.
   * @return the material list
   */
  string_array const & getMaterialList() const {return m_materialList;}

  /**
   * @brief Get the name of the constiutive in the element region.
   * @tparam CONSTITUTIVE_TYPE the type of the constitutive model
   * @return the string_array of the constitutive names
   */
  template< typename CONSTITUTIVE_TYPE >
  string_array getConstitutiveNames() const;


  ///@}

  /**
   * @name Functor-based loops over subregions
   */
  ///@{


/**
 * @brief Apply a lambda to all subregions.
 * @param lambda the functor to be applied
 */
  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda ) const
  {
    forElementSubRegions< CellElementSubRegion, FaceElementSubRegion, WellElementSubRegion, EmbeddedSurfaceSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

/**
 * @copydoc forElementSubRegions( LAMBDA && lambda ) const
 */
  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda )
  {
    forElementSubRegions< CellElementSubRegion, FaceElementSubRegion, WellElementSubRegion, EmbeddedSurfaceSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

/**
 * @brief Apply LAMBDA to the subregions with the specific subregion types
 *        listed in the template.
 * @param lambda the functor to be applied
 */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda ) const
  {
    this->getGroup( viewKeyStruct::elementSubRegions() ).forSubGroups< SUBREGIONTYPE, SUBREGIONTYPES... >( std::forward< LAMBDA >( lambda ) );
  }

/**
 * @copydoc forElementSubRegions( LAMBDA && lambda ) const
 */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda )
  {
    this->getGroup( viewKeyStruct::elementSubRegions() ).forSubGroups< SUBREGIONTYPE, SUBREGIONTYPES... >( std::forward< LAMBDA >( lambda ) );
  }

/**
 * @brief Apply LAMBDA to the subregions, loop using subregion indices.
 * @tparam LAMBDA type of functor to call
 * @param lambda the functor to be applied to all subregions
 */
  template< typename LAMBDA >
  void forElementSubRegionsIndex( LAMBDA && lambda ) const
  {
    forElementSubRegionsIndex< CellElementSubRegion, FaceElementSubRegion, WellElementSubRegion, EmbeddedSurfaceSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

/**
 * @copydoc forElementSubRegionsIndex( LAMBDA && lambda ) const
 */
  template< typename LAMBDA >
  void forElementSubRegionsIndex( LAMBDA && lambda )
  {
    forElementSubRegionsIndex< CellElementSubRegion, FaceElementSubRegion, WellElementSubRegion, EmbeddedSurfaceSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

/**
 * @brief Apply LAMBDA to the subregions with the specific subregion types
 *        listed in the template, loop using subregion indices.
 * @tparam SUBREGIONTYPE the first subregion type
 * @tparam SUBREGIONTYPES a variadic list of subregion types
 * @tparam LAMBDA type of functor to call
 * @param lambda the functor to be applied
 */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegionsIndex( LAMBDA && lambda ) const
  {
    for( localIndex esr=0; esr<this->numSubRegions(); ++esr )
    {
      ElementSubRegionBase const & subRegion = this->getSubRegion( esr );
      applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto const & castedSubRegion )
      {
        lambda( esr, castedSubRegion );
      } );
    }
  }

/**
 * @copydoc forElementSubRegionsIndex( LAMBDA && lambda ) const
 */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegionsIndex( LAMBDA && lambda )
  {
    for( localIndex esr=0; esr<this->numSubRegions(); ++esr )
    {
      ElementSubRegionBase & subRegion = this->getSubRegion( esr );
      applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto & castedSubRegion )
      {
        lambda( esr, castedSubRegion );
      } );
    }
  }

  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {
    /// @return String key for the material list
    static constexpr char const * materialListString() { return "materialList"; }
    /// @return String key for the element subregions
    static constexpr char const * elementSubRegions() { return "elementSubRegions"; }
  };

private:

  ElementRegionBase & operator=( const ElementRegionBase & rhs );

  /// List of materials for the element region
  string_array m_materialList;

  /// Name of the numerical method
  string m_numericalMethod;

};



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////


/**
 * @brief Get the names of all constitutive models of a specific type
 *
 * @tparam CONSTITUTIVE_TYPE type of constitutive model
 * return string array with the names of the constitutive models
 */
template< typename CONSTITUTIVE_TYPE >
string_array ElementRegionBase::getConstitutiveNames() const
{
  string_array rval;
  for( string const & matName : m_materialList )
  {
    if( this->getSubRegion( 0 ).getConstitutiveModels().hasGroup< CONSTITUTIVE_TYPE >( matName ) )
    {
      rval.emplace_back( matName );
    }
  }
  return rval;
}

}



#endif /* GEOSX_MESH_ELEMENTREGIONBASE_HPP */
