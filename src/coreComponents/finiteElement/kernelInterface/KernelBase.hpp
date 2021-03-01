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
 * @file KernelBase.hpp
 */

#ifndef GEOSX_FINITEELEMENT_KERNELBASE_HPP_
#define GEOSX_FINITEELEMENT_KERNELBASE_HPP_

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/MeshLevel.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "LvArray/src/jitti/Cache.hpp"

#if defined(__APPLE__)
/// Use camp::tuple to hold constructor params.
#define CONSTRUCTOR_PARAM_OPTION 2
#else
/// Use std::tuple to hold constructor params.
#define CONSTRUCTOR_PARAM_OPTION 1
#endif

#if CONSTRUCTOR_PARAM_OPTION==1
namespace std
{
namespace detail
{
/**
 * @brief Implementation of std::make_from_tuple()
 * @tparam T
 * @tparam Tuple
 * @tparam I
 * @tparam T
 * @tparam Tuple
 * @tparam I
 * @param t
 * @param
 * @return
 */
template< class T, class Tuple, std::size_t... I >
constexpr T make_from_tuple_impl( Tuple && t, std::index_sequence< I... > )
{
  return T( std::get< I >( std::forward< Tuple >( t ))... );
}
} // namespace detail

/**
 * @brief Implementation of std::make_from_tuple()
 * @tparam T
 * @tparam Tuple
 * @tparam T
 * @tparam Tuple
 * @param t
 * @return
 */
template< class T, class Tuple >
constexpr T make_from_tuple( Tuple && t )
{
  return detail::make_from_tuple_impl< T >( std::forward< Tuple >( t ),
                                            std::make_index_sequence< std::tuple_size< std::remove_reference_t< Tuple > >::value >{} );
}

}
#elif CONSTRUCTOR_PARAM_OPTION==2
  #include "camp/camp.hpp"
#endif

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
namespace finiteElement
{

/**
 * @class KernelBase
 * @brief Define the base interface for finite element kernels.
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 * @tparam CONSTITUTIVE_TYPE The type of constitutive model present in the
 *                           subregion.
 * @tparam NUM_TEST_SUPPORT_POINTS_PER_ELEM The number of test space support
 *                                          points per element.
 * @tparam NUM_TRIAL_SUPPORT_POINTS_PER_ELEM The number of trial space support
 *                                           points per element.
 * @tparam NUM_DOF_PER_TEST_SP The number of DOF per test support point.
 * @tparam NUM_DOF_PER_TRIAL_SP The number of DOF per trial support point.
 *
 * ### General KernelBase Description
 *
 * KernelBase defines an interface for implementing finite element kernels
 * that will be callable by the family of kernel launch functions. Specific
 * physics kernels may or may not derive from KernelBase, but must follow
 * the same interface in order to be callable from the generic launching
 * functions.
 *
 * The template parameters of KernelBase should be duplicated as part of the
 * interface, EXCEPT for @p NUM_DOF_PER_TEST_SP and @p NUM_DOF_PER_TRIAL_SP.
 * These values should be set internally by the physics solver since each
 * physics discretization will have a constant intrinsic value for these
 * quantities. For example, when solving or the heat equation with scalar
 * temperature as the primary variable at the support point, these will have
 * a value of 1. In contrast, when solving a solid mechanics problem, with
 * vector displacement as the primary variable at the support point, these
 * will have a value of 3. Note that the interface provided by
 * geosx::finiteElement::RegionBasedKernelApplication will construct a
 * kernel assuming only the first 4 template arguments.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE,
          int NUM_DOF_PER_TEST_SP,
          int NUM_DOF_PER_TRIAL_SP >
class KernelBase
{
public:
  /// Compile time value for the number of test function support points per
  /// element.
  static constexpr int numTestSupportPointsPerElem  = FE_TYPE::numNodes;

  /// Compile time value for the number of trial function support points per
  /// element.
  static constexpr int numTrialSupportPointsPerElem = FE_TYPE::numNodes;

  /// Compile time value for the number of degrees of freedom per test function
  /// support point.
  static constexpr int numDofPerTestSupportPoint    = NUM_DOF_PER_TEST_SP;

  /// Compile time value for the number of degrees of freedom per trial
  /// function support point.
  static constexpr int numDofPerTrialSupportPoint   = NUM_DOF_PER_TRIAL_SP;

  /// Compile time value for the number of quadrature points per element.
  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

  /**
   * @brief Constructor
   * @param elementSubRegion Reference to the SUBREGION_TYPE(class template
   *                         parameter) object.
   * @param inputConstitutiveType The constitutive object.
   * @param finiteElementSpace Placeholder for the finite element space object,
   *                           which currently doesn't do much.
   */
  KernelBase( SUBREGION_TYPE const & elementSubRegion,
              FE_TYPE const & finiteElementSpace,
              CONSTITUTIVE_TYPE & inputConstitutiveType ):
    m_elemsToNodes( elementSubRegion.nodeList().toViewConst() ),
    m_elemGhostRank( elementSubRegion.ghostRank() ),
    m_constitutiveUpdate( inputConstitutiveType.createKernelUpdates() ),
    m_finiteElementSpace( finiteElementSpace )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables allocated on the stack.
   *
   * ### ImplicitKernelBase::StackVariables Description
   *
   * Contains variables that will be allocated on the stack of the main kernel.
   * This will typically consist of local arrays to hold data mapped from the
   * global data arrays, and/or local storage for the residual and jacobian
   * contributions.
   */
  struct StackVariables
  {};

  /**
   * @brief Performs the setup phase for the kernel.
   * @tparam STACK_VARIABLE_TYPE The type of StackVariable that holds the stack
   *                             variables. This is most likely a defined in a
   *                             type that derives from KernelBase.
   * @param k The element index.
   * @param stack The StackVariable object that hold the stack variables.
   *
   * ### KernelBase::setup() Description
   *
   * The operations typically found in setup are thing such as the collection
   * of global data into local stack storage.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( stack );
  }

  /**
   * @brief Performs a state update at a quadrature point.
   * @tparam STACK_VARIABLE_TYPE The type of StackVariable that holds the stack
   *                             variables. This is most likely a defined in a
   *                             type that derives from KernelBase.
   * @param k The element index.
   * @param q The quadrature point index.
   * @param stack The StackVariable object that hold the stack variables.
   *
   * ### KernelBase::quadraturePointKernel() Description
   *
   * The operations found here are the mapping from the support points to the
   * quadrature point, calculation of gradients, etc. From this data the
   * state of the constitutive model is updated if required by the physics
   * package.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( q );
    GEOSX_UNUSED_VAR( stack );
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @tparam STACK_VARIABLE_TYPE The type of StackVariable that holds the stack
   *                             variables. This is most likely a defined in a
   *                             type that derives from KernelBase.
   * @param k The element index.
   * @param stack The StackVariable object that hold the stack variables.
   * @return The maximum contribution to the residual.
   *
   * ### KernelBase::complete() Description
   *
   * The operations typically found in complete are the mapping of the local
   * Jacobian and Residual into the global Jacobian and Residual.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( stack );
    return 0;
  }


  /**
   * @brief Kernel Launcher.
   * @tparam POLICY The RAJA policy to use for the launch.
   * @tparam NUM_QUADRATURE_POINTS The number of quadrature points per element.
   * @tparam KERNEL_TYPE The type of Kernel to execute.
   * @param numElems The number of elements to process in this launch.
   * @param kernelComponent The instantiation of KERNEL_TYPE to execute.
   * @return The maximum residual contribution.
   *
   * This is a generic launching function for all of the finite element kernels
   * that follow the interface set by KernelBase.
   */
  //START_kernelLauncher
  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;

    // Define a RAJA reduction variable to get the maximum residual contribution.
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxResidual( 0 );

    forAll< POLICY >( numElems,
                      [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );
      for( integer q=0; q<numQuadraturePointsPerElem; ++q )
      {
        kernelComponent.quadraturePointKernel( k, q, stack );
      }
      maxResidual.max( kernelComponent.complete( k, stack ) );
    } );
    return maxResidual.get();
  }
  //END_kernelLauncher

protected:
  /// The element to nodes map.
  traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > const m_elemsToNodes;

  /// The element ghost rank array.
  arrayView1d< integer const > const m_elemGhostRank;

  /// The constitutive update object used to update the constitutive state,
  /// and extract constitutive data.
  typename CONSTITUTIVE_TYPE::KernelWrapper const m_constitutiveUpdate;

  /// The finite element space/discretization object for the element type in
  /// the SUBREGION_TYPE.
  FE_TYPE const & m_finiteElementSpace;
};


//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

jitti::CompilationInfo getCompilationInfo();

template< typename POLICY,
          typename SUB_REGION,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE,
          template< typename,
                    typename,
                    typename > class KERNEL_TEMPLATE,
          typename KERNEL_CONSTRUCTOR_PARAMS >
real64 fooBar( NodeManager & nodeManager,
               EdgeManager & edgeManager,
               FaceManager & faceManager,
               KERNEL_CONSTRUCTOR_PARAMS const & kernelConstructorParamsTuple,
               ElementSubRegionBase & subRegionBase,
               localIndex const numElems,
               constitutive::ConstitutiveBase & constitutiveBase,
               FiniteElementBase const & feBase )
{
  // Define an alias for the kernel type for easy use.
  using KERNEL_TYPE = KERNEL_TEMPLATE< SUB_REGION, CONSTITUTIVE_TYPE, FE_TYPE >;

  SUB_REGION & elementSubRegion = dynamicCast< SUB_REGION & >( subRegionBase );
  FE_TYPE const & finiteElement = dynamicCast< FE_TYPE const & >( feBase );
  CONSTITUTIVE_TYPE & castedConstitutiveRelation = dynamicCast< CONSTITUTIVE_TYPE & >( constitutiveBase );

  // 1) Combine the tuple containing the physics kernel specific constructor parameters with
  // the parameters common to all physics kernels that use this interface,
  // 2) Instantiate the kernel.
  // note: have two options, using std::tuple and camp::tuple. Due to a bug in the OSX
  // implementation of std::tuple_cat, we must use camp on OSX. In the future, we should
  // only use one option...most likely camp since we can easily fix bugs.
#if CONSTRUCTOR_PARAM_OPTION==1
  auto temp = std::forward_as_tuple( nodeManager,
                                     edgeManager,
                                     faceManager,
                                     elementSubRegion,
                                     finiteElement,
                                     castedConstitutiveRelation );

  auto fullKernelComponentConstructorArgs = std::tuple_cat( temp,
                                                            kernelConstructorParamsTuple );

  KERNEL_TYPE kernelComponent = std::make_from_tuple< KERNEL_TYPE >( fullKernelComponentConstructorArgs );

#elif CONSTRUCTOR_PARAM_OPTION==2
  auto temp = camp::forward_as_tuple( nodeManager,
                                      edgeManager,
                                      faceManager,
                                      elementSubRegion,
                                      finiteElement,
                                      castedConstitutiveRelation );
  auto fullKernelComponentConstructorArgs = camp::tuple_cat_pair( temp,
                                                                  kernelConstructorParamsTuple );
  KERNEL_TYPE kernelComponent = camp::make_from_tuple< KERNEL_TYPE >( fullKernelComponentConstructorArgs );

#endif

  // Call the kernelLaunch function, and store the maximum contribution to the residual.
  return KERNEL_TYPE::template kernelLaunch< POLICY, KERNEL_TYPE >( numElems, kernelComponent );
}

//START_regionBasedKernelApplication
/**
 * @brief Performs a loop over specific regions (by type and name) and calls
 *        a kernel launch on the subregions with compile time knowledge of
 *        sub-loop bounds such as number of nodes and quadrature points per
 *        element.
 * @tparam POLICY The RAJA launch policy to pass to the kernel launch.
 * @tparam CONSTITUTIVE_BASE The common base class for constitutive
 *                           pass-thru/dispatch which gives the kernel launch
 *                           compile time knowledge of the constitutive model.
 *                           This is achieved through a call to the
 *                           ConstitutivePassThru function which
 *                           should have a specialization for CONSTITUTIVE_BASE
 *                           implemented in order to perform the compile time
 *                           dispatch.
 * @tparam REGION_TYPE The type of region to loop over. TODO make this a
 *                     parameter pack?
 * @tparam KERNEL_TEMPLATE The type of template for the physics kernel, which
 *                         conforms to the interface specified by KernelBase.
 * @tparam KERNEL_CONSTRUCTOR_PARAMS The template parameter pack to hold the
 *                                   parameter pack (i.e. custom) parameters
 *                                   that are sent to the constructor for the
 *                                   @p KERNEL_TEMPLATE.
 * @param mesh The MeshLevel object.
 * @param targetRegions The names of the target regions(of type @p REGION_TYPE)
 *                      to apply the @p KERNEL_TEMPLATE.
 * @param finiteElementName The name of the finite element.
 * @param constitutiveNames The names of the constitutive models present in the
 *                          Region.
 * @param kernelConstructorParams The parameter list for corresponding to the
 *                                parameter @p KERNEL_CONSTRUCTOR_PARAMS that
 *                                are passed to the @p KERNEL_TEMPLATE
 *                                constructor.
 * @return The maximum contribution to the residual, which may be used to scale
 *         the residual.
 *
 * Loops over all regions Applies/Launches a kernel specified by the @p KERNEL_TEMPLATE through
 * #::geosx::finiteElement::KernelBase::kernelLaunch().
 */
template< typename POLICY,
          typename CONSTITUTIVE_BASE,
          typename REGION_TYPE,
          template< typename SUBREGION_TYPE,
                    typename CONSTITUTIVE_TYPE,
                    typename FE_TYPE > class KERNEL_TEMPLATE,
          typename ... KERNEL_CONSTRUCTOR_PARAMS >
static
real64 regionBasedKernelApplication( MeshLevel & mesh,
                                     arrayView1d< string const > const & targetRegions,
                                     string const & finiteElementName,
                                     arrayView1d< string const > const & constitutiveNames,
                                     string const & kernelName,
                                     KERNEL_CONSTRUCTOR_PARAMS && ... kernelConstructorParams )
{
  GEOSX_MARK_FUNCTION;

  // Create a tuple that contains the kernelConstructorParams, as the lambda does not properly catch the parameter pack
  // until c++20
#if CONSTRUCTOR_PARAM_OPTION==1
  std::tuple< KERNEL_CONSTRUCTOR_PARAMS &... > kernelConstructorParamsTuple = std::forward_as_tuple( kernelConstructorParams ... );
#elif CONSTRUCTOR_PARAM_OPTION==2
  camp::tuple< KERNEL_CONSTRUCTOR_PARAMS &... > kernelConstructorParamsTuple = camp::forward_as_tuple( kernelConstructorParams ... );
#endif

  using JIT_FUNCTION_TYPE = real64 (*)( NodeManager &,
                                        EdgeManager &,
                                        FaceManager &,
                                        decltype( kernelConstructorParamsTuple ) const &,
                                        ElementSubRegionBase &,
                                        localIndex const,
                                        constitutive::ConstitutiveBase &,
                                        FiniteElementBase const & );

  jitti::CompilationInfo info = getCompilationInfo();
  static jitti::Cache< JIT_FUNCTION_TYPE > s_jitCache( info.compilationTime, "./lib/jitti" );

  // save the maximum residual contribution for scaling residuals for convergence criteria.
  real64 maxResidualContribution = 0;

  NodeManager & nodeManager = mesh.getNodeManager();
  EdgeManager & edgeManager = mesh.getEdgeManager();
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager & elementRegionManager = mesh.getElemManager();

  // Loop over all sub-regions in regiongs of type REGION_TYPE, that are listed in the targetRegions array.
  elementRegionManager.forElementSubRegions< REGION_TYPE >( targetRegions,
                                                            [&constitutiveNames,
                                                             &maxResidualContribution,
                                                             &nodeManager,
                                                             &edgeManager,
                                                             &faceManager,
                                                             &kernelConstructorParamsTuple,
                                                             &finiteElementName,
                                                             &kernelName,
                                                             &info]
                                                              ( localIndex const targetRegionIndex, auto & elementSubRegion )
  {
    localIndex const numElems = elementSubRegion.size();

    // Create an alias for the type of subregion we are in, which is now known at compile time.
    typedef TYPEOFREF( elementSubRegion ) SUBREGIONTYPE;

    // Get the constitutive model...and allocate a null constitutive model if required.
    constitutive::ConstitutiveBase * constitutiveRelation = nullptr;
    constitutive::NullModel * nullConstitutiveModel = nullptr;
    if( targetRegionIndex <= constitutiveNames.size()-1 )
    {
      constitutiveRelation = &elementSubRegion.template getConstitutiveModel( constitutiveNames[targetRegionIndex] );
    }
    else
    {
      nullConstitutiveModel = &elementSubRegion.template registerGroup< constitutive::NullModel >( "nullModelGroup" );
      constitutiveRelation = nullConstitutiveModel;
    }

    FiniteElementBase &
    subRegionFE = elementSubRegion.template getReference< FiniteElementBase >( finiteElementName );

    info.templateParams = LvArray::system::demangleType< POLICY >() + ", " +
                          LvArray::system::demangleType< SUBREGIONTYPE >() + ", " +
                          LvArray::system::demangleType( *constitutiveRelation ) + ", " +
                          LvArray::system::demangleType( subRegionFE ) + ", " +
                          kernelName + ", " +
                          LvArray::system::demangleType( kernelConstructorParamsTuple );

    auto const & jitFunction = s_jitCache.getOrLoadOrCompile( info );

    maxResidualContribution =
      std::max( maxResidualContribution,
                jitFunction( nodeManager,
                             edgeManager,
                             faceManager,
                             kernelConstructorParamsTuple,
                             elementSubRegion,
                             numElems,
                             *constitutiveRelation,
                             subRegionFE ) );

    // Remove the null constitutive model (not required, but cleaner)
    if( nullConstitutiveModel )
    {
      elementSubRegion.deregisterGroup( "nullModelGroup" );
    }

  } );

  return maxResidualContribution;
}
//END_regionBasedKernelApplication

} // namespace finiteElement
} // namespace geosx

#include "physicsSolvers/solidMechanics/SolidMechanicsSmallStrainExplicitNewmarkKernel.hpp"

#endif /* GEOSX_FINITEELEMENT_KERNELBASE_HPP_ */
