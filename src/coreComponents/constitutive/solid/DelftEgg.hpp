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
 *  @file DelftEgg.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_DELFTEGG_HPP
#define GEOSX_CONSTITUTIVE_SOLID_DELFTEGG_HPP

#include "SolidBase.hpp"
#include "linearAlgebra/interfaces/BlasLapackLA.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class DelftEggUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class DelftEggUpdates : public SolidBaseUpdates
{
public:
  /**
   * @brief Constructor
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] stress The ArrayView holding the stress data for each quadrature point.
   */
  DelftEggUpdates( arrayView1d< real64 const > const & shearModulus,
                   arrayView2d< real64 const > const & refPressure,
                   arrayView2d< real64 const > const & refStrainVol,
                   arrayView1d< real64 const > const & recompressionIndex,
                   arrayView1d< real64 const > const & virginCompressionIndex,
                   arrayView1d< real64 const > const & cslSlope,
                   arrayView1d< real64 const > const & shapeParameter,
                   arrayView2d< real64 > const & newPreConsolidationPressure,
                   arrayView2d< real64 > const & oldPreConsolidationPressure,
                   arrayView3d< real64, solid::STRESS_USD > const & newStress,
                   arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                   arrayView3d< real64, solid::STRESS_USD > const & stress):
    SolidBaseUpdates( stress ),
    m_shearModulus( shearModulus ),
    m_refPressure( refPressure ),
    m_refStrainVol( refStrainVol ),
    m_recompressionIndex( recompressionIndex ),
    m_virginCompressionIndex( virginCompressionIndex ),
    m_cslSlope( cslSlope ),
    m_shapeParameter( shapeParameter ),
    m_oldPreConsolidationPressure( oldPreConsolidationPressure ),
    m_newPreConsolidationPressure( newPreConsolidationPressure ),
    m_newStress( newStress ),
    m_oldStress( oldStress )
  {}

  /// Default copy constructor
  DelftEggUpdates( DelftEggUpdates const & ) = default;

  /// Default move constructor
  DelftEggUpdates( DelftEggUpdates && ) = default;

  /// Deleted default constructor
  DelftEggUpdates() = delete;

  /// Deleted copy assignment operator
  DelftEggUpdates & operator=( DelftEggUpdates const & ) = delete;

  /// Deleted move assignment operator
  DelftEggUpdates & operator=( DelftEggUpdates && ) =  delete;

  GEOSX_HOST_DEVICE
  virtual void SmallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( & strainIncrement )[6],
                                  real64 ( & stress )[6],
                                  real64 ( & stiffness )[6][6] ) override final;
  
  
  GEOSX_HOST_DEVICE
  virtual void SaveConvergedState( localIndex const k,
                                   localIndex const q ) override final;
    
private:

  /// A reference to the ArrayView holding the shear modulus for each element.
  arrayView1d< real64 const > const m_shearModulus;
  
  /// A reference to the ArrayView holding the reference pressure (p0) for each element.
  arrayView2d< real64 const > const m_refPressure;

    /// A reference to the ArrayView holding the shear modulus for each element.
  arrayView2d< real64 const > const m_refStrainVol;

  /// A reference to the ArrayView holding the recompression index for each element.
  arrayView1d< real64 const > const m_recompressionIndex;
  
  /// A reference to the ArrayView holding the virgin compression index for each element.
  arrayView1d< real64 const > const m_virginCompressionIndex;
  
  /// A reference to the ArrayView holding the slope of the critical state line for each element.
  arrayView1d< real64 const > const m_cslSlope;
  
  /// A reference to the ArrayView holding the shape parameter for each integration point
  arrayView1d< real64 const > const m_shapeParameter;

  /// A reference to the ArrayView holding the old preconsolidation pressure for each integration point
  arrayView2d< real64 > const m_oldPreConsolidationPressure;

  /// A reference to the ArrayView holding the new preconsolidation pressure for each integration point
  arrayView2d< real64 > const m_newPreConsolidationPressure;
  
  /// A reference to the ArrayView holding the new stress for each integration point
  arrayView3d< real64, solid::STRESS_USD > const m_newStress;
  
  /// A reference to the ArrayView holding the old stress for each integration point
  arrayView3d< real64, solid::STRESS_USD > const m_oldStress;
};

///////////////////////// new proposal /////////////////////////////////////////

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DelftEggUpdates::SmallStrainUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const ( & strainIncrement )[6],
                                         real64 ( & stress )[6],
                                         real64 ( & stiffness )[6][6] )
{
  real64 const oldPc  = m_oldPreConsolidationPressure[k][q]; //pre-consolidation pressure
  real64 const mu     = m_shearModulus[k];
  real64 const p0     = m_refPressure[k][q];

  real64 const eps_v0 = m_refStrainVol[k][q];
  real64 const M      = m_cslSlope[k];
  real64 const Cr     = m_recompressionIndex[k];
  real64 const Cc     = m_virginCompressionIndex[k];
  real64 const alpha  = m_shapeParameter[k];

  real64        pc    = oldPc; 
  real64 bulkModulus  = -p0/Cr;

  real64 sqrt23 = std::sqrt(2./3.);
  array1d< real64 > identity(6);
    
  for(localIndex i=0; i<3; ++i)
  {
    identity[i] = 1.0;
    identity[i+3] = 0.0;
  }

  // two-invariant decomposition of old stress in P-Q space (mean & deviatoric stress)

  real64 oldP = 0.0;
  real64 oldQ = 0.0;
  for(localIndex i=0; i<6; ++i)
  {
    stress[i] = m_oldStress[k][q][i];
  }

  for(localIndex i=0; i<3; ++i)
  {
    oldP += stress[i];
  }
  oldP /= 3.0;

  if (std::abs(oldP) < 1e-15)
  {
    oldP=p0;
    stress[0] = p0; stress[1] = p0; stress[2] = p0;
  }

  array1d< real64 > oldDeviator(6);  // array allocation
  for(localIndex i=0; i<3; ++i)
  {
    oldDeviator[i] = stress[i]-oldP;
    oldDeviator[i+3] = stress[i+3];
  }
  
  for(localIndex i=0; i<3; ++i)
  {
    oldQ += oldDeviator[i]*oldDeviator[i];
    oldQ += 2. * oldDeviator[i+3]*oldDeviator[i+3];
  }

  oldQ = std::sqrt(oldQ); //+ 1e-15; // ;

    for(localIndex i=0; i<6; ++i)
  {
    oldDeviator[i] /= (oldQ+ 1e-15); // normalized deviatoric direction, "nhat" from previous step
                                    // also perturbed to avoid divide by zero when Q=0
  }
  
  oldQ *= std::sqrt(3./2.);

  // Recover elastic strains from the previous step, based on stress from the previous step
  // [Note: in order to minimize data transfer, we are not storing and passing elastic strains] 

  real64 oldElasticStrainVol = std::log(oldP/p0) * Cr * (-1.0) + eps_v0 ;
  real64 oldElasticStrainDev = oldQ/3./mu;  

  // Now recover the old strain tensor from the strain invariants. 
  // Note that we need the deviatoric direction (n-hat) from the previous step.

  array1d< real64 > strainElasticTrial(6);
  array1d< real64 > oldStrainElastic(6);
  real64 strainElasticTrialVol = 0.0 ;
  
   for(localIndex i=0; i<6; ++i)
  {
    oldStrainElastic[i] = oldDeviator[i] * sqrt23 * oldElasticStrainDev + 1./3. * oldElasticStrainVol * identity[i];
    strainElasticTrial[i] = oldStrainElastic[i] + strainIncrement[i];
  }

     for(localIndex i=0; i<3; ++i)
  {
    oldStrainElastic[i+3] *= 2.0; // Voigt form
    strainElasticTrialVol += strainElasticTrial[i];
  }

  // elastic predictor
  // newP= oldP * exp(-1/Cr* strainIncrementVol)
  // newQ = oldQ + 3 * mu * strainIncrementDev

  real64 eps_v_trial = strainElasticTrialVol; 
  real64 temp = eps_v_trial/3.;
  array1d< real64 > deviator(6);  // array allocation
  for(localIndex i=0; i<3; ++i)
  {
     deviator[i] = strainElasticTrial[i]-temp;
     deviator[i+3] = strainElasticTrial[i+3];
   }

   real64 eps_s_trial = 0.0; //deviatoric part of strain increment
   for(localIndex i=0; i<3; ++i)
   {
     eps_s_trial += deviator[i]*deviator[i];
     eps_s_trial += deviator[i+3]*deviator[i+3]/2. ; //because of Voigt form of strain increment vector
   }

     eps_s_trial  = std::sqrt(eps_s_trial);// + 1e-15; // 

    for(localIndex i=0; i<6; ++i)
  {
    deviator[i] /= (eps_s_trial + 1e-15); // normalized deviatoric direction, "nhat" from previous step
                                            // also perturbed to avoid divide by zero when Q=0;
  }
  // Calculate the normalized deviatoric direction, "nhat"
  
    for(localIndex i=0; i<3; ++i)
  {
    deviator[i+3] /= 2. ; // because of Voigt form 
  }
   eps_s_trial *= std::sqrt(2./3.);

  real64 trialP = p0 * std::exp(-1./Cr* (eps_v_trial-eps_v0));
  real64 trialQ = 3. * mu * eps_s_trial;

  for(localIndex i=0; i<6; ++i)
  {
    stress[i] = trialP * identity[i] + trialQ * sqrt23 *deviator[i];
  }

  // set stiffness to elastic predictor
  
   bulkModulus = -trialP/Cr ; 
  real64 lame = bulkModulus - 2./3. * mu;

  for(localIndex i=0; i<6; ++i)
  {
    for(localIndex j=0; j<6; ++j)
    {
      stiffness[i][j] = 0;
    }
  }
  
  stiffness[0][0] = lame + 2.*mu;
  stiffness[0][1] = lame;
  stiffness[0][2] = lame;

  stiffness[1][0] = lame;
  stiffness[1][1] = lame + 2.*mu;
  stiffness[1][2] = lame;

  stiffness[2][0] = lame;
  stiffness[2][1] = lame;
  stiffness[2][2] = lame + 2.*mu;

  stiffness[3][3] = mu;
  stiffness[4][4] = mu;
  stiffness[5][5] = mu;

  // check yield function F <= 0
  
  real64 yield = trialQ*trialQ/(M*M)- alpha*alpha*trialP *(2*alpha/(alpha+1)*pc-trialP)+alpha*alpha*(alpha-1)/(alpha+1)* pc*pc;
  
  if(yield > 1e-9) // plasticity branch
  {
    std::cout << "plastic " <<  "\n " << std::endl;
    // the return mapping can in general be written as a newton iteration.
    // here we have a linear problem, so the algorithm will converge in one
    // iteration, but this is a template for more general models with either
    // nonlinear hardening or curved yield surfaces.
    
    // for GPU we can simplify the DP model to avoid Newton, but Cam-Clay and
    // others will need it.
    
    array1d< real64 > solution(3), residual(3), delta(3);
    array2d< real64 > jacobian(3,3), jacobianInv(3,3);
    
    solution[0] = eps_v_trial; // initial guess for elastic volumetric strain
    solution[1] = eps_s_trial; // initial guess for elastic deviatoric strain
    solution[2] = 0;      // initial guess for plastic multiplier
    
    real64 norm,normZero = 1e30;
    jacobian.setValues< serialPolicy >( 0 );    

    for(localIndex iter=0; iter<10; ++iter) // could be fixed at one iter
    {

      trialP = p0 * std::exp(-1./Cr* (solution[0] - eps_v0));
      trialQ = 3. * mu * solution[1];
      bulkModulus = -trialP/Cr;
      pc = oldPc * std::exp(-1./(Cc-Cr)*(eps_v_trial-solution[0]));

      yield = trialQ*trialQ/(M*M)- alpha*alpha*trialP *(2.*alpha/(alpha+1.)*pc-trialP)+alpha*alpha*(alpha-1.)/(alpha+1.)* pc*pc;
      
      // derivatives of yield surface
      real64 alphaTerm = 2. * alpha*alpha*alpha / (alpha+1.); 
      real64 df_dp = -alphaTerm * pc + 2. * alpha * trialP;
      real64 df_dq = 2. * trialQ /(M*M); 
      real64 df_dpc = 2. * alpha*alpha*(alpha-1.) /(alpha+1.) * pc - alphaTerm * trialP ;
      real64 dpc_dve = -1./(Cc-Cr) * pc;

      real64 df_dp_dve = 2. * alpha * alpha * bulkModulus - alphaTerm * dpc_dve;
      real64 df_dq_dse = 2. /(M*M) * 3. * mu; 
      //real64 df_dpc_dve = -alphaTerm * bulkModulus + 2*alpha*alpha*(alpha-1) /(alpha+1) * dpc_dve;

      // assemble residual system
      
      residual[0] = solution[0] - eps_v_trial + solution[2]*df_dp; // strainElasticDev - strainElasticTrialDev + dlambda*dG/dPQ = 0
      residual[1] = solution[1] - eps_s_trial + solution[2]*df_dq;       // strainElasticVol - strainElasticTrialVol + dlambda*dG/dQ = 0
      residual[2] = yield;    // F = 0
      
      // check for convergence (can be avoided for linear model)
      
      norm = LvArray::tensorOps::l2Norm<3>(residual);
      
     std::cout << "iter= " << iter << ", norm =  " << norm << std::endl;
      
      if(iter==0)
      {
        normZero = norm;
      }
      
      if(norm < 1e-8*(normZero+1)) 
      {
        break;
      }
      
      // solve Newton system
      
      jacobian(0,0) = 1. + solution[2] * df_dp_dve;
      jacobian(0,2) = df_dp;
      jacobian(1,1) = 1. + solution[2]*df_dq_dse;
      jacobian(1,2) = df_dq;
      jacobian(2,0) = bulkModulus * df_dp - dpc_dve * df_dpc;
      jacobian(2,1) = 3.0 * mu * df_dq;
      jacobian(2,2) = 0.0;
      
      LvArray::tensorOps::invert<3>(jacobianInv,jacobian);
      LvArray::tensorOps::AijBj<3,3>(delta,jacobianInv,residual);
     
      for(localIndex i=0; i<3; ++i)
      {
        solution[i] -= delta[i];
      }
    } //end of iteration
    
    // construct stress = P*eye + sqrt(2/3)*Q*nhat

     //array1d< real64 > deviator2(3)   [no need to construct elastic strain?]
    for(localIndex i=0; i<6; ++i)
    {
      stress[i] = trialP * identity[i] + trialQ * sqrt23 *deviator[i];
    }
    
    // construct consistent tangent operator
    
    for(localIndex i=0; i<6; ++i)
    {
      for(localIndex j=0; j<6; ++j)
      {
        stiffness[i][j] = 0;
      }
    }

    array2d< real64 > BB(2,2);
    BB.setValues< serialPolicy >( 0 );

    real64 dpc_dve = -1./(Cc-Cr) * pc;
    real64 a1= 1. + solution[2]*dpc_dve; 
    real64 a2 = trialP * dpc_dve;  

    bulkModulus = -trialP/Cr; 

    BB[0][0] = bulkModulus*(a1*jacobianInv[0][0]+a2*jacobianInv[0][2]);
    BB[0][1] =bulkModulus*jacobianInv[0][1];
    BB[1][0] =3. * mu*(a1*jacobianInv[1][0]+a2*jacobianInv[1][2]);
    BB[1][1] = 3. * mu*jacobianInv[1][1];

    real64 c1; 
    
    if(eps_s_trial<1e-10) // confirm eps_s_trial != 0
    {
      c1 = 2. * mu;
    }else{
      c1 = 2. * trialQ/(3. * eps_s_trial);
    }

    real64 c2 = BB[0][0] - c1/3.;
    real64 c3 = sqrt23 * BB[0][1];
    real64 c4 = sqrt23 * BB[1][0];
    real64 c5 = 2./3. * BB[1][1] - c1;
    
    for(localIndex i=0; i<3; ++i)
    {
      stiffness[i][i] = c1;
      stiffness[i+3][i+3] = 0.5*c1;
    }
    
    for(localIndex i=0; i<6; ++i)
    {
      for(localIndex j=0; j<6; ++j)
      {
        stiffness[i][j] +=   c2 * identity[i] * identity[j]
                           + c3 * identity[i] * deviator[j]
                           + c4 * deviator[i] * identity[j]
                           + c5 * deviator[i] * deviator[j];
      }
    }
    
  } // end plastic branch
  
  // remember history variables before returning
  
  m_newPreConsolidationPressure[k][q] = pc;
  
  for(localIndex i=0; i<6; ++i)
  {
    m_newStress[k][q][i] = stress[i];
  }
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DelftEggUpdates::SaveConvergedState( localIndex const k,
                                               localIndex const q )
{
  m_oldPreConsolidationPressure[k][q] = m_newPreConsolidationPressure[k][q];
  
  for(localIndex i=0; i<6; ++i)
  {
    m_oldStress[k][q][i] = m_newStress[k][q][i];
  }
}


/**
 * @class DelftEgg
 *
 * Drucker-Prager material model.
 */
class DelftEgg : public SolidBase
{
public:

  /// @typedef Alias for DelftEggUpdates
  using KernelWrapper = DelftEggUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  DelftEgg( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~DelftEgg() override;

  virtual void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;
  
  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "DelftEgg";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string CatalogName() { return m_catalogNameString; }

  virtual string GetCatalogName() override { return CatalogName(); }

  ///@}

  /**
   * @struct Set of "char const *" and keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {

    /// string/key for default shear modulus (mu)
    static constexpr auto defaultShearModulusString = "defaultShearModulus";

    /// string/key for default reference pressure (p0)
    static constexpr auto defaultRefPressureString = "defaultRefPressure";
    
    /// string/key for default reference volumetric strain (eps_v0)
    static constexpr auto defaultRefStrainVolString = "defaultRefStrainVol";
    
    /// string/key for default recompression index (Cr)
    static constexpr auto defaultRecompressionIndexString = "defaultRecompressionIndex";
    
    /// string/key for default virgin compression index (Cc)
    static constexpr auto defaultVirginCompressionIndexString = "defaultVirginCompressionIndex";

    /// string/key for default slope of the critical state line (M)
    static constexpr auto defaultCslSlopeString = "defaultCslSlope";

    /// string/key for default shape parameter for yield surface (alpha)
    static constexpr auto defaultShapeParameterString = "defaultShapeParameter";

    /// string/key for default preconsolidation pressure (pc)
    static constexpr auto defaultPreConsolidationPressureString = "defaultPreConsolidationPressure";
        
    /// string/key for shear modulus
    static constexpr auto shearModulusString = "ShearModulus";
    
    /// string/key for reference pressure
    static constexpr auto refPressureString  = "RefPressure";
    
    /// string/key for reference volumetric strain
    static constexpr auto refStrainVolString  = "RefStrainVol";
    
    /// string/key for recompression index
    static constexpr auto recompressionIndexString  = "RecompressionIndex";
    
    /// string/key for virgin compression index
    static constexpr auto virginCompressionIndexString  = "VirginCompressionIndex";
    
    /// string/key for slope of the critical state line
    static constexpr auto cslSlopeString  = "CslSlope";
    
    /// string/key for shape parameter for the yield surface 
    static constexpr auto shapeParameterString  = "ShapeParameter";

    /// string/key for new preconsolidation pressure
    static constexpr auto newPreConsolidationPressureString  = "NewPreConsolidationPressure";

     /// string/key for old preconsolidation pressure
    static constexpr auto oldPreConsolidationPressureString  = "OldPreConsolidationPressure";

    /// string/key for new stress
    static constexpr auto newStressString  = "NewStress";
    
    /// string/key for old stress
    static constexpr auto oldStressString  = "OldStress";
  };

  /**
   * @brief Create a instantiation of the DelftEggUpdate class that refers to the data in this.
   * @return An instantiation of DelftEggUpdate.
   */
  DelftEggUpdates createKernelWrapper()
  {
    return DelftEggUpdates( m_shearModulus,
                                 m_refPressure,
                                 m_refStrainVol,
                                 m_recompressionIndex,
                                 m_virginCompressionIndex,
                                 m_cslSlope,
                                 m_shapeParameter,
                                 m_newPreConsolidationPressure,
                                 m_oldPreConsolidationPressure,
                                 m_newStress,
                                 m_oldStress,
                                 m_stress ); // TODO: m_stress is redundant with m_newStress
  }

protected:
  virtual void PostProcessInput() override;

private:
  //TODO: maybe define some helper structs for defaults, arrays, arrayViews, etc.
  
  /// Material parameter: The default value of the shear modulus
  real64 m_defaultShearModulus;
  
  /// Material parameter: The default value of reference pressure
  real64 m_defaultRefPressure;
  
  /// Material parameter: The default value of reference volumetric strain
  real64 m_defaultRefStrainVol;
  
  /// Material parameter: The default value of the recompression index
  real64 m_defaultRecompressionIndex;
  
  /// Material parameter: The default value of the virgin compresion index
  real64 m_defaultVirginCompressionIndex;

  /// Material parameter: The default value of slope of the critical state line
  real64 m_defaultCslSlope;

  /// Material parameter: The default value of the shape parameter of yield surface
  real64 m_defaultShapeParameter;

  /// Material parameter: The default value of preconsolidation pressure
  real64 m_defaultPreConsolidationPressure;

  /// Material parameter: The shear modulus for each element
  array1d< real64 > m_shearModulus;
  
  /// Material parameter: The reference pressure for each quadrature point
  array2d< real64 > m_refPressure;
  
  /// Material parameter: The reference volumetric strain for each quadrature point
  array2d< real64 > m_refStrainVol;
  
  /// Material parameter: The compression index for each element
  array1d< real64 > m_recompressionIndex;
  
  /// History variable: The virgin compresion index for each element
  array1d< real64 > m_virginCompressionIndex;

    /// History variable: The slope of the critical state line for each element
  array1d< real64 > m_cslSlope;

    /// History variable: The shape parameter of yield surface for each element
  array1d< real64 > m_shapeParameter;

    /// History variable: The current preconsolidation pressure for each quadrature point
  array2d< real64 > m_newPreConsolidationPressure;
  
  /// History variable: The previous preconsolidation pressure for each quadrature point
  array2d< real64 > m_oldPreConsolidationPressure;
  
  /// History variable: The current stress for each quadrature point
  array3d< real64, solid::STRESS_PERMUTATION > m_newStress; //TODO: redundant storage with base m_stress
  
  /// History variable: The previous stress for each quadrature point.
  array3d< real64, solid::STRESS_PERMUTATION > m_oldStress;
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DELFTEGG_HPP_ */