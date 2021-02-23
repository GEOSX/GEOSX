#
# Specify all headers
#
set( constitutive_headers
     ConstitutiveManager.hpp
     ConstitutiveBase.hpp
     ConstitutivePassThruHandler.hpp
     ExponentialRelation.hpp
     NullModel.hpp
     capillaryPressure/CapillaryPressureBase.hpp
     capillaryPressure/capillaryPressureSelector.hpp
     capillaryPressure/BrooksCoreyCapillaryPressure.hpp
     capillaryPressure/VanGenuchtenCapillaryPressure.hpp
     contact/ContactRelationBase.hpp
     contact/Coulomb.hpp
     fluid/MultiPhaseMultiComponentFluid.hpp     
     fluid/PVTFunctions/BrineCO2DensityFunction.hpp
     fluid/PVTFunctions/BrineViscosityFunction.hpp
     fluid/PVTFunctions/CO2SolubilityFunction.hpp
     fluid/PVTFunctions/FenghourCO2ViscosityFunction.hpp
     fluid/PVTFunctions/BrineEnthalpyFunction.hpp
     fluid/PVTFunctions/BrineInternalEnergyFunction.hpp
     fluid/PVTFunctions/CO2EnthalpyFunction.hpp
     fluid/PVTFunctions/CO2InternalEnergyFunction.hpp
     fluid/PVTFunctions/FlashModelBase.hpp
     fluid/PVTFunctions/PVTFunctionBase.hpp
     fluid/PVTFunctions/SpanWagnerCO2DensityFunction.hpp
     fluid/PVTFunctions/UtilityFunctions.hpp
     fluid/SingleFluidBase.hpp
     fluid/singleFluidSelector.hpp
     fluid/MultiFluidBase.hpp
     fluid/multiFluidSelector.hpp
     fluid/MultiFluidUtils.hpp
     fluid/CompressibleSinglePhaseFluid.hpp
     fluid/SlurryFluidBase.hpp
     fluid/slurryFluidSelector.hpp
     fluid/ProppantSlurryFluid.hpp
     fluid/ParticleFluidBase.hpp
     fluid/particleFluidSelector.hpp
     solid/Damage.hpp
     solid/DamageVolDev.hpp
     solid/DamageSpectral.hpp
     solid/DamageSpectralUtilities.hpp
     solid/DruckerPrager.hpp
     solid/DruckerPragerExtended.hpp
     solid/ElasticIsotropic.hpp
     solid/ElasticTransverseIsotropic.hpp
     solid/PoroElastic.hpp
     solid/PoreVolumeCompressibleSolid.hpp
     solid/PropertyConversions.hpp
     solid/SolidBase.hpp
     solid/SolidModelDiscretizationOps.hpp
     solid/SolidModelDiscretizationOpsFullyAnisotroipic.hpp
     solid/SolidModelDiscretizationOpsIsotropic.hpp
     solid/SolidModelDiscretizationOpsTransverseIsotropic.hpp
     relativePermeability/RelativePermeabilityBase.hpp
     relativePermeability/BrooksCoreyRelativePermeability.hpp
     relativePermeability/BrooksCoreyBakerRelativePermeability.hpp
     relativePermeability/relativePermeabilitySelector.hpp
     relativePermeability/VanGenuchtenBakerRelativePermeability.hpp
     relativePermeability/TableRelativePermeability.hpp
     relativePermeability/RelativePermeabilityInterpolators.hpp
     capillaryPressure/CapillaryPressureBase.hpp
     capillaryPressure/BrooksCoreyCapillaryPressure.hpp
     capillaryPressure/VanGenuchtenCapillaryPressure.hpp
   )
#
# Specify all sources
#
set( constitutive_sources
     ConstitutiveManager.cpp
     ConstitutiveBase.cpp
     NullModel.cpp
     capillaryPressure/CapillaryPressureBase.cpp
     capillaryPressure/BrooksCoreyCapillaryPressure.cpp
     capillaryPressure/VanGenuchtenCapillaryPressure.cpp
     contact/ContactRelationBase.cpp
     contact/Coulomb.cpp
     fluid/CompressibleSinglePhaseFluid.cpp
     fluid/MultiPhaseMultiComponentFluid.cpp     
     fluid/PVTFunctions/BrineCO2DensityFunction.cpp
     fluid/PVTFunctions/BrineViscosityFunction.cpp 
     fluid/PVTFunctions/CO2SolubilityFunction.cpp
     fluid/PVTFunctions/FenghourCO2ViscosityFunction.cpp
     fluid/PVTFunctions/SpanWagnerCO2DensityFunction.cpp
     fluid/PVTFunctions/BrineEnthalpyFunction.cpp
     fluid/PVTFunctions/BrineInternalEnergyFunction.cpp
     fluid/PVTFunctions/CO2EnthalpyFunction.cpp
     fluid/PVTFunctions/CO2InternalEnergyFunction.cpp    
     fluid/PVTFunctions/UtilityFunctions.cpp
     fluid/SingleFluidBase.cpp
     fluid/MultiFluidBase.cpp
     fluid/CompressibleSinglePhaseFluid.cpp
     fluid/SlurryFluidBase.cpp
     fluid/ProppantSlurryFluid.cpp
     fluid/ParticleFluidBase.cpp
     fluid/ParticleFluid.cpp               
     solid/Damage.cpp
     solid/DamageVolDev.cpp
     solid/DamageSpectral.cpp
     solid/DruckerPrager.cpp
     solid/DruckerPragerExtended.cpp
     solid/ElasticIsotropic.cpp
     solid/ElasticTransverseIsotropic.cpp
     solid/PoroElastic.cpp
     solid/PoreVolumeCompressibleSolid.cpp
     solid/SolidBase.cpp
     relativePermeability/RelativePermeabilityBase.cpp
     relativePermeability/BrooksCoreyRelativePermeability.cpp
     relativePermeability/BrooksCoreyBakerRelativePermeability.cpp
     relativePermeability/VanGenuchtenBakerRelativePermeability.cpp
     relativePermeability/TableRelativePermeability.cpp     
     capillaryPressure/CapillaryPressureBase.cpp
     capillaryPressure/BrooksCoreyCapillaryPressure.cpp
     capillaryPressure/VanGenuchtenCapillaryPressure.cpp
   )

if( BUILD_OBJ_LIBS)
  set( dependencyList dataRepository)
else()
  set( dependencyList common )
endif()



if( ENABLE_PVTPackage )
    set( constitutive_headers
         ${constitutive_headers}
         fluid/MultiFluidPVTPackageWrapper.hpp
         fluid/CompositionalMultiphaseFluid.hpp
         fluid/BlackOilFluid.hpp
       )
    set( constitutive_sources
         ${constitutive_sources}
         fluid/MultiFluidPVTPackageWrapper.cpp
         fluid/CompositionalMultiphaseFluid.cpp
         fluid/BlackOilFluid.cpp
       )
    set( dependencyList ${dependencyList} PVTPackage )
endif()

if ( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  constitutive
                 SOURCES               ${constitutive_sources}
                 HEADERS               ${constitutive_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )
               
target_include_directories( constitutive PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX constitutive )
