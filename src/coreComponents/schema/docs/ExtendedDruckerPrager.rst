

============================ ====== ======== =============================================================================== 
Name                         Type   Default  Description                                                                     
============================ ====== ======== =============================================================================== 
defaultBulkModulus           real64 -1       Elastic Bulk Modulus Parameter                                                  
defaultCohesion              real64 0        Initial cohesion                                                                
defaultDilationRatio         real64 1        Dilation ratio [0,1] (ratio = tan dilationAngle / tan frictionAngle)            
defaultGrainDensity          real64 required Default grain density. It's used to set the default value of the grain density. 
defaultHardening             real64 0        Hardening parameter (hardening rate is faster for smaller values)               
defaultInitialFrictionAngle  real64 30       Initial friction angle (degrees)                                                
defaultPoissonRatio          real64 -1       Poisson's ratio                                                                 
defaultReferencePorosity     real64 required Default value of the reference porosity                                         
defaultResidualFrictionAngle real64 30       Residual friction angle (degrees)                                               
defaultShearModulus          real64 -1       Elastic Shear Modulus Parameter                                                 
defaultYoungsModulus         real64 -1       Elastic Young's Modulus.                                                        
grainBulkModulus             real64 required Grain bulk modulus                                                              
name                         string required A name is required for any non-unique nodes                                     
============================ ====== ======== =============================================================================== 


