

============================ ====== ======== ==================================================================== 
Name                         Type   Default  Description                                                          
============================ ====== ======== ==================================================================== 
BiotCoefficient              real64 1        Biot's coefficient                                                   
compressibility              real64 0        Pore volume compressibilty                                           
defaultBulkModulus           real64 -1       Elastic Bulk Modulus Parameter                                       
defaultCohesion              real64 0        Initial cohesion                                                     
defaultDensity               real64 required Default Material Density                                             
defaultDilationRatio         real64 1        Dilation ratio [0,1] (ratio = tan dilationAngle / tan frictionAngle) 
defaultHardening             real64 0        Hardening parameter (hardening rate is faster for smaller values)    
defaultInitialFrictionAngle  real64 30       Initial friction angle (degrees)                                     
defaultPoissonRatio          real64 -1       Poisson's ratio                                                      
defaultResidualFrictionAngle real64 30       Residual friction angle (degrees)                                    
defaultShearModulus          real64 -1       Elastic Shear Modulus Parameter                                      
defaultYoungModulus          real64 -1       Elastic Young's Modulus.                                             
name                         string required A name is required for any non-unique nodes                          
referencePressure            real64 0        ReferencePressure                                                    
============================ ====== ======== ==================================================================== 


