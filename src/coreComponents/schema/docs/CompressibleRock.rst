

======================== ============== ======== =============================================================================== 
Name                     Type           Default  Description                                                                     
======================== ============== ======== =============================================================================== 
compressibility          real64         required Solid compressibility                                                           
defaultGrainDensity      real64         0        Default grain density. It's used to set the default value of the grain density. 
defaultReferencePorosity real64         required Default value of the reference porosity                                         
grainBulkModulus         real64         0        Grain bulk modulus                                                              
grainDensity             real64_array2d {{0}}    Grain density                                                                   
name                     string         required A name is required for any non-unique nodes                                     
referencePressure        real64         required Reference pressure for solid compressibility                                    
======================== ============== ======== =============================================================================== 


