

======================== ========================== ======== =================================================================================================================================================== 
Name                     Type                       Default  Description                                                                                                                                         
======================== ========================== ======== =================================================================================================================================================== 
control                  geosx_WellControls_Control required | Well control. Valid options:                                                                                                                        
                                                             | * BHP                                                                                                                                               
                                                             | * phaseVolRate                                                                                                                                      
                                                             | * totalVolRate                                                                                                                                      
injectionStream          real64_array               {-1}     Global component densities for the injection stream                                                                                                 
name                     string                     required A name is required for any non-unique nodes                                                                                                         
referenceElevation       real64                     required Reference elevation where BHP control is enforced                                                                                                   
surfacePressure          real64                     0        Surface pressure used to compute volumetric rates when surface conditions are used                                                                  
surfaceTemperature       real64                     0        Surface temperature used to compute volumetric rates when surface conditions are used                                                               
targetBHP                real64                     -1       Target bottom-hole pressure                                                                                                                         
targetBHPTableName       string                              Name of the BHP table when the rate is a time dependent function                                                                                    
targetPhaseName          string                              Name of the target phase                                                                                                                            
targetPhaseRate          real64                     0        Target phase volumetric rate                                                                                                                        
targetPhaseRateTableName string                              Name of the phase rate table when the rate is a time dependent function                                                                             
targetTotalRate          real64                     0        Target total volumetric rate                                                                                                                        
targetTotalRateTableName string                              Name of the rate table when the rate is a time dependent function                                                                                   
type                     geosx_WellControls_Type    required | Well type. Valid options:                                                                                                                           
                                                             | * producer                                                                                                                                          
                                                             | * injector                                                                                                                                          
useSurfaceConditions     integer                    0        | Flag to specify whether rates are checked at surface or reservoir conditions.                                                                       
                                                             | Equal to 1 for surface conditions, and to 0 for reservoir conditions                                                                                
======================== ========================== ======== =================================================================================================================================================== 


