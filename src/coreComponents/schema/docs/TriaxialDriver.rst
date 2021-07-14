

============= ======= ======== ===================================================================== 
Name          Type    Default  Description                                                           
============= ======= ======== ===================================================================== 
axialControl  string  required Function controlling axial stress or strain (depending on test mode)  
baseline      path    none     Baseline file                                                         
initialStress real64  required Initial stress (scalar used to set an isotropic stress state)         
logLevel      integer 0        Log level                                                             
material      string  required Solid material to test                                                
mode          string  required Test mode [stressControl, strainControl, mixedControl]                
name          string  required A name is required for any non-unique nodes                           
output        string  none     Output file                                                           
radialControl string  required Function controlling radial stress or strain (depending on test mode) 
steps         integer required Number of load steps to take                                          
============= ======= ======== ===================================================================== 


