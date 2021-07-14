

<<<<<<< HEAD
=================== ======= ======== =============================================================== 
Name                Type    Default  Description                                                     
=================== ======= ======== =============================================================== 
material            string  required Solid material to test                                          
mode                string  required Test mode [triaxial, volumetric, oedometer]                     
name                string  required A name is required for any non-unique nodes                     
output              string  required Output file                                                     
steps               integer required Number of load steps to take                                    
strainFunction      string  required Function controlling strain loading (role depends on test mode) 
stressFunction      string  required Function controlling stress loading (role depends on test mode) 
useNumericalTangent integer 0        Flag to use numerical tangent                                   
=================== ======= ======== =============================================================== 
=======
============== ======= ======== =============================================================== 
Name           Type    Default  Description                                                     
============== ======= ======== =============================================================== 
baseline       path    none     Baseline file                                                   
logLevel       integer 0        Log level                                                       
material       string  required Solid material to test                                          
mode           string  required Test mode [triaxial, volumetric, oedometer]                     
name           string  required A name is required for any non-unique nodes                     
output         string  none     Output file                                                     
steps          integer required Number of load steps to take                                    
strainFunction string  required Function controlling strain loading (role depends on test mode) 
stressFunction string  required Function controlling stress loading (role depends on test mode) 
============== ======= ======== =============================================================== 
>>>>>>> develop


