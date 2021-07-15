

=============== ============ ======== =========================================================================================== 
Name            Type         Default  Description                                                                                 
=============== ============ ======== =========================================================================================== 
fieldName       string       required The name of the (packable) field associated with the specified object to retrieve data from 
name            string       required A name is required for any non-unique nodes                                                 
objectPath      string       required The name of the object from which to retrieve field values.                                 
onlyOnSetChange localIndex   0        Whether or not to only collect when the collected sets of indices change in any way.        
setNames        string_array {}       The set(s) for which to retrieve data.                                                      
=============== ============ ======== =========================================================================================== 


