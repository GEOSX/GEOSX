import numpy as np

# Rotate a vector in local coodinate of an inclined borehole to the global coordinate
def vectorRotation(x,y,z,phi_x,phi_z):
    rotx = np.array([[np.cos(phi_x), np.sin(phi_x),0.],[-np.sin(phi_x), np.cos(phi_x),0.],[0.,0.,1.]])
    rotz = np.array([[np.cos(phi_z),0., np.sin(phi_z)],[0.,1.,0.],[-np.sin(phi_z),0., np.cos(phi_z)]])
    
    localCoord = np.array([x,y,z])
    return np.dot( rotz, np.dot( rotx,localCoord ) )

phi_x = -90./180.*3.1416
phi_z = 45./180.*3.1416

start_point = (0, 0, 0)

# Rotate from local coordinates of a derived wellbore to global coordinates
end = vectorRotation( 4, 0, 0, phi_x , phi_z )
end_point = (end[0], end[1], end[2])
print(end)

OpenDatabase("localhost:"+"."+"/siloFiles/plot_* database", 0)

fieldName = "Omega_Solid_MaterialFields"
AddPlot("Pseudocolor", fieldName+"/stress_11", 1, 1)
DrawPlots()
SetQueryFloatFormat("%g")
Query("Lineout", end_point=end_point, num_samples=50, start_point=start_point, use_sampling=0, vars=(fieldName+"/stress_11"))
SetActiveWindow(2)
SetTimeSliderState(TimeSliderGetNStates()-1)
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = "."
SaveWindowAtts.fileName = "stress_11"
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.CURVE  
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
DeleteWindow()
SetActiveWindow(1)
SetActivePlots(0)
DeleteActivePlots()


AddPlot("Pseudocolor", fieldName+"/stress_12", 1, 1)
DrawPlots()
SetQueryFloatFormat("%g")
Query("Lineout", end_point=end_point, num_samples=50, start_point=start_point, use_sampling=0, vars=(fieldName+"/stress_12"))
SetActiveWindow(2)
SetTimeSliderState(TimeSliderGetNStates()-1)
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = "."
SaveWindowAtts.fileName = "stress_12"
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.CURVE  
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
DeleteWindow()
SetActiveWindow(1)
SetActivePlots(0)
DeleteActivePlots()


AddPlot("Pseudocolor", fieldName+"/stress_13", 1, 1)
DrawPlots()
SetQueryFloatFormat("%g")
Query("Lineout", end_point=end_point, num_samples=50, start_point=start_point, use_sampling=0, vars=(fieldName+"/stress_13"))
SetActiveWindow(2)
SetTimeSliderState(TimeSliderGetNStates()-1)
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = "."
SaveWindowAtts.fileName = "stress_13"
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.CURVE  
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
DeleteWindow()
SetActiveWindow(1)
SetActivePlots(0)
DeleteActivePlots()


AddPlot("Pseudocolor", fieldName+"/stress_22", 1, 1)
DrawPlots()
SetQueryFloatFormat("%g")
Query("Lineout", end_point=end_point, num_samples=50, start_point=start_point, use_sampling=0, vars=(fieldName+"/stress_22"))
SetActiveWindow(2)
SetTimeSliderState(TimeSliderGetNStates()-1)
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = "."
SaveWindowAtts.fileName = "stress_22"
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.CURVE  
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
DeleteWindow()
SetActiveWindow(1)
SetActivePlots(0)
DeleteActivePlots()

# Extract sig11
AddPlot("Pseudocolor", fieldName+"/stress_23", 1, 1)
DrawPlots()
SetQueryFloatFormat("%g")
Query("Lineout", end_point=end_point, num_samples=50, start_point=start_point, use_sampling=0, vars=(fieldName+"/stress_23"))
SetActiveWindow(2)
SetTimeSliderState(TimeSliderGetNStates()-1)
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = "."
SaveWindowAtts.fileName = "stress_23"
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.CURVE  
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
DeleteWindow()
SetActiveWindow(1)
SetActivePlots(0)
DeleteActivePlots()


AddPlot("Pseudocolor", fieldName+"/stress_33", 1, 1)
DrawPlots()
SetQueryFloatFormat("%g")
Query("Lineout", end_point=end_point, num_samples=50, start_point=start_point, use_sampling=0, vars=(fieldName+"/stress_33"))
SetActiveWindow(2)
SetTimeSliderState(TimeSliderGetNStates()-1)


SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = "."
SaveWindowAtts.fileName = "stress_33"
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.CURVE  
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
DeleteWindow()
SetActiveWindow(1)
SetActivePlots(0)
DeleteActivePlots()

fieldName = "Omega_ElementFields/pressure"
AddPlot("Pseudocolor", fieldName, 1, 1)
DrawPlots()
SetQueryFloatFormat("%g")
Query("Lineout", end_point=end_point, num_samples=50, start_point=start_point, use_sampling=0, vars=(fieldName))
SetActiveWindow(2)
SetTimeSliderState(TimeSliderGetNStates()-1)
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = "."
SaveWindowAtts.fileName = "pressure"
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.CURVE  
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
DeleteWindow()
SetActiveWindow(1)
SetActivePlots(0)
DeleteActivePlots()
