'''
Created on 8/02/2021

@author: macpro
'''

import numpy as np
from mesh import *
from acquisition import *
import pygeosx
import sys
#from mpi4py import MPI
from segy import *
from print import *

    

def shot_simul(shot_list, dt):
    """
    Parameters
    ----------

    shot_list : list
        A list containing sets of Source and ReceiverSet objects

    dt : float
        Time step for the solver
    """
      
    problem = initialize()      
    do_shots(problem, shot_list, dt)
    
    
   
def initialize():
    """ Grouping of pygeox initializations
    
    Return
    ------
    problem :
        The pygeosx ProblemManager
        
    Notes
    -----
    Need to give MPI rank at this point for initialization. 
    Conflict with first initialization to get the list of shots
    """
    
    problem = pygeosx.reinit(sys.argv)
    #Should be : 
    #    problem = pygeosx.initialize()(rank, sys.argv), 
    #Once the seismic acquisition from segy file is implemented
    
    pygeosx.apply_initial_conditions()
    
    return problem



def do_shots(problem, shot_list, dt):
    """ Given a GEOSX problem, a list of shots, and a time step, 
        solve wave eqn with different configurations 
        
    Parameters
    ----------   
    problem :
        The pygeosx ProblemManager
        
    shot_list : list
        A list containing sets of Source and ReceiverSet objects

    dt : float
        Time step for the solver
        
    Notes
    -----    
    Export pressure to segy depending on output flag defined in XML 
    """

    #Get Acoustic group
    acoustic_group  = problem.get_group("Solvers/acousticSolver")
    
    
    #Get Wrappers
    src_pos_geosx   = acoustic_group.get_wrapper("sourceCoordinates").value()
    src_pos_geosx.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    
    rcv_pos_geosx   = acoustic_group.get_wrapper("receiverCoordinates").value()
    rcv_pos_geosx.set_access_level(pygeosx.pylvarray.RESIZEABLE)
    
    pressure_geosx  = acoustic_group.get_wrapper("pressureNp1AtReceivers").value()
    
    pressure_nm1 = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/pressure_nm1").value()
    pressure_nm1.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    
    pressure_n = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/pressure_n").value()
    pressure_n.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    
    pressure_np1 = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/pressure_np1").value()
    pressure_np1.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    
    outputSismoTrace = acoustic_group.get_wrapper("outputSismoTrace").value()
    
    dt_geosx        = problem.get_wrapper("Events/solverApplications/forceDt").value()
    maxT            = problem.get_wrapper("Events/maxTime").value()
    cycle_freq      = problem.get_wrapper("Events/python/cycleFrequency").value()
    cycle           = problem.get_wrapper("Events/python/lastCycle").value()
    cflFactor       = problem.get_wrapper("Solvers/acousticSolver/cflFactor").value()[0]
    curr_time       = problem.get_wrapper("Events/time").value()
    
    
    cycle_freq[0] = 1
    dt_geosx[0]   = dt       
    maxCycle      = int(maxT[0]/dt)
    maxT[0]       = (maxCycle+1) * dt 
    dt_cycle      = 1 
    
    pressure_at_receivers = np.zeros((maxCycle+1, shot_list[0].getReceiverSet().getNumberOfReceivers())) 
    nb_shot = len(shot_list)
    
    ishot = 0
    print_shot_config(shot_list, ishot)
    
    #Set first source and receivers positions in GEOSX        
    src_pos          = shot_list[ishot].getSource().getCoord()       
    rcv_pos_list     = shot_list[ishot].getReceiverSet().getSetCoord() 
    
    src_pos_geosx.to_numpy()[0] = src_pos 
    rcv_pos_geosx.resize(len(rcv_pos_list))    
    rcv_pos_geosx.to_numpy()[:] = rcv_pos_list[:]
    
    #Update shot flag
    shot_list[ishot].flagUpdate("In Progress")
    print_flag(shot_list)
    
    while (np.array([shot.getFlag() for shot in shot_list]) == "Done").all() != True and pygeosx.run() != pygeosx.COMPLETED:        
        #Save pressure
        if cycle[0] < (ishot+1) * maxCycle:
            pressure_at_receivers[cycle[0] - ishot * maxCycle, :] = pressure_geosx.to_numpy()[:]
            
        else:
            pressure_at_receivers[maxCycle, :] = pressure_geosx.to_numpy()[:]
            print_pressure(pressure_at_receivers, ishot)
           
            #Segy export and flag update
            if outputSismoTrace == 1 :
                export_to_segy(pressure_at_receivers, shot_list[0].getReceiverSet().getSetCoord(), ishot, dt_cycle)
                
            shot_list[ishot].flagUpdate("Done")
            
            #Reset time/pressure to 0
            curr_time[0]    = 0.0
            pressure_nm1.to_numpy()[:] = 0.0
            pressure_n.to_numpy()[:]   = 0.0
            pressure_np1.to_numpy()[:] = 0.0
           
            #Increment shot
            ishot += 1             
            if ishot < nb_shot:
                print_shot_config(shot_list, ishot)

                #Set new receivers and source positions in GEOSX
                src_pos          = shot_list[ishot].getSource().getCoord()       
                rcv_pos_list     = shot_list[ishot].getReceiverSet().getSetCoord() 
    
                src_pos_geosx.to_numpy()[0] = src_pos 
                rcv_pos_geosx.resize(len(rcv_pos_list))    
                rcv_pos_geosx.to_numpy()[:] = rcv_pos_list[:]
                
                #Update shot flag
                shot_list[ishot].flagUpdate("In Progress") 
                
            print_flag(shot_list) 




def ricker(maxT, dt, f0):
    """ Source function
    
    Parameters
    ----------
    maxT : float
        The max time for simulation
    
    dt : float
        The time step for simulation
    
    f0 : float
        Intensity
    
    Return
    ------
    fi :
        np array containing source value at all timestep
    """
        
    T0 = 1.0/f0;
    fi = np.zeros(int(maxT/dt))
    
    for t in range(int(maxT/dt)):
        t0 = dt*t
        if t0 <= -0.9*T0 or t0 >= 2.9*T0:
            fi[t] = 0.0;
        else:
            tmp      = f0*t0-1.0
            f0tm1_2  = 2*(tmp*np.pi)*(tmp*np.pi)
            gaussian = m.exp( -f0tm1_2 )
            fi[t]    = -(t0-1)*gaussian
               
    return fi



#Calculate dt using order of space discretization method, Wave velocity, and the radius of the included sphere in element 
def calculDt(mesh):
    """Calcul timestep for simulation based on mesh parameters
    
    Parameters
    ----------
    mesh : object
        A Mesh python object derived from GEOSX mesh 
        
    Return 
    ------
    dt : float
        Time step for simulation
    """
    
    if mesh.getOrd()==1:
        nx = 1
        ny = 2
        nz = 4
    elif mesh.getOrd()==3:
        nx = 3
        ny = 12
        nz = 48
    elif mesh.getOrd()==5:
        nx = 5
        ny = 30
        nz = 180
    

    h    = np.linalg.norm(mesh.getElem_List()[0].getNode_List()[0].getCoords() - mesh.getElem_List()[0].getNode_List()[nx].getCoords())/2
    Vmax = mesh.getElem_List()[0].getSpeed()
    V    = mesh.getElem_List()[0].getVolume()
    
    for elem in (mesh.getElem_List()):
        
        if elem.getVolume()<V: 
        
            xradius = np.linalg.norm(elem.getNode_List()[0].getCoords() - elem.getNode_List()[nx].getCoords())/2
            yradius = np.linalg.norm(elem.getNode_List()[0].getCoords() - elem.getNode_List()[ny].getCoords())/2
            zradius = np.linalg.norm(elem.getNode_List()[0].getCoords() - elem.getNode_List()[nz].getCoords())/2
            
            if xradius < h:
                h = xradius
            if yradius < h:
                h = yradius
            if zradius < h:
                h = zradius
            
        if Vmax < elem.getSpeed():
            Vmax = elem.getSpeed()
            
        
    dt = h/(Vmax*mesh.ord)           
            
    return dt
    

            
