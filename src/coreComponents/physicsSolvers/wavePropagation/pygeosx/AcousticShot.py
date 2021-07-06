'''
Created on 8/02/2021

@author: macpro
'''

import numpy as np
import pygeosx
from math import ceil

from segyManager import *
from fileManager import rootPath
import os
from print import *


def acoustic_shots(rank, problem, shot_list, tracePath):
    """ Given a GEOSX problem, a list of shots, and a time step,
        solve wave eqn with different configurations

    Parameters
    ----------
    rank : int
        MPI rank

    problem :
        The pygeosx ProblemManager

    shot_list : list
        A list containing sets of Source and ReceiverSet objects

    tracePath : string
        Directory where to output .sgy files

    Notes
    -----
    Export pressure to segy depending on output flag defined in XML
    """

    #Get Acoustic group
    acoustic_group  = problem.get_group("Solvers/acousticSolver")


    #Get Wrappers
    src_pos_geosx  = acoustic_group.get_wrapper("sourceCoordinates").value()
    src_pos_geosx.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    rcv_pos_geosx  = acoustic_group.get_wrapper("receiverCoordinates").value()
    rcv_pos_geosx.set_access_level(pygeosx.pylvarray.RESIZEABLE)

    pressure_geosx = acoustic_group.get_wrapper("pressureNp1AtReceivers").value()

    pressure_nm1 = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/pressure_nm1").value()
    pressure_nm1.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    pressure_n = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/pressure_n").value()
    pressure_n.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    pressure_np1 = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/pressure_np1").value()
    pressure_np1.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    nsamples = acoustic_group.get_wrapper("nSampleSismoTrace").value()
    outputSismoTrace = acoustic_group.get_wrapper("outputSismoTrace").value()
    indexSismoTrace = acoustic_group.get_wrapper("indexSismoTrace").value()

    dt_geosx   = problem.get_wrapper("Events/solverApplications/forceDt").value()
    dt_sismo   = acoustic_group.get_wrapper("dtSismoTrace").value()
    maxT       = problem.get_wrapper("Events/maxTime").value()
    cycle_freq = problem.get_wrapper("Events/python/cycleFrequency").value()
    cycle      = problem.get_wrapper("Events/python/lastCycle").value()
    cycle2     = problem.get_wrapper("Events/cycle").value()
    curr_time  = problem.get_wrapper("Events/time").value()

    dt            = shot_list[0].getSource().getTimeStep()
    dt_geosx[0]   = dt
    maxCycle      = int(ceil(maxT[0]/dt))
    cycle_freq[0] = maxCycle
    maxT[0]       = (maxCycle+1)*dt
    nsamples[0]   = int(maxT[0]/dt_sismo)

    pressure_at_receivers = np.zeros((nsamples[0], shot_list[0].getReceiverSet().getNumberOfReceivers()))
    nb_shot = len(shot_list)

    if outputSismoTrace == 1 :
        if rank==0:
            if os.path.exists(os.path.join(rootPath,"outputSismoTrace/")):
                pass
            else:
                os.mkdir(os.path.join(rootPath,"outputSismoTrace/"))
            if os.path.exists(tracePath):
                pass
            else:
                os.mkdir(tracePath)
            create_segy(shot_list, "pressure", nsamples[0], tracePath)


    ishot = 0
    if rank==0:
        print_shot_config(shot_list, ishot)

    #Set first source and receivers positions in GEOSX
    src_pos      = shot_list[ishot].getSource().getCoord()
    rcv_pos_list = shot_list[ishot].getReceiverSet().getSetCoord()

    src_pos_geosx.to_numpy()[0] = src_pos
    rcv_pos_geosx.resize(len(rcv_pos_list))
    rcv_pos_geosx.to_numpy()[:] = rcv_pos_list[:]

    #Update shot flag
    shot_list[ishot].flagUpdate("In Progress")
    if rank==0:
        print_flag(shot_list)

    while (np.array([shot.getFlag() for shot in shot_list]) == "Done").all() != True and pygeosx.run() != pygeosx.COMPLETED:
        #Save pressure
        if cycle[0] !=0 :
            pressure_at_receivers[:, :] = pressure_geosx.to_numpy().transpose()

            #Segy export and flag update
            if outputSismoTrace == 1 :
                segyFile = os.path.join(tracePath, "pressure_Shot"+ str(ishot) + ".sgy")
                export_to_segy(pressure_at_receivers,
                               shot_list[ishot].getReceiverSet().getSetCoord(),
                               segyFile)

            shot_list[ishot].flagUpdate("Done")

            #Reset time to -dt and pressure to 0
            curr_time[0]               = -dt
            cycle_freq[0]              = maxCycle + 1
            #cycle2[0]                  = 0
            pressure_nm1.to_numpy()[:] = 0.0
            pressure_n.to_numpy()[:]   = 0.0
            pressure_np1.to_numpy()[:] = 0.0
            indexSismoTrace[0]         = 0

            #Increment shot
            ishot += 1
            if ishot < nb_shot:
                if rank==0:
                    print_shot_config(shot_list, ishot)

                #Set new receivers and source positions in GEOSX
                src_pos          = shot_list[ishot].getSource().getCoord()
                rcv_pos_list     = shot_list[ishot].getReceiverSet().getSetCoord()

                src_pos_geosx.to_numpy()[0] = src_pos
                rcv_pos_geosx.resize(len(rcv_pos_list))
                rcv_pos_geosx.to_numpy()[:] = rcv_pos_list[:]

                #Update shot flag
                shot_list[ishot].flagUpdate("In Progress")

            if rank==0:
                print_flag(shot_list)
