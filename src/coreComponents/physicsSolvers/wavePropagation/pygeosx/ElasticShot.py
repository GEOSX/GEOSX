import numpy as np
import pygeosx

from segyManager import *
from fileManager import rootPath
import os
from print import *


def elastic_shots(rank, problem, shot_list, tracePath):
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
    Export displacement values to different segy based on direction of displacement (x, y, z), if output flag defined in XML == 1
    """

    #Get Acoustic group
    elastic_group  = problem.get_group("Solvers/elasticSolver")

    #Get Wrappers
    src_pos_geosx   = elastic_group.get_wrapper("sourceCoordinates").value()
    src_pos_geosx.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    rcv_pos_geosx   = elastic_group.get_wrapper("receiverCoordinates").value()
    rcv_pos_geosx.set_access_level(pygeosx.pylvarray.RESIZEABLE)


    displacement_geosx  = elastic_group.get_wrapper("displacementNp1AtReceivers").value()


    displacementx_nm1 = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/displacementx_nm1").value()
    displacementx_nm1.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    displacementy_nm1 = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/displacementy_nm1").value()
    displacementy_nm1.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    displacementz_nm1 = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/displacementz_nm1").value()
    displacementz_nm1.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    displacementx_n = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/displacementx_n").value()
    displacementx_n.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    displacementy_n = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/displacementy_n").value()
    displacementy_n.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    displacementz_n = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/displacementz_n").value()
    displacementz_n.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    displacementx_np1 = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/displacementx_np1").value()
    displacementx_np1.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    displacementy_np1 = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/displacementy_np1").value()
    displacementy_np1.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    displacementz_np1 = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/displacementz_np1").value()
    displacementz_np1.set_access_level(pygeosx.pylvarray.MODIFIABLE)


    outputSismoTrace = elastic_group.get_wrapper("outputSismoTrace").value()

    dt_geosx   = problem.get_wrapper("Events/solverApplications/forceDt").value()
    maxT       = problem.get_wrapper("Events/maxTime").value()
    cycle_freq = problem.get_wrapper("Events/python/cycleFrequency").value()
    cycle      = problem.get_wrapper("Events/python/lastCycle").value()
    curr_time  = problem.get_wrapper("Events/time").value()

    dt            = shot_list[0].getSource().getTimeStep()
    dt_geosx[0]   = dt
    maxCycle      = int(maxT[0]/dt)
    cycle_freq[0] = 1
    maxT[0]       = (maxCycle+1) * dt
    dt_cycle      = 1 #to be changed

    displacementx_at_receivers = np.zeros((maxCycle+1, shot_list[0].getReceiverSet().getNumberOfReceivers()))
    displacementy_at_receivers = np.zeros((maxCycle+1, shot_list[0].getReceiverSet().getNumberOfReceivers()))
    displacementz_at_receivers = np.zeros((maxCycle+1, shot_list[0].getReceiverSet().getNumberOfReceivers()))

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
            create_segy(shot_list, "displacementx", maxCycle+1, tracePath)
            create_segy(shot_list, "displacementy", maxCycle+1, tracePath)
            create_segy(shot_list, "displacementz", maxCycle+1, tracePath)


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
        if cycle[0] < (ishot+1) * maxCycle:
            displacementx_at_receivers[cycle[0] - ishot * maxCycle, :] = displacement_geosx.to_numpy()[:, 0]
            displacementy_at_receivers[cycle[0] - ishot * maxCycle, :] = displacement_geosx.to_numpy()[:, 1]
            displacementz_at_receivers[cycle[0] - ishot * maxCycle, :] = displacement_geosx.to_numpy()[:, 2]

        else:
            displacementx_at_receivers[maxCycle, :] = displacement_geosx.to_numpy()[:, 0]
            displacementy_at_receivers[maxCycle, :] = displacement_geosx.to_numpy()[:, 1]
            displacementz_at_receivers[maxCycle, :] = displacement_geosx.to_numpy()[:, 2]

            #Segy export and flag update
            if outputSismoTrace == 1 :
                segyFilex = os.path.join(tracePath, "displacementx_Shot"+ str(ishot) + ".sgy")
                export_to_segy(displacementx_at_receivers,
                               shot_list[ishot].getReceiverSet().getSetCoord(),
                               segyFilex)

                segyFiley = os.path.join(tracePath, "displacementy_Shot"+ str(ishot) + ".sgy")
                export_to_segy(displacementy_at_receivers,
                               shot_list[ishot].getReceiverSet().getSetCoord(),
                               segyFiley)

                segyFilez = os.path.join(tracePath, "displacementz_Shot"+ str(ishot) + ".sgy")
                export_to_segy(displacementz_at_receivers,
                               shot_list[ishot].getReceiverSet().getSetCoord(),
                               segyFilez)


            shot_list[ishot].flagUpdate("Done")

            #Reset time to -dt and pressure to 0
            curr_time[0] = -dt
            displacementx_nm1.to_numpy()[:] = 0.0
            displacementx_n.to_numpy()[:]   = 0.0
            displacementx_np1.to_numpy()[:] = 0.0

            displacementy_nm1.to_numpy()[:] = 0.0
            displacementy_n.to_numpy()[:]   = 0.0
            displacementy_np1.to_numpy()[:] = 0.0

            displacementz_nm1.to_numpy()[:] = 0.0
            displacementz_n.to_numpy()[:]   = 0.0
            displacementz_np1.to_numpy()[:] = 0.0

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
