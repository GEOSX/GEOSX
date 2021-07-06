from fileManager import *
import multiprocessing as mp
from dask.distributed import Client, as_completed
from dask_jobqueue import SLURMCluster


def multiProcessing(shot_list, GEOSX, xml,
                    processes = 1,
                    cluster = False,
                    mpiranks = 1,
                    nodes = 1,
                    nodeName = None,
                    x = 1, y = 1, z = 1):

    """Cut shot list into desired amount of processes and run each process in parrallel
       using the Multiprocessing library

    Parameters
    ----------
    shot_list : list of Shot object
        list of shot configuration

    GEOSX : path
        location of GEOSX directory

    xml : xml file
        Time step

    processes : int
        Number of processes

    mpirank : int
        Number of MPI processes for parrallel solving

    x : int
        Number of MPI partition along x axe (domain partition)

    y : int
        Number of MPI partition along x axe (domain partition)

    z : int
        Number of MPI partition along x axe (domain partition)
    """


    p = []
    nb_shot_m1 = len(shot_list)
    ind = 0

    #Loop over the process launch
    for i in range(processes):
        nb_shot = int(nb_shot_m1/(processes - i))

        shot_file = exportShotList(i, shot_list[ind:ind + nb_shot])

        p.append( mp.Process(target = callProcess,
                             args   = (shot_file, GEOSX, xml, cluster, mpiranks, nodes, nodeName, x, y, z)) )

        ind = ind + nb_shot
        nb_shot_m1 = nb_shot_m1 - nb_shot

        #Start process
        p[i].start()


    for i in range(processes):
        p[i].join()

    #remove_shots_files()
def callProcess(shot_file, GEOSX, xml, cluster, mpiranks, nodes, nodeName, x, y, z):

    """Sub function of multriProcessing(). Call of different processes

    Parameters
    ----------
    shot_list : list of Shot object
        list of shot configuration

    GEOSX : path
        location of GEOSX directory

    xml : xml file
        Time step

    mpirank : int
        Number of MPI processes for parrallel solving

    x : int
        Number of MPI partition along x axe (domain partition)

    y : int
        Number of MPI partition along x axe (domain partition)

    z : int
        Number of MPI partition along x axe (domain partition)
    """

    wavePropagationPath = os.path.join(GEOSX, "GEOSX/src/coreComponents/physicsSolvers/wavePropagation/")
    pygeosxPath = os.path.join(wavePropagationPath, "pygeosx/")

    tracePath = os.path.join(rootPath, "outputSismoTrace/")
    traceProcPath = os.path.join(tracePath, "Process" + str(mp.current_process()._identity)[1] + "/")

    if cluster == False:
        os.system("mpirun -np " + str(mpiranks) + " python " + pygeosxPath + "geosxProcess.py -i " + \
                  xml + " -x " + str(x) + " -y " + str(y) + " -z " + str(z) + " " + shot_file + " " + traceProcPath)
    else:
        fname = "process"+str(mp.current_process()._identity)[1]+".sh"
        fsh = open(fname, 'w+')
        fsh.write("#!/usr/bin/bash \n")
        fsh.write("#SBATCH --job-name=pygeosx_seismic_acquisition \n")
        fsh.write("#SBATCH -n " + str(mpiranks) + " \n")
        fsh.write("#SBATCH -N " + str(nodes) + " \n")
        if nodeName != None:
            fsh.write("#SBATCH -C " + nodeName + " \n")
        fsh.write("cd /beegfs/jbesset/codes/GEOSX/build-test_module-release/lib/PYGEOSX/bin/ \n")
        fsh.write("source activate \n")
        fsh.write("cd /beegfs/jbesset/codes/GEOSX/src/coreComponents/physicsSolvers/wavePropagation/pygeosx/ \n")
        fsh.write("module load physics/geosx_deps \n")
        fsh.write("module load physics/pygeosx \n")
        fsh.write("mpirun -np " + str(mpiranks) + " --map-by node python " + pygeosxPath + "/geosxProcess.py -i " + \
                  xml + " -x " + str(x) + " -y " + str(y) + " -z " + str(z) + " " + shot_file + " " + traceProcPath)
        fsh.close()
        os.system("/usr/bin/sbatch " + fname )



def daskProcessing(shot_list, GEOSX, xml,
                   memory = None,
                   cores = 1,
                   processes = 1,
                   nodes = 1,
                   x = 1, y = 1, z = 1,
                   cluster_type = "local"):

    """Cut shot list into desired amount of processes and run each process in parrallel
       each on a single node on a cluster using the Dask library

    Parameters
    ----------
    shot_list : list of Shot object
        list of shot configuration

    GEOSX : string
        location of GEOSX directory

    xml : string
        Time step

    processes : int
        Number of processes

    cores : int
        Number of cores for parrallel solving

    processes : int
        Number of processes per node

    nodes : int
        Number of nodes to use

    x : int
        Number of MPI partition along x axe (domain partition)

    y : int
        Number of MPI partition along x axe (domain partition)

    z : int
        Number of MPI partition along x axe (domain partition)

    cluster_type : sting
        Type of cluster (SLURM)
    """

    futures = []
    nb_shot_m1 = len(shot_list)
    ind = 0

    wavePropagationPath = os.path.join(GEOSX, "GEOSX/src/coreComponents/physicsSolvers/wavePropagation/")
    pygeosxPath = os.path.join(wavePropagationPath, "pygeosx/")

    if cluster_type == "SLURM":
        cluster = SLURMcluster(cores = cores, memory = memory, processes = processes, job_cpu = cores)
        cluster.scale(nodes)
    elif cluster_type == "local":
        cluster = None

    client = Client(cluster)

    for i in range(nodes * processes):
        nb_shot = int(nb_shot_m1/(nodes*processes - i))

        shot_file = exportShotList(i, shot_list[ind:ind + nb_shot])

        tracePath = os.path.join(rootPath, "outputSismoTrace/")
        if os.path.exists(tracePath):
            pass
        else:
            os.mkdir(tracePath)

        traceProcPath = os.path.join(tracePath, "traceProc" + str(i) + "/")
        if os.path.exists(traceProcPath):
            pass
        else:
            os.mkdir(traceProcPath)

        fsh = open("bash"+str(i)+".sh", 'w+')
        fsh.write("#!/bin/bash \n")
        fsh.write("cd " + GEOSX + "/GEOSX/build-test_module-release/lib/PYGEOSX/bin/ \n")
        if cluster_type == "SLURM":
            fsh.write("module load physics/pygeosx \n")
            fsh.write("module load physics/geosx_deps \n")
        fsh.write("source activate \n")
        fsh.write("cd " + pygeosxPath + " \n")
        if cores == 1:
            fsh.write("python geosxProcess.py -i " + xml + " " + shot_file + " " + traceProcPath)
        else:
            fsh.write("mpirun -np cores --oversubscribe python geosxProcess.py -i " + xml + " -x " + str(x) + " -y " + str(y) + " -z " + str(z) + " " + shot_file + " " + traceProcPath)
        fsh.close()

        futures.append(client.submit(os.system, "bash bash"+str(i)+".sh", key = "proc"+str(i)))
        ind = ind + nb_shot
        nb_shot_m1 = nb_shot_m1 - nb_shot

    for future, result in as_completed(futures, with_results=True):
        print(future, result)

    remove_shots_files()
    remove_bash_files()
    remove_dask_files()
    client.close()



def firstInit(GEOSX, xml,
              mpiranks = None,
              x = 1, y = 1, z = 1):

    """Function to initialize GEOSX and get optimal time step dt calculate for mesh, and also
       the max time of simulation in XML

    Parameters
    ----------

    GEOSX : string
        location of GEOSX directory

    xml : string
        xml file

    mpiranks : int
        Number of MPI processes for parrallel solving

    x : int
        Number of MPI partition along x axe (domain partition)

    y : int
        Number of MPI partition along x axe (domain partition)

    z : int
        Number of MPI partition along x axe (domain partition)


    Return
    ------
    maxT : float
        Max time for GEOSX simulation

    dt : float
        Time step
    """

    wavePropagationPath = os.path.join(GEOSX, "GEOSX/src/coreComponents/physicsSolvers/wavePropagation/")
    pygeosxPath = os.path.join(wavePropagationPath, "pygeosx/")

    if mpiranks == None:
        os.system("python " + pygeosxPath + "firstInit.py -i " + xml)
    else:
        os.system("mpirun -np " + str(mpiranks) + " python firstInit.py -i " + xml + " -x " + str(x) + " -y " + str(y) + " -z " + str(z))

    maxT, dt = readInitVariable()

    return maxT, dt
