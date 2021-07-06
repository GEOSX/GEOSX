import ast

import argparse
from shot import *
from source import *
from receiver import *

import os

rootPath = os.path.abspath(os.getcwd())


def exportShotList(proc, shot_list):
    """Export shot configuration corresponding to one process to .txt file

    Parameters
    ----------
    proc : string
        Name of the process

    shot_list : list of Shot object
        Part of the full shot list
    """

    path = os.path.join(rootPath, "shots_lists")
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)
    shot_file = os.path.join(path, "shot_list" + str(proc) + ".txt")
    f = open(shot_file, 'w+')
    dt = shot_list[0].getSource().getTimeStep()
    wavelet = shot_list[0].getSource().getFunction()
    f.write(str(dt) + '\n')
    f.write(str(wavelet.tolist()) +'\n')
    rcvlist = []
    for shot in (shot_list):
        src = shot.getSource().getCoord().tolist()
        rcvs = shot.getReceiverSet().getSetCoord()
        for rcv in rcvs:
            rcvlist.append(rcv.tolist())

        f.write("{}\v{}\n".format(str(src), str(rcvlist)))
        rcvlist = []

    f.close()
    return shot_file

def readShotList(shot_file):
    """Read .txt file containing a list of shot configuration

    Parameters
    ----------
    shot_file : string
        .txt file containing a shot list


    Return
    ------
    shot_list : list of Shot object
        List of shots configuration
    """

    f = open(shot_file, 'r')
    dt = float(f.readline())
    wavelet = ast.literal_eval(f.readline())
    shot_list = []
    for line in f.readlines():
        src = ast.literal_eval(line.splitlines()[0])
        rcvs = ast.literal_eval(line.splitlines()[1])
        shot_list.append(Shot(Source(src, wavelet, dt), ReceiverSet([Receiver(rcv) for rcv in(rcvs)])))

    f.close()
    return shot_list

def exportInitVariable(maxT, dt):
    """Export initial parameters maxTime and time step dt

    Parameters
    ----------
    maxT : float
        Max time of GEOSX simulation

    dt : float
        Time step
    """

    f = open(rootPath + "/init_variable.txt", 'w+')
    f.write(str(maxT) + "\n")
    f.write(str(dt) + "\n")

    f.close()

def readInitVariable():
    """Get maxTime and time step dt from .txt file

    Return
    ------
    maxT : float
        Max time of GEOSX simulation

    dt : float
        Time step
    """

    f = open(rootPath + "/init_variable.txt", 'r')
    maxT = float(f.readline())
    dt = float(f.readline())

    f.close()
    os.remove(os.path.join(rootPath, "init_variable.txt"))
    return maxT, dt


def parse_args():
    """Get parameters from python command and put them into args variable

    Return
    ------
    args : list
        Variable containing command parameters
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--xml", help = "input xml file", required = True)
    parser.add_argument("--geosx", help = "GEOSX folder location", required = True)
    args = parser.parse_args()

    return args


def remove_shots_files():
    """Remove temporary shot_files that were created for each processes"""

    for root, dir, files in os.walk(os.path.join(rootPath, "shots_lists")):
        for file in files:
            os.remove(os.path.join(root, file))

    os.rmdir(os.path.join(rootPath, "shots_lists"))

def remove_bash_files():
    """Remove temporary bash_files that were created for each processes"""

    bashFiles = glob.glob(os.path.join(rootPath, "bash*"))
    for file in bashFiles:
        os.remove(file)

def remove_dask_files():
    """Remove temporary dask_files that were created for each processes"""

    for root, dir, files in os.walk(os.path.join(rootPath, "dask-worker-space")):
        for file in files:
            os.remove(os.path.join(root, file))

    os.rmdir(os.path.join(rootPath, "dask-worker-space"))
