import pygeosx
from mpi4py import MPI
import sys

from AcousticShot import *
from ElasticShot import *
from fileManager import *

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def main():
    if len(sys.argv) == 5:
        shot_file = sys.argv[3]
        tracePath = sys.argv[4]

        shot_list = readShotList(shot_file)

        shot_simul(rank, sys.argv[0:3], shot_list, tracePath)

    elif len(sys.argv) == 7:
        shot_file = sys.argv[5]
        tracePath = sys.argv[6]

        shot_list = readShotList(shot_file)

        shot_simul(rank, sys.argv[0:5], shot_list, tracePath)

    elif len(sys.argv) == 9:
        shot_file = sys.argv[7]
        tracePath = sys.argv[8]

        shot_list = readShotList(shot_file)

        shot_simul(rank, sys.argv[0:7], shot_list, tracePath)

    elif len(sys.argv) == 11:
        shot_file = sys.argv[9]
        tracePath = sys.argv[10]

        shot_list = readShotList(shot_file)

        shot_simul(rank, sys.argv[0:9], shot_list, tracePath)



def shot_simul(rank, args, shot_list, tracePath):
    """
    Parameters
    ----------
    rank : int
        MPI rank

    xml : string
        XML file

    shot_list : list
        A list containing sets of Source and ReceiverSet objects

    tracePath : string
        Directory for .sgy files output
    """

    problem = pygeosx.initialize(rank, args)
    pygeosx.apply_initial_conditions()

    solver = str(problem.get_group("Solvers").groups()[0]).split("::")[1][:-2]

    if solver == "AcousticWaveEquationSEM":
        acoustic_shots(rank, problem, shot_list, tracePath)
    elif solver == "ElasticWaveEquationSEM":
        elastic_shots(rank, problem, shot_list, tracePath)


if __name__ == "__main__":
    main()
