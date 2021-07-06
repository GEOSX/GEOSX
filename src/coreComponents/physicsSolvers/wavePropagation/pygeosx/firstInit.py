import pygeosx
from mpi4py import MPI
import sys

from mesh import calculDt, domainBoundary
from fileManager import exportInitVariable
from print import print_group

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def main():

    problem = pygeosx.initialize(rank, sys.argv)
    pygeosx.apply_initial_conditions()

    if rank==0:
        maxT = problem.get_wrapper("Events/maxTime").value()[0]
        dt = calculDt(problem)

        exportInitVariable(maxT, dt)


if __name__ == "__main__":
    main()
