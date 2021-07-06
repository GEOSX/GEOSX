'''
Created on 9/02/2021

@author: macpro
'''

from scheduler import firstInit, multiProcessing, daskProcessing
from seismicAcquisition import equispaced_acquisition
from source import ricker
from fileManager import parse_args

args = parse_args()

def main():

    #maxT, dt = firstInit(args.geosx, args.xml)
    maxT = 2.0
    dt = 0.002
    frequency = 5.0
    wavelet   = ricker(maxT, dt, frequency)


    shot_list = equispaced_acquisition([[0,2000],[0,2000],[0,2000]],
                                       wavelet,
                                       dt,
                                       start_source_pos    = [751, 1001],
                                       end_source_pos      = [1501, 1001],
                                       start_receivers_pos = [[21, 1001]],
                                       end_receivers_pos   = [[1981, 1001]],
                                       number_of_sources   = 2,
                                       number_of_receivers = [100],
                                       source_depth = 101,
                                       receivers_depth = 51
                                       )

    multiProcessing(shot_list, args.geosx, args.xml, mpiranks=4, x=2, y=2)
    #daskProcessing(shot_list, args.geosx, args.xml)

if __name__ == "__main__":
    main()
