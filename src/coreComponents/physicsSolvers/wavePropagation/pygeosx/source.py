'''
Created on 10/02/2021

@author: macpro
'''
import numpy as np
import math as m


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



class Source:
    """A class representing a point source

    Attributes

    ----------
    coords :
        Coordinates of the source
    f :
        Source function (Ricker)

    domainFlag :
        0 the source is not in the domain,
        1 the source is in the domain
    """
    def __init__(self, pos, wavelet, dt):
        """Constructor for the point source

        Parameters
        ----------
        meshBox :
            Numpy array containing min/max boundary coordinates of the domain

        pos :
            Coordinates for the point source

        wavelet :
            Source function (Ricker)
        """

        self.coords = np.array([pos[0], pos[1], pos[2]])
        self.f = wavelet #Ricker
        self.dt = dt

    def __repr__(self):
        return '('+str(self.coords[0])+','+str(self.coords[1])+','+str(self.coords[2])+')'

    def getCoord(self):
        return self.coords

    def getDomainFlag(self):
        return self.domainFlag

    def getFunction(self):
        return self.f

    def getTimeStep(self):
        return self.dt

    def setTimeStep(self, dt):
        self.dt = dt


    def x(self):
        return self.coords[0]

    def y(self):
        return self.coords[1]

    def z(self):
        return self.coords[2]


class SourceSet:
    """A class representing a point source

    Attributes
    ----------
    source_list :
        List of Source object
    n :
        Number of Source
    """

    def __init__(self, source_list=None):
       if source_list is None:
           self.source_list = []
           self.n = 0
       else:
           self.source_list = source_list
           self.n = len(source_list) #Number of sources


    def append(self, source):
       self.source_list.append(source)
