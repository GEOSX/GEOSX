'''
Created on 10/02/2021

@author: macpro
'''
import numpy as np

class Receiver:
    """A class representing a receiver

    Attributes
    ----------
    coords :
        Coordinates of the source
    """

    def __init__(self, pos):
        """Constructor for the receiver

        Parameters
        ----------
        pos :
            Coordinates for the receiver
        """

        self.coords = np.array([pos[0], pos[1], pos[2]])

    def __repr__(self):
        return '('+str(self.coords[0])+','+str(self.coords[1])+','+str(self.coords[2])+') \n'

    def getCoord(self):
        return self.coords

    def setReceiverPos(self, coord, value):
        if coord==0:
            self.coords[0] = value
        elif coord==1:
            self.coords[1] = value
        elif coord==2:
            self.coords[2] = value

    #Transpose the receiver position
    def linearMove(self, coord, value):
        if coord==0:
            self.coords[0] += value
        elif coord==1:
            self.coords[1] += value
        elif coord==2:
            self.coords[2] += value

    def x(self):
        return self.coords[0]

    def y(self):
        return self.coords[1]

    def z(self):
        return self.coords[2]



class ReceiverSet:
    """A class representing a set receiver

    Attributes
    ----------
    receiver_list :
        List of Receiver

    n :
        Number of Receiver
    """

    def __init__(self, receiver_list = None):
        """Constructor for the point source

        Parameters
        ----------
        receiver_list :
            List of Receiver
        """
        if receiver_list is None:
            self.receiver_list = []
            self.n = 0
        else:
            self.receiver_list = receiver_list
            self.n = len(receiver_list)

    def __repr__(self):
        if self.n >=10:
            return str(self.receiver_list[0:4])[:-1] + '...' + '\n' + str(self.receiver_list[-4:])[1:]
        else:
            return str(self.receiver_list)

    def getInsideDomain(self, meshBoundaries):
        list=[]

        for r in self.receiver_list:
            if meshBoundaries[0][0] <= r.coords[0] and r.coords[0] <= meshBoundaries[0][1] \
            and meshBoundaries[1][0] <= r.coords[1] and r.coords[1] <= meshBoundaries[1][1] \
            and meshBoundaries[2][0] <= r.coords[2] and r.coords[2] <= meshBoundaries[2][1]:
                list.append(r)

        final_receiver_set = ReceiverSet(list)
        return final_receiver_set

    def getList(self):
        return self.receiver_list

    def getNumberOfReceivers(self):
        return self.n

    def getSetCoord(self):
        listRcv = []
        for i in range(self.n):
            listRcv.append(np.array([self.receiver_list[i].coords[0], self.receiver_list[i].coords[1], self.receiver_list[i].coords[2]]))
        return listRcv

    def append(self, ReceiverSet):
        for i in range(ReceiverSet.getNumberOfReceivers()):
            self.receiver_list.append(ReceiverSet.getList()[i])
        self.n += ReceiverSet.getNumberOfReceivers()


    #Transpose all receivers position
    def linearSetMove(self, coord, value):
        for i in range(self.n):
            self.receiver_list[i].linearMove(coord, value)
