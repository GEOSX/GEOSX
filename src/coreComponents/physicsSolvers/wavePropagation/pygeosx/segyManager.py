import segyio
import os

def create_segy(shot_list, physicalName, nsamples, tracePath):

    if os.path.exists(tracePath):
        pass
    else:
        os.mkdir(tracePath)

    for ishot in range(len(shot_list)):
        rcvCoord = shot_list[ishot].getReceiverSet().getSetCoord()
        srcCoord = shot_list[ishot].getSource().getCoord()
        spec = segyio.spec()
        ilines  = [x[0] for x in rcvCoord]
        xlines  = [x[1] for x in rcvCoord]
        spec.ilines = list(set(ilines))
        spec.xlines = list(set(xlines))
        spec.samples = list(range(nsamples))
        spec.sorting = 2
        spec.format  = 1

        fileName = fileName = physicalName + "_Shot"+ str(ishot) + ".sgy"
        with segyio.create(os.path.join(tracePath, fileName), spec) as f:
            for i in range(len(rcvCoord)):
                f.header[i] = {segyio.su.scalco : -100,
                               segyio.su.scalel : -100,
                               segyio.su.sx : int(srcCoord[0]*100),
                               segyio.su.sy : int(srcCoord[1]*100),
                               segyio.su.sdepth : int(srcCoord[2]*100),
                               segyio.su.gx : int(rcvCoord[i][0]*100),
                               segyio.su.gy : int(rcvCoord[i][1]*100),
                               segyio.su.gelev : int(rcvCoord[i][2]*100)}
                f.trace[i] = [0.0] * nsamples
            f.bin.update(tsort = segyio.TraceSortingFormat.INLINE_SORTING)


def export_to_segy(physicalValue, rcvCoord, segyFile):
    """Export the pressure value calculated by GEOSX to a segy file

    Parameters
    ----------
    physicalVAlue :
        Numpy array with physical values at all receivers coordinates

    rcvCoord :
        Receiver coordinates

    segyFile :
        File to which we export the values
    """

    print(physicalValue[:,51])
    with segyio.open(segyFile, 'r+', ignore_geometry=True) as f:
        for i in range(len(rcvCoord)):
            if any(physicalValue[1:,i])==True:
                f.trace[i] = physicalValue[:, i]


def export_for_acquisition(shot_list, acq_name):
    segyPath = os.path.join(os.path.abspath(os.getcwd()),"segyAcquisition/")
    if os.path.exists(segyPath):
        pass
    else:
        os.mkdir(segyPath)

    if os.path.exists(os.path.join(segyPath, acq_name)):
        pass
    else:
        os.mkdir(os.path.join(segyath, acq_name))

    for j in range(len(shot_list)):
        srcCoord = shot_list[j].getSource().getCoord()
        rcvCoord = shot_list[j].getReceiverSet().getSetCoord()

        spec = segyio.spec()
        ilines  = [x[0] for x in rcvCoord]
        xlines  = [x[1] for x in rcvCoord]
        spec.ilines = list(set(ilines))
        spec.xlines = list(set(xlines))
        spec.samples = [0]
        spec.sorting = 2
        spec.format  = 1

        with segyio.create(os.path.join(segyPath, acq_name) + "/acquisition_" + str(j) + ".sgy", spec) as f:
            for i in range(len(rcvCoord)):
                f.header[i] = {segyio.su.scalco : -100,
                               segyio.su.scalel : -100,
                               segyio.su.sx : int(srcCoord[0]*100),
                               segyio.su.sy : int(srcCoord[1]*100),
                               segyio.su.sdepth : int(srcCoord[2]*100),
                               segyio.su.gx : int(rcvCoord[i][0]*100),
                               segyio.su.gy : int(rcvCoord[i][1]*100),
                               segyio.su.gelev : int(rcvCoord[i][2]*100)}
                f.bin.update(tsort = segyio.TraceSortingFormat.INLINE_SORTING)
