import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
import xml.etree.ElementTree as ElementTree
from mpmath import *
import math
import argparse

class Sneddon:

    def __init__(self, mechanicalParameters, length, pressure):
        K = mechanicalParameters["bulkModulus"]
        G = mechanicalParameters["shearModulus"]

        E = (9 * K * G) / (3*K+G)
        nu = E / (2 * G) - 1
        #
        print("young modulus = ", E, " Pa")
        print("poisson ratio = ", nu)
        print("fracture length = ", 2*length, " m")

        self.scaling = ( 4 * (1 - nu**2) ) * pressure / E;
        self.halfLength = length;

    def computeAperture(self, x):
        return self.scaling * ( self.halfLength**2  - x**2 )**0.5;


def getMechanicalParametersFromXML( xmlFilePath ):
    tree = ElementTree.parse(xmlFilePath)

    param = tree.find('Constitutive/ElasticIsotropic')

    mechanicalParameters = dict.fromkeys(["bulkModulus", "shearModulus"])

    mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
    mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))

    return mechanicalParameters

def getFracturePressureFromXML( xmlFilePath ):
    tree = ElementTree.parse(xmlFilePath)

    param = tree.findall('FieldSpecifications/FieldSpecification')

    found_traction = False
    for elem in param:
        if elem.get("fieldName") == "fractureTraction" and elem.get("component") == "0":
            pressure = float(elem.get("scale"))*(-1)
            found_traction = True

        if found_traction: break

    return pressure


def getFractureLengthFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)
    boundedPlane = tree.find('Geometry/BoundedPlane')
    dimensions = boundedPlane.get("dimensions")
    dimensions = [float(i) for i in dimensions[1:-1].split(",")]
    length = dimensions[0] / 2

    origin = boundedPlane.get("origin")
    origin = [float(i) for i in origin[1:-1].split(",")]
    return length, origin[0]

def main(filesPaths):
    # File path
    hdf5File1Path = filesPaths[0]
    xmlFilePath   = filesPaths[1]

    # Read HDF5
    hf = h5py.File(hdf5File1Path, 'r')
    jump = hf.get('displacementJump')
    jump = np.array(jump)
    aperture = jump[0,:,0]

    hf = h5py.File(hdf5File1Path, 'r')
    x = hf.get('displacementJump elementCenter')
    x = x[0,:,0]

    # Filter out extra entries in the hdf5 file. It is just to make the plot look nicer
    voidIndexes = np.asarray( np.where(x == 0) )
    if voidIndexes.size !=0:
        lastValue = voidIndexes[0][0]
        aperture = aperture[0:lastValue]
        x = x[0:lastValue]

    # Extract info from XML
    mechanicalParameters = getMechanicalParametersFromXML(xmlFilePath)
    appliedPressure = getFracturePressureFromXML(xmlFilePath)

    # Get length of the fracture
    length, origin = getFractureLengthFromXML(xmlFilePath)

    x = x - origin
    # Initialize Sneddon's analytical solution
    sneddonAnalyticalSolution = Sneddon(mechanicalParameters, length, appliedPressure)

    # Plot analytical (continuous line) and numerical (markers) aperture solution
    x_analytical = np.linspace(-length, length, 101, endpoint=True)
    aperture_analytical = np.empty(len(x_analytical))

    apertureMaxAnalytical =  sneddonAnalyticalSolution.computeAperture( 0 )

    print("max aperture = ", apertureMaxAnalytical, " m")

    cmap = plt.get_cmap("tab10")
    i=0
    for xCell in x_analytical:
        aperture_analytical[i] = sneddonAnalyticalSolution.computeAperture( xCell )
        i += 1
    plt.plot(x_analytical, aperture_analytical, color=cmap(-1), label='analytical solution')
    plt.plot(x, aperture, 'o', color=cmap(2), label='numerical solution')

    plt.grid()
    plt.xlabel('length [m]')
    plt.ylabel('aperture [m]')
    plt.legend(bbox_to_anchor=(0.5, 0.2), loc='center', borderaxespad=0.)
    plt.show()

def parseArguments():
    parser = argparse.ArgumentParser(description="Sneddon comparison with analytical solution")

    parser.add_argument("-dp",
                        "--datapath",
                        type=str,
                        default="",
                        help="path to the hdf5 files with the data.")

    parser.add_argument("-xp",
                        "--xmlpath",
                        type=str, default="",
                        help="path to xml file.")

    args, unknown_args = parser.parse_known_args()

    filesPaths = []
    if args.datapath != "":
        # get path to the hdf5 files
        filesPaths.append(args.datapath + "displacementJump_history.hdf5")
    else:
        print("The path to the hdf5 files must be specified. Use -dp to specify it.")

    if args.xmlpath != "":
        filesPaths.append(args.xmlpath + "Sneddon-Validation.xml")
    else:
        print("The path to the xml file must be specified. Use -xp to specify it.")

    return filesPaths


if __name__ == "__main__":
    filesPaths = parseArguments()
    main(filesPaths)

