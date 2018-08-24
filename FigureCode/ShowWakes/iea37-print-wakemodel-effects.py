import numpy as np
import sys
import yaml
import matplotlib.pyplot as plt # For Debug
from disk_grid_regular_count import disk_grid_regular_count # For validation
from disk_grid_regular import disk_grid_regular # for validation
from matplotlib import cm                   # For plotting colormap
from collections import OrderedDict         # For plotting colormap
# Written by Nicholas F. Baker to validate the IEA37 simplified Bastankhah Gaussian Wake Model
# BYU FLOW lab
# Completed 16 June 2018

def WindFrame(turbineX, turbineY, windDirectionDeg):
    """ Calculates the locations of each turbine in the wind direction reference frame """

    # convert from meteorological polar system (CW, 0 deg.=N) to standard polar system (CCW, 0 deg.=E)
    windDirectionDeg = 270. - windDirectionDeg
    if windDirectionDeg < 0.:
        windDirectionDeg += 360.
    windDirectionRad = np.pi*windDirectionDeg/180.0    # inflow wind direction in radians

    # convert to downwind(x)-crosswind(y) coordinates
    turbineXw = turbineX*np.cos(-windDirectionRad)-turbineY*np.sin(-windDirectionRad)
    turbineYw = turbineX*np.sin(-windDirectionRad)+turbineY*np.cos(-windDirectionRad)

    return turbineXw, turbineYw

def getTurbLocYAML(sFileName):
    turbineX = np.array([])
    turbineY = np.array([])

    # Read in the .yaml file
    with open(sFileName, 'r') as f:
        doc = yaml.load(f)

    # rip the x- and y-coordinates
    turbineX = np.asarray(doc['definitions']['position']['items']['xc']) # Convert from <list> to <ndarray>
    turbineY = np.asarray(doc['definitions']['position']['items']['yc'])
    # rip the expected AEP, used for comparison
    # AEP = doc['definitions']['plant_energy']['properties']['annual_energy_production']['default']

    return turbineX, turbineY#, AEP

def getWindFreqYAML(sFileName):
    windFreq = np.array([])

    # Read in the .yaml file
    with open(sFileName, 'r') as f:
        doc = yaml.load(f)

    # rip wind frequency distribution array
    windFreq = np.asarray(doc['definitions']['wind']['properties']['probability']['default']) # Convert from <list> to <ndarray>

    return windFreq

def getWindSpeedYAML(sFileName):

    # Read in the .yaml file
    with open(sFileName, 'r') as f:
        doc = yaml.load(f)

    # rip wind speed distribution array
    windSpeed = float(doc['definitions']['wind']['properties']['speed']['default']) # Convert from <list> to <float>

    return windSpeed

def getWindDirYAML(sFileName):
    windFreq = np.array([])

    # Read in the .yaml file
    with open(sFileName, 'r') as f:
        doc = yaml.load(f)

    # rip wind sampling bins array
    windFreq = np.asarray(doc['definitions']['wind']['properties']['direction']['bins']) # Convert from <list> to <ndarray>

    return windFreq

# --- BELOW FOR DEBUGGING --- #

def GWakeGrid(turbineXw, turbineYw, gridPointsXw, gridPointsYw):
    """ Returns total loss from wakes at each grid point"""
    # Equations and values explained in <iea37-wakemodel.pdf>
    nTurbines = len(turbineXw)      # Get the number of turbines
    nPoints = len(gridPointsXw)     # For the number for points we have

    CT = 4.0*1./3.*(1.0-1./3.)  # constant thrust coefficient
    k = 0.0324555   # constant turbulence

    D = 130.  # IEA37 3.35MW onshore reference turbine rotor diameter

    loss = np.zeros(nPoints) # store a loss for each grid point
    
    for i in range(nPoints):                # For every point we wish to evaluate at
        loss_array = np.zeros(nTurbines)    # Calculate the loss contribution from each turbine at that point
        for j in range(nTurbines):          # Loop through each turbine
            x = gridPointsXw[i]-turbineXw[j]
            y = gridPointsYw[i]-turbineYw[j]
            if x > 0.:  # Simplified Bastankhah Gaussian wake model, applied to downstream turbines
                sigma = k*(x)+D/np.sqrt(8.)
                loss_array[j] = (1.-np.sqrt(1.-CT/(8.*sigma**2/D**2)))*np.exp(-0.5*(y/sigma)**2)
            else:
                loss_array[j] = 0.
        # total wake loss, sqrt of sum of sqrs
        loss[i] = np.sqrt(np.sum(loss_array**2))

    return loss

def calcWindSpeeds(turbineX, turbineY, windFreq, windSpeed, windDir, testDirNum, gridPoints):
    """ Get the wind speeds at each point in the farm, for a single direction on windrose """
    testWindDir = np.array([1], dtype=float)
    testWindFreq = np.array([1], dtype=float)
    gridPointsX = gridPoints[:,0]
    gridPointsY = gridPoints[:,1]

    testWindDir[0] = windDir[testDirNum]
    testWindFreq[0] = windFreq[testDirNum]

    # turbines and evaluation points in wind frame coordinates
    turbineXw, turbineYw = WindFrame(turbineX, turbineY, testWindDir)
    gridPointsXw, gridPointsYw = WindFrame(gridPointsX, gridPointsY, testWindDir)
    
    # Evaluate loss due to wakes at all points
    loss = GWakeGrid(turbineXw, turbineYw, gridPointsXw, gridPointsYw)  # wake losses
    effWindSpeed = windSpeed*(1.-loss)
    effWindSpeed = effWindSpeed.reshape((-1, 1))
    # Store in 3rd column of speedMatrix
    speedMatrix = np.append(gridPoints, effWindSpeed, axis=1)

    return speedMatrix


def printWindSpeeds(turbineX, turbineY, windFreq, windSpeed, windDir, testDirNum, numTestDivs):
    # Make an array of all the points to test
    fieldRad = 1300.
    center = np.array([0.0, 0.0])
    ng = disk_grid_regular_count(numTestDivs, fieldRad, center) # Figures out the total number of points given the subdivisions (n)
    gridPoints = disk_grid_regular(numTestDivs, fieldRad, center, ng) # returns the points for the given number of subdivisions (n)
    gridPoints = np.transpose(gridPoints)

    # Calculate windspeeds at all points in array (speedMatrix), with third column = windspeeds.
    speedMatrix = calcWindSpeeds(turbineX, turbineY, windFreq, windSpeed, windDir, testDirNum, gridPoints)

    # Print Turbine Lcoations
    nNumTurb = 16  # 16 Turbine Farm
    rtrDiam = 130  # Diameter of NREL 3.35 MW turbine rotor
    colorNum = 0 # Should be 0-5 for different colors
    plotFarmSpeeds(turbineX, turbineY, rtrDiam, nNumTurb, fieldRad, colorNum, speedMatrix)

    return

def plotFarmSpeeds(fXcoords, fYcoords, rtrDiam, nNumTurb, fieldRad, cNum, speedMat):
    #color = (74./255., 145./255., 200./255.)  # For nice MATLAB blue
    color = np.empty(6, dtype='string')
    color[0] = 'blue'
    color[1] = 'red'
    color[2] = 'green'
    color[3] = 'yellow'
    color[4] = 'magenta'
    color[5] = 'cyan'

    plt.figure(1)
    # Print the windspeeds
    plt.scatter(speedMat[:, 0], speedMat[:, 1], c=speedMat[:, 2], marker='.', cmap='winter')

    # Print the Turbines
    for i in range(nNumTurb):
        circ_opt = plt.Circle((fXcoords[i]*1., fYcoords[i]*1.), rtrDiam/2., facecolor=color[cNum], edgecolor=color[cNum], alpha=0.5)
        plt.gca().add_patch(circ_opt)
    circ_outer = plt.Circle((0, 0), fieldRad, linestyle='dashed', edgecolor='black', facecolor='None', label='Boundaries')
    plt.gca().add_patch(circ_outer)

    plt.axis('equal')
    plt.axis('off')
    plt.title("16 Farm")
    plt.savefig("GaussianWakeField.pdf", transparent=True)
    plt.show()
    plt.gcf().clear()
# --- Above for debugging --- #

if __name__ == "__main__":
    """ Used for demonstration """
    turbineX = np.array([])
    turbineY = np.array([])

    # For Python .yaml capability, in the terminal type "pip install pyyaml".
    # An example command line syntax to run this file is "python iea37-aepcalc.py iea37-ex16.yaml iea37-windrose.yaml"

    # Read necessary values from .yaml files
    turbineX, turbineY = getTurbLocYAML(sys.argv[1])                # Get turbine locations from .yaml file
    windDir = getWindDirYAML(sys.argv[2])                           # Get the array wind sampling bins
    windFreq = getWindFreqYAML(sys.argv[2])                         # Get wind frequency distribution from .yaml file
    windSpeed = getWindSpeedYAML(sys.argv[2])                       # Get the wind speed from the .yaml file

    # --- Debug to print the different windspeeds on the farm --- #
    numTestPoints = 70  # How many points to sample across the farm
    testDirNum = 3 # 0-15, which wind direction bucket to test
    printWindSpeeds(turbineX, turbineY, windFreq, windSpeed, windDir, testDirNum,  numTestPoints)
    # --- End of Tests --- #
