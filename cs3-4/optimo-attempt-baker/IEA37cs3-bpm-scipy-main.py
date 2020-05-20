#-- Necessary Headers --#
from __future__ import print_function   # For Python 3 compatibility
from scipy.interpolate import interp1d  # To create our splines
from scipy import optimize  # To create our splines
from os import path
import scipy.spatial.distance as ssd
import matplotlib.pyplot as plt
import numpy as np
import sys
import yaml                             # For reading .yaml files
import time
# For the AEP calculation code
import baker_cs34_functions as Iea37sb
import iea37_aepcalc as iea37aepC
from scipy.special import binom         # "Combination", for deterimining unique turbine pairs
from math import radians as DegToRad    # For converting degrees to radians
from math import log as ln  # For natural logrithm

"""
This is the "main" function to run an optimization for IEA37's cs3 using:
    Boundry constraint: Boundary partition method, (full boundary, if clauses)
    Optimizer: scipy.minimize()
"""
if __name__ == "__main__":
    scaledAEP = 1#e5
    scaledTC = 1#e3
    strCase = 'cs4'  # Which case study we're doing. 'cs3' or 'cs4'
    numGridLines = 10                   # How many gridlines we'll use for the splining
    nNumRegions = 5                     # Number of reigons we're using (cs4 = 5, cs3 = 1)

    #- Rip the boundary coordinates from the .yaml file -# 
    fn = "../startup-files/iea37-boundary-" + strCase + ".yaml"
    [coordList3a, coordList3b, coordList4a, coordList4b, coordList4c] = Iea37sb.getTurbAtrbtCs4YAML(fn)
    clsdBP3a = Iea37sb.closeBndryList(coordList3a)    # Duplicate the first coordinate for a closed boundary
    coordList3b = np.roll(coordList3b, 1)             # Shift our points to the right so rightmost vertex is zero
    clsdBP3b = Iea37sb.closeBndryList(coordList3b)    # Duplicate the first coordinate for a closed boundary
    coordList4a = np.roll(coordList4a, -3)            # Shift our points to the left so rightmost vertex is zero
    clsdBP4a = Iea37sb.closeBndryList(coordList4a)    # Duplicate the first coordinate for a closed boundary
    clsdBP4b = Iea37sb.closeBndryList(coordList4b)    # Duplicate the first coordinate for a closed boundary
    clsdBP4c = Iea37sb.closeBndryList(coordList4c)    # Duplicate the first coordinate for a closed boundary

    #-- Spline all Boundaries --#
    vertexList3a = [0, 6, 8, 9, 18]
    numSides3a = len(vertexList3a) - 1      # The number of sides for our original coordinate system.
    splineList3a, segCoordList3a = Iea37sb.makeCs3BndrySplines(vertexList3a, clsdBP3a, numGridLines)
    vertexList3b = [0, 1, 2, 3, 8]       # Hard code the vertices (though this could be done algorithmically)
    numSides3b = len(vertexList3b) - 1
    splineList3b, segCoordList3b = Iea37sb.makeCs3BndrySplines(vertexList3b, clsdBP3b, numGridLines)
    vertexList4a = [0, 1, 2, 3, 6]       # Hard code the vertices (though this could be done algorithmically)
    numSides4a = len(vertexList4a) - 1      # The number of sides for our original coordinate system.
    splineList4a, segCoordList4a = Iea37sb.makeCs3BndrySplines(vertexList4a, clsdBP4a, numGridLines)
    vertexList4b = [0, 1, 2, 3]       # Hard code the vertices (though this could be done algorithmically)
    numSides4b = len(vertexList4b) - 1      # The number of sides for our original coordinate system.
    splineList4b, segCoordList4b = Iea37sb.makeCs3BndrySplines(vertexList4b, clsdBP4b, numGridLines)
    vertexList4c = [0, 1, 4, 5]       # Hard code the vertices (though this could be done algorithmically)
    numSides4c = len(vertexList4c) - 1      # The number of sides for our original coordinate system.
    splineList4c, segCoordList4c = Iea37sb.makeCs3BndrySplines(vertexList4c, clsdBP4c, numGridLines)

    #clsdBP = Iea37sb.closeBndryList(bndryPts)   # Duplicate the 1st coord for a closed boundary
    #- Load the turbine and windrose atributes -#
    fname_turb = "../startup-files/iea37-10mw.yaml"
    fname_wr = "../startup-files/iea37-windrose-cs3.yaml"
    wind_dir, wind_dir_freq, wind_speeds, wind_speed_probs, num_speed_bins, min_speed, max_speed = iea37aepC.getWindRoseYAML(
        fname_wr)
    turb_ci, turb_co, rated_ws, rated_pwr, turb_diam = iea37aepC.getTurbAtrbtYAML(
        fname_turb)
    turb_diam = turb_diam/scaledTC

    #- Some necessary variables -#
    numSides = [numSides3a, numSides3b, numSides4a, numSides4b, numSides4c]  # Collate all the sides
    fMinTurbDist = (turb_diam * 2)
    #- args in the correct format for optimization -#
    Args = dict([('wind_dir_freq', wind_dir_freq), \
                ('wind_speeds', wind_speeds), \
                ('wind_speed_probs', wind_speed_probs), \
                ('wind_dir', wind_dir), \
                ('turb_diam', turb_diam), \
                ('turb_ci', turb_ci), \
                ('turb_co', turb_co), \
                ('rated_ws', rated_ws), \
                ('rated_pwr', rated_pwr), \
                ('fAEPscale', scaledAEP), \
                ('fTCscale', scaledTC),
                ('fMinTurbDist', fMinTurbDist)])
    # Make the matrix to pass in to the function
    splineMatDict = {'3a':splineList3a, '3b':splineList3b, '4a':splineList4a, '4b':splineList4b, '4c':splineList4c}
    coordsCornersDict = {'3a':clsdBP3a[vertexList3a[0:4]], '3b':clsdBP3b[vertexList3b[0:4]], '4a':clsdBP4a[vertexList4a[0:4]], '4b':clsdBP4b[vertexList4b[0:3]], '4c':clsdBP4c[vertexList4c[0:3]]}

    # Make an array of the number of turbines in each region
    nRegionNumTurbs = np.zeros(nNumRegions, dtype=np.int32)     # Our number of turbines per each region
    for i in range(nNumRegions):
        nRegionNumTurbs[i] = Iea37sb.cs34Regions().getNumTurbs(Iea37sb.cs34Regions().getRegionName(i))

    numRestarts = 3                     # Number of restarts we're doing
    print("Running region: " + strCase)
    print("Running: " + str(numRestarts) + " restarts.")

    #- Initialize variables to hold results -#
    listAEP = np.zeros(numRestarts)     # An array holding our AEP values
    listTurbLocs = np.zeros((numRestarts, (np.sum(nRegionNumTurbs)*2)))
    bestResult = np.zeros(2)            # Index, AEP number
    timeArray = np.zeros(numRestarts)   # To hold timing information

    #- Loop for every restart -#
    for cntr in range(numRestarts):
        print("Restart #" + str(cntr+1) + "/" +  str(numRestarts) + " (index " + str(cntr) + ")")
        timeStart = time.time() # Start the clock
        #-- Use our pregenerated turbine locations --#
        x0l = []                        # Initialize our turbine <coord> list
        for i in range(nNumRegions):    # Loop through our regions
            PreStarts = np.loadtxt('./results/randostarts-' + Iea37sb.cs34Regions().getRegionName(i) + '-200.csv', delimiter=',')
            x0l.extend(Iea37sb.makeArrayCoord(PreStarts[cntr]))
        x0  = Iea37sb.makeCoordListArray(x0l)
        x0s = Iea37sb.makeArrayCoord(x0)

        #- Get our turbine list ready for processing -#
        #x0 = Iea37sb.makeCoordArray(x0s)#/ Args['fTCscale']         # Get a random turbine placement and scale it
        startAEP = Iea37sb.optimoFun(x0, Args)
        print("Start AEP = " + str(startAEP*scaledAEP))#*Args['fAEPscale']))

        cons = ({'type': 'ineq', 'fun': lambda x:  Iea37sb.checkBndryConsCs4(x, nRegionNumTurbs, splineMatDict, coordsCornersDict)}, {
                'type': 'ineq', 'fun': lambda x: Iea37sb.checkTurbSpacingCs4(x, nRegionNumTurbs, fMinTurbDist)})
        res = optimize.minimize(Iea37sb.optimoFun, x0, args=Args, method='SLSQP',
                                constraints=cons, options={'disp': True, 'maxiter': 1000}, tol=1e-7)

        #-- Save our results (scaling back up to normal) --#
        listAEP[cntr] =  (Iea37sb.optimoFun(res.x* Args['fTCscale'], Args) * Args['fAEPscale']) #(Iea37sb.optimoFun(res.x, Args))
        listTurbLocs[cntr] = res.x * Args['fTCscale']
        timeEnd = time.time()
        timeArray[cntr] = timeEnd - timeStart
        print("End AEP = " + str(listAEP[cntr]))
        print("Elapesed time: " + str(timeArray[cntr]))
        print()

    # Find the best result and save it    
    for j in range(numRestarts):
        if (bestResult[1] > listAEP[j]):  # If our new AEP is better (Remember negative switches)
            bestResult[1] = listAEP[j]    # Save it
            bestResult[0] = j             # And the index of which run we're on
    print("Best run: " + str(int(bestResult[0])+1))
    print("Best AEP: " + str(listAEP[int(bestResult[0])]))  # Print the best one, have to make sure the index is an (int)
    percentImprovement = ((listAEP[int(bestResult[0])] - startAEP)/startAEP) *100
    print("Improvement: " + str(percentImprovement))

    bestTurbs =  Iea37sb.makeArrayCoord(listTurbLocs[int(bestResult[0])]) #/scaledTC)
    print("Best result was index: " + str(int(bestResult[0])) )
    # If we've already saved this kind of run, give it a new name
    for i in range(100):
        if(path.exists('./results/turblocs-bpm-' + str(numRestarts) + 'run-' + strCase + '-(' + str(i) + ').csv') == False):
            np.savetxt('./results/turblocs-bpm-' + str(numRestarts) + 'run-' + strCase + '-(' + str(i) + ').csv', listTurbLocs, delimiter=',')
            break
    
    for i in range(100):
        if(path.exists('./results/turblocs-bpm-' + str(numRestarts) + 'run-' + strCase + '-(' + str(i) + ')-time.csv') == False):
            np.savetxt('./results/turblocs-bpm-' + str(numRestarts) + 'run-' + strCase + '-(' + str(i) + ')-time.csv', timeArray, delimiter=',')
            break

    #--- Print it all out and save the overlay ---#
    x0sTemp = []                        # Initialize our turbine <coord> list
    for i in range(nNumRegions):    # Loop through our regions
        PreStarts = np.loadtxt('./results/randostarts-' + Iea37sb.cs34Regions().getRegionName(i) + '-200.csv', delimiter=',')
        x0sTemp.extend(Iea37sb.makeArrayCoord(PreStarts[int(bestResult[0])]))
    x0sStart = Iea37sb.makeArrayCoord(Iea37sb.makeCoordListArray(x0sTemp))
    plt.figure("Start point")
    plt.figure(figsize=(20,10))
    plt.hold = True
    #- Plot all the boundaries -#
    Iea37sb.printBndryCoord(clsdBP3a, 5, 1)
    Iea37sb.printBndryCoord(clsdBP3b, 5, 1)
    Iea37sb.printBndryCoord(clsdBP4a, 5, 1)
    Iea37sb.printBndryCoord(clsdBP4b, 5, 1)
    Iea37sb.printBndryCoord(clsdBP4c, 5, 1)
    plt.plot(x0sStart.x*scaledTC, x0sStart.y*scaledTC, marker='o', color='black', linestyle='', markersize=7)
    plt.axis('scaled')                      # Trim the white space
    plt.axis('off')                         # Turn off the framing

    plt.figure("End point")
    plt.hold = True
    #- Plot all the boundaries -#
    Iea37sb.printBndryCoord(clsdBP3a, 5, 1)
    Iea37sb.printBndryCoord(clsdBP3b, 5, 1)
    Iea37sb.printBndryCoord(clsdBP4a, 5, 1)
    Iea37sb.printBndryCoord(clsdBP4b, 5, 1)
    Iea37sb.printBndryCoord(clsdBP4c, 5, 1)
    plt.plot(bestTurbs.x, bestTurbs.y, marker='o', color='black', linestyle='', markersize=7)
    plt.axis('scaled')                      # Trim the white space
    plt.axis('off')                         # Turn off the framing

    plt.figure("Overlay")
    plt.hold = True
    #- Plot all the boundaries -#
    Iea37sb.printBndryCoord(clsdBP3a, 5, 1)
    Iea37sb.printBndryCoord(clsdBP3b, 5, 1)
    Iea37sb.printBndryCoord(clsdBP4a, 5, 1)
    Iea37sb.printBndryCoord(clsdBP4b, 5, 1)
    Iea37sb.printBndryCoord(clsdBP4c, 5, 1)
    plt.plot(x0sStart.x*scaledTC, x0sStart.y*scaledTC, marker='o', color=Iea37sb.getPltClrs().getColor(5), linestyle='', markersize=7)
    plt.plot(bestTurbs.x, bestTurbs.y, marker='o', color=Iea37sb.getPltClrs().getColor(1), linestyle='', markersize=7)
    plt.axis('scaled')                      # Trim the white space
    plt.axis('off')                         # Turn off the framing
    plt.savefig('./results/' + strCase + '(' + str(i) + ').pdf')

    plt.show()