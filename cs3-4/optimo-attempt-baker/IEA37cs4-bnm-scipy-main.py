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
    Boundry constraint: Boundary normal method, (partial boundary, convex)
    Optimizer: scipy.minimize()
"""
if __name__ == "__main__":
    numTurbs = 25
    scaledAEP = 1#e5
    scaledTC = 1#e3

    #- Load the boundary (with scaling) -#
    fn = "../startup-files/iea37-boundary-cs3.yaml"
    #bndryPts = Iea37sb.getTurbAtrbtCs3YAML(fn)  # Normal read
    #- Scaled read -# 
    tempPtsCoord = Iea37sb.getTurbAtrbtCs3YAML(fn)                  # Read in as <coord>
    tempPtsArray = Iea37sb.makeCoordArray(tempPtsCoord) / scaledTC  # Convert to an array to scale 
    bndryPts = Iea37sb.makeArrayCoord(tempPtsArray)                 # Convert back to <coord> type

    clsdBP = Iea37sb.closeBndryList(bndryPts)   # Duplicate the 1st coord for a closed boundary
    cncvBP = Iea37sb.makeSimpleCs3Bndry(clsdBP) # Make the boundary concave
    BndryNormals = Iea37sb.bndryNormals(cncvBP)

    #- Load the turbine and windrose atributes -#
    fname_turb = "../startup-files/iea37-10mw.yaml"
    fname_wr = "../startup-files/iea37-windrose-cs3.yaml"
    wind_dir, wind_dir_freq, wind_speeds, wind_speed_probs, num_speed_bins, min_speed, max_speed = iea37aepC.getWindRoseYAML(
        fname_wr)
    turb_ci, turb_co, rated_ws, rated_pwr, turb_diam = iea37aepC.getTurbAtrbtYAML(
        fname_turb)
    turb_diam = turb_diam/scaledTC 


    #- Some display variables -#
    numLinspace = 10
    numGridLines = 10                   # How many gridlines we'll use for the visualization
    vertexList = [0, 6, 8, 9, 18]       # Hard code the vertices (though this could be done algorithmically)
    numSides = len(vertexList) - 1
    coordsCorners = bndryPts[vertexList[0:4]] # Just the "corners" we've selected.
    fMinTurbDist = (turb_diam * 2)
#    fMinTurbDistScaled = fMinTurbDist / scaledTC
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

    # Spline up the boundary
    [splineList, segCoordList] = Iea37sb.makeCs3BndrySplines(vertexList, clsdBP, numGridLines)
    # Load the premade restarts
    PreStarts = np.loadtxt('./results/randostarts-3a-200.csv', delimiter=',')
    PreStarts = PreStarts / scaledTC

    numRestarts = 1                     # Number of restarts we're doing
    print("Running: " + str(numRestarts) + " restarts.")

    #- Initialize variables to hold results -#
    listAEP = np.zeros(numRestarts)     # An array holding our AEP values
    listTurbLocs = np.zeros((numRestarts, (numTurbs*2)))
    bestResult = np.zeros(2)            # Index, AEP number
    timeArray = np.zeros(numRestarts)   # To hold timing information

    #- Loop for every restart -#
    for cntr in range(numRestarts):
        print("Restart #" + str(cntr+1) + "/" +  str(numRestarts) + " (index " + str(cntr) + ")")
        timeStart = time.time() # Start the clock
        #-- Use our pregenerated turbine locations --#
        x0s = Iea37sb.makeArrayCoord(PreStarts[cntr])
        
        #- Get our turbine list ready for processing -#
        x0 = Iea37sb.makeCoordArray(x0s)#/ Args['fTCscale']         # Get a random turbine placement and scale it
        startAEP = Iea37sb.optimoFun(x0, Args)
        print("Start AEP = " + str(startAEP*scaledAEP))#*Args['fAEPscale']))

        cons = ({'type': 'ineq', 'fun': lambda x:  Iea37sb.calcDistNorms(x, coordsCorners, BndryNormals)}, {
                'type': 'ineq', 'fun': lambda x: Iea37sb.checkTurbSpacing(x, fMinTurbDist)})
        res = optimize.minimize(Iea37sb.optimoFun, x0, args=Args, method='SLSQP',
                                constraints=cons, options={'disp': True, 'maxiter': 1000}, tol=1e-10)

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
    print("Best run: " + str(int(bestResult[0])))
    print("Best AEP: " + str(listAEP[int(bestResult[0])]))  # Print the best one, have to make sure the index is an (int)
    percentImprovement = ((listAEP[int(bestResult[0])] - startAEP)/startAEP) *100
    print("Improvement: " + str(percentImprovement))

    bestTurbs =  Iea37sb.makeArrayCoord(listTurbLocs[int(bestResult[0])]) #/scaledTC)

    # If we've already saved this kind of run, give it a new name
    # for i in range(100):
    #     if(path.exists('./results/turblocs-bnm-' + str(numRestarts) + 'run(' + str(i) + ').csv') == False):
    #         np.savetxt('./results/turblocs-bnm-' + str(numRestarts) + 'run(' + str(i) + ').csv', listTurbLocs, delimiter=',')
    #         break
    
    # for i in range(100):
    #     if(path.exists('./results/turblocs-bnm-' + str(numRestarts) + 'run(' + str(i) + ')-time.csv') == False):
    #         np.savetxt('./results/turblocs-bnm-' + str(numRestarts) + 'run(' + str(i) + ')-time.csv', timeArray, delimiter=',')
    #         break

    x0Start = Iea37sb.makeArrayCoord(PreStarts[int(bestResult[0])])
    plt.figure("Start point")
    plt.hold = True
    plt.plot(clsdBP.x*scaledTC, clsdBP.y*scaledTC, color=Iea37sb.getPltClrs().getColor(5), linewidth=1)
    plt.plot(x0Start.x*scaledTC, x0Start.y*scaledTC, marker='o', color='black', linestyle='', markersize=7)
    plt.axis('scaled')                      # Trim the white space
    plt.axis('off')                         # Turn off the framing

    plt.figure("End point")
    plt.hold = True
    plt.plot(clsdBP.x*scaledTC, clsdBP.y*scaledTC, color=Iea37sb.getPltClrs().getColor(5), linewidth=1)
    plt.plot(bestTurbs.x, bestTurbs.y, marker='o', color='black', linestyle='', markersize=7)
    plt.axis('scaled')                      # Trim the white space
    plt.axis('off')                         # Turn off the framing

    plt.figure("Overlay")
    plt.hold = True
    plt.plot(clsdBP.x*scaledTC, clsdBP.y*scaledTC, color=Iea37sb.getPltClrs().getColor(5), linewidth=1)
    plt.plot(x0Start.x*scaledTC, x0Start.y*scaledTC, marker='o', color=Iea37sb.getPltClrs().getColor(5), linestyle='', markersize=7)
    plt.plot(bestTurbs.x, bestTurbs.y, marker='o', color=Iea37sb.getPltClrs().getColor(1), linestyle='', markersize=7)
    plt.axis('scaled')                      # Trim the white space
    plt.axis('off')                         # Turn off the framing


    plt.show()