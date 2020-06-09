#- Load all previously written and tested functions -#
from os import path
import scipy.spatial.distance as ssd
import matplotlib.pyplot as plt
import numpy as np
import random
import sys
import yaml                             # For reading .yaml files
import time
# For the AEP calculation code
from scipy.interpolate import interp1d  # To create our splines
import baker_cs34_functions as Iea37sb
import iea37_aepcalc as iea37aepC
from baker_cs34_SNOPT_functions import *
from scipy.special import binom         # "Combination", for deterimining unique turbine pairs
from math import radians as DegToRad    # For converting degrees to radians
from math import log as ln  # For natural logrithm

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
    print(coordList4a)
    clsdBP3b = Iea37sb.closeBndryList(coordList3b)    # Duplicate the first coordinate for a closed boundary
    coordList4a = np.roll(coordList4a, -3)            # Shift our points to the left so rightmost vertex is zero
    print(coordList4a)
    clsdBP4a = Iea37sb.closeBndryList(coordList4a)    # Duplicate the first coordinate for a closed boundary
    clsdBP4b = Iea37sb.closeBndryList(coordList4b)    # Duplicate the first coordinate for a closed boundary
    clsdBP4c = Iea37sb.closeBndryList(coordList4c)    # Duplicate the first coordinate for a closed boundary

    


    # #-- Spline all Boundaries --#
    # vertexList3a = [0, 6, 8, 9, 18]
    # numSides3a = len(vertexList3a) - 1      # The number of sides for our original coordinate system.
    # splineList3a, segCoordList3a = Iea37sb.makeCs3BndrySplines(vertexList3a, clsdBP3a, numGridLines)
    # vertexList3b = [0, 1, 2, 3, 8]       # Hard code the vertices (though this could be done algorithmically)
    # numSides3b = len(vertexList3b) - 1
    # splineList3b, segCoordList3b = Iea37sb.makeCs3BndrySplines(vertexList3b, clsdBP3b, numGridLines)
    # vertexList4a = [0, 1, 2, 3, 6]       # Hard code the vertices (though this could be done algorithmically)
    # numSides4a = len(vertexList4a) - 1      # The number of sides for our original coordinate system.
    # splineList4a, segCoordList4a = Iea37sb.makeCs3BndrySplines(vertexList4a, clsdBP4a, numGridLines)
    # vertexList4b = [0, 1, 2, 3]       # Hard code the vertices (though this could be done algorithmically)
    # numSides4b = len(vertexList4b) - 1      # The number of sides for our original coordinate system.
    # splineList4b, segCoordList4b = Iea37sb.makeCs3BndrySplines(vertexList4b, clsdBP4b, numGridLines)
    # vertexList4c = [0, 1, 4, 5]       # Hard code the vertices (though this could be done algorithmically)
    # numSides4c = len(vertexList4c) - 1      # The number of sides for our original coordinate system.
    # splineList4c, segCoordList4c = Iea37sb.makeCs3BndrySplines(vertexList4c, clsdBP4c, numGridLines)