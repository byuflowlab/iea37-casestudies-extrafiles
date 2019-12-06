#-- Necessary Headers --#
from __future__ import print_function   # For Python 3 compatibility
from scipy.interpolate import interp1d  # To create our splines
import matplotlib.pyplot as plt
import numpy as np
import sys
import yaml                             # For reading .yaml files
from math import radians as DegToRad    # For converting degrees to radians
from math import log as ln  # For natural logrithm

# Structured datatype for holding coordinate pair
coordinate = np.dtype([('x', 'f8'), ('y', 'f8')])  

# -- Helper Functions --#
# Read in boundary (originally written in MATLAB, needed to be translated)
def getTurbAtrbtYAML(file_name):
    '''Retreive boundary coordinates from the <.yaml> file'''
    # Read in the .yaml file
    with open(file_name, 'r') as f:
        ptList = yaml.safe_load(f)['boundaries']['IIIa']

    # (Convert from <list> to <coordinate> array)
    # Determine how many points we have for the boundary
    numCoords = len(ptList)
    # Initialize the point list array
    coordList = np.recarray(numCoords, coordinate)
    for i in range(numCoords):
        coordList[i].x = float(ptList[i][0])
        coordList[i].y = float(ptList[i][1])

    return coordList


def findNewPtOnLine(pt0, pt1, dist):
    # Given (x0,y0) and (x1,y1), finds (x2,y2) a distance dist from (x0,y0)
    # along line defined by (x0,y0) and (x1,y1).
    # https://math.stackexchange.com/questions/175896/finding-a-point-along-a-line-a-certain-distance-away-from-another-point
    newPt = np.recarray(1, coordinate)
    distTot = abs(coordDist(pt0, pt1))
    t = dist/distTot
    newPt.x = ((1-t)*pt0.x) + t*pt1.x
    newPt.y = ((1-t)*pt0.y) + t*pt1.y

    return newPt


def sliceBoundary(totCoordList, numDivs):
    # Given a set of <coordinate> on a line (or boundary), and number of divisions,
    # This will return a set of pts that have sliced that line <numDivs> many times.
    # TERMINOLOGY: 'Segment' slice of boundary line, 'Leg' means portion of boundary we're on
    tol = 1e-6  # Tolerance for reaching the end of the line
    # Initialize slice list (extra spot for starting point)
    segCoordList = np.recarray(numDivs+1, coordinate)
    rsArcLen = getArcLength(totCoordList)    # Get total edge length
    # The length of each subdivided segment.
    rsSegmentLen = rsArcLen/numDivs

    # Traverse along the edge, marking where a segment ends
    distNeeded = rsSegmentLen  # Distance remaining till next marker
    distRemain = distNeeded   # Distance remaining until segment ends
    segCoordList[0] = totCoordList[0]  # Don't forget the first point.
    cntrLeg = 0               # Start at the first leg
    #-- Traverse and place --#
    leftPt = totCoordList[0]
    rightPt = totCoordList[1]
    distRemLeg = abs(coordDist(leftPt, rightPt))
    # Loop through all portions, marking as we go
    for cntrSeg in range(numDivs):
        # If the remainder of our segment is smaller than we need
        while(distRemLeg < distRemain):
            if (abs(distRemain - distRemLeg) > tol):
                distRemain = distRemain - distRemLeg    # Take out that much distance
                cntrLeg = cntrLeg + 1                   # Move to the next leg
                leftPt = totCoordList[cntrLeg]          # Move our left point
                rightPt = totCoordList[cntrLeg+1]
                distRemLeg = abs(coordDist(leftPt, rightPt))
            else:                                       # If rounding error makes us go past our line
                break

        segCoordList[cntrSeg+1] = findNewPtOnLine(leftPt, rightPt, distRemain)
        leftPt = segCoordList[cntrSeg+1]
        distRemLeg = distRemLeg - distRemain       # Take out how much we used
        distRemain = distNeeded                    # Reset how much we need

    return segCoordList


def coordDist(pt0, pt1):
    # Returns the euclidean distance between two <coordinate> points
    xDiff = pt0.x - pt1.x
    yDiff = pt0.y - pt1.y
    return np.sqrt(xDiff**2 + yDiff**2)


def getArcLength(coordList):
    # Returns the total distance of a given <np.ndarray> of <coordinate> type
    numCoordPairs = len(coordList)-1
    totDist = 0                   # Initialize to zero
    for i in range(numCoordPairs):  # Go through all the point pairs
        totDist = totDist + coordDist(coordList[i], coordList[i+1])

    return totDist


#--- Code for visualizing farm area ---#
def closeBndryList(bndryPts):
    # Appends the first element to the end of an <np.ndarray> of type <coordinate> for closed boundary #
    # Determine how many points we have for the boundary
    numCoords = len(bndryPts)
    # Initialize the point list array
    coordList = np.recarray((numCoords+1), coordinate)
    for i in range(numCoords):
        coordList[i].x = float(bndryPts[i].x)
        coordList[i].y = float(bndryPts[i].y)

    coordList[i+1].x = bndryPts[0].x
    coordList[i+1].y = bndryPts[0].y

    return coordList


def printBoundary(bndryPts, size):
    #-- Print the windfarm boundary. bndryPts must be <np.ndarray> of type <coordinate>
    plt.plot(bndryPts.x, bndryPts.y)
    #plt.xlim(coordList.x.min(), coordList.x.max()) # scales the x-axis to only include the boundary
    plt.axis('scaled')                      # Trim the white space
    plt.axis('off')                         # Turn off the framing
    plt.gcf().set_size_inches(size.x, size.y)   # Make it big and readable


def printBoundaryArray(bndryPtsX, bndryPtsY, size):
    #-- Print the windfarm boundary. bndryPts must be <np.ndarray> of type <coordinate>
    plt.plot(bndryPtsX, bndryPtsY)
    plt.axis('scaled')                      # Trim the white space
    plt.axis('off')                         # Turn off the framing
    plt.gcf().set_size_inches(size.x, size.y)   # Make it big and readable


def printVerticies(coordList, vertList, colorName):
    plt.hold = True
    plt.plot(coordList[vertList].x,
             coordList[vertList].y, 'o', color=colorName)
