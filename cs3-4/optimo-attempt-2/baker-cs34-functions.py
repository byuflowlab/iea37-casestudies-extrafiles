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
def getTurbAtrbtCs4YAML(file_name):
    '''Retreive boundary coordinates from the <.yaml> file'''
    # Read in the .yaml file
    with open(file_name, 'r') as f:
        bndrs = yaml.safe_load(f)['boundaries']
        
    ptList3a = bndrs['IIIa']
    ptList3b = bndrs['IIIb']
    ptList4a = bndrs['IVa']
    ptList4b = bndrs['IVb']
    ptList4c = bndrs['IVc']

    # (Convert from <list> to <coordinate> array)
    # Determine how many points we have for the boundary
    numCoords3a = len(ptList3a)
    numCoords3b = len(ptList3b)
    numCoords4a = len(ptList4a)
    numCoords4b = len(ptList4b)
    numCoords4c = len(ptList4c)
    # Initialize the point list array
    coordList3a = np.recarray(numCoords3a, coordinate)
    coordList3b = np.recarray(numCoords3b, coordinate)
    coordList4a = np.recarray(numCoords4a, coordinate)
    coordList4b = np.recarray(numCoords4b, coordinate)
    coordList4c = np.recarray(numCoords4c, coordinate)

    for i in range(numCoords3a):
        coordList3a[i].x = float(ptList3a[i][0])
        coordList3a[i].y = float(ptList3a[i][1])
    for i in range(numCoords3b):
        coordList3b[i].x = float(ptList3b[i][0])
        coordList3b[i].y = float(ptList3b[i][1])
    for i in range(numCoords4a):
        coordList4a[i].x = float(ptList4a[i][0])
        coordList4a[i].y = float(ptList4a[i][1])
    for i in range(numCoords4b):
        coordList4b[i].x = float(ptList4b[i][0])
        coordList4b[i].y = float(ptList4b[i][1])
    for i in range(numCoords4c):
        coordList4c[i].x = float(ptList4c[i][0])
        coordList4c[i].y = float(ptList4c[i][1])

    return coordList3a, coordList3b, coordList4a, coordList4b, coordList4c


def getTurbAtrbtCs3YAML(file_name):
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


def sliceBoundary(totCoordList, numGridLines):
    # Given a set of <coordinate> on a line (or boundary), and number of grid lines,
    # This will return a set of pts that have sliced that boundary line <numGridLines> many times.
    # TERMINOLOGY: 'Segment' slice of boundary line, 'Leg' means portion of boundary we're on
    tol = 1e-6  # Tolerance for reaching the end of the line
    numDivs = numGridLines-1    # If we have 5 grid lines, that gives us 4 sectors (or divisions)
    # Initialize slice list (extra spot for starting point)
    segCoordList = np.recarray(numGridLines, coordinate)
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

def computeVertLines(splineList, segCoordList, numGridLines, numLinspace, numSides):
    # Make our frame-of reference indices
    rgt = 0
    btm = 1
    lft = 2
    top = 3  
    # Get our vertical lines (in the translated space)
    gridLinspace = np.linspace(0, 1, numLinspace, endpoint=True)
    vertLineArrayTrs = np.recarray([numGridLines, numLinspace], coordinate)
    horizLineArrayTrs = np.recarray([numLinspace, numGridLines], coordinate)
    # This will hold our percentages
    perimDist = np.empty([numSides, len(segCoordList[top])], np.float64)

    # Invert the list so it's not counter-clockwise
    segCoordList[btm] = segCoordList[btm][::-1]
    segCoordList[rgt] = segCoordList[rgt][::-1]

    #-- Do the vertical lines --#
    # x-distribution weight for our top boundary
    topDenom = segCoordList[top][len(
        segCoordList[top])-1].x - segCoordList[top][0].x
    btmDenom = segCoordList[btm][len(
        segCoordList[btm])-1].x - segCoordList[btm][0].x
    for j in range(numLinspace):
        topNum = segCoordList[top][j].x - segCoordList[top][0].x
        btmNum = segCoordList[btm][j].x - segCoordList[btm][0].x
        # Get the percentage along the path
        perimDist[top][j] = (topNum / topDenom)
        # Get the percentage along the path
        perimDist[btm][j] = (btmNum / btmDenom)

    yStartCoord = splineList[btm](segCoordList[btm].x)
    yEndCoord = splineList[top](segCoordList[top].x)
    for i in range(numGridLines):
        for j in range(numLinspace):
            vertLineArrayTrs[i][j].y = (
                (yEndCoord[i] - yStartCoord[i]) * gridLinspace[j]) + yStartCoord[i]
            xEndCoord = segCoordList[rgt][j].x
            xStartCoord = segCoordList[lft][j].x
            # "Where the magic happens" according to Eduardo
            weightTerm = (gridLinspace[j] * perimDist[top][i]) + \
                ((1-gridLinspace[j]) * perimDist[btm][i])
            vertLineArrayTrs[i][j].x = (
                (xEndCoord - xStartCoord) * weightTerm) + xStartCoord
    
    # Do the horizontal lines
    # x-distribution weight for our top and btm boundaries
    rgtDenom = segCoordList[rgt][len(
        segCoordList[rgt])-1].x - segCoordList[rgt][0].x
    lftDenom = segCoordList[lft][len(
        segCoordList[lft])-1].x - segCoordList[lft][0].x
    for j in range(numGridLines):
        rgtNum = segCoordList[rgt][j].x - segCoordList[rgt][0].x
        lftNum = segCoordList[lft][j].x - segCoordList[lft][0].x
        # Get the percentage along the path
        perimDist[rgt][j] = (rgtNum / rgtDenom)
        # Get the percentage along the path
        perimDist[lft][j] = (lftNum / lftDenom)
        
    for i in range(numLinspace):
        xStartCoord = segCoordList[lft][i].x
        xEndCoord = segCoordList[rgt][i].x
        for j in range(numGridLines):
            horizLineArrayTrs[i][j].x = (
                (xEndCoord - xStartCoord) * gridLinspace[j]) + xStartCoord
            yEndCoord = splineList[top](segCoordList[top][j].x)
            yStartCoord = splineList[btm](segCoordList[btm][j].x)
            # "Where the magic happens" according to Eduardo
            weightTerm = (gridLinspace[j] * perimDist[lft][i]) + \
                ((1-gridLinspace[j]) * perimDist[lft][i])
            horizLineArrayTrs[i][j].y = (
                (yEndCoord - yStartCoord) * weightTerm) + yStartCoord

    return horizLineArrayTrs, vertLineArrayTrs

