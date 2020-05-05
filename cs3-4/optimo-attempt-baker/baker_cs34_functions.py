#-- Necessary Headers --#
from __future__ import print_function   # For Python 3 compatibility
from scipy.interpolate import interp1d  # To create our splines
from scipy import optimize  # To create our splines
import matplotlib.pyplot as plt
import numpy as np
import sys
import yaml                             # For reading .yaml files
# For the AEP calculation code
import iea37_aepcalc as iea37aepC       # All the given cs3/cs4 functions
from scipy.special import binom         # "Combination", for deterimining unique turbine pairs
from math import radians as DegToRad    # For converting degrees to radians
from math import log as ln              # For natural logrithm

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

#--- Code for making splined boundaries ---#
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
    # If we have 5 grid lines, that gives us 4 sectors (or divisions)
    numDivs = numGridLines-1
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


def makeCs3BndrySplines(vertexList, clsdBP, numGridLines):
    #-- Spline the boundary --#
    numSides = len(vertexList) - 1      # The number of sides for our original coordinate system. Usually (4) to Euclidean, but could be any number)
    splineList = np.empty(numSides, interp1d)                  # Init. array IOT save the Splines for each "side"
    buf = np.zeros((numSides, numGridLines, 2))                # Used to initalize the recarray to zeros
    segCoordList = np.recarray([numSides, numGridLines], dtype=coordinate, buf=buf)

    #- Create the splines for each side (<numSides> many)-#
    for i in range(numSides):
        BndPts = clsdBP[vertexList[i]:(vertexList[i+1]+1)]       # Extract the points for the "edge" we want
        segCoordList[i] = sliceBoundary(BndPts, numGridLines)    # Reparameterize the boundry to be defined by <numGridLines> many points
        splineList[i] = interp1d(segCoordList[i].x, segCoordList[i].y, kind='linear')   # Make the spline using NumPy's <interp1d>

    return splineList, segCoordList

class getPltClrs(object):
    # Made into a MATLAB function 02.Dec.18, converted to Python 26.Feb.20
    # For easy and pretty color plotting
    def getColor(self, argument):
        """Dispatch method"""
        method_name = 'color_' + str(argument)
        # Get the method from 'self'. Default to a lambda.
        method = getattr(self, method_name, lambda: "Invalid month")
        # Call the method as we return it
        return method()

    def color_1(self):
        return [0.6350, 0.0780, 0.1840]  # Red

    def color_2(self):
        return [0.9290, 0.6940, 0.1250]  # Yellow

    def color_3(self):
        return [0.4940, 0.1840, 0.5560]  # Purple

    def color_4(self):
        return [0.4660, 0.6740, 0.1880]  # Green

    def color_5(self):
        return [0, 0.4470, 0.7410]  # Blue

    def color_6(self):
        return [0.8500, 0.3250, 0.0980]  # Orange

    def color_7(self):
        return [0.7500, 0.7500, 0.0000]  # Puke yellow

    def color_8(self):
        return [0.0000, 0.4078, 0.3412]  # loyolagreen

    def color_9(self):
        return [1.0000, 0.8000, 1.0000]  # Pink

    def color_10(self):
        return [0.7843, 0.7843, 0.7843]  # loyolagray

    def color_11(self):
        return [0.6000, 0.4000, 0.2000]  # Brown

    def color_12(self):
        return [0.3010, 0.7450, 0.9330]  # SkyBlue
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


def printBoundary(bndryPts):
    #-- Print the windfarm boundary. bndryPts must be <np.ndarray> of type <coordinate>
    plt.plot(bndryPts.x, bndryPts.y)
    plt.axis('scaled')                      # Trim the white space
    plt.axis('off')                         # Turn off the framing


def printBoundaryClr(bndryPts, colorNum):
    #-- Print the windfarm boundary. bndryPts must be <np.ndarray> of type <coordinate>
    plt.plot(bndryPts.x, bndryPts.y, color=getPltClrs().getColor(colorNum))
    #plt.xlim(coordList.x.min(), coordList.x.max()) # scales the x-axis to only include the boundary
    plt.axis('scaled')                      # Trim the white space
    plt.axis('off')                         # Turn off the framing


def printBoundaryArray(bndryPtsX, bndryPtsY, colorNum):
    #-- Print the windfarm boundary. bndryPts must be <np.ndarray> of type <coordinate>
    plt.plot(bndryPtsX, bndryPtsY, color=getPltClrs().getColor(colorNum))
    plt.axis('scaled')                      # Trim the white space
    plt.axis('off')                         # Turn off the framing


def printVerticies(coordList, vertList, colorNum):
    plt.hold = True
    plt.plot(coordList[vertList].x,
             coordList[vertList].y, 'o', color=getPltClrs().getColor(colorNum))


def printTurbines(coordList, colorName, turbRadius, bShowIndx=False):
    plt.hold = True
    plt.scatter(coordList.x, coordList.y, s=turbRadius, color=colorName)
    if bShowIndx:
        for i in range(len(coordList)):
            # Plot the index number offset enough so the turb circle doesn't cover it
            plt.text(coordList[i].x+turbRadius,
                     coordList[i].y+turbRadius, str(i))

#---- Functions for coordinate data manipulation --#
def makeCoordArray(x0s):
    # Takes <coordinate> input of [(x1,y1), (x2,y2),...]
    # and gives [x1, x2, ... xn, y1 , y2, ...yn]
    x0 = np.concatenate((x0s.x, x0s.y))
    return x0


def makeArrayCoord(x0):
    # Takes array of values [x1, x2, ... xn, y1 , y2, ...yn]
    # and puts them into type <coordinate> for [(x1,y1), (x2,y2),...]
    # For optimzer use.
    nNumRtrs = int(len(x0)/2)

    x0s = np.recarray(nNumRtrs, coordinate)
    x0s.x = x0[0:nNumRtrs]
    x0s.y = x0[nNumRtrs:len(x0)]

    return x0s


def makeArrayMatrix(x0):
    # Takes array of values [x1, x2, ... xn, y1 , y2, ...yn]
    # and makes them [[x1, y1], [x2, y1], ...]
    nNumRtrs = int(len(x0)/2)

    x0m = np.zeros((nNumRtrs, 2))
    x0m[:, 0] = x0[0:nNumRtrs]
    x0m[:, 1] = x0[nNumRtrs:len(x0)]
    return x0m


def makeCoordMatrix(x0s):
    # Takes <coordinate> type of values [(x1,y1), (x2,y2),...]
    # and gives values in matrix form [[x1, y1], [x2, y1], ...]
    nNumRtrs = len(x0s)
    x0m = np.zeros((nNumRtrs, 2))
    x0m[:, 0] = x0s.x
    x0m[:, 1] = x0s.y
    return x0m


def makeMatrixCoord(x0m):
    # Takes ripped .yaml coordinates (in matrix form [[x1, y1], [x2, y2], ...]
    # and puts them into our struct format.
    nNumRtrs = len(x0m)
    x0s = np.recarray(nNumRtrs, coordinate)

    x0s.x = x0m[:, 0]
    x0s.y = x0m[:, 1]

    return x0s


#---- Turbine Spacing constraint function ----#
def checkTurbSpacing(x0, fMinTurbDist):
    #-- Returns an array of the distance between every pair of coordinates
    x0s = makeArrayCoord(x0)
    #-- turbCoords should be of <coordinate> type.
    nNumTurbs = len(x0s)         # Our number of turbines
    # Number of unique turbine pairs > C(numTurbs, 2) = numTurbs! / (2*(numTurbs-2)!).
    nNumPairs = int(binom(nNumTurbs, 2))
    # Array holding the dist. between each pair
    fTurbSpace = np.zeros(nNumPairs)
    bSpacing = np.ones(nNumPairs)  # False means pair is too close, True means they're ok

    nCntr = 0  # Logs where on the list we are
    for i in range(nNumTurbs):          # For every turbine
        for j in range(i):              # Check the space between pairs we haven't calculated
            fTurbSpace[nCntr] = coordDist(x0s[i], x0s[j])
            if ((fTurbSpace[nCntr] - fMinTurbDist) < 0):
                bSpacing[nCntr] = False
            nCntr = nCntr + 1

    # Constrain that the turbines are less than 2 diams apart
    fSpaceConst = fTurbSpace - fMinTurbDist
    return fSpaceConst#, bSpacing  # Negative if ok, positive if too close


#-- Specific for the scipy optimizaiton --#
def optimoFun(x0, args):
    # Assume coordinates are coming in scaled, so upscale accordingly
    # Rip and parse the data we need
    newCoords = makeArrayMatrix(x0) * args['fTCscale'] # Make sure to scale them up

    # Calculate the AEP, remembering to scale the coordinates
    dirAEP = iea37aepC.calcAEPcs3(newCoords, args['wind_dir_freq'], args['wind_speeds'], args['wind_speed_probs'],
                        args['wind_dir'], args['turb_diam'], args['turb_ci'], args['turb_co'], args['rated_ws'], args['rated_pwr'])
    scaledAEP = dirAEP / args['fAEPscale']
    return -np.sum(scaledAEP)

#-- Makes random start locations for cs3 --#
def iea37cs3randomstarts(numTurbs, splineList, coordsCorners, turb_diam):
    #-- Initialize our array --#
    buf = np.zeros((numTurbs, 2))
    turbRandoList = np.recarray((numTurbs), dtype=coordinate, buf=buf)
    minTurbDist = 2*turb_diam

    #-- Get the x-values --#
    xmin = coordsCorners[2].x   # Our minimum x-value
    xmax = coordsCorners[0].x   # our maximum x-value
    for i in range(numTurbs):
        turbRandoList[i].x = np.random.uniform(xmin, xmax)

    #-- Get the y-values --#
    #- Determine the upper and lower splines to use for the given x -#
    # Fake for-loop here for proximity checking.
    i = 0
    while i < numTurbs:
        ymin, ymax = getUpDwnYvals(
            turbRandoList[i].x, splineList, coordsCorners)
        # Get a random number in our bounds
        turbRandoList[i].y = np.random.uniform(ymin, ymax)
        # Check it doesn't conflict with nearby turbines
        for j in range(i):  # Check only the ones we've place so far
            # If this turbine has a proximity conflict
            if (coordDist(turbRandoList[i], turbRandoList[j]) < minTurbDist):
                turbRandoList[i].x = np.random.uniform(
                    xmin, xmax)  # Give it a new x-val
                i = i-1  # Redo the y-val too
                break  # Stop checking for conflicts and redo the y-values
        i = i+1

    return turbRandoList

#---- Boundary Partition Method constraint functions ----#
#-- Returns the max and min y-vals for a given x-val in the cs3 boundary --#
def getUpDwnYvals(xCoord, splineList, coordsCorners):
    # Given there are 4 splines (with 0&1 below, 2&3 above),
    # returns the indecies of splines the given x-coord falls between
    
    #-- Upper. If it's out of bounds
    if (xCoord < coordsCorners[2].x):
        ymax = coordsCorners[2].y # Give it the y-value of our leftmost point
    # If it's to the left of the right upper spline
    elif (xCoord < coordsCorners[3].x):
        ymax = splineList[2](xCoord)  # Make it the left upper spline
    # If it's to the left of the rightmost point
    elif (xCoord < coordsCorners[0].x):
        ymax = splineList[3](xCoord)  # Make it the left upper spline
    else:
        ymax = coordsCorners[0].y # Give it the y-value of our rightmost point

    #-- Lower. If it's to the left of the right lower spline
    if (xCoord < coordsCorners[2].x):
        ymin = coordsCorners[2].y # Give it the y-value of our leftmost point
    elif (xCoord < coordsCorners[1].x):
        ymin = splineList[1](xCoord) # Make it the left upper spline
    elif (xCoord < coordsCorners[0].x):
        ymin = splineList[0](xCoord)
    else:
        ymin = coordsCorners[0].y # Give it the y-value of our rightmost point
    
    return ymin,ymax

#-- Returns neg numbers for out of bounds coordinates, positive if it's in bounds --#
def checkBndryCons(x0, splineList, coordsCorners):
    x0s = makeArrayCoord(x0)
    numTurbs = int(len(x0s))

    # Check to make sure our poits are in
    bndryCons = np.zeros(numTurbs*4)   # four values (two x and two y) for every turbine
    xmin = coordsCorners[2].x   # Our minimum x-value
    xmax = coordsCorners[0].x   # our maximum x-value
    
    # For every turbine
    for i in range(numTurbs):
        #- Check x-vals
        bndryCons[4*i] = (xmax - x0s[i].x)   # Positive good, neg bad
        bndryCons[4*i+1] = (x0s[i].x - xmin) # pos good, neg bad
        
        #- Check y-vals
        ymin,ymax = getUpDwnYvals(x0s[i].x, splineList, coordsCorners)
        bndryCons[4*i+2] = (ymax - x0s[i].y)
        bndryCons[4*i+3] = (x0s[i].y - ymin)
            
    return bndryCons

#---- Boundary Normal Method constraint functions ----#
#-- Functions for checking if turbines are inside a concave boundary --#
def bndryNormals(bndryList):
    # Rewritten and adapted from Jared Thomas' code on 25.Mar.20
    # Number of verticies in our boundary (minus the repeat for closure)
    nVerts = len(bndryList)-1
    unit_normals = np.zeros([nVerts, 2])  # For our unit Normals
    unit_normals_coords = np.recarray(nVerts, coordinate)

    # determine if point is inside or outside of each face, and distance from each face
    for j in range(nVerts):
        # calculate the unit normal vector of the current face (taking points CCW)
        if j < nVerts:  # all but the set of point that close the shape
            normal = np.array([bndryList[j+1].y - bndryList[j].y,
                               -(bndryList[j+1].x-bndryList[j].x)])
            unit_normals[j] = normal/np.linalg.norm(normal)
        else:   # the set of points that close the shape
            normal = np.array([bndryList[0].y-bndryList[j].y,
                               -(bndryList[0].x-bndryList[j].x)])
            unit_normals[j] = normal/np.linalg.norm(normal)

    #Convert to our <coordinate> data type
    for i in range(nVerts):
        unit_normals_coords[i].x = unit_normals[i][0]
        unit_normals_coords[i].y = unit_normals[i][1]

    return unit_normals_coords


def calcDistNorms(x0, vertices, unit_normals):
    # Rewritten and adapted from Jared Thomas' code on 25.Mar.20
    # print points.shape, vertices.shape, unit_normals.shape
    points = makeArrayCoord(x0)
    nPoints = points.shape[0]
    nVertices = vertices.shape[0] - 1
    # initialize array to hold distances from each point to each face
    face_distance = np.zeros([nPoints, nVertices])
    # init bool array indicating whether a pt is in the hull or not
    inside = np.zeros(nPoints)
    # Temp vector from turbine to boundary face
    pa = np.zeros(2)
    d_vec = np.zeros(2)               # Temp vector
    #-- Convert from <coordinate> --#
    unit_norms = makeCoordMatrix(unit_normals)
    turbCoords = makeCoordMatrix(points)
    bndryCoords = makeCoordMatrix(vertices)

    for i in range(nPoints):          # loop through pts and find dist to each face
        # determine if pt is in or out of each face, and dist from each face
        for j in range(nVertices):
            # define the vector from the point of interest to the first point of the face
            pa = [[bndryCoords[j, 0]-turbCoords[i, 0],
                   bndryCoords[j, 1]-turbCoords[i, 1]]]
            # find perpendicular distance from point to current surface (vector projection)
            d_vec = np.vdot(pa, unit_norms[j])*unit_norms[j]
            # calculate the sign of perpendicular distance from point to current face (- is inside, + is outside)
            face_distance[i, j] = np.vdot(d_vec, unit_norms[j])
        # check if the point is inside the convex hull by checking the sign of the distance
        if np.all(face_distance[i] <= 0):
            inside[i] = True
    return face_distance#, inside


def line(p1, p2):
    # Makes a line from the given <coordinate> points
    A = (p1.y - p2.y)
    B = (p2.x - p1.x)
    C = (p1.x*p2.y - p2.x*p1.y)
    return A, B, -C

def intersection(L1, L2):
    # Finds intersection of the given lines
    D = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x, y
    else:
        return False

#- Construct the simplified (concave) cs3 boundary -#
def makeSimpleCs3Bndry(clsdBP):
    # Given the boundary points for cs3, this simplifies them into a non-concave shape
    numVertices = 4  # Our new closed boundary. Will be made of only 4 values, to stay convex
    # Our new closed boundary. Will be made of only 4 values, to stay convex
    newVertices = np.recarray(numVertices, coordinate)

    #-- Make our line equations --#
    lineRgt = line(clsdBP[0], clsdBP[1])  # - Right side line equation -#
    lineBtm = line(clsdBP[6], clsdBP[8])  # - Bottom side line equation -#
    lineLft = line(clsdBP[8], clsdBP[9])  # - Left side line equation -#
    #- Top side line equation -#
    crdTemp = np.recarray(1, coordinate)
    crdTemp.x = 0
    crdTemp.y = clsdBP[12].y
    # point [12] is the min y-val of that curve
    lineTop = line(crdTemp[0], clsdBP[12])

    #-- Figure out intersection points --#
    [newVertices[0].x, newVertices[0].y] = intersection(lineRgt, lineTop)
    [newVertices[1].x, newVertices[1].y] = intersection(lineRgt, lineBtm)
    newVertices[2] = clsdBP[8]
    [newVertices[3].x, newVertices[3].y] = intersection(lineLft, lineTop)
    newVertices = closeBndryList(newVertices)  # Close up the loop for plotting

    return newVertices