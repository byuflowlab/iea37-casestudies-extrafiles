#- Load all previously written and tested functions -#
import random
from baker_cs34_SNOPT_functions import *
from scipy.interpolate import interp1d  # To create our splines
import baker_cs34_functions as Iea37sb
import iea37_aepcalc as iea37aepC
import numpy as np

#-- Setting up objFun --#
def cs3posObjFun(posDict):
    #- Target function (AEP calculation) -#

    # Take our turbine array [x1 ... y1...] and make a matrix
    x0m = Iea37sb.makeArrayMatrix(posDict['aTurbCoords'])
    #nTurbs = len(x0m)  # get the number of turbines
    #nPairs = int(binom(nTurbs, 2))  # Number of unique turbine pairs
    funcs = {}
    AEP = g(x0m)         # Our tricky global function
    totAEP = np.sum(AEP)
    funcs['obj'] = - (totAEP)  # / fScaleFactorAEP) # Negative to minimize

    #- Prep data for constraints -#
    #fScaleFactorTurbLoc = posDict['fTCscale']
    x0s = Iea37sb.makeArrayCoord(posDict['aTurbCoords'])
    [fNormVals, __] = Iea37sb.calcDistNorms(
        x0s, dictParams['cncvVerts'], dictParams['cncvNorms'])
    funcs['bndry'] = fNormVals.flatten()
    fTurbSpace, __ = Iea37sb.checkTurbSpacing(x0s, dictParams['minTurbSpace'])
    funcs['spacing'] = fTurbSpace

    #- Constraints (Pairwise distance [C(numTurbs, 2)], and all boundary checks [4* numTurbs]) -#
    #nBndryCnstr = nTurbs * 4 # each turbine (25) has one constraint for each boundary (4)
    #nSpaceCnstr = nPairs
    #conval = [0]*(nBndryCnstr + nSpaceCnstr)
    # List all the boundary and spacing constraints
    #conval = np.concatenate(cnstrSpacing, nBndryCnstr)
    #funcs['con'] = conval

    fail = False

    return funcs, fail

if __name__ == "__main__":

    #--- Load boundary, turb attributes, and windrose data ---#
        #-- Load the Boundary --#
    fn = "iea37-boundary-cs3.yaml"
    bndryPts = Iea37sb.getTurbAtrbtCs3YAML(fn)      # Pull the boundary vertices
    clsdBP = Iea37sb.closeBndryList(bndryPts)       # repeat the first so it's 'closed'
    cncvVerts = Iea37sb.makeSimpleCs3Bndry(clsdBP)  # Make the simplified Concave shape
    # Calculate the normals for the concave shape
    cncvNorms = Iea37sb.bndryNormals(cncvVerts)
    #- Load the turbine and windrose atributes -#
    fname_turb = "iea37-10mw.yaml"
    fname_wr = "iea37-windrose-cs3.yaml"
    wind_dir, wind_dir_freq, wind_speeds, wind_speed_probs, num_speed_bins, min_speed, max_speed = iea37aepC.getWindRoseYAML(
        fname_wr)
    turb_ci, turb_co, rated_ws, rated_pwr, turb_diam = iea37aepC.getTurbAtrbtYAML(
        fname_turb)
    fAEPscale = 1.0
    fTCscale = 1.0
    #- Make a dictionary for variable passing -#
    dictParams = makeParamsDict(wind_dir_freq, wind_speeds, wind_speed_probs, wind_dir, turb_diam,
                                turb_ci, turb_co, rated_ws, rated_pwr, fAEPscale, fTCscale, cncvNorms, cncvVerts)

    #-- Make boundary splines and random turbine locations--#
    #-- Spline the boundary --#
    # Hard code the vertices (though this could be done algorithmically)
    numTurbs =25
    vertexList = [0, 6, 8, 9, 18]
    # The number of sides for our original coordinate system. Usually (4) to Euclidean, but could be any number)
    numSides = len(vertexList) - 1
    # Init. array IOT save the Splines for each "side"
    splineList = np.empty(numSides, interp1d)
    # Used to initalize the recarray to zeros
    numGridLines = 10
    buf = np.zeros((numSides, numGridLines, 2))
    segCoordList = np.recarray([numSides, numGridLines], dtype=coordinate, buf=buf)

    #- Create the splines for each side (<numSides> many)-#
    for i in range(numSides):
        # Extract the points for the "edge" we want
        BndPts = clsdBP[vertexList[i]:(vertexList[i+1]+1)]
        # Reparameterize the boundry to be defined by <numGridLines> many points
        segCoordList[i] = Iea37sb.sliceBoundary(BndPts, numGridLines)
        # Make the spline using NumPy's <interp1d>
        splineList[i] = interp1d(
            segCoordList[i].x, segCoordList[i].y, kind='linear')

    #-- Make random turbine locations --#
    # Make an array of just the four "vertex" points we're using
    vertexPts = bndryPts[vertexList[0:4]]
    turbRandoList = Iea37sb.iea37cs3randomstarts(
        numTurbs, splineList, vertexPts, turb_diam)
    # Use the random locations calculated previously
    randoMatrix = Iea37sb.makeCoordMatrix(turbRandoList)


    #--- Running the optimization ---#
    #- Constants --#
    nNumTurbs = len(turbRandoList)    # Number of turbines we're passed
    nNumPairs = int(binom(nNumTurbs, 2))  # Number of unique turbine pairs
    g = f(dictParams)

    #-- Setup optimization --#
    optProb = Optimization('CaseStudy3', cs3posObjFun)
    optProb.addObj('obj')
    #-- Setup variables to alter (turbine locations) --#
    x0rando = Iea37sb.makeCoordArray(turbRandoList)   # Make our coordinate list an array
    # two (2) values for each turbine (x&y)
    optProb.addVarGroup('aTurbCoords', 2*nNumTurbs, type='c', value=x0rando)

    #-- Boundary constraints (setup to stay positive) --#
    # 4 boundaries to chek for each turbine
    optProb.addConGroup('bndry', 4*nNumTurbs, lower=0.0, upper=None)
    #-- Spacing constraint (setup to stay positive) --#
    optProb.addConGroup('spacing', nNumPairs, lower=0.0, upper=None)

    #-- Run the actual optimization --#
    opt = SNOPT()
    opt.setOption('Scale option', 1)
    opt.setOption('Iterations limit', 1000000)
    opt.setOption('Major optimality tolerance', 1.e-5)
    opt.setOption('Major feasibility tolerance', 1.e-6)
    #opt.setOption('Print file', 'print_file.out')
    #opt.setOption('Summary file','summary_file.out')

    #print(optProb)

    sol = opt(optProb)
    print(sol)
