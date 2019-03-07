import numpy as np
import sys
import csv
import matplotlib.pyplot as plt

#RipTurbineLocations
def getTurbLoc(sXFile, sYFile):

    turbineX = np.array([],[])
    turbineY = np.array([],[])

    turbineX = np.loadtxt(open(sXFile, "rb"), delimiter=',')
    turbineY = np.loadtxt(open(sYFile, "rb"), delimiter=',')

    return turbineX, turbineY

def getAEPlist(sAEPFile):
    AEPlist = np.array([])                                      # Declare var
    AEPlist = np.loadtxt(open(sAEPFile, "rb"), delimiter='\n')
    return AEPlist

def getTopAEP(fAEPlist, nNumTop):
    # Pulls the indices of the highest AEPs (unsorted)
    bestInd = np.argpartition(fAEPlist, -nNumTop)[-nNumTop:]
    # Sort them in the right order
    bestInd = bestInd[np.argsort(fAEPlist[bestInd])]
    return bestInd

def plotCircleFarm(fXcoords, fYcoords, rtrDiam, nNumTurb, fieldRad, nBest, sType, color):
    #color = (74./255., 145./255., 200./255.)  # For nice MATLAB blue
    """print everything here"""
    #opt_name = "bidirectional_seclearccond"

    plt.figure(1)
    for i in range(nNumTurb):
        circ_opt = plt.Circle(
            (fXcoords[i]*1., fYcoords[i]*1.), rtrDiam/2., facecolor=color, edgecolor=color, alpha=0.5)
        plt.gca().add_patch(circ_opt)
    circ_outer = plt.Circle((0, 0), fieldRad, linestyle='dashed',
                            edgecolor='black', facecolor='None', label='Boundaries')
    plt.gca().add_patch(circ_outer)
    plt.axis('equal')
    plt.title("" + sType+"_"+str(nBest+1))
    if (nBest == 5):
        plt.savefig("" + sType + "-overlay-" + str(nBest+1) + ".pdf", transparent=True)
    plt.savefig("" +sType+ "_" +str(nBest+1)+ ".pdf", transparent=True)
    plt.show()
    plt.gcf().clear()

## Main
if __name__ == "__main__":
    RoseType = 4 # 4 = quad, 3 = tri, 2 = bi-off
    nNumTurb = 16 # 16 Turbine Farm
    rtrDiam = 130 # Diameter of NREL 3.35 MW turbine rotor
    fieldRad = 1300.
    nNumTopPlots = 6  # Denotes number of top AEP layouts to plot

    # Figure out which wind rose to do
    if (RoseType == 4):
        sAEPFile = 'AEP-quad.csv'
        sXcoords = 'x-quad.csv'
        sYcoords = 'y-quad.csv'
        sType = 'quad'
    elif (RoseType == 3):
        sAEPFile = 'AEP-tri.csv'
        sXcoords = 'x-tri.csv'
        sYcoords = 'y-tri.csv'
        sType = 'tri'
    elif (RoseType == 2):
        sAEPFile = 'AEP-bioff.csv'
        sXcoords = 'x-bioff.csv'
        sYcoords = 'y-bioff.csv'
        sType = 'bioff'

    color = np.empty(6, dtype='string')
    color[0] = 'blue'
    color[1] = 'red'
    color[2] = 'green'
    color[3] = 'yellow'
    color[4] = 'magenta'
    color[5] = 'cyan'

    # Rop turb locations and AEP values from files
    fTurbX, fTurbY = getTurbLoc(sXcoords, sYcoords) # Rip X/Y coords
    fAEPlist = getAEPlist(sAEPFile)                 # Rip AEPs
    
    # Organize by AEP
    bestInd = getTopAEP(fAEPlist, nNumTopPlots)          # Gets top AEP indices
    print fAEPlist[bestInd]
    # Plot top farms
    for i in range(nNumTopPlots):
        plotCircleFarm(fTurbX[bestInd[i]], fTurbY[bestInd[i]], rtrDiam, nNumTurb, fieldRad, i, sType, color[i])
