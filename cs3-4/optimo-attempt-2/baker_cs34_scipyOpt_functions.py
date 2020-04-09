import numpy as np
from baker_cs34_functions_sandbox import makeCoordMatrix
import iea37_aepcalc as iea37aepC

def optimoStripArgTuples(args):
    # Takes the passed array and strips items into their proper variables.
    cntr = 0
    wind_dir_freq_len = int(args[0])
    cntr = cntr + 1
    wind_dir_freq = args[cntr:(wind_dir_freq_len + cntr)]
    cntr = cntr + wind_dir_freq_len
    wind_speeds_len = int(args[cntr])
    cntr = cntr + 1
    wind_speeds = args[cntr:(wind_speeds_len + cntr)]
    cntr = cntr + wind_speeds_len
    wind_speed_prob_len = int(args[cntr])
    cntr = cntr + 1
    wind_speed_prob_wid = int(args[cntr])
    cntr = cntr + 1
    wind_speed_prob_size = wind_speed_prob_len * wind_speed_prob_wid
    wind_speed_probs = args[cntr:(wind_speed_prob_size + cntr)]
    wind_speed_probs = wind_speed_probs.reshape(
        wind_speed_prob_len, wind_speed_prob_wid)
    cntr = cntr + wind_speed_prob_size
    wind_dir_len = int(args[cntr])
    cntr = cntr + 1
    wind_dir = args[cntr:(wind_dir_len + cntr)]
    cntr = cntr + wind_dir_len
    turb_diam = float(args[cntr])
    cntr = cntr + 1
    turb_ci = float(args[cntr])
    cntr = cntr + 1
    turb_co = float(args[cntr])
    cntr = cntr + 1
    rated_ws = float(args[cntr])
    cntr = cntr + 1
    rated_pwr = float(args[cntr])
    cntr = cntr + 1
    scaledAEP = float(args[cntr])
    cntr = cntr + 1
    scaledTC = float(args[cntr])

    return [wind_dir_freq, wind_speeds, wind_speed_probs, wind_dir,
            turb_diam, turb_ci, turb_co, rated_ws, rated_pwr, scaledAEP, scaledTC]

def optimoMakeArgTuples(wind_dir_freq, wind_speeds, wind_speed_probs, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr, scaledAEP, scaledTC):
    # Takes the list of passed variables and concatenates them into one long array with size markers for matricies
    wind_speed_probs_size = len(wind_speed_probs) * len(wind_speed_probs[0])
    argLen = len(wind_dir_freq) \
        + len(wind_speeds) \
        + wind_speed_probs_size \
        + len(wind_dir) \
        + 4 \
        + 8  # 4 For lengths of first four entries
    # 7 For single element items <turb_diam, turb_ci, turb_co, rated_ws, rated_pwr, scaledAEP, scaledTurbCoords>

    args = np.zeros(argLen)
    cntr = 0  # Where we are right now.

    args[0] = len(wind_dir_freq)
    cntr = cntr + 1
    args[cntr:(len(wind_dir_freq)+cntr)] = wind_dir_freq
    cntr = cntr + len(wind_dir_freq)
    args[cntr] = len(wind_speeds)
    cntr = cntr + 1
    args[cntr:len(wind_speeds)+cntr] = wind_speeds
    cntr = cntr + len(wind_speeds)
    args[cntr] = len(wind_speed_probs)
    cntr = cntr + 1
    args[cntr] = len(wind_speed_probs[0])
    cntr = cntr + 1
    args[cntr:wind_speed_probs_size +
         cntr] = np.asarray(wind_speed_probs).reshape(-1)
    cntr = cntr + wind_speed_probs_size
    args[cntr] = len(wind_dir)
    cntr = cntr + 1
    args[cntr:len(wind_dir)+cntr] = wind_dir
    cntr = cntr + len(wind_dir)
    args[cntr] = turb_diam
    cntr = cntr + 1
    args[cntr] = turb_ci
    cntr = cntr + 1
    args[cntr] = turb_co
    cntr = cntr + 1
    args[cntr] = rated_ws
    cntr = cntr + 1
    args[cntr] = rated_pwr
    cntr = cntr + 1
    args[cntr] = scaledAEP
    cntr = cntr + 1
    args[cntr] = scaledTC

    return args

#-- Specific for the optimizaiton --#
def optimoFun(x0, args):
    # Rip and parse the data we need
    [wind_freq, wind_speeds, wind_speed_probs, wind_dir, turb_diam,
        turb_ci, turb_co, rated_ws, rated_pwr, scaledAEP, scaledTC] = optimoStripArgTuples(args)
    newCoords = makeCoordMatrix(x0)

    # Calculate the AEP, remembering to scale the coordinates
    dirAEP = iea37aepC.calcAEPcs3((newCoords * scaledTC), wind_freq, wind_speeds, wind_speed_probs,
                                  wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)
    scaledAEP = dirAEP / scaledAEP
    return -np.sum(scaledAEP)
    #return np.sum(dirAEP)


if __name__ == "__main__":
#    numTurbs = 25

    # #- Load the boundary -#
#    fn = "iea37-boundary-cs3.yaml"
#    bndryPts = getTurbAtrbtCs3YAML(fn)
    # Duplicate the 1st coord for a closed boundary
#    clsdBP = closeBndryList(bndryPts)
    # #- Load the turbine and windrose atributes -#
#    fname_turb = "iea37-10mw.yaml"
    # fname_wr = "iea37-windrose-cs3.yaml"
    # wind_dir, wind_dir_freq, wind_speeds, wind_speed_probs, num_speed_bins, min_speed, max_speed = iea37aepC.getWindRoseYAML(
    #     fname_wr)
#    turb_ci, turb_co, rated_ws, rated_pwr, turb_diam = iea37aepC.getTurbAtrbtYAML(
#        fname_turb)

    # #- Some display variables -#
    # displaySize = np.recarray(1, coordinate)
    # displaySize.x = 5
    # displaySize.y = 5
    # numLinspace = 10
    # numGridLines = 10                   # How many gridlines we'll use for the visualization
    # vertexList = [0, 6, 8, 9, 18]       # Hard code the vertices (though this could be done algorithmically)
    # numSides = len(vertexList) - 1
    # scaledAEP = 1e5
    # scaledTC = 1e3
    # #- args in the correct format for optimization -#
    # Args = optimoMakeArgTuples(wind_dir_freq, wind_speeds, wind_speed_probs,
    #                            wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr, scaledAEP, scaledTC)

    # # Spline up the boundary
    # [splineList, segCoordList] = makeCs3BndrySplines(vertexList, clsdBP, numGridLines)

    # #-- Use the example layout --#
#    fn = "iea37-ex-opt3.yaml"
#    turb_coords, fname_turb, fname_wr = iea37aepC.getTurbLocYAML(fn)
#    x0s = makeFirstCoordStruct(turb_coords)
    # x0s = makeCoordStruct(x0a)
#    x0a = makeCoordArray(x0s)
    # # printTurbines(x0s, getPltClrs().getColor(1), turb_diam/2)
    # # startAEP = optimoFun(x0a, Args)
    # numRestarts = 2                     # Number of restarts we're doing
    # listAEP = np.zeros(numRestarts)     # An array holding our AEP values
    # listTurbLocs = np.zeros((numRestarts, (numTurbs*2)))
    # bestResult = np.zeros(2)            # Index, AEP number

#     for cntr in range(numRestarts):
#         #-- Make some random starting turbine places --#
#         # turbRandoList = iea37cs3randomstarts(
#         #     numTurbs, splineList, vertexList, bndryPts, turb_diam)
#         #- Get our turbine list ready for processing -#
#         x0 = makeCoordArray(x0s)/scaledTC         # Get a random turbine placement and scale it
#         # startAEP = optimoFun(x0, Args)
#         # print(startAEP)

#         cons = ({'type': 'ineq', 'fun': lambda x:  checkBndryCons(x, splineList, vertexList, bndryPts, scaledTC)}, {
#                 'type': 'ineq', 'fun': lambda x: checkTurbSpacing(x, turb_diam, scaledTC)})
#         res = optimize.minimize(optimoFun, x0, args=Args, method='SLSQP',
#                                 constraints=cons, options={'disp': True, 'maxiter': 1000})
# #tol= 1e-6
#         #-- Save our results --#
#         listAEP[cntr] = optimoFun(res.x, Args)
#         listTurbLocs[cntr] = res.x

#     for j in range(numRestarts):
#         if (bestResult[1] > listAEP[j]):  # If our new AEP is better (Remember negative switches)
#             bestResult[1] = listAEP[j]    # Save it
#             bestResult[0] = j             # And the index of which run we're on

#     listTurbLocs = listTurbLocs * scaledTC
#     print("Best AEP:")
#     print(listAEP[int(bestResult[0])])  # Print the best one, have to make sure the index is an (int)
#    print(optimoFun(x0a, Args))
#     # save to csv file
#     np.savetxt('turblocs-10run-ex.csv', listTurbLocs[int(bestResult[0])], delimiter=',') # Save the best setup

    # for i in range(numSides):
    #     plt.hold = True
    #     printBoundaryArray(segCoordList[i].x, splineList[i](
    #         segCoordList[i].x), 5)

    # plt.axis('scaled')                      # Trim the white space
    # plt.axis('off')                         # Turn off the framing
    # plt.show()
