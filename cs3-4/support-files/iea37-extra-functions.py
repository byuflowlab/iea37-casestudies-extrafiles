from __future__ import print_function   # For Python 3 compatibility
import numpy as np
import sys
import yaml                             # For reading .yaml files
from math import radians as DegToRad    # For converting degrees to radians
from math import log as ln              # For natural logrithm


def calcWeibull(x, L, k):
    """Given a set of <Lambda>, <k> Weibull variables,
        returns the function evaluation at <x>
    """
    if x < 0:
        return 0
    else:
        fFirstPart = (k/L) * ((x/L)**(k-1))
        fSecondPart = np.exp(-((x/L)**k))
        return fFirstPart * fSecondPart


def calcWeibullInt(x, L, k):
    """Given a set of <Lambda>, <k> Weibull variables,
        returns the integral of the Weibull function at <x>
    """
    if x < 0:
        return 0
    else:
        return 1 - np.exp(-((x/L)**k))
    return


def binWindSpeeds(speed_weib, num_speed_bins, min_speed, max_speed):
    """Given an array of [Lambda, k] Weibull variables in <speed_weib>,
        creates <num_speed_bins> many probabilities and speeds
    """
    bin_divs = np.linspace(min_speed, max_speed, num=(
        num_speed_bins+1))
    num_dir_bins = len(speed_weib)
    # Matrix of proabilities for each speed bin and direcitonal bin
    wind_speeds = np.zeros((num_dir_bins, num_speed_bins))
    wind_speed_probs = np.zeros((num_dir_bins, num_speed_bins))
    #print(bin_divs)
    #print(num_dir_bins)
    # For each bin, find the relevant values
    for i in range(num_dir_bins):
        for j in range(num_speed_bins):
            # The area under the curve is the diff of the integral evaluated
            # at start and end of bin.
            wind_speed_probs[i][j] = \
                calcWeibullInt(bin_divs[j+1], speed_weib[i][0], speed_weib[i][1]) - \
                calcWeibullInt(bin_divs[j], speed_weib[i][0], speed_weib[i][1])
            # Calculations for finding wind speed of 50% frequency in bin
            halfIntegral = wind_speed_probs[i][j]/2
            RHS = halfIntegral + \
                calcWeibullInt(bin_divs[j], speed_weib[i][0], speed_weib[i][1])
            # x = Lambda * (kth root of (-ln(1 - (1/2)Integral)))
            wind_speeds[i][j] = round(speed_weib[i][0] *
                                      ((- ln(1-RHS)) ** (1/speed_weib[i][1])), 2)
        # Normalize to 6 decimal places
        wind_speed_probs[i] = np.around(wind_speed_probs[i], decimals=6)
        wind_speed_probs[i] = wind_speed_probs[i] / sum(wind_speed_probs[i])

    return wind_speeds, wind_speed_probs


if __name__ == "__main__":

    #-- Bin the wind speeds from continuous Weibull --#
    # 60 bins:
    #speed_weib = [  [9.02, 2.10],  [8.85, 2.09],  [8.72, 2.09],  [8.63, 2.09],  [8.55, 2.10],
    #                [8.49, 2.12],  [8.44, 2.15],  [8.41, 2.18],  [8.42, 2.22],  [8.47, 2.26],
    #                [8.57, 2.29],  [8.75, 2.32],  [8.99, 2.35],  [9.27, 2.39],  [9.57, 2.43],
    #                [9.88, 2.48], [10.17, 2.54], [10.44, 2.60], [10.68, 2.65], [10.88, 2.71],
    #               [11.02, 2.75], [11.11, 2.78], [11.16, 2.8],  [11.17, 2.80], [11.18, 2.79],
    #               [11.18, 2.77], [11.18, 2.75], [11.20, 2.72], [11.23, 2.69], [11.25, 2.65], 
    #               [11.28, 2.62], [11.32, 2.57], [11.34, 2.53], [11.36, 2.49], [11.38, 2.44], 
    #               [11.37, 2.40], [11.35, 2.36], [11.31, 2.33], [11.26, 2.30], [11.21, 2.28],
    #               [11.16, 2.26], [11.12, 2.26], [11.10, 2.25], [11.08, 2.26], [11.08, 2.26], 
    #               [11.09, 2.27], [11.11, 2.28], [11.13, 2.29], [11.14, 2.29], [11.13, 2.30], 
    #               [11.09, 2.29], [11.02, 2.29], [10.90, 2.28], [10.74, 2.26], [10.55, 2.23],
    #               [10.32, 2.21], [10.06, 2.18], [ 9.78, 2.15],  [9.50, 2.12],  [9.24, 2.11]]
    # 20 bins:
    speed_weib = [  [ 9.02, 2.1],   [8.63, 2.09],  [8.44, 2.15],  [8.47, 2.26],  [8.99, 2.35],
                    [ 9.88, 2.48], [10.68, 2.65], [11.11, 2.78], [11.18, 2.79], [11.20, 2.72],
                    [11.28, 2.62], [11.36, 2.49], [11.35, 2.36], [11.21, 2.28], [11.10, 2.25],
                    [11.09, 2.27], [11.14, 2.29], [11.02, 2.29], [10.55, 2.23],  [9.78, 2.15]]
    num_speed_bins = len(speed_weib)
    min_speed = 0.0
    max_speed = 25.0
    wind_speeds, wind_speed_probs = binWindSpeeds(speed_weib, num_speed_bins, min_speed,
                                                  max_speed)
    wind_speed_probs = np.around(wind_speed_probs, decimals=10)
    wind_speed_probs = wind_speed_probs / np.sum(wind_speed_probs, axis=1)
    #print(wind_speed_probs)
    #print(np.sum(wind_speed_probs, axis=1))
    np.savetxt('windspeeds.txt', wind_speeds,
               fmt='%.2f', delimiter=', ')
    np.savetxt('windspeedprobs.txt', wind_speed_probs,
               fmt='%.10f', delimiter=', ')
