import numpy as np
import sys
import csv
import yaml


def getTurbLocCSV(sFileName):
    turbineX = np.array([])
    turbineY = np.array([])
    csvFile = csv.reader(open(sFileName, "rb"))
    for row in csvFile:
            turbineX = np.append(turbineX, float(row[0]))
            turbineY = np.append(turbineY, float(row[1]))
    return turbineX, turbineY

def getTurbLocYAML(sFileName):

    turbineX = np.array([])
    turbineY = np.array([])

    # Read in the .yaml file
    with open(sFileName, 'r') as f:
        doc = yaml.load(f)
    
    # rip the x-coordinates
    turbineX = doc['definitions']['position']['items']['xc']
    turbineY = doc['definitions']['position']['items']['yc']

    return turbineX, turbineY


if __name__ == "__main__":
    """for testing during development"""
    turbineX = np.array([])
    turbineY = np.array([])

    # Reads turbine locations from a .csv file
    # example command line syntax is "python AEPcal.py iea37-optocs-example16.yaml"
    turbineX, turbineY = getTurbLocYAML(sys.argv[1])
    print turbineX
    print turbineY
