# This converts the turbine locations in csv format to .yaml files
from os import path
import sys
sys.path.append('..')
import yaml
import numpy as np
import baker_cs34_functions as Iea37sb

csv_file = 'baker-cs3-bpm-snopt.csv'
TurbLocsA = np.loadtxt(csv_file, delimiter=',')
TurbsLocsM = Iea37sb.makeArrayMatrix(TurbLocsA)

result_file = 'baker-cs4-bpm-snopt.yaml'

print(TurbsLocsM)

coords_yaml = [{'coords': TurbsLocsM}]

with open(result_file, 'w') as file:
    documents = yaml.dump(coords_yaml, file)


