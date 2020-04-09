from __future__ import print_function   # For Python 3 compatibility
import numpy as np
import sys
import yaml                             # For reading .yaml files
from math import radians as DegToRad    # For converting degrees to radians
from math import log as ln              # For natural logrithm

def getTurbLocYAML(file_name):
    """ Retrieve turbine locations and auxiliary file names from <.yaml> file.

    Auxiliary (reference) files supply wind rose and turbine attributes.
    """
    # Read in the .yaml file
    with open(file_name, 'r') as f:
        defs = yaml.safe_load(f)['definitions']

    # Rip the (x,y) coordinates (Convert from <list> to <ndarray>)
    turb_coords = np.asarray(defs['position']['items'])

    # Rip the expected AEP, used for comparison
    # AEP = defs['plant_energy']['properties']
    #           ['annual_energy_production']['default']

    # Read the auxiliary filenames for the windrose and the turbine attributes
    ref_list_turbs = defs['wind_plant']['properties']['turbine']['items']
    ref_list_wr = (defs['plant_energy']['properties']
                       ['wind_resource']['properties']['items'])

    # Iterate through all listed references until we find the one we want
    # The one we want is the first reference not internal to the document
    # Note: internal references use '#' as the first character
    fname_turb = next(ref['$ref']
                      for ref in ref_list_turbs if ref['$ref'][0] != '#')
    fname_wr = next(ref['$ref']
                    for ref in ref_list_wr if ref['$ref'][0] != '#')

    # Return turbine (x,y) locations, and the filenames for the others .yamls
    return turb_coords, fname_turb, fname_wr


def getWindRoseWeibYAML(file_name):
    """Retrieve wind rose data (bins, freqs, speeds) from <.yaml> file."""
    # Read in the .yaml file
    with open(file_name, 'r') as f:
        props = yaml.safe_load(f)['definitions']['wind_inflow']['properties']

    # Rip wind directional bins, their frequency, and the windspeed parameters for each bin
    # (Convert from <list> to <ndarray>)
    wind_dir = np.asarray(props['direction']['bins'])
    wind_dir_freq = np.asarray(props['direction']['frequency'])
    # (Convert from <list> to <float>)
    wind_speeds = np.asarray(props['speed']['bins'])
    wind_speed_probs = np.asarray(props['speed']['frequency'])
    # Get default number of windspeed bins per direction
    num_speed_bins = wind_speeds.shape[0]
    min_speed = props['speed']['minimum']
    max_speed = props['speed']['maximum']

    return wind_dir, wind_dir_freq, wind_speeds, wind_speed_probs, num_speed_bins, min_speed, max_speed


def getTurbAtrbtYAML(file_name):
    '''Retreive turbine attributes from the <.yaml> file'''
    # Read in the .yaml file
    with open(file_name, 'r') as f:
        defs = yaml.safe_load(f)['definitions']
        ops = defs['operating_mode']
        turb = defs['wind_turbine']
        rotor = defs['rotor']

    # Rip the turbine attributes
    # (Convert from <list> to <float>)
    turb_ci = float(ops['cut_in_wind_speed']['default'])
    turb_co = float(ops['cut_out_wind_speed']['default'])
    rated_ws = float(ops['rated_wind_speed']['default'])
    rated_pwr = float(turb['rated_power']['maximum'])
    turb_diam = float(rotor['diameter']['default'])

    return turb_ci, turb_co, rated_ws, rated_pwr, turb_diam


if __name__ == "__main__":
    """Used for demonstration.

    An example command line syntax to run this file is:

        python iea37-aepcalc.py iea37-ex-opt3.yaml

    For Python .yaml capability, in the terminal type "pip install pyyaml".
    """
    #-- Read necessary values from .yaml files --#
    # Get turbine locations and auxiliary <.yaml> filenames
    # turb_coords, fname_turb, fname_wr = getTurbLocYAML(sys.argv[1])
    # fname_wr = "iea37-windrose-cs4.yaml"
    fname_wr = "iea37-windrose-cs4.yaml"
    # Get the array wind sampling bins, frequency at each bin, and wind speed
    wind_dir, wind_dir_freq, wind_speeds, wind_speed_probs, num_speed_bins, min_speed, max_speed = getWindRoseWeibYAML(
        fname_wr)
    # Pull the needed turbine attributes from file
    # turb_ci, turb_co, rated_ws, rated_pwr, turb_diam = getTurbAtrbtYAML(
    #     fname_turb)
    np.set_printoptions(precision=20)
    # print(sum(wind_dir_freq))
    # wind_dir_freq = np.around(wind_dir_freq, decimals=5)
    # wind_dir_freq = wind_dir_freq/ sum(wind_dir_freq) # Normalize
    # wind_dir_freq = np.around(wind_dir_freq, decimals=5)
    # print(sum(wind_dir_freq))
    # print(wind_dir_freq)
    
    print()
    # print(sum(wind_speed_probs))
    # print(sum(sum(wind_speed_probs)))
    wind_speed_probs = np.around(wind_speed_probs, decimals=10)
    # wind_dir_freq = wind_dir_freq/ sum(wind_dir_freq) # Normalize
    # wind_dir_freq = np.around(wind_dir_freq, decimals=5)
    # print(sum(wind_dir_freq))
    # print(wind_dir_freq)
    np.savetxt('wind_speed_probs.txt', wind_speed_probs,
               fmt='%.10f', delimiter=', ')

    # wind_speed_probs = np.around(wind_speed_probs, decimals=10)
    # wind_speed_probs = wind_speed_probs / np.sum(wind_speed_probs, axis=1)
    #print(wind_speed_probs)
    #print(np.sum(wind_speed_probs, axis=1))