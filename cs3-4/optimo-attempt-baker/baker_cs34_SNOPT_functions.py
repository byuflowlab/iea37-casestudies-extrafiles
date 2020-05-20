from pyoptsparse import Optimization, SNOPT, pyOpt_solution, NSGA2
from iea37_aepcalc import calcAEPcs3

def f(ParamsDict):
    # To follow Dr. Ning's rubric
    return lambda x0: calcAEPcs3(x0, ParamsDict['wind_dir_freq'],
                                           ParamsDict['wind_speeds'],
                                           ParamsDict['wind_speed_probs'],
                                           ParamsDict['wind_dir'],
                                           ParamsDict['turb_diam'],
                                           ParamsDict['turb_ci'],
                                           ParamsDict['turb_co'],
                                           ParamsDict['rated_ws'],
                                           ParamsDict['rated_pwr'])


def makeParamsDict(wind_dir_freq, wind_speeds, wind_speed_probs,
                   wind_dir, turb_diam, turb_ci, turb_co, rated_ws,
                   rated_pwr, fAEPscale, fTCscale, cncvNorms, cncvVerts):
    # A helper function to make a dictionary of all the necessary constants for AEP calculation
    ieq37Dict = dict([('wind_dir_freq', wind_dir_freq),
                      ('wind_speeds', wind_speeds),
                      ('wind_speed_probs', wind_speed_probs),
                      ('wind_dir', wind_dir),
                      ('turb_diam', turb_diam),
                      ('minTurbSpace', 2*turb_diam),
                      ('turb_ci', turb_ci),
                      ('turb_co', turb_co),
                      ('rated_ws', rated_ws),
                      ('rated_pwr', rated_pwr),
                      ('fAEPscale', fAEPscale),
                      ('fTCscale', fTCscale),
                      ('cncvNorms', cncvNorms),
                      ('cncvVerts', cncvVerts)])

    return ieq37Dict
