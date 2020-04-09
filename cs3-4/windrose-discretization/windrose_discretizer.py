from ruamel.yaml import YAML as _yaml
import xarray as _xr
import numpy as _np
import matplotlib.pyplot as _plt


#%% load original wind resource
with open("iea37-windrose-cs3.yaml") as f:
    _yaml.allow_duplicate_keys = True
    wind_res_orig = _yaml(typ='safe').load(f)

#%% extract data from structure
props_orig = wind_res_orig['definitions']['wind_inflow']['properties']
directions_orig = props_orig['direction']['bins']
direction_pmf_orig = props_orig['direction']['frequency']
speeds_orig = props_orig['speed']['bins']
speed_cpmf_orig = props_orig['speed']['frequency']

#%% create xarray dataset
ds_orig = _xr.Dataset(
    data_vars={'direction_pmf': ('direction', direction_pmf_orig),
               'speed_cpmf': (('direction', 'speed'), speed_cpmf_orig)},
    coords={'direction': directions_orig, 'speed': speeds_orig}
)

#%% create cyclical dataset
dirs = ds_orig.coords['direction']
dirs_cyc = _xr.concat([dirs, dirs[0]], 'direction')
ds_cyc = ds_orig.sel({'direction': dirs_cyc})
dirs_cyc[-1] += 360.0
ds_cyc.coords['direction'] = dirs_cyc

#%% interpolate & normalize
ds = ds_cyc.interp(direction=_np.linspace(0, 359, 360), method='linear')
ds['direction_pmf'] /= ds.direction_pmf.sum('direction')

#%% graphical output for visual inspection
_plt.figure()
(ds_orig.speed_cpmf * ds_orig.direction_pmf).plot()
_plt.figure()
(ds.speed_cpmf * ds.direction_pmf).plot()

#%% create output data structure
props = props_orig
props['direction']['bins'] = ds.direction.values.tolist()
props['direction']['frequency'] = ds.direction_pmf.values.tolist()
props['speed']['bins'] = ds.speed.values.tolist()
props['speed']['frequency'] = ds.speed_cpmf.values.tolist()
wind_res = wind_res_orig
wind_res['definitions']['wind_inflow']['properties'] = props

#%% modify description
wind_res['description'] = (
    "Wind resource conditions for a modified IEA37 WFLO case study 4.")

#%% write out interpolated wind resource
with open("iea37-windrose-cs4.yaml", "w") as f:
    _yaml(typ='safe').dump(wind_res, f)

#%% also save as a netcdf4 file
ds.to_netcdf('iea37-windrose-cs4.nc')
