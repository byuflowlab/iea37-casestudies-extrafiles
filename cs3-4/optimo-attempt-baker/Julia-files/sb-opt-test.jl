# This is a testfile to ensure the splined_boundary() function works.
import YAML
using Printf
using Parameters
include("baker_cs34_functions.jl")
include("../aepcalc-baker1.jl")

### Complete, functional ###
function getBndryCs4YAML(file_name)
    """Retreive boundary coordinates from the <.yaml> file"""

    # Read in the .yaml file
    f = YAML.load(open(file_name))
    bndrs = f["boundaries"]
    # Initialize some variables
    nRegions = 5;
    nPts = fill(0, 1, nRegions)
    ptList = []
    # Go through all our regions to get the boundary points
    for cntr in 1:nRegions
        ptList = push!(ptList, bndrs[getCs34NameYAML(cntr)])
        nPts[cntr] = floor(length(bndrs[getCs34NameYAML(cntr)]))
    end

    # Reorder points from 3b and 4a to be CCW from NE corner
    ptList[2] = circshift(ptList[2],1)
    ptList[3] = circshift(ptList[3],-3)
    # Change from CW -> CCW
    for i in 1:nRegions
        ptList[i] = reverse(ptList[i])
    end

    # Initialize our arrays for the coordinates
    x_boundary_coords = [ Float64[] for i in 1:nRegions ]
    y_boundary_coords = [ Float64[] for i in 1:nRegions ]
    # Read in all the coordinates
    for i in 1:nRegions             # Looping through all regions
        for j in 1:nPts[i]          # Looping through all point sin this region
            x_boundary_coords[i] = push!(x_boundary_coords[i], ptList[i][j][1])
            y_boundary_coords[i] = push!(y_boundary_coords[i], ptList[i][j][2])
        end
    end

    return x_boundary_coords, y_boundary_coords
end

function writeCs4TurbLocs(x_coords, y_coords, file_name, fname_turb, fname_wr, binned_AEP)
    """Write the passed turbine locations to a file in the iea37 format """
    # Combine coordinates

    # Header info
    Data = Dict()
    Data["title"] = "IEA Wind Task 37 case study 4"
    Data["definitions"] = "Result file of Nick Baker using SNOPT and BPM"

    defs = Data["definitions"]
    # Turbine info
    wp = defs["wind_plant"]
    wp["description"] = "specific plant design including turbine selection and placement'"
    wp["properties"]["type"] = "array"
    # wp["properties"]["items{2}"] = fname_turb

    # # Turbine Locations
    # Data.definitions.position.description = 'an array of x-coordinates [x0, x1, ...] and y-coordinates [y0, y1, ...] of wind turbine positions in cartesian coordinates'
    # Data.definitions.position.units = 'm'
    # # num_coords = length(turb_coords);
    # # x_coords = turb_coords(1:(num_coords/2));
    # # y_coords = turb_coords((num_coords/2 +1):num_coords);
    # # Data.definitions.position.items.xc = round(x_coords, 6);
    # # Data.definitions.position.items.yc = round(y_coords, 6);

    # # Wake model
    # Data.definitions.plant_energy.description = 'energy production from simplified Bastankhah Gaussian wake model'
    # Data.definitions.plant_energy.properties.wake_model.description = 'wake model used to calculate AEP'
    # Data.definitions.plant_energy.properties.wake_model.items{1}.ref = '"iea37-aepcalc.py"'

    # # WindRose
    # Data.definitions.plant_energy.properties.wind_resource.description = 'specific wind resource used to calculate AEP'
    # Data.definitions.plant_energy.properties.wind_resource.description = 'MWh'
    # Data.definitions.plant_energy.properties.wind_resource.items{1}.ref = fname_wr

    # # AEP
    # Data.definitions.plant_energy.properties.annual_energy_production.description = 'binned and total (default) annual energy production for a wind plant given a layout and binned wind rose'
    # Data.definitions.plant_energy.properties.annual_energy_production.units = 'MWh'
    # Data.definitions.plant_energy.properties.annual_energy_production.binned = round(binned_AEP,digits=6)
    # Data.definitions.plant_energy.properties.annual_energy_production.default = round(sum(binned_AEP),digits=6)

    YAML.write_file(file_name, Data)
end

#--- Read in the data ---#
# scaledAEP = 1#e5
# scaledTC = 1#e3
# strCase = "cs4"  # Which case study we're doing. 'cs3' or 'cs4'
# numGridLines = 10                   # How many gridlines we'll use for the splining
# nRegions = 5                     # Number of reigons we're using (cs4 = 5, cs3 = 1)

# #- Rip the boundary coordinates from the .yaml file -# 
# bnry_file_name = "../../startup-files/iea37-boundary-" * strCase * ".yaml"
# x_bndry_coords, y_bndry_coords = getBndryCs4YAML(bnry_file_name)
# x_bndry_coords_clsd, y_bndry_coords_clsd = closeBndryLists(x_bndry_coords, y_bndry_coords)

file_dir = "../../startup-files/"
file_name = "iea37-ex-opt4.yaml"
file_name = string(file_dir,file_name)
# Get turbine locations
turb_locs = getTurbLocYAML(file_name)
turb_coords = turb_locs[1]
fname_turb = turb_locs[2]
fname_turb = string(file_dir,fname_turb)
fname_wr = turb_locs[3]
fname_wr = string(file_dir,fname_wr)

# Get turbine attributes
turb_data = getTurbAtrbtYAML(fname_turb)
turb_ci = turb_data[1]
turb_co = turb_data[2]
rated_ws = turb_data[3]
rated_pwr = turb_data[4]
turb_diam = turb_data[5]

# Get windrose info
wr_data = getWindRoseYAML(fname_wr)
wind_dir = wr_data[1]
wind_dir_freq = wr_data[2]
wind_speeds = wr_data[3]
wind_speed_probs = wr_data[4]
num_speed_bins = wr_data[5]
min_speed = wr_data[6]
max_speed = wr_data[7]

AEP = calcAEPcs3(turb_coords, wind_dir_freq, wind_speeds, wind_speed_probs, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)

print(AEP, "\n")
print("AEP = ", sum(AEP), "\n")

writeCs4TurbLocs(1, 2, "test_julia.yaml", fname_turb, fname_wr, AEP)

# for cntr in 1:nRegions
#     println(x_bndry_coords_clsd[cntr])
# end
# println()
# for cntr in 1:nRegions
#     println(y_bndry_coords_clsd[cntr])
# end


# Augment the data so it's CCW and from NW "corner"