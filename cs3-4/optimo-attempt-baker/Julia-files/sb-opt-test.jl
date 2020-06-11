# This is a testfile to ensure the splined_boundary() function works.
import YAML
using PyPlot
using Printf
using Parameters
using FlowFarm; const FF = FlowFarm
include("baker_cs34_functions.jl")
include("../aepcalc-baker1.jl")

function writeCs4TurbLocs(x_coords, y_coords, binned_AEP, file_name)
    """Write the passed turbine locations to a file in the iea37 format """
    #set_bigfloat_precision(5)
    # Combine coordinates

    # Header info
    Data = Dict()
    # Data["Turbines"] = x_coords

    Data["AEP"] = Dict()
    AEP = Data["AEP"]
    #AEP["binned"] = binned_AEP
    AEP["default"] = big(sum(binned_AEP))

    # Data["title"] = "IEA Wind Task 37 case study 4"
    # Data["description"] = "Result file of Nick Baker using SNOPT and BPM"

    # Data["definitions"] = Dict()
    # defs = Data["definitions"]
    # # Turbine info
    # defs["wind_plant"] = Dict()
    # wp = defs["wind_plant"]
    # wp["description"] = "specific plant design including turbine selection and placement'"
    # wp["properties"] = Dict()
    # wp["properties"]["type"] = "array"
    # wp["properties"]["items"] = fname_turb

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
    return Data
end

function checkBndryConsCs4(turbine_x, turbine_y, turbs_per_region, bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies)
end

########################################## MAIN ################################

#--- Read in the data ---#
scaledAEP = 1#e5
scaledTC = 1#e3
strCase = "cs4"  # Which case study we're doing. 'cs3' or 'cs4'
numGridLines = 10                   # How many gridlines we'll use for the splining
nRegions = 5                     # Number of reigons we're using (cs4 = 5, cs3 = 1)

#- Rip the boundary coordinates from the .yaml file -# 
bnry_file_name = "../../startup-files/iea37-boundary-" * strCase * ".yaml"
bndry_x, bndry_y = getBndryCs4YAML(bnry_file_name)
bndry_x_clsd, bndry_y_clsd = closeBndryLists(bndry_x, bndry_y)

#--- Read in Turbine data to calculate AEP ---#
file_dir = "../../startup-files/"
file_name_orig = "iea37-ex-opt4.yaml"
file_name = string(file_dir,file_name_orig)
# Get turbine locations
turbine_x, turbine_y, fname_turb_orig, fname_wr_orig  = FF.get_turb_loc_YAML(file_name)
fname_turb = string(file_dir,fname_turb_orig)
fname_wr = string(file_dir,fname_wr_orig)

# Get turbine attributes
turb_ci, turb_co, rated_ws, rated_pwr, rotor_diameter, turb_height = FF.get_turb_atrbt_YAML(fname_turb)

# Get windrose info
#wind_dir, wind_speeds, wind_dir_freq, ti = FF.get_wind_rose_YAML(fname_wr)
wr_data = getWindRoseYAML(fname_wr)
wind_dir = wr_data[1]
wind_dir_freq = wr_data[2]
wind_speeds = wr_data[3]
wind_speed_probs = wr_data[4]
num_speed_bins = wr_data[5]
min_speed = wr_data[6]
max_speed = wr_data[7]

# Plot the boundaries
for cntr in 1:nRegions
    plot(bndry_x_clsd[cntr], bndry_y_clsd[cntr])
    num_bndry_bpts = length(bndry_x_clsd[cntr])
    for i in 1:num_bndry_bpts
        plt.gcf().gca().add_artist(plt.Circle((bndry_x_clsd[cntr][i],bndry_y_clsd[cntr][i]), rotor_diameter/2.0, fill=true,color="black"))
        plt.text(bndry_x_clsd[cntr][i]+rotor_diameter,bndry_y_clsd[cntr][i]+rotor_diameter, string(i))
    end
end

# Plot the turbines
# for i = 1:length(turbine_x)
#     plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter/2.0, fill=true,color="black"))
# end
axis("square")
axis("off")
plt.show()


# AEP = calcAEPcs3(turb_coords, wind_dir_freq, wind_speeds, wind_speed_probs, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)

# print(AEP, "\n")
# print("AEP = ", sum(AEP), "\n")

#--- Print Boundary coordinates (Debug) ---#
# for cntr in 1:nRegions
#     println(x_bndry_coords_clsd[cntr])
# end
# println()
# for cntr in 1:nRegions
#     println(y_bndry_coords_clsd[cntr])
# end



#--- YAML reading and writing ---#
# f = YAML.load(open(file_name))
# #println(f)

# g = writeCs4TurbLocs(turb_coords, 2,  AEP, "test_julia.yaml")
# println(g)