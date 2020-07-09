# This is a testfile to ensure the splined_boundary() function works.
import YAML
using PyPlot
using Printf
using Parameters
using FlowFarm; const FF = FlowFarm
using FLOWMath; const FM = FLOWMath
include("baker_cs34_functions.jl")
include("../aepcalc-baker1.jl")

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

# For checking boundary constraints
turbs_per_region =                                                                                                 
bndry_corner_indcies = [ Int64[] for i in 1:nRegions ]

for cntr in 1:nRegions
    bndry_corner_indcies[cntr] =
        append!(bndry_corner_indcies[cntr],getCs34VertList(getCs34Name(cntr)))
    turbs_per_region[cntr] = floor(getCs34NumTurbs(getCs34Name(cntr)))
end

bndry_cons = splined_boundary_discreet_regions(turbine_x, turbine_y, bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies, turbs_per_region)

println(bndry_cons)
# for cntr in 1:nRegions
#     println(bndry_cons[cntr])
# end

# Plot the boundaries
# for cntr in 1:nRegions
#     plot(bndry_x_clsd[cntr], bndry_y_clsd[cntr])
#     num_bndry_bpts = length(bndry_x_clsd[cntr])
#     # for i in 1:num_bndry_bpts
#     #     plt.gcf().gca().add_artist(plt.Circle((bndry_x_clsd[cntr][i],bndry_y_clsd[cntr][i]), rotor_diameter/2.0, fill=true,color="black"))
#     #     plt.text(bndry_x_clsd[cntr][i]+rotor_diameter,bndry_y_clsd[cntr][i]+rotor_diameter, string(i))
#     # end
#     # for j in 1:length(bndry_corner_indcies[cntr])
#     #     plt.gcf().gca().add_artist(plt.Circle((bndry_x_clsd[cntr][bndry_corner_indcies[cntr][j]],bndry_y_clsd[cntr][bndry_corner_indcies[cntr][j]]), rotor_diameter/2.0, fill=true,color="red"))
#     #     plt.text(bndry_x_clsd[cntr][bndry_corner_indcies[cntr][j]]+rotor_diameter,bndry_y_clsd[cntr][bndry_corner_indcies[cntr][j]]+rotor_diameter, string(j))
#     # end
# end

# # Plot the turbines
# for i = 1:length(turbine_x)
#     plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter/2.0, fill=true,color="black"))
#     plt.text(turbine_x[i]+rotor_diameter,turbine_y[i]+rotor_diameter, string(i))
# end

# axis("square")
# axis("off")
# plt.show()


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