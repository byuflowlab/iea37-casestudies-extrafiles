# This is a testfile to ensure the splined_boundary() function works.
import YAML
using Printf
using Parameters
include("baker_cs34_functions.jl")

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



#--- Read in the data ---#
scaledAEP = 1#e5
scaledTC = 1#e3
strCase = "cs4"  # Which case study we're doing. 'cs3' or 'cs4'
numGridLines = 10                   # How many gridlines we'll use for the splining
nRegions = 5                     # Number of reigons we're using (cs4 = 5, cs3 = 1)

#- Rip the boundary coordinates from the .yaml file -# 
bnry_file_name = "../../startup-files/iea37-boundary-" * strCase * ".yaml"
x_bndry_coords, y_bndry_coords = getBndryCs4YAML(bnry_file_name)
x_bndry_coords_clsd, y_bndry_coords_clsd = closeBndryLists(x_bndry_coords, y_bndry_coords)

# for cntr in 1:nRegions
#     println(x_bndry_coords_clsd[cntr])
# end
# println()
# for cntr in 1:nRegions
#     println(y_bndry_coords_clsd[cntr])
# end


# Augment the data so it's CCW and from NW "corner"