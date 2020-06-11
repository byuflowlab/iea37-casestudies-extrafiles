#### TO DO ####
# [ ] Get read in of boundary written
# [ ] Get splining of boundary with flowmath written
# [ ] Get checkBoundary written

import YAML
using Printf
using Parameters

#-- Begin functions specific to the cs3/cs4 case studies and their regions --#
function getCs34NameYAML(indx::Int)
    name = "no such region"

    if(indx == 1)
        name = "IIIa"
    elseif(indx == 2)
        name = "IIIb"
    elseif(indx == 3)
        name = "IVa"
    elseif(indx == 4)
        name = "IVb"
    elseif(indx == 5)
        name = "IVc"
    end
    
    return name
end

function getCs34Name(indx::Int)
    name = "no such region"

    if(indx == 1)
        name = "3a"
    elseif(indx == 2)
        name = "3b"
    elseif(indx == 3)
        name = "4a"
    elseif(indx == 4)
        name = "4b"
    elseif(indx == 5)
        name = "4c"
    end
    
    return name
end

function getCs34NumTurbs(sReg::String)
    numTurbs = zeros(1)

    if(sReg == "cs3")
        numTurbs = 25
    elseif(sReg == "cs4")
        numTurbs = 81
    elseif(sReg == "3a")
        numTurbs = 31
    elseif(sReg == "3b")
        numTurbs = 11
    elseif(sReg == "4a")
        numTurbs = 16
    elseif(sReg == "4b")
        numTurbs = 14
    elseif(sReg == "4c")
        numTurbs = 9
    end
    
    return numTurbs
end

function getCs34VertList(sReg::String)
    if(sReg == "cs3")
        vertList = [1, 10, 11, 13, 19]
    elseif(sReg == "3a")
        vertList = [1, 10, 11, 13, 19]
    elseif(sReg == "3b")
        vertList = [1, 6, 7, 8, 9]
    elseif(sReg == "4a")
        vertList = [1, 4, 5, 6, 7]
    elseif(sReg == "4b")
        vertList = [1, 2, 3, 4]
    elseif(sReg == "4c")
        vertList = [1, 2, 5, 6]
    end
    
    return vertList
end

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
    ptList[1] = circshift(ptList[1],-1)
    ptList[3] = circshift(ptList[3],-4)
    ptList[4] = circshift(ptList[4],-1)
    ptList[5] = circshift(ptList[5],-1)
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
#-- End specific cs3/cs4 functions --#

### Complete and functional ###
function getTurbLocYAML(file_name)
    ### Retrieve turbine locations and auxiliary file names from <.yaml> file.
    ### Auxiliary (reference) files supply wind rose and turbine attributes.

    # Read in the .yaml file
    f = YAML.load(open(file_name))
    defs = f["definitions"]

    # Rip the (x,y) coordinates (Convert from <list> to <ndarray>)
    turb_coords = defs["position"]["items"]

    # Rip the expected AEP, used for comparison
    AEP = defs["plant_energy"]["properties"]["annual_energy_production"]["default"]

    # Read the auxiliary filenames for the windrose and the turbine attributes (first one)
    fname_turb = string.(values(defs["wind_plant"]["properties"]["turbine"]["items"][1]))[1]
    fname_wr = string.(values(defs["plant_energy"]["properties"]["wind_resource"]["properties"]["items"][1]))[1]

    # Return turbine (x,y) locations, and the filenames for the others .yamls
    return turb_coords, fname_turb, fname_wr
end
### Complete and functional ###
function getTurbAtrbtYAML(file_name)
    ###Retreive turbine attributes from the <.yaml> file###

    # Read in the .yaml file
    f = YAML.load(open(file_name))
    defs = f["definitions"]
    ops = defs["operating_mode"]
    turb = defs["wind_turbine"]
    rotor = defs["rotor"]

    # Rip the turbine attributes
    # (Convert from <list> to <float>)
    turb_ci = float(ops["cut_in_wind_speed"]["default"])
    turb_co = float(ops["cut_out_wind_speed"]["default"])
    rated_ws = float(ops["rated_wind_speed"]["default"])
    rated_pwr = float(turb["rated_power"]["maximum"])
    turb_diam = float(rotor["diameter"]["default"])

    return turb_ci, turb_co, rated_ws, rated_pwr, turb_diam
end
### Complete and functional ###
function getWindRoseYAML(file_name)
    ### Retrieve wind rose data (bins, freqs, speeds) from <.yaml> file.

    # Read in the .yaml file
    f = YAML.load(open(file_name))
    props = f["definitions"]["wind_inflow"]["properties"]

    # Rip wind directional bins, their frequency, and the windspeed parameters for each bin
    # (Convert from <list> to <ndarray>)
    wind_dir = props["direction"]["bins"]
    wind_dir_freq = props["direction"]["frequency"]
    # (Convert from <list> to <float>)
    wind_speeds = props["speed"]["bins"]
    wind_speed_probs = props["speed"]["frequency"]
    # Get default number of windspeed bins per direction
    num_speed_bins = length(wind_speeds)
    min_speed = props["speed"]["minimum"]
    max_speed = props["speed"]["maximum"]

    return wind_dir, wind_dir_freq, wind_speeds, wind_speed_probs, num_speed_bins, min_speed, max_speed
end

function closeBndryLists(bndryPts_x, bndryPts_y)
    """ Appends 1st element to the end of each array for a closed bndry """
    # Determine how many regions and points per region were passed
    nRegions = length(bndryPts_x)
    # Append the initial points to the end
    for i in 1:nRegions
        bndryPts_x[i] = push!(bndryPts_x[i], bndryPts_x[i][1])
        bndryPts_y[i] = push!(bndryPts_y[i], bndryPts_y[i][1])
    end

    return bndryPts_x, bndryPts_y
end