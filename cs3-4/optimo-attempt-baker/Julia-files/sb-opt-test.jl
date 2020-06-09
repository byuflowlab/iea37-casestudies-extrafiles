# This is a testfile to ensure the splined_boundary() function works.
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
        vertList = [1, 7, 9, 10, 19]
    elseif(sReg == "3a")
        vertList = [1, 7, 9, 10, 19]
    elseif(sReg == "3b")
        vertList = [1, 2, 3, 4, 9]
    elseif(sReg == "4a")
        vertList = [1, 2, 3, 4, 7]
    elseif(sReg == "4b")
        vertList = [1, 2, 3, 4]
    elseif(sReg == "4c")
        vertList = [1, 2, 5, 6]
    end
    
    return vertList
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
### Incomplete, untested ###
function getBndryCs4YAML(file_name)
    ###Retreive boundary coordinates from the <.yaml> file

    # Read in the .yaml file
    f = YAML.load(open(file_name))
    bndrs = f["boundaries"]
    
    nRegions = 5;
    nPts = fill(0, 1, nRegions)
    ptList = []
    for cntr in 1:nRegions
        #println(cntr)
        #println(length(bndrs[getCs34NameYAML(cntr)]))
        ptList = push!(ptList, bndrs[getCs34NameYAML(cntr)])
        nPts[cntr] = floor(length(bndrs[getCs34NameYAML(cntr)]))
    end

    for i in 1:nRegions
        println(ptList[i])
    end

    println()

    x_boundary_coords = fill(0.0, 1, nRegions)
    y_boundary_coords = fill(0.0, 1, nRegions)
    for i in 1:nRegions
        for j in 1:nPts[i]
            x_boundary_coords[i] = pushlast!(x_boundary_coords[i], ptList[i][j][1])
            y_boundary_coords[i] = pushlast!(x_boundary_coords[i], ptList[i][j][2])
        end
    end

    println(x_boundary_coords)
    println(y_boundary_coords)
end


#--- Read in the data ---#
scaledAEP = 1#e5
scaledTC = 1#e3
strCase = "cs4"  # Which case study we're doing. 'cs3' or 'cs4'
numGridLines = 10                   # How many gridlines we'll use for the splining
nNumRegions = 5                     # Number of reigons we're using (cs4 = 5, cs3 = 1)

#- Rip the boundary coordinates from the .yaml file -# 
bnry_file_name = "../../startup-files/iea37-boundary-" * strCase * ".yaml"
getBndryCs4YAML(bnry_file_name)

# Augment the data so it's CCW and from NW "corner"