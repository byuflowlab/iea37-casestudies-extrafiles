import YAML
using Printf
using Parameters

# Define our coordinate structure
struct Coordinate
    xCoord::Float64
    yCoord::Float64
    Coordinate(xCoord,yCoord) = new(xCoord,yCoord)
    Base.zero(::Type{Coordinate}) = Coordinate(0,0)
end

### Complete, tested ###
function WindFrame(turb_coords, wind_dir_deg)
    ###Convert map coordinates to downwind/crosswind coordinates.
    num_turb = length(turb_coords)
    frame_coords = zeros(Coordinate, num_turb)
    # Convert from meteorological polar system (CW, 0 deg.=N)
    # to standard polar system (CCW, 0 deg.=W)
    # Shift so North comes "along" x-axis, from left to right.
    wind_dir_deg = 270.0 - wind_dir_deg
    # Convert inflow wind direction from degrees to radians
    wind_dir_rad = deg2rad(wind_dir_deg)

    # Constants to use below
    cos_dir = cos(-wind_dir_rad)
    sin_dir = sin(-wind_dir_rad)
    # Convert to downwind(x) & crosswind(y) coordinates
    for n  = 1:num_turb 
        xCoord = (turb_coords[n][1] * cos_dir) - (turb_coords[n][2] * sin_dir)
        yCoord = (turb_coords[n][1] * sin_dir) + (turb_coords[n][2] * cos_dir)
        frame_coords[n] = Coordinate(xCoord, yCoord)
    end
    return frame_coords
end

### Complete, tested ###
function DirPower(frame_coords, dir_loss, wind_speed, turb_ci, turb_co, rated_ws, rated_pwr)
    num_turb = length(frame_coords)
    # Effective windspeed is freestream multiplied by wake deficits
    ws_eff = wind_speed*(1.0 .- dir_loss)

    # By default, the turbine's power output is zero
    turb_pwr = zeros(Float64, num_turb)
    # Check to see if turbine produces power for experienced wind speed
    for n  = 1:num_turb
        if ws_eff[n] > turb_co
            continue
        elseif ws_eff[n] > rated_ws
            turb_pwr[n] = rated_pwr
        elseif ws_eff[n] > turb_ci
            turb_pwr[n] = rated_pwr * (((ws_eff[n]-turb_ci)/ (rated_ws-turb_ci))^3)
        end
    end
    # Sum the power from all turbines for this direction
    return sum(turb_pwr)
end

### Complete, tested ###
function GaussianWake(frame_coords, turb_diam)
    ###Return each turbine's total loss due to wake from upstream turbines
    # Equations and values explained in <iea37-wakemodel.pdf>
    num_turb = length(frame_coords)
    CT = 4*(1/3).*(1-(1/3))         # Constant thrust coefficient
    k = 0.0324555                   # Constant, relating to a turbulence intensity of 0.075
    loss = zeros(num_turb)          # Array holding the wake deficit seen at each turbine

    for i= 1:num_turb               # Looking at each turb (Primary)
        loss_array = zeros(num_turb)                              # Calculate the loss from all others 
        for j = 1:num_turb                                        # Looking at all other turbs (Target)
            x = frame_coords[i].xCoord - frame_coords[j].xCoord   # Calculate the x-dist
            y = frame_coords[i].yCoord - frame_coords[j].yCoord   # And the y-offset
            if x > 0                                              # If Primary is downwind of the Target
                sigma = k*x + turb_diam/sqrt(8)                   # Calculate the wake loss
                # Simplified Bastankhah Gaussian wake model
                exponent = -0.5 * (y/sigma)^2
                radical = 1 - CT/(8*sigma^2 / turb_diam^2)
                loss_array[j] = (1-sqrt(radical)) * exp(exponent)
            # Note that if the Target is upstream, loss is defaulted to zero
            end
        # Total wake losses from all upstream turbs, using sqrt of sum of sqrs
        end
        loss[i] = sqrt(sum(loss_array.^2))
    end

    return loss
end

### Complete, tested ###
function calcAEPcs3(turb_coords, wind_freq, wind_speeds, wind_speed_probs, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)
    ### Calculate the wind farm AEP.
    num_dir_bins = length(wind_freq)       # Number of bins used for our windrose
    num_speed_bins = length(wind_speeds)   # Number of wind speed bins

    # Power produced by the wind farm from each wind direction
    pwr_prod_dir = zeros(num_dir_bins)
    #  Power produced by the wind farm at a given windspeed
    pwr_prod_ws = zeros(Float64, num_dir_bins, num_speed_bins)

    # For each directional bin
    for i = 1:num_dir_bins
        # For each wind speed bin
        # Shift coordinate frame of reference to downwind/crosswind
        frame_coords = WindFrame(turb_coords, wind_dir[i])
        # Use the Simplified Bastankhah Gaussian wake model for wake deficits
        dir_loss = GaussianWake(frame_coords, turb_diam)

        for j = 1:num_speed_bins
            # Find the farm's power for the current direction and speed,
            # multiplied by the probability that the speed will occur
            pwr_prod_ws[i,j] = DirPower(frame_coords, dir_loss, wind_speeds[j], turb_ci, turb_co, rated_ws, rated_pwr) * wind_speed_probs[i][j]
        end
        pwr_prod_dir[i] = sum(pwr_prod_ws[i,:]) * wind_freq[i]
    end

    #  Convert power to AEP
    hrs_per_year = 365.0 * 24.0
    AEP = hrs_per_year * pwr_prod_dir
    AEP /= 1.E6  # Convert to MWh

    return AEP
end

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

# file_name = "iea37-ex-opt3.yaml"
# # Get turbine locations
# turb_locs = getTurbLocYAML(file_name)
# turb_coords = turb_locs[1]
# fname_turb = turb_locs[2]
# fname_wr = turb_locs[3]

# # Get turbine attributes
# turb_data = getTurbAtrbtYAML(fname_turb)
# turb_ci = turb_data[1]
# turb_co = turb_data[2]
# rated_ws = turb_data[3]
# rated_pwr = turb_data[4]
# turb_diam = turb_data[5]

# # Get windrose info
# wr_data = getWindRoseYAML(fname_wr)
# wind_dir = wr_data[1]
# wind_dir_freq = wr_data[2]
# wind_speeds = wr_data[3]
# wind_speed_probs = wr_data[4]
# num_speed_bins = wr_data[5]
# min_speed = wr_data[6]
# max_speed = wr_data[7]

# AEP = calcAEPcs3(turb_coords, wind_dir_freq, wind_speeds, wind_speed_probs, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)

# print(AEP, "\n")
# print("AEP = ", sum(AEP), "\n")