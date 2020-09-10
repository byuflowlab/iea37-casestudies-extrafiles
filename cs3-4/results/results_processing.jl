# This file is for processing the yaml results files for IEA CS4

# import required packages
using FlowFarm; const ff = FlowFarm
using PyPlot

# set file names to read
submissions_directory = "../submissions/iea-37-cs4-shared-files/"
startupfiles_directory = "../startup-files/"
extension = ".yaml"
files = [
        "opt_layout_BYU",
        "ifpen_layout",
        "opt_layout_DavidBieniek",
        "opt_layout_SebastianSanchez",
        "nrel_optimal_results",
        "opt_layout-Tilli-ADREMOG", 
        "opt_layout-Quaeghebeur-ADREMOG+PG", 
        "opt_layout_nickR"
        ]

# set numper of turbines
nTurbines = 81

# set base AEP
AEPT = 3446.535439743956 # Gwh
AEPI = 2851.096412523507 # GWh

# calculate base case wake loss
WLI = 100.0*(1 - AEPI/AEPT)

# initialize arrays to store AEP and wake loss
AEPO = zeros(length(files), 2)
WLO = zeros(length(files), 2)

# initialize array to store x and y locations
locations = zeros((length(files), 2*nTurbines))

# initialize array to store number of turbines in each regions
region_count = zeros(Int64, length(files), 5)

# get boundary coordinates
include("boundary_normals_calculator.jl")
include("../optimo-attempt-baker/Julia-files/baker_cs34_functions.jl")
strCase = "cs4"  # Which case study we're doing. 'cs3' or 'cs4'
bnry_file_name = "../startup-files/iea37-boundary-" * strCase * ".yaml"
bndry_x, bndry_y = getBndryCs4YAML(bnry_file_name)
boundary_vertices = [[bndry_x[1][:] bndry_y[1][:]], [bndry_x[2][:] bndry_y[2][:]], [bndry_x[3][:] bndry_y[3][:]], [bndry_x[4][:] bndry_y[4][:]], [bndry_x[5][:] bndry_y[5][:]]]
boundary_normals = [boundary_normals_calculator(boundary_vertices[1]), boundary_normals_calculator(boundary_vertices[2]), boundary_normals_calculator(boundary_vertices[3]), boundary_normals_calculator(boundary_vertices[4]), boundary_normals_calculator(boundary_vertices[5])]

# loop through results
for i = 1:length(files)

    println("Load ", files[i])

    ################## read yaml ##################################
    # get initial turbine x and y locations
    layout_file_name = string(submissions_directory, files[i], extension)
    turbine_x, turbine_y, fname_turb, fname_wr, AEP_provided = ff.get_turb_loc_YAML(layout_file_name)

    # calculate the number of turbines
    nturbines = length(turbine_x)

    # set turbine base heights
    turbine_z = zeros(nturbines)

    # set turbine yaw values
    turbine_yaw = zeros(nturbines)

    # set turbine design parameters
    turbine_file_name = string(startupfiles_directory,fname_turb)
    turb_ci, turb_co, rated_ws, rated_pwr, turb_diam, turb_hub_height = ff.get_turb_atrbt_YAML(turbine_file_name)

    rotor_diameter = zeros(nturbines) .+ turb_diam # m
    hub_height = zeros(nturbines) .+ turb_hub_height   # m
    cut_in_speed = zeros(nturbines) .+ turb_ci  # m/s
    cut_out_speed = zeros(nturbines) .+ turb_co  # m/s
    rated_speed = zeros(nturbines) .+ rated_ws # m/s
    rated_power = zeros(nturbines) .+ rated_pwr # W
    generator_efficiency = zeros(nturbines) .+ 1.0

    # rotor swept area sample points (normalized by rotor radius)
    rotor_points_y = [0.0]
    rotor_points_z = [0.0]

    # set flow parameters
    windrose_file_name = string(startupfiles_directory,fname_wr)
    winddirections, windspeeds, windprobabilities, ambient_ti = ff.get_wind_rose_YAML(windrose_file_name)

    nstates = length(winddirections)
    winddirections *= pi/180.0

    air_density = 1.1716  # kg/m^3
    shearexponent = 0.15
    ambient_tis = zeros(nstates) .+ ambient_ti
    measurementheight = zeros(nstates) .+ turb_hub_height

    # initialize power model
    power_model = ff.PowerModelPowerCurveCubic()
    power_models = Vector{typeof(power_model)}(undef, nturbines)
    for i = 1:nturbines
        power_models[i] = power_model
    end

    # load thrust curve
    ct = 4.0*(1.0/3.0)*(1.0 - 1.0/3.0)

    # initialize thurst model
    ct_model = ff.ThrustModelConstantCt(ct)
    ct_models = Vector{typeof(ct_model)}(undef, nturbines)
    for i = 1:nturbines
        ct_models[i] = ct_model
    end

    # initialize wind shear model
    wind_shear_model = ff.PowerLawWindShear(shearexponent)

    # get sorted indecies 
    sorted_turbine_index = sortperm(turbine_x)

    # initialize the wind resource definition
    windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_tis, wind_shear_model)

    # set up wake and related models
    k = 0.0324555
    wakedeficitmodel = ff.GaussSimple(k)

    wakedeflectionmodel = ff.JiminezYawDeflection()
    wakecombinationmodel = ff.SumOfSquaresFreestreamSuperposition()
    localtimodel = ff.LocalTIModelNoLocalTI()

    # initialize model set
    model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)

    # recalculate AEP
    AEP_calculated = 1e-6*ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set,
                rotor_sample_points_y=rotor_points_y, rotor_sample_points_z=rotor_points_z, hours_per_year=365.0*24.0)

    # store AEP
    AEPO[i, 1] = deepcopy(AEP_provided)
    AEPO[i, 2] = deepcopy(AEP_calculated)

    # store wake loss percentage
    WLO[i, 1] = 100.0*(1 - deepcopy(AEP_provided)/(AEPT*1e3))
    WLO[i, 2] = 100.0*(1 - deepcopy(AEP_calculated)/(AEPT*1e3))

    # store x and y locations
    locations[i, 1:nTurbines] = deepcopy(turbine_x)
    locations[i, nTurbines+1:end] = deepcopy(turbine_y)

    # plot locations with name and AEP overlaid

    # add wind farm boundary to plot
    clf()
    for i = 1:length(boundary_vertices)
        plt.gcf().gca().plot([boundary_vertices[i][:,1];boundary_vertices[i][1,1]],[boundary_vertices[i][:,2];boundary_vertices[i][1,2]], color="C2",zorder=0)
    end

    # add final turbine layout to plot
    for j = 1:nTurbines
        plt.gcf().gca().add_artist(plt.Circle((turbine_x[j],turbine_y[j]), rotor_diameter[1]/2.0, fill=false,color="C0",zorder=10))
    end

    # set up plot window
    axis("square")
    xlim(-500, 11000)
    ylim(-500, 13000)
    title("A$i", fontsize=25)
    xticks(range(0,stop=11000,step=2500))
    yticks(range(0,stop=13000,step=2500))
    xlabel("x (m)", fontsize=13)
    ylabel("y (m)", fontsize=13)

    # save the figure
    savefig("A$(i)_" * files[i] * ".pdf", transparent=true)

    # find the nearest boundary for each turbine
    closed_boundary_vertices = copy(boundary_vertices)
    nearest_region = zeros(Int64, nturbines)
    for k = 1:length(boundary_vertices)
        closed_boundary_vertices[k] = [closed_boundary_vertices[k]; closed_boundary_vertices[k][1,1] closed_boundary_vertices[k][1,2]]
    end
    for m = 1:nturbines
        nearest_region_distance = 1.0e30
        for k = 1:length(boundary_vertices)
            # get vector from turbine to the first vertex in first face
            turbine_to_first_facepoint = closed_boundary_vertices[k][1, :] - [locations[i,m]; locations[i,nturbines+m]]
            for j = 1:length(boundary_vertices[k][:,1])
                # define the vector from the turbine to the second point of the face
                turbine_to_second_facepoint = closed_boundary_vertices[k][j+1, :] - [locations[i,m]; locations[i,nturbines+m]]
                # find perpendicular distance from turbine to current face (vector projection)
                boundary_vector = closed_boundary_vertices[k][j+1, :] - closed_boundary_vertices[k][j, :]
                # check if perpendicular distance is the shortest
                if sum(boundary_vector .* -turbine_to_first_facepoint) > 0 && sum(boundary_vector .* turbine_to_second_facepoint) > 0
                    # perpendicular distance from turbine to face
                    turbine_to_face_distance = abs(sum(turbine_to_first_facepoint .* boundary_normals[k][j,:]))
                # check if distance to first facepoint is shortest
                elseif sum(boundary_vector .* -turbine_to_first_facepoint) <= 0
                    # distance from turbine to first facepoint
                    turbine_to_face_distance = sqrt(sum(turbine_to_first_facepoint.^2))
                # distance to second facepoint is shortest
                else
                    # distance from turbine to second facepoint
                    turbine_to_face_distance = sqrt(sum(turbine_to_second_facepoint.^2))
                end
                if turbine_to_face_distance < nearest_region_distance
                    nearest_region_distance = turbine_to_face_distance
                    nearest_region[m] = k
                end
                # reset for next face iteration
                turbine_to_first_facepoint = turbine_to_second_facepoint        # (for efficiency, so we don't have to recalculate for the same vertex twice)
            end
        end
    end

    for j = 1:5
        region_count[i,j] = sum(nearest_region .== j)
    end

end