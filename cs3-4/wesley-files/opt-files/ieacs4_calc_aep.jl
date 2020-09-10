function calc_aep4(yaml_opt_file)
    # calculates AEP and state directional AEPs from the ieacs4 YAML opt file

    # set model
    # include("./model_sets/model_set_7_ieacs4.jl")

    # get initial turbine x and y locations
    layout_file_name = string(yaml_opt_file)
    println(layout_file_name)
    turbine_x2, turbine_y2, fname_turb, fname_wr, AEP_provided = ff.get_turb_loc_YAML(layout_file_name)

    # calculate the number of turbines
    nturbines = length(turbine_x2)

    # set turbine base heights
    turbine_z = zeros(nturbines)

    # set turbine yaw values
    turbine_yaw = zeros(nturbines)

    # set turbine design parameters
    turbine_file_name = string("../../startup-files/",fname_turb)
    println(turbine_file_name)
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
    windrose_file_name = string("../../startup-files/",fname_wr)
    println(windrose_file_name)
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
    sorted_turbine_index2 = sortperm(turbine_x2)

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


    # create x and y coordinate vectors
    base = YAML.load(open(yaml_opt_file))
    turb_loc = base["definitions"]["position"]["items"]
    println(turb_loc)
    turbine_x = zeros(length(turb_loc))
    turbine_y = zeros(length(turb_loc))
    for i = 1:length(turbine_x)
        turbine_x[i] = turb_loc[i][1]
        turbine_y[i] = turb_loc[i][2]
    end

    # calculate AEP
    AEP = ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set,
                rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z)

    return AEP
end

function calc_dir_aep(yaml_opt_file)
    # calculates AEP from the ieacs4 YAML opt file

    # set model
    include("./model_sets/model_set_7_ieacs4.jl")

    # create x and y coordinate vectors
    base = YAML.load(open(yaml_opt_file))
    turb_loc = base["definitions"]["position"]["items"]
    turbine_x = zeros(length(turb_loc))
    turbine_y = zeros(length(turb_loc))
    for i = 1:length(turbine_x)
        turbine_x[i] = turb_loc[i][1]
        turbine_y[i] = turb_loc[i][2]
    end

    # calculate state and directional AEPs
    state_aeps = ff.calculate_state_aeps(turbine_x, turbine_y, turbine_z, rotor_diameter,
                    hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                    cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set;
                    rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0], hours_per_year=365.0*24.0)
    dir_aep = zeros(360)
    for i in 1:360
        for j in 1:20
            dir_aep[i] += state_aeps[(i-1)*20 + j]
        end
    end

    return dir_aep
end

# function create_wind_farm_plot(yaml_opt_file)
#     # creates a plot of the wind farm with turbine locations

#     # create x and y coordinate vectors
#     base = YAML.load(open(yaml_opt_file))
#     turb_loc = base["definitions"]["position"]["items"]
#     turbine_x = zeros(length(turb_loc))
#     turbine_y = zeros(length(turb_loc))
#     for i = 1:length(turbine_x)
#         turbine_x[i] = turb_loc[i][1]
#         turbine_y[i] = turb_loc[i][2]
#     end
filename = "../../submissions/iea-37-cs4-shared-files/opt_layout_BYU.yaml"
    