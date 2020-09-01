using Distributed, DelimitedFiles, YAML, ClusterManagers

addprocs(SlurmManager(parse(Int, ENV["SLURM_NTASKS"])-1))
@everywhere import FlowFarm; const ff = FlowFarm

@everywhere function aep_wrapper(turbine_x,turbine_y,params)

    params.turbine_z
    params.rotor_diameter
    params.hub_height
    params.turbine_yaw
    params.ct_models
    params.generator_efficiency
    params.cut_in_speed
    params.cut_out_speed
    params.rated_speed
    params.rated_power
    params.windresource
    params.power_models
    params.model_set
    params.rotor_points_y
    params.rotor_points_z

    AEP = ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set,
                rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z)
end

# set model
@everywhere include("./model_sets/model_set_7_ieacs4.jl")

# get layout number
layout_number = Base.parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
# layout_number = 1
println("Initial layout number: ", layout_number)

# read in x locations
# open("x_history-$layout_number.txt", "a") do io
#     writedlm(io, [rand(81)*10000 rand(81)*10000])
# end
x_history = readdlm("../results/x_history-$layout_number.txt", '\t', Float64, '\n')

# check if there are the correct number of turbine coordinates
if mod(length(x_history[:,1]),nturbines)!==0
    error("incorrect number of turbine locations")
end

struct params_struct1{}
    model_set
    rotor_points_y
    rotor_points_z
    turbine_z
    ambient_ti
    rotor_diameter
    nturbines
    hub_height
    turbine_yaw
    ct_models
    generator_efficiency
    cut_in_speed
    cut_out_speed
    rated_speed
    rated_power
    windresource
    power_models
end

params_full = params_struct1(model_set, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
    rotor_diameter, nturbines, hub_height, turbine_yaw, 
    ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
    windresource, power_models)

n_iterations = Int(length(x_history[:,1])/nturbines)
aep_history = zeros(n_iterations)
for i = 1:n_iterations
    aep_history[i] = 1e-6*aep_wrapper(x_history[(i-1)*nturbines+1:i*nturbines,1],x_history[(i-1)*nturbines+1:i*nturbines,2],params_full)
end

include("write_opt_log_YAML.jl")
filename = "../results/iea37-byu-log4-$layout_number.yaml"
write_opt_log_YAML(filename, aep_history; baseyaml=string(@__DIR__, "/default_cs4_log.yaml"), title="", titledescription="", 
gradient_based=true, algorithm_name="", program_language="Julia", total_optimizations=1, total_wall_time=[], 
units="s", aepunits="MWh")