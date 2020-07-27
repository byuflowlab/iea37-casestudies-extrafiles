using Distributed
using ClusterManagers
using Snopt
using DelimitedFiles 
using PyPlot
import ForwardDiff

using CSV
using DataFrames

using BenchmarkTools


# addprocs(SlurmManager(parse(Int, ENV["SLURM_NTASKS"])-1))
# addprocs(3)
# using YAML
@everywhere import FlowFarm; const ff = FlowFarm
#using ClusterManagers

#allocatedcpus = run("echo $SLURM_JOB_NODELIST")
#addprocs_slurm(allocatedcpus)
#println(allocatedcpus)

# const IN_SLURM = "SLURM_JOBID" in keys(ENV)

# IN_SLURM && using ClusterManagers

# if IN_SLURM
#     pids = addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]),dir=pwd())
#     print("\n")
# else
#     pids = addprocs()
# end
# @sync println(pids)
# set up boundary constraint wrapper function
@everywhere function boundary_wrapper(x, params)
    # include relevant globals
    params.boundary_vertices_a
    params.boundary_normals_a
    params.boundary_vertices_b
    params.boundary_normals_b
    params.boundary_vertices_c
    params.boundary_normals_c
    params.boundary_vertices_d
    params.boundary_normals_d
    params.boundary_vertices_e
    params.boundary_normals_e
    
    # get number of turbines
    nturbines = Int(length(x)/2)
    
    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    boundcon_a = ff.ray_trace_boundary(boundary_vertices_a, boundary_normals_a, turbine_x[1:31], turbine_y[1:31])
    boundcon_b = ff.ray_trace_boundary(boundary_vertices_b, boundary_normals_b, turbine_x[32:42], turbine_y[32:42])
    boundcon_c = ff.ray_trace_boundary(boundary_vertices_c, boundary_normals_c, turbine_x[43:58], turbine_y[43:58])
    boundcon_d = ff.ray_trace_boundary(boundary_vertices_d, boundary_normals_d, turbine_x[59:72], turbine_y[59:72])
    boundcon_e = ff.ray_trace_boundary(boundary_vertices_e, boundary_normals_e, turbine_x[73:81], turbine_y[73:81])

    # get and return boundary distances
    return [boundcon_a; boundcon_b; boundcon_c; boundcon_d; boundcon_e]
end

# set up spacing constraint wrapper function
@everywhere function spacing_wrapper(x, params)
    # include relevant globals
    params.rotor_diameter

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    spacecon_a = 2.0*rotor_diameter[1] .- ff.turbine_spacing(turbine_x[1:31], turbine_y[1:31])
    spacecon_b = 2.0*rotor_diameter[1] .- ff.turbine_spacing(turbine_x[32:42], turbine_y[32:42])
    spacecon_c = 2.0*rotor_diameter[1] .- ff.turbine_spacing(turbine_x[43:58], turbine_y[43:58])
    spacecon_d = 2.0*rotor_diameter[1] .- ff.turbine_spacing(turbine_x[59:72], turbine_y[59:72])
    spacecon_e = 2.0*rotor_diameter[1] .- ff.turbine_spacing(turbine_x[73:81], turbine_y[73:81])

    # get and return spacing distances
    return [spacecon_a; spacecon_b; spacecon_c; spacecon_d; spacecon_e]
end

# set up objective wrapper function
@everywhere function aep_wrapper(x, params)
    # include relevant globals
    turbine_z = params.turbine_z
    rotor_diameter = params.rotor_diameter
    hub_height = params.hub_height
    turbine_yaw = params.turbine_yaw
    ct_models = params.ct_models
    generator_efficiency = params.generator_efficiency
    cut_in_speed = params.cut_in_speed
    cut_out_speed = params.cut_out_speed
    rated_speed = params.rated_speed
    rated_power = params.rated_power
    windresource = params.windresource
    power_models = params.power_models
    model_set = params.model_set
    rotor_points_y = params.rotor_points_y
    rotor_points_z = params.rotor_points_z
    obj_scale = params.obj_scale

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines] 
    turbine_y = x[nturbines+1:end]

    # calculate AEP
    AEP = obj_scale*ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set,
                rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z, hours_per_year=365.0*24.0)

    # return the objective as an array
    return [AEP]
end

# set up optimization problem wrapper function
function wind_farm_opt(x, params)
    it = params.it

    # calculate spacing constraint value and jacobian
    spacing_con = spacing_wrapper(x)
    ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

    # calculate boundary constraint and jacobian
    boundary_con = boundary_wrapper(x)
    db_dx = ForwardDiff.jacobian(boundary_wrapper, x)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
    c = [spacing_con; boundary_con]
    dcdx = [ds_dx; db_dx]

    # calculate the objective function and jacobian (negative sign in order to maximize AEP)
    AEP = -aep_wrapper(x)[1]
    dAEP_dx = -ForwardDiff.jacobian(aep_wrapper,x)

    it[1] += 1
    params.funcalls_AEP[it[1]] = it[1]
    params.iter_AEP[it[1]] = -AEP

    # set fail flag to false
    fail = false

    # return objective, constraint, and jacobian values
    return AEP, c, dAEP_dx, dcdx, fail
end

# set globals for iteration history
iter_AEP = zeros(Float64, 10000)
funcalls_AEP = zeros(Float64, 10000)

# import model set with wind farm and related details
include("./model_sets/model_set_7_ieacs4.jl")

# scale objective to be between 0 and 1
obj_scale = 1E-12

# set wind farm boundary parameters
boundary_vertices_a = [10363.8 6490.3; 9449.7 1602.2; 9387.0 1056.6; 9365.1 625.5; 9360.8 360.2; 9361.5 126.9; 9361.3 137.1; 7997.6 1457.9; 6098.3 3297.5;
    8450.3 6455.3; 8505.4 6422.3; 9133.0 6127.4; 9332.8 6072.6; 9544.2 6087.1; 9739.0 6171.2; 9894.9 6316.9; 10071.8 6552.5; 10106.9 6611.1]
boundary_normals_a = [0.9829601758936983 -0.1838186405319916; 0.9934614633172962 -0.11416795042154541; 0.9987121579438882 -0.050734855622757584; 
    0.9998686751666075 -0.01620593781838486; 0.9999954987444023 0.0030004151269687495; -0.9998078216567232 -0.019604074934516894; -0.6957179389375846 -0.718315076718037; 
    -0.6957275377423737 -0.7183057797532565; -0.8019887481131871 0.5973391397688945; 0.5138086803485797 0.8579047965820281; 0.4252760929807897 0.905063668886888; 
    0.2645057513093967 0.9643841078762402; -0.0684295708121141 0.9976559496331737; -0.39636379138742883 0.9180935381958544; -0.6828023205475376 0.7306031693435896; 
    -0.7996740386176392 0.6004343694034798; -0.8578802011411015 0.5138497450520954; 0.42552559023380465 0.9049463918134445]
boundary_vertices_b = [5588.4 3791.3; 4670.7 4680.2; 7274.9 7940.8; 7369.9 7896.2; 7455.1 7784.3; 7606.5 7713.0; 7638.9 7708.4; 8297.1 7398.9]
boundary_normals_b = [-0.6957460043611584 -0.7182878931288504; -0.7813688797257963 0.6240694462926818; 0.4249708760634733 0.9052070230051488; 0.7956275395848184 0.6057861159305391; 
    0.4260560153872896 0.9046967844268629; 0.14056568619461773 0.9900713549359138; 0.4255255464063141 0.9049464124220882; 0.7996806883794807 -0.6004255129763556]
boundary_vertices_c = [3267.1 10100.6; 4164.1 9586.6; 5749.8 9068.6; 6054.7 8925.3; 1468.5 7781.7; 107.4 9100.0]
boundary_normals_c = [0.49718026396417986 0.8676472699919642; 0.31052117525343714 0.9505664625470563; 0.42535384615162936 0.9050271297392228; 0.24194817066179167 -0.9702891747893577; 
    -0.6957228969594285 -0.7183102746351193; -0.30189947425802094 0.9533397649540959]
boundary_vertices_d = [6764.9 8399.7; 4176.8 5158.6; 2047.8 7220.7]
boundary_normals_d = [0.7814306689309158 -0.6239920749930895; -0.6957310325444781 -0.7183023947855072; -0.24248239299288069 0.9701558066045093]
boundary_vertices_e = [8953.7 11901.5; 7048.3 9531.5; 6127.7 9962.7; 4578.1 10464.9; 4524.1 10498.7]
boundary_normals_e = [0.7793586677376737 -0.6265780613955122; -0.4241667101838764 -0.9055841219742026; -0.30829751674447764 -0.9512899879475178; -0.5305632140423848 -0.847645371546978; -0.3019099610801309 0.9533364439695956]

# set globals for use in wrapper functions
struct params_struct{MS, AF, F, ACTM, WR, APM, AF2, AI}
    model_set::MS
    rotor_points_y::AF
    rotor_points_z::AF
    turbine_z::AF
    rotor_diameter::AF
    obj_scale::F
    hub_height::AF
    turbine_yaw::AF
    ct_models::ACTM
    generator_efficiency::AF
    cut_in_speed::AF
    cut_out_speed::AF
    rated_speed::AF
    rated_power::AF
    windresource::WR
    power_models::APM

    boundary_vertices_a::AF2
    boundary_normals_a::AF2
    boundary_vertices_b::AF2
    boundary_normals_b::AF2
    boundary_vertices_c::AF2
    boundary_normals_c::AF2
    boundary_vertices_d::AF2
    boundary_normals_d::AF2
    boundary_vertices_e::AF2
    boundary_normals_e::AF2
    iter_AEP::AF
    funcalls_AEP::AF
    it::AI
end

paramemeters = params_struct(model_set, rotor_points_y, rotor_points_z, turbine_z, 
    rotor_diameter, obj_scale, hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed, 
    cut_out_speed, rated_speed, rated_power, windresource, power_models, boundary_vertices_a, boundary_normals_a, boundary_vertices_b, boundary_normals_b,
    boundary_vertices_c, boundary_normals_c, boundary_vertices_d, boundary_normals_d, boundary_vertices_e,
    boundary_normals_e, iter_AEP, funcalls_AEP, [0])

# initialize design variable array
x = [copy(turbine_x);copy(turbine_y)]

# get number of design variables
n_designvariables = length(x)

# get number of constraints
n_spacingconstraints = length(spacing_wrapper(x, paramemeters))
n_boundaryconstraints = length(boundary_wrapper(x, paramemeters))
n_constraints = n_spacingconstraints + n_boundaryconstraints

# set general lower and upper bounds for design variables
lb = ones(n_designvariables) * -Inf
ub = ones(n_designvariables) * Inf

# set lower and upper bounds for constraints
lb_g = ones(n_constraints) * -Inf
ub_g = zeros(n_constraints)

# generate wrapper function surrogates
spacing_wrapper(x) = spacing_wrapper(x, paramemeters)
aep_wrapper(x) = aep_wrapper(x, paramemeters)
boundary_wrapper(x) = boundary_wrapper(x, paramemeters)
wind_farm_opt(x) = wind_farm_opt(x, paramemeters)

c1 = aep_wrapper(x)
c2 = ForwardDiff.jacobian(aep_wrapper, x)

# # start benchmarking
println()
println("Benchmarking aep_wrapper")
println("using ", nworkers(), " worker(s)")
println("nturbines: ", nturbines)
println("nstates: ", nstates)
println("AEP value (MWh): ", aep_wrapper(x)[1]*1E5)
t1 = @belapsed aep_wrapper(x) 
println("aep call time (s): ", t1)
t2 = @belapsed ForwardDiff.jacobian(aep_wrapper, x)
println("aep jac call time (s): ", t2)
# start benchmarking
# println()
# println("Benchmarking aep_wrapper")
# t1 = @benchmark aep_wrapper(x) samples=1
# println("min: ", minimum(t1))
# println("max: ", maximum(t1))
# println("mean: ", mean(t1))
# println("median: ", median(t1))
# println("workers: ", nworkers())

# println()
# println("Benchmarking spacing_wrapper")
# t2 = @benchmark spacing_wrapper(x) samples=1
# println("min: ", minimum(t2))
# println("max: ", maximum(t2))
# println("mean: ", mean(t2))
# println("median: ", median(t2))
# println("workers: ", nworkers())

# println()
# println("Benchmarking boundary_wrapper")
# t3 = @benchmark boundary_wrapper(x) samples=1
# println("min: ", minimum(t3))
# println("max: ", maximum(t3))
# println("mean: ", mean(t3))
# println("median: ", median(t3))
# println("workers: ", nworkers())
