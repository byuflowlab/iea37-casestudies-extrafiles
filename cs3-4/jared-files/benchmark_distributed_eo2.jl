using Distributed
using ClusterManagers
using Snopt
using DelimitedFiles 
import ForwardDiff
using BenchmarkTools

# addprocs(SlurmManager(parse(Int, ENV["SLURM_NTASKS"])-1))
# addprocs(2)
# const IN_SLURM = "SLURM_JOBID" in keys(ENV)
@everywhere import FlowFarm; const ff = FlowFarm
# # IN_SLURM && using ClusterManagers

# # if IN_SLURM
# #     pids = addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]),dir=pwd(), tunneling=true)
# #     print("\n")
# # else
# #     pids = addprocs()
# # end
# # @sync println(pids)

# # using ClusterManagers

# if IN_SLURM
#     addprocs(SlurmManager(parse(Int, ENV["SLURM_NTASKS"])-1))
# end

# set up boundary constraint wrapper function
function boundary_wrapper(x, params)
    # include relevant globals
    params.boundary_center
    params.boundary_radius

    # get number of turbines
    nturbines = Int(length(x)/2)
    
    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return boundary distances
    return ff.circle_boundary(boundary_center, boundary_radius, turbine_x, turbine_y)
end

# set up spacing constraint wrapper function
function spacing_wrapper(x, params)
    # include relevant globals
    params.rotor_diameter

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return spacing distances
    return 2.0*rotor_diameter[1] .- ff.turbine_spacing(turbine_x,turbine_y)
end

# set up objective wrapper function
@everywhere function aep_wrapper(x, params)

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines] 
    turbine_y = x[nturbines+1:end]

    # calculate AEP
    AEP = obj_scale*ff.calculate_aep(turbine_x, turbine_y, params.turbine_z, params.rotor_diameter,
    params.hub_height, params.turbine_yaw, params.ct_models, params.generator_efficiency, params.cut_in_speed,
    params.cut_out_speed, params.rated_speed, params.rated_power, params.windresource, params.power_models, params.model_set,
                rotor_sample_points_y=params.rotor_points_y,rotor_sample_points_z=params.rotor_points_z)
    
    # return the objective as an array
    return [AEP]
end

# set up optimization problem wrapper function
function wind_farm_opt(x)

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

    # set fail flag to false
    fail = false

    # return objective, constraint, and jacobian values
    return AEP, c, dAEP_dx, dcdx, fail
end

# set globals for use in wrapper functions
struct params_struct2{MS, AF, F, I, ACTM, WR, APM}
    model_set::MS
    rotor_points_y::AF
    rotor_points_z::AF
    turbine_z::AF
    rotor_diameter::AF
    boundary_center::AF
    boundary_radius::F
    obj_scale::I
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
end

# import model set with wind farm and related details
@everywhere include("./model_sets/model_set_6.jl")
nstates = length(windresource.wind_directions)
# scale objective to be between 0 and 1
obj_scale = 1E-11

# set wind farm boundary parameters
boundary_center = [0.0,0.0]
boundary_radius = 300.0

# initialize struct for opt params
params = params_struct2(model_set, rotor_points_y, rotor_points_z, turbine_z, 
    rotor_diameter, boundary_center, boundary_radius, obj_scale, hub_height, turbine_yaw, 
    ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
    windresource, power_models)

# initialize design variable array
x = [copy(turbine_x);copy(turbine_y)]

# generate wrapper function surrogates
spacing_wrapper(x) = spacing_wrapper(x, params)
aep_wrapper(x) = aep_wrapper(x, params)
boundary_wrapper(x) = boundary_wrapper(x, params)
obj_func(x) = wind_farm_opt(x)

# make sure code for benchmarking has been pre-compiled
println("pre-compiling functions to be timed")
c1 = aep_wrapper(x)
c2 = ForwardDiff.jacobian(aep_wrapper, x)
# c2 = spacing_wrapper(x)
# c3 = boundary_wrapper(x)

# run and time optimization
println("nturbines: ", nturbines)
println("nstates: ", nstates)
println("rotor diameter: ", rotor_diameter[1])
println("starting AEP value (MWh): ", aep_wrapper(x)[1]*1E5)

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

# println("min: ", minimum(t1))
# println("max: ", maximum(t1))
# println("mean: ", mean(t1))
# println("median: ", median(t1))
# println("workers: ", nworkers())

# println()
# println("Benchmarking spacing_wrapper")
# t2 = @elapsed spacing_wrapper(x)
# # t2 = @time spacing_wrapper(x)
# println(t2)
# # # t2 = @benchmark spacing_wrapper(x) samples=1
# # # println("min: ", minimum(t2))
# # # println("max: ", maximum(t2))
# # # println("mean: ", mean(t2))
# # # println("median: ", median(t2))
# # # println("workers: ", nworkers())

# println()
# println("Benchmarking boundary_wrapper")
# t3 = @elapsed boundary_wrapper(x)
# println(t3)
# t3 = @time boundary_wrapper(x)
# # t3 = @benchmark boundary_wrapper(x) samples=1
# # println("min: ", minimum(t3))
# # println("max: ", maximum(t3))
# # println("mean: ", mean(t3))
# # println("median: ", median(t3))
# # println("workers: ", nworkers())
