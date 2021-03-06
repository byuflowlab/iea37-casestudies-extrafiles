using Distributed
using ClusterManagers
using Snopt
using DelimitedFiles 
using PyPlot
import ForwardDiff

using CSV
using DataFrames

using BenchmarkTools

addprocs(SlurmManager(parse(Int, ENV["SLURM_NTASKS"])-1))
@everywhere import FlowFarm; const ff = FlowFarm


function nondiscrete_boundary_wrapper(x, params)
    # include relevant globals
    params.boundary_vertices_nondiscrete
    params.boundary_normals_nondiscrete

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return boundary distances
    return ff.ray_trace_boundary(boundary_vertices_nondiscrete, boundary_normals_nondiscrete, turbine_x, turbine_y)
end

# set up boundary constraint wrapper function
function discrete_boundary_wrapper(x, params)
    # include relevant globals
    params.boundary_vertices
    params.boundary_normals

    # get number of turbines
    nturbines = Int(length(x)/2)
    
    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return boundary distances
    return ff.ray_trace_boundary(boundary_vertices, boundary_normals, turbine_x, turbine_y, discrete = true)
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
    # include relevant globals
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
    params.obj_scale
    global μ

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines] 
    turbine_y = x[nturbines+1:end]

    # calculate AEP
    AEP = obj_scale*ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set,
                rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z)
    
    # calculate discrete region penalty
    penalty = -μ*sum(max.(0, discrete_boundary_wrapper(x, params)).^2)
    
    # calculate discrete region penalty
    penalty = -μ*sum(max.(0, discrete_boundary_wrapper(x, params)).^2)

    # return the objective as an array
    return [AEP + penalty]
end

# set up optimization problem wrapper function
function wind_farm_opt(x, params)
    it = params.it

    # calculate spacing constraint value and jacobian
    spacing_con = spacing_wrapper(x)
    ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

    # calculate boundary constraint and jacobian
    boundary_con = nondiscrete_boundary_wrapper(x)
    db_dx = ForwardDiff.jacobian(nondiscrete_boundary_wrapper, x)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
    c = [spacing_con; boundary_con]
    dcdx = [ds_dx; db_dx]

    # calculate the objective function and jacobian (negative sign in order to maximize AEP)
    AEP = -aep_wrapper(x)[1]
    dAEP_dx = -ForwardDiff.jacobian(aep_wrapper,x)
    # global funcalls_AEP
    # append!(funcalls_AEP, -AEP)

    it[1] += 1
    params.funcalls_AEP[it[1]] = it[1]
    params.iter_AEP[it[1]] = -AEP

    # set fail flag to false
    fail = false

    # return objective, constraint, and jacobian values
    return AEP, c, dAEP_dx, dcdx, fail
end

# import model set with wind farm and related details
@everywhere include("./model_sets/model_set_7_ieacs4.jl")

# scale objective to be between 0 and 1
obj_scale = 1E-11

# set wind farm boundary parameters
include("boundary_normals_calculator.jl")
include("../../optimo-attempt-baker/Julia-files/baker_cs34_functions.jl")
strCase = "cs4"  # Which case study we're doing. 'cs3' or 'cs4'
bnry_file_name = "../../startup-files/iea37-boundary-" * strCase * ".yaml"
bndry_x, bndry_y = getBndryCs4YAML(bnry_file_name)
boundary_vertices = [[bndry_x[1][:] bndry_y[1][:]], [bndry_x[2][:] bndry_y[2][:]], [bndry_x[3][:] bndry_y[3][:]], [bndry_x[4][:] bndry_y[4][:]], [bndry_x[5][:] bndry_y[5][:]]]
boundary_normals = [boundary_normals_calculator(boundary_vertices[1]), boundary_normals_calculator(boundary_vertices[2]), boundary_normals_calculator(boundary_vertices[3]), boundary_normals_calculator(boundary_vertices[4]), boundary_normals_calculator(boundary_vertices[5])]
# bndry_x_clsd, bndry_y_clsd = closeBndryLists(bndry_x, bndry_y)
boundary_vertices_nondiscrete = [10363.8 6490.3; 9449.7 1602.2; 9387.0 1056.6; 9365.1 625.5; 9360.8 360.2; 9361.5 126.9; 9361.3 137.1; 7997.6 1457.9; 6098.3 3297.5; 
    5588.4 3791.3; 4670.7 4680.2; 4176.8 5158.6; 2047.8 7220.7; 1468.5 7781.7; 107.4 9100.0; 3267.1 10100.6; 4524.1 10498.7; 8953.7 11901.5; 7048.3 9531.5;
    6764.9 8399.7; 7274.9 7940.8; 7369.9 7896.2; 7455.1 7784.3; 7606.5 7713.0; 7638.9 7708.4; 8297.1 7398.9; 8450.3 6455.3; 8505.4 6422.3; 9133.0 6127.4; 
    9332.8 6072.6; 9544.2 6087.1; 9739.0 6171.2; 9894.9 6316.9; 10071.8 6552.5; 10106.9 6611.1]
boundary_normals_nondiscrete = boundary_normals_calculator(boundary_vertices_nondiscrete)

# set globals for iteration history
iter_AEP = zeros(Float64, 10000)
funcalls_AEP = zeros(Float64, 10000)

# set globals for use in wrapper functions
struct params_struct{}
    model_set
    rotor_points_y
    rotor_points_z
    turbine_z
    ambient_ti
    rotor_diameter
    boundary_vertices
    boundary_normals
    boundary_vertices_nondiscrete
    boundary_normals_nondiscrete
    obj_scale
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
    iter_AEP
    funcalls_AEP
    it
end

params = params_struct(model_set, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
    rotor_diameter, boundary_vertices, boundary_normals, boundary_vertices_nondiscrete, boundary_normals_nondiscrete, obj_scale, hub_height, turbine_yaw, 
    ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
    windresource, power_models, iter_AEP, funcalls_AEP, [0])

# initialize design variable array
x = [copy(turbine_x);copy(turbine_y)]

# penalty parameters
global μ
μ = 0.0
ρ = 5.0

# report initial objective value
println("Nturbines: ", nturbines)
println("Rotor diameter: ", rotor_diameter[1])
println("Starting AEP value (GWh): ", aep_wrapper(x, params)[1]*1e-9/obj_scale)
# println("Directional AEP at start: ", dir_aep.*1E-6)

t1 = time()
for i in 1:10
    println(i)
    aep_wrapper(x, params)[1]*1e-9/obj_scale
end
t2 = time()
at = (t2-t1)/10.0
act = at/7200.0
println("average time: ", at)
println("fcal time: ", act)

# add initial turbine location to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C0"))
end

# set general lower and upper bounds for design variables
lb = zeros(length(x)) .+ minimum(boundary_vertices_nondiscrete)
ub = zeros(length(x)) .+ maximum(boundary_vertices_nondiscrete)

# set up options for SNOPT
options = Dict{String, Any}()
options["Derivative option"] = 1
options["Verify level"] = 3
options["Major optimality tolerance"] = 1e-2 #1e-5
options["Major iteration limit"] = 1e6
options["Summary file"] = "summary-ieacs4-dist.out"
options["Print file"] = "print-ieacs4-dist.out"

# generate wrapper function surrogates
spacing_wrapper(x) = spacing_wrapper(x, params)
aep_wrapper(x) = aep_wrapper(x, params)
nondiscrete_boundary_wrapper(x) = nondiscrete_boundary_wrapper(x, params)
discrete_boundary_wrapper(x) = discrete_boundary_wrapper(x, params)
wind_farm_opt(x) = wind_farm_opt(x, params)

# run and time optimization
println
t1 = time()
iter = 1
xopt_intermediate = zeros(nturbines, 2)
while in(1,discrete_boundary_wrapper(x) .> 1e-4) && iter < 20
    println(iter)
    global x, μ, iter, xopt_intermediate, xopt
    xopt, fopt, info = snopt(wind_farm_opt, x, lb, ub, options)
    # xopt = x .+ 1.0
    x = xopt
    xopt_intermediate = cat(dims=3, xopt_intermediate, [xopt[1:nturbines] xopt[nturbines+1:end]]) 
    if μ == 0.0
        μ += 1
    else
        μ = μ*ρ
    end
    iter += 1
end
t2 = time()
clkt = t2-t1

# print optimization results
println("Finished in : ", clkt, " (s)")
# println("info: ", info)
println("end objective value: ", aep_wrapper(xopt))

# extract final turbine locations
turbine_x = copy(xopt[1:nturbines])
turbine_y = copy(xopt[nturbines+1:end])

# add final turbine locations to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C1", linestyle="--")) 
end

# add wind farm boundary to plot
plt.gcf().gca().plot([boundary_vertices[1][:,1];boundary_vertices[1][1,1]],[boundary_vertices[1][:,2];boundary_vertices[1][1,2]], color="C2")
plt.gcf().gca().plot([boundary_vertices[2][:,1];boundary_vertices[2][1,1]],[boundary_vertices[2][:,2];boundary_vertices[2][1,2]], color="C2")
plt.gcf().gca().plot([boundary_vertices[3][:,1];boundary_vertices[3][1,1]],[boundary_vertices[3][:,2];boundary_vertices[3][1,2]], color="C2")
plt.gcf().gca().plot([boundary_vertices[4][:,1];boundary_vertices[4][1,1]],[boundary_vertices[4][:,2];boundary_vertices[4][1,2]], color="C2")
plt.gcf().gca().plot([boundary_vertices[5][:,1];boundary_vertices[5][1,1]],[boundary_vertices[5][:,2];boundary_vertices[5][1,2]], color="C2")

# set up and show plot
axis("square")
xlim(0, 11000)
ylim(-500, 13000)
plt.show()
savefig("opt_plot")

# write results to csv files
dataforcsv_xopt = DataFrame(xopt = xopt)
dataforcsv_funceval = DataFrame(function_value = funcalls_AEP)
CSV.write("cs4_functionvalue_log.csv", dataforcsv_funceval)
CSV.write("cs4_xopt.csv", dataforcsv_xopt)