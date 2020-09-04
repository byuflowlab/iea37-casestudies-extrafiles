using Snopt
using DelimitedFiles 
using PyPlot
using LazySets
import ForwardDiff
using Distributed
using ClusterManagers
using CSV
using DataFrames
using BenchmarkTools

addprocs(SlurmManager(parse(Int, ENV["SLURM_NTASKS"])-1))
@everywhere import FlowFarm; const ff = FlowFarm

# set up discrete boundary constraint wrapper function
function discrete_boundary_wrapper(x, params)
    # include relevant globals
    params.boundary_vertices
    params.boundary_normals
    params.nearest_region

    # get number of turbines
    nturbines = Int(length(x)/2)
    
    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    display(nearest_region)
    println(nearest_region)
    display(nearest_region.==1)
    println(nearest_region.==1)
    
    turbine_x_region_1 = turbine_x[nearest_region.==1]
    turbine_y_region_1 = turbine_y[nearest_region.==1]
    turbine_x_region_2 = turbine_x[nearest_region.==2]
    turbine_y_region_2 = turbine_y[nearest_region.==2]
    turbine_x_region_3 = turbine_x[nearest_region.==3]
    turbine_y_region_3 = turbine_y[nearest_region.==3]
    turbine_x_region_4 = turbine_x[nearest_region.==4]
    turbine_y_region_4 = turbine_y[nearest_region.==4]
    turbine_x_region_5 = turbine_x[nearest_region.==5]
    turbine_y_region_5 = turbine_y[nearest_region.==5]

    if in(1,nearest_region)
        boundcon_region_1 = ff.ray_trace_boundary(boundary_vertices[1], boundary_normals[1], turbine_x_region_1, turbine_y_region_1)
    else
        boundcon_region_1 = []
    end
    if in(2,nearest_region)
        boundcon_region_2 = ff.ray_trace_boundary(boundary_vertices[2], boundary_normals[2], turbine_x_region_2, turbine_y_region_2)
    else
        boundcon_region_2 = []
    end
    if in(3,nearest_region)
        boundcon_region_3 = ff.ray_trace_boundary(boundary_vertices[3], boundary_normals[3], turbine_x_region_3, turbine_y_region_3)
    else
        boundcon_region_3 = []
    end
    if in(4,nearest_region)
        boundcon_region_4 = ff.ray_trace_boundary(boundary_vertices[4], boundary_normals[4], turbine_x_region_4, turbine_y_region_4)
    else
        boundcon_region_4 = []
    end
    if in(5,nearest_region)
        boundcon_region_5 = ff.ray_trace_boundary(boundary_vertices[5], boundary_normals[5], turbine_x_region_5, turbine_y_region_5)
    else
        boundcon_region_5 = []
    end

    # get and return boundary distances
    return [boundcon_region_1; boundcon_region_2; boundcon_region_3; boundcon_region_4; boundcon_region_5]
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

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines] 
    turbine_y = x[nturbines+1:end]

    # calculate AEP
    AEP = obj_scale*ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set,
                rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z,hours_per_year=365.0*24.0)
    
    # return the objective as an array
    return [AEP]
end

# set up optimization problem wrapper function
function wind_farm_opt_discrete(x, layout_number, params)
    it = params.it

    # calculate spacing constraint value and jacobian
    spacing_con = spacing_wrapper(x)
    ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

    # calculate boundary constraint and jacobian
    boundary_con = discrete_boundary_wrapper(x)
    db_dx = ForwardDiff.jacobian(discrete_boundary_wrapper, x)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
    c = [spacing_con; boundary_con]
    dcdx = [ds_dx; db_dx]

    # calculate the objective function and jacobian (negative sign in order to maximize AEP)
    AEP = -aep_wrapper(x)[1]
    dAEP_dx = -ForwardDiff.jacobian(aep_wrapper,x)

    # add to x log
    it[1] += 1
    params.funcalls_AEP_WEC[it[1]] = -AEP*1e-6/obj_scale
    open("../results/x_history-"*lpad(layout_number,3,"0")*".txt", "a") do io
        writedlm(io, [x[1:nturbines] x[nturbines+1:end]])
    end

    # set fail flag to false
    fail = false

    # return objective, constraint, and jacobian values
    return AEP, c, dAEP_dx, dcdx, fail
end

# get slurm variables
# layout_number = Base.parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
# layout_number = 1
layout_number = Base.parse(Int, Base.ARGS[1])
println("Initial layout number: ", layout_number)

# import model set with full wind rose
@everywhere include("./model_sets/model_set_7_ieacs4.jl")

# set boundary global parameters
include("boundary_normals_calculator.jl")
include("../../optimo-attempt-baker/Julia-files/baker_cs34_functions.jl")
strCase = "cs4"  # Which case study we're doing. 'cs3' or 'cs4'
bnry_file_name = "../../startup-files/iea37-boundary-" * strCase * ".yaml"
bndry_x, bndry_y = getBndryCs4YAML(bnry_file_name)
boundary_vertices = [[bndry_x[1][:] bndry_y[1][:]], [bndry_x[2][:] bndry_y[2][:]], [bndry_x[3][:] bndry_y[3][:]], [bndry_x[4][:] bndry_y[4][:]], [bndry_x[5][:] bndry_y[5][:]]]
boundary_normals = [boundary_normals_calculator(boundary_vertices[1]), boundary_normals_calculator(boundary_vertices[2]), boundary_normals_calculator(boundary_vertices[3]), boundary_normals_calculator(boundary_vertices[4]), boundary_normals_calculator(boundary_vertices[5])]
boundary_vertices_nondiscrete = [10363.8 6490.3; 9449.7 1602.2; 9387.0 1056.6; 9365.1 625.5; 9360.8 360.2; 9361.5 126.9; 9361.3 137.1; 7997.6 1457.9; 6098.3 3297.5; 
    5588.4 3791.3; 4670.7 4680.2; 4176.8 5158.6; 2047.8 7220.7; 1468.5 7781.7; 107.4 9100.0; 3267.1 10100.6; 4524.1 10498.7; 8953.7 11901.5; 7048.3 9531.5;
    6764.9 8399.7; 7274.9 7940.8; 7369.9 7896.2; 7455.1 7784.3; 7606.5 7713.0; 7638.9 7708.4; 8297.1 7398.9; 8450.3 6455.3; 8505.4 6422.3; 9133.0 6127.4; 
    9332.8 6072.6; 9544.2 6087.1; 9739.0 6171.2; 9894.9 6316.9; 10071.8 6552.5; 10106.9 6611.1]
boundary_normals_nondiscrete = boundary_normals_calculator(boundary_vertices_nondiscrete)

# import turbine locations from previous optimization
intermediate_yaml = YAML.load(open("../results/iea37-byu-opt4-intermediate-"*lpad(layout_number,3,"0")*".yaml"))
nturbines = length(intermediate_yaml["definitions"]["position"]["items"])
x = zeros(nturbines*2)
for i = 1:2, j = 1:nturbines
    x[(i-1)*nturbines+j] = intermediate_yaml["definitions"]["position"]["items"][j][i]
end

# set up for WEC optimization
wec_steps = 6
wec_max = 3.0
wec_end = 1.0
wec_values = collect(LinRange(wec_max, wec_end, wec_steps))
println(wec_values)
noptimizations = length(wec_values) + 2

# set other globals for iteration history
nearest_region = readdlm("nearest_region_ieacs4-"*lpad(layout_number,3,"0")*".txt", '\t', Int, '\n')
println("closest regions: ", nearest_region)
obj_scale = 1E-12
funcalls_AEP_WEC = zeros(Float64, 50000*8)

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
    nturbines
    nearest_region
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
    funcalls_AEP_WEC
    it
end

params_full = params_struct(model_set, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
    rotor_diameter, boundary_vertices, boundary_normals, boundary_vertices_nondiscrete, 
    boundary_normals_nondiscrete, nturbines, nearest_region, obj_scale, hub_height, turbine_yaw, 
    ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
    windresource, power_models, funcalls_AEP_WEC, [0])

# initialize xopt array
nturbines = params_full.nturbines
noptimizations = wec_steps + 2
xopt_all = zeros(2*nturbines,noptimizations)

# report initial objective value
println("Nturbines: ", nturbines)
println("Rotor diameter: ", rotor_diameter[1])
println("Starting AEP value (GWh): ", aep_wrapper(x, params_full)[1]*1e-9/obj_scale)
println()

t1 = time()
for i in 1:10
    println(i)
    aep_wrapper(x, params_full)[1]*1e-9/obj_scale
end
t2 = time()
at = (t2-t1)/10.0
act = at/360.0
println("average time: ", at)
println("fcal time: ", act)
println()

# set general lower and upper bounds for design variables
lb = zeros(length(x)) .+ minimum(boundary_vertices_nondiscrete)
ub = zeros(length(x)) .+ maximum(boundary_vertices_nondiscrete)

# set up options for SNOPT
tol = Base.parse(Float64, Base.ARGS[2])
options = Dict{String, Any}()
options["Derivative option"] = 1
options["Verify level"] = 1
options["Major optimality tolerance"] = tol
options["Major iteration limit"] = 1e6
options["Summary file"] = "snopt-output/summary-ieacs4-WEC-"*lpad(layout_number,3,"0")*"-discrete-fullwindrose.out"
options["Print file"] = "snopt-output/print-ieacs4-WEC-"*lpad(layout_number,3,"0")*"-discrete-fullwindrose.out"
println()
println("Objective scaling factor: ", params_full.obj_scale)
println("Major optimality tolerance: ", tol)
println()

# generate wrapper function surrogates
spacing_wrapper(x) = spacing_wrapper(x, params_full)
aep_wrapper(x) = aep_wrapper(x, params_full)
discrete_boundary_wrapper(x) = discrete_boundary_wrapper(x, params_full)
wind_farm_opt_discrete(x) = wind_farm_opt_discrete(x, layout_number, params_full)

# start time again for discrete boundary optimization
t1t = time()

# run optimization with discrete regions
println()
println("x input into snopt: ", x[1:10])
xopt, fopt, info = snopt(wind_farm_opt_discrete, x, lb, ub, options)
        # xopt = deepcopy(x).+100
        # fopt = 50.0
        # info = []
println("xopt output after snopt: ", xopt)
println()
xopt = deepcopy(xopt)

# stop time after discrete boundaries optimization
t2t = time()

# calculate total time
clkt = (t2t - t1t)

# print optimization results
println("Finished in : ", clkt, " (s)")
println("info: ", info)
println("end objective value (snopt output): ", -fopt)
println("initial locations ", x[1:10])
println("optimal locations (WEC = " * "$(wec_values[1])" * ") ", xopt[1:10])

final_objective = aep_wrapper(xopt, params_full)[1]
final_AEP = final_objective*1e-9/obj_scale
println("end objective value: ", final_objective)
println("Ending AEP value (GWh): ", final_AEP)
println()

# add turbine locations after discrete optimization to plot
clf()
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((xopt[:][i],xopt[:][nturbines+i]), rotor_diameter[1]/2.0, fill=false,color="C4", linestyle="--")) 
end

# add wind farm boundary to plot
plt.gcf().gca().plot([boundary_vertices[1][:,1];boundary_vertices[1][1,1]],[boundary_vertices[1][:,2];boundary_vertices[1][1,2]], color="C2")
plt.gcf().gca().plot([boundary_vertices[2][:,1];boundary_vertices[2][1,1]],[boundary_vertices[2][:,2];boundary_vertices[2][1,2]], color="C2")
plt.gcf().gca().plot([boundary_vertices[3][:,1];boundary_vertices[3][1,1]],[boundary_vertices[3][:,2];boundary_vertices[3][1,2]], color="C2")
plt.gcf().gca().plot([boundary_vertices[4][:,1];boundary_vertices[4][1,1]],[boundary_vertices[4][:,2];boundary_vertices[4][1,2]], color="C2")
plt.gcf().gca().plot([boundary_vertices[5][:,1];boundary_vertices[5][1,1]],[boundary_vertices[5][:,2];boundary_vertices[5][1,2]], color="C2")

# set up plot window
axis("square")
xlim(0, 11000)
ylim(-500, 13000)

# save current figure
savefig("../results/opt_plot-"*lpad(layout_number,3,"0")*"-5")

# extract final turbine locations
turbine_x = copy(xopt[1:nturbines])
turbine_y = copy(xopt[nturbines+1:end])

# calculate state and directional AEPs
state_aeps = 1e-6*ff.calculate_state_aeps(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set;
                rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0], hours_per_year=365.0*24.0)
dir_aep = zeros(360)
for i in 1:360
    for j in 1:20
        dir_aep[i] += state_aeps[(i-1)*20 + j]
    end
end

# write results to yaml files
ff.write_turb_loc_YAML("../results/iea37-byu-opt4-"*lpad(layout_number,3,"0")*".yaml",turbine_x,turbine_y,
    title="IEA Wind Task 37 case study 4, BYU's  optimal layout",
    titledescription="BYU's optimal layout for the 81 turbine wind plant model for IEA Task 37 case study 4",
    turbinefile="iea37-10mw.yaml",
    locunits="m",
    wakemodelused="iea37-aepcalc.py",
    windresourcefile="iea37-windrose-cs4.yaml",
    aeptotal=final_AEP*1e3,
    aepdirs=dir_aep,
    aepunits="MWh",
    baseyaml="default_cs4.yaml")