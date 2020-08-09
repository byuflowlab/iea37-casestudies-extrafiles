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

# set up nondiscrete boundary constraint wrapper function
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

# # set up boundary constraint wrapper function
# function discrete_boundary_wrapper(x, params)
#     # include relevant globals
#     params.boundary_vertices
#     params.boundary_normals

#     # get number of turbines
#     nturbines = Int(length(x)/2)
    
#     # extract x and y locations of turbines from design variables vector
#     turbine_x = x[1:nturbines]
#     turbine_y = x[nturbines+1:end]

#     # get and return boundary distances
#     return ff.ray_trace_boundary(boundary_vertices, boundary_normals, turbine_x, turbine_y, discrete = true)
# end

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
                rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z)
    
    # return the objective as an array
    return [AEP]
end

# set up optimization problem wrapper function
function wind_farm_opt_nondiscrete(x, params)
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

    it[1] += 1
    params.funcalls_AEP[it[1]] = it[1]
    params.iter_AEP[it[1]] = -AEP

    # set fail flag to false
    fail = false

    # return objective, constraint, and jacobian values
    return AEP, c, dAEP_dx, dcdx, fail
end

# import model set with wind farm and related details
@everywhere include("./model_sets/model_set_7_ieacs4_reduced_wind_rose.jl")

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
savefig("../results/opt_plot1_test")

# set general lower and upper bounds for design variables
lb = zeros(length(x)) .+ minimum(boundary_vertices_nondiscrete)
ub = zeros(length(x)) .+ maximum(boundary_vertices_nondiscrete)

# set up options for SNOPT
options = Dict{String, Any}()
options["Derivative option"] = 1
options["Verify level"] = 3
options["Major optimality tolerance"] = 1e-5
options["Major iteration limit"] = 1e6
options["Summary file"] = "summary-ieacs4-WEC-discrete_test.out"
options["Print file"] = "print-ieacs4-WEC-discrete_test.out"

# generate wrapper function surrogates
spacing_wrapper(x) = spacing_wrapper(x, params)
aep_wrapper(x) = aep_wrapper(x, params)
nondiscrete_boundary_wrapper(x) = nondiscrete_boundary_wrapper(x, params)
wind_farm_opt_nondiscrete(x) = wind_farm_opt_nondiscrete(x, params)

# set up for WEC optimization
wec_steps = 6
wec_max = 3.0
wec_end = 1.0
wec_values = collect(LinRange(wec_max, wec_end, wec_steps))
println(wec_values)
info = fill("",wec_steps)

# start optimization timer
t1t = time()

# first, run optimization with nondiscrete boundaries and WEC=3
params.model_set.wake_deficit_model.wec_factor[1] = wec_values[1]
# xopt_nondiscrete, fopt_nondiscrete, info_nondiscrete = snopt(wind_farm_opt_nondiscrete, x, lb, ub, options)
xopt_nondiscrete = convert(Matrix,DataFrame!(CSV.File("xopt_ieacs4_WEC_discrete.csv")))[:,1]
fopt_nondiscrete = aep_wrapper(xopt_nondiscrete)[1]
info_nondiscrete = []
x = xopt_nondiscrete

# time after nondiscrete boundaries optimization
t2t = time()

# print nondiscrete boundary constraint optimization results
println("info: ", info_nondiscrete)
println("end objective value (nondiscrete): ", aep_wrapper(xopt_nondiscrete)[1])
println("locations ", x[1:5])
println("locations opt (nondiscrete) ", xopt_nondiscrete[1:5])

# add turbine locations after nondiscrete optimization to plot
clf()
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((xopt_nondiscrete[i],xopt_nondiscrete[nturbines+i]), rotor_diameter[1]/2.0, fill=false,color="C3", linestyle="--")) 
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
savefig("../results/opt_plot2_test")

# find the nearest boundary for each turbine
nearest_region = zeros(Int64, nturbines)
closed_boundary_vertices = copy(boundary_vertices)
for k = 1:length(boundary_vertices)
    closed_boundary_vertices[k] = [closed_boundary_vertices[k]; closed_boundary_vertices[k][1,1] closed_boundary_vertices[k][1,2]]
end
global nearest_region, closed_boundary_vertices
for i = 1:nturbines
    global nearest_region, closed_boundary_vertices, boundary_vertices
    nearest_region_distance = 1.0e30
    for k = 1:length(boundary_vertices)
        # get vector from turbine to the first vertex in first face
        turbine_to_first_facepoint = closed_boundary_vertices[k][1, :] - [xopt_nondiscrete[1]; xopt_nondiscrete[nturbines+1]]
        for j = 1:length(boundary_vertices[k][:,1])
            # define the vector from the turbine to the second point of the face
            turbine_to_second_facepoint = closed_boundary_vertices[k][j+1, :] - [xopt_nondiscrete[i]; xopt_nondiscrete[nturbines+i]]
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
                nearest_region[i] = k
            end
            # reset for next face iteration
            turbine_to_first_facepoint = turbine_to_second_facepoint        # (for efficiency, so we don't have to recalculate for the same vertex twice)
        end
    end
end

println()
println("closest regions: ", nearest_region)
println()

# set up discrete boundary constraint wrapper function
function discrete_boundary_wrapper(x, params)
    # include relevant globals
    params.boundary_vertices
    params.boundary_normals
    global nearest_region


    # get number of turbines
    nturbines = Int(length(x)/2)
    
    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

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

# set up optimization problem wrapper function
function wind_farm_opt_discrete(x, params)
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

    it[1] += 1
    params.funcalls_AEP[it[1]] = it[1]
    params.iter_AEP[it[1]] = -AEP

    # set fail flag to false
    fail = false

    # return objective, constraint, and jacobian values
    return AEP, c, dAEP_dx, dcdx, fail
end

# set up function wrapper surrogates for discrete boundary problem
discrete_boundary_wrapper(x) = discrete_boundary_wrapper(x, params)
wind_farm_opt_discrete(x) = wind_farm_opt_discrete(x, params)

# start time again for discrete boundary optimization
t3t = time()

# run optimization with discrete regions and WEC=3
println("input x to discrete optimization: ", x)
xopt_discrete, fopt_discrete, info_discrete = snopt(wind_farm_opt_discrete, x, lb, ub, options)
x = xopt_discrete

# stop time after discrete boundaries optimization
t4t = time()

# add turbine locations after discrete optimization to plot
clf()
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((xopt_discrete[i],xopt_discrete[nturbines+i]), rotor_diameter[1]/2.0, fill=false,color="C4", linestyle="--")) 
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
savefig("../results/opt_plot3_test")



t5t = 0
t6t = 0
t7t = 0
t8t = 0
xopt = xopt_discrete
fopt = fopt_discrete
info = info_discrete



# calculate total time
clkt = (t2t - t1t) + (t4t - t3t) + (t6t - t5t) + (t8t - t7t)

# print optimization results
println("Finished in : ", clkt, " (s)")
println("info: ", info)
println("end objective value: ", aep_wrapper(xopt))

# extract final turbine locations
turbine_x = copy(xopt[1:nturbines])
turbine_y = copy(xopt[nturbines+1:end])

# add final turbine locations to plot
clf()
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C1")) 
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
savefig("../results/opt_plot5_test")

# write results to csv files
dataforcsv_xopt = DataFrame(xopt = xopt)
dataforcsv_funceval = DataFrame(function_value = funcalls_AEP)
CSV.write("functionvalue_log_ieacs4_WEC_discrete_test.csv", dataforcsv_funceval)
CSV.write("xopt_ieacs4_WEC_discrete_test.csv", dataforcsv_xopt)