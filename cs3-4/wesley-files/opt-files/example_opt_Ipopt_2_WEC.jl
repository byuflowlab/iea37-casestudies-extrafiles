using Ipopt
using DelimitedFiles 
using PyPlot
using LazySets
import ForwardDiff
using Distributed
using ClusterManagers
using CSV
using DataFrames

### uses a Julia interface to the Ipopt nonlinear solver

# set up for distributed processing
addprocs(SlurmManager(parse(Int, ENV["SLURM_NTASKS"])-1))
@everywhere import FlowFarm; const ff = FlowFarm

# set up nondiscrete boundary constraint wrapper function
function nondiscrete_boundary_wrapper(x, params)
    # include relevant globals
    params.boundary_hull_vertices
    params.boundary_hull_normals

    # get number of turbines
    nturbines = Int(length(x)/2)
    
    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return boundary distances
    return ff.convex_boundary(boundary_hull_vertices, boundary_hull_normals, turbine_x, turbine_y)
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
#     return ff.ray_trace_boundary(boundary_vertices, boundary_normals, turbine_x, turbine_y, discrete=true)
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

# objective function
function obj(x) 
    return -aep_wrapper(x)[1]
end

# initial constraint function (convex hull boundary)
function con_nondiscrete(x, g)
    # calculate spacing constraint value
    spacing_con = spacing_wrapper(x)

    # calculate boundary constraint
    boundary_con = nondiscrete_boundary_wrapper(x)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
    g[:] = [spacing_con; boundary_con]
end

# constraint function (discrete boundaries)
function con_discrete(x, g)
    # calculate spacing constraint value
    spacing_con = spacing_wrapper(x)

    # calculate boundary constraint
    boundary_con = discrete_boundary_wrapper(x)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
    g[:] = [spacing_con; boundary_con]
end

# objective gradient function
function obj_grad(x, grad_f)
    grad_f[:] = -ForwardDiff.jacobian(aep_wrapper,x)
end

# initial constraint gradients function (convex hull boundary)
function con_grad_nondiscrete(x, mode, rows, cols, values)
    if mode == :Structure
        # report the sparcity structure of the jacobian
        for i = 1:prob_nondiscrete.m
            rows[(i-1)*prob_nondiscrete.n+1:i*prob_nondiscrete.n] = i*ones(Int,prob_nondiscrete.n)
            cols[(i-1)*prob_nondiscrete.n+1:i*prob_nondiscrete.n] = 1:prob_nondiscrete.n
        end
    else
        # calculate spacing constraint jacobian
        ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

        # calculate boundary constraint jacobian
        db_dx = ForwardDiff.jacobian(nondiscrete_boundary_wrapper, x)

        # combine constaint jacobians into overall constaint jacobian arrays
        for i = 1:prob_nondiscrete.m
            for j = 1:prob_nondiscrete.n
                values[(i-1)*prob_nondiscrete.n+j] = [ds_dx; db_dx][i, j]
            end
        end
    end
end

# constraint gradients function (discrete regions)
function con_grad_discrete(x, mode, rows, cols, values)
    if mode == :Structure
        # report the sparcity structure of the jacobian
        for i = 1:prob_discrete.m
            rows[(i-1)*prob_discrete.n+1:i*prob_discrete.n] = i*ones(Int,prob_discrete.n)
            cols[(i-1)*prob_discrete.n+1:i*prob_discrete.n] = 1:prob_discrete.n
        end
    else
        # calculate spacing constraint jacobian
        ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

        # calculate boundary constraint jacobian
        db_dx = ForwardDiff.jacobian(discrete_boundary_wrapper, x)

        # combine constaint jacobians into overall constaint jacobian arrays
        for i = 1:prob_discrete.m
            for j = 1:prob_discrete.n
                values[(i-1)*prob_discrete.n+j] = [ds_dx; db_dx][i, j]
            end
        end
    end
end

# import model set with wind farm and related details
@everywhere include("./model_sets/model_set_6.jl")

# scale objective to be between 0 and 1
obj_scale = 1E-11

# set wind farm boundary parameters
boundary_vertices_utah = ([0 0; 1 0; 1 .75; .75 .75; .75 1; 0 1] .+ [-1.15 -.5]) .* 350
boundary_normals_utah = [0 1.0; -1 0; 0 -1; -1 0; 0 -1; 1 0]
boundary_vertices_zigzag = ([0 1.5; .25 1.25; .75 1.25; .75 .875; 1.125 .5; .5 0; 1.75 0; 1.25 .25; 1.5 .5; 1 .75; 1.25 1; 1 1.5] .+ [.15 -.5]) .* 250
boundary_normals_zigzag = [0.7071067811865475 0.7071067811865475; 0.0 1.0; 1.0 0.0; 0.7071067811865475 0.7071067811865475; 0.6246950475544243 -0.7808688094430304; 0.0 1.0; -0.4472135954999579 -0.8944271909999159; -0.7071067811865475 0.7071067811865475; -0.4472135954999579 -0.8944271909999159; -0.7071067811865475 0.7071067811865475; -0.8944271909999159 -0.4472135954999579; 0.0 -1.0]
boundary_vertices = [boundary_vertices_utah, boundary_vertices_zigzag]
boundary_normals = [boundary_normals_utah, boundary_normals_zigzag]
global combined_boundary_vertices
combined_boundary_vertices = zeros(0,2)
for m = 1:length(boundary_vertices)
    global combined_boundary_vertices
    combined_boundary_vertices = [combined_boundary_vertices; boundary_vertices[:][m]]
end

# shift turbine starting positions
turbine_x .*= 1.0
turbine_x .+= -200.0
turbine_y .*= 1.0
turbine_y .+= 0.0

# get the convex hull of wind farm boundary
nvertices = length(combined_boundary_vertices[:,1])
v = fill(Float64[], nvertices)
for i = 1:nvertices
    v[i] = combined_boundary_vertices[i,:]
end
boundary_hull_vertices_array = convex_hull(v)
boundary_hull_vertices = zeros(length(boundary_hull_vertices_array), 2)
for i = 1:length(boundary_hull_vertices_array)
    boundary_hull_vertices[i,:] = [boundary_hull_vertices_array[i][1] boundary_hull_vertices_array[i][2]]
end
boundary_hull_vertices = reverse(boundary_hull_vertices, dims=1)
include("boundary_normals_calculator.jl")
boundary_hull_normals = boundary_normals_calculator(boundary_hull_vertices)

# set globals for use in wrapper functions
struct params_struct2{}
    model_set
    rotor_points_y
    rotor_points_z
    turbine_z
    ambient_ti
    rotor_diameter
    boundary_vertices
    boundary_normals
    boundary_hull_vertices
    boundary_hull_normals
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
end

params = params_struct2(model_set, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
    rotor_diameter, boundary_vertices, boundary_normals, boundary_hull_vertices, boundary_hull_normals, obj_scale, hub_height, turbine_yaw, 
    ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
    windresource, power_models)

# initialize design variable array
x = [copy(turbine_x);copy(turbine_y)]
global x

# report initial objective value
println("starting objective value: ", aep_wrapper(x, params)[1])

# add initial turbine location to plot
clf()
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C0"))
end

# add boundaries to plot
plt.gcf().gca().plot([boundary_vertices[1][:,1];boundary_vertices[1][1,1]],[boundary_vertices[1][:,2];boundary_vertices[1][1,2]], color="C2")
plt.gcf().gca().plot([boundary_vertices[2][:,1];boundary_vertices[2][1,1]],[boundary_vertices[2][:,2];boundary_vertices[2][1,2]], color="C2")
plt.gcf().gca().plot([boundary_hull_vertices[:,1];boundary_hull_vertices[1,1]],[boundary_hull_vertices[:,2];boundary_hull_vertices[1,2]], color="C3")

# set up plot window
axis("square")
# xlim(minimum(combined_boundary_vertices) - (maximum(combined_boundary_vertices)-minimum(combined_boundary_vertices))/5, maximum(combined_boundary_vertices) + (maximum(combined_boundary_vertices)-minimum(combined_boundary_vertices))/5)
# ylim(minimum(combined_boundary_vertices) - (maximum(combined_boundary_vertices)-minimum(combined_boundary_vertices))/5, maximum(combined_boundary_vertices) + (maximum(combined_boundary_vertices)-minimum(combined_boundary_vertices))/5)
xlim(-600, 600)
ylim(-400, 400)

# save current figure
savefig("../results/opt_plot1")

# get number of design variables
n_designvariables = length(x)

# get number of constraints for initial optimization (convex hull boundary)
function numberofspacingconstraints(nturb)
    # calculates number of spacing constraints needed for given number of turbines
    ncon = 0
    for i = 1:nturb-1; ncon += i; end
    return ncon
end
n_spacingconstraints = numberofspacingconstraints(nturbines)
n_boundaryconstraints_nondiscrete = length(nondiscrete_boundary_wrapper(x, params))
n_constraints_nondiscrete = n_spacingconstraints + n_boundaryconstraints_nondiscrete

# set general lower and upper bounds for design variables
lb = ones(n_designvariables) * -Inf
ub = ones(n_designvariables) * Inf

# set lower and upper bounds for constraints (convex hull boundary)
lb_g_nondiscrete = ones(n_constraints_nondiscrete) * -Inf
ub_g_nondiscrete = zeros(n_constraints_nondiscrete)

# generate wrapper function surrogates
spacing_wrapper(x) = spacing_wrapper(x, params)
aep_wrapper(x) = aep_wrapper(x, params)
nondiscrete_boundary_wrapper(x) = nondiscrete_boundary_wrapper(x, params)

# create the problems
prob_nondiscrete = createProblem(n_designvariables, lb, ub, n_constraints_nondiscrete, lb_g_nondiscrete, ub_g_nondiscrete, n_designvariables*n_constraints_nondiscrete, 0,
    obj, con_nondiscrete, obj_grad, con_grad_nondiscrete)
addOption(prob_nondiscrete, "hessian_approximation", "limited-memory")

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
prob_nondiscrete.x = x
status_nondiscrete = solveProblem(prob_nondiscrete)
xopt_nondiscrete = prob_nondiscrete.x
fopt_nondiscrete = prob_nondiscrete.obj_val
info_nondiscrete = Ipopt.ApplicationReturnStatus[status_nondiscrete]

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

# add boundaries to plot
plt.gcf().gca().plot([boundary_vertices[1][:,1];boundary_vertices[1][1,1]],[boundary_vertices[1][:,2];boundary_vertices[1][1,2]], color="C2")
plt.gcf().gca().plot([boundary_vertices[2][:,1];boundary_vertices[2][1,1]],[boundary_vertices[2][:,2];boundary_vertices[2][1,2]], color="C2")
plt.gcf().gca().plot([boundary_hull_vertices[:,1];boundary_hull_vertices[1,1]],[boundary_hull_vertices[:,2];boundary_hull_vertices[1,2]], color="C3")

# set up plot window
axis("square")
xlim(-600, 600)
ylim(-400, 400)

# save current figure
savefig("../results/opt_plot2")

# find the nearest boundary for each turbine
nearest_region = zeros(Int64, nturbines)
size(boundary_hull_vertices)
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
println(size_boundary_vertices)

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

    boundcon_region_1 = ff.ray_trace_boundary(boundary_vertices[1], boundary_normals[1], turbine_x_region_1, turbine_y_region_1)
    boundcon_region_2 = ff.ray_trace_boundary(boundary_vertices[2], boundary_normals[2], turbine_x_region_2, turbine_y_region_2)

    # get and return boundary distances
    return [boundcon_region_1; boundcon_region_2]
end

# generate discrete regions wrapper surrogate
discrete_boundary_wrapper(x) = discrete_boundary_wrapper(x, params)

# get number of constraints for regular optimization (discrete boundaries)
n_boundaryconstraints_discrete = length(discrete_boundary_wrapper(x, params))
n_constraints_discrete = n_spacingconstraints + n_boundaryconstraints_discrete

# set lower and upper bounds for constraints (discrete boundaries)
lb_g_discrete = ones(n_constraints_discrete) * -Inf
ub_g_discrete = zeros(n_constraints_discrete)

# set up discrete regions optimization problem
prob_discrete = createProblem(n_designvariables, lb, ub, n_constraints_discrete, lb_g_discrete, ub_g_discrete, n_designvariables*n_constraints_discrete, 0,
obj, con_discrete, obj_grad, con_grad_discrete)
addOption(prob_discrete, "hessian_approximation", "limited-memory")
prob_discrete.x = prob_nondiscrete.x

# start time again for discrete boundary optimization
t3t = time()

# run optimization with discrete regions and WEC=3
status_discrete = solveProblem(prob_discrete)
xopt_discrete = prob_discrete.x
fopt_discrete = prob_discrete.obj_val
info_discrete = Ipopt.ApplicationReturnStatus[status_discrete]

# time after discrete boundaries optimization
t4t = time()

# add turbine locations after discrete optimization to plot
clf()
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((xopt_discrete[i],xopt_discrete[nturbines+i]), rotor_diameter[1]/2.0, fill=false,color="C4", linestyle="--")) 
end

# add boundaries to plot
plt.gcf().gca().plot([boundary_vertices[1][:,1];boundary_vertices[1][1,1]],[boundary_vertices[1][:,2];boundary_vertices[1][1,2]], color="C2")
plt.gcf().gca().plot([boundary_vertices[2][:,1];boundary_vertices[2][1,1]],[boundary_vertices[2][:,2];boundary_vertices[2][1,2]], color="C2")
plt.gcf().gca().plot([boundary_hull_vertices[:,1];boundary_hull_vertices[1,1]],[boundary_hull_vertices[:,2];boundary_hull_vertices[1,2]], color="C3")

# set up plot window
axis("square")
xlim(-600, 600)
ylim(-400, 400)

# save current figure
savefig("../results/opt_plot3")

# start time again for WEC optimization
t5t = time()

# optimization with decreasing WEC values
x = prob_discrete.x
for i in 1:length(wec_values)
    global xopt
    println("Running with WEC = ", wec_values[i])
    params.model_set.wake_deficit_model.wec_factor[1] = wec_values[i]
    
    println(prob_discrete.x)
    x_initial = deepcopy(prob_discrete.x)
    t1 = time()
    status = solveProblem(prob_discrete)
    t2 = time()
    clk = t2-t1
    xopt = prob_discrete.x
    fopt = prob_discrete.obj_val
    info = Ipopt.ApplicationReturnStatus[status]

    # print optimization results
    println("Finished in : ", clk, " (s)")
    println("info: ", info)
    println("end objective value: ", -fopt)
    println("locations ", x_initial[1:5])
    println("locations opt ", xopt[1:5])
end

# stop time for WEC optimization
t6t = time()

# stop time
clkt = (t2t - t1t) + (t4t - t3t) + (t6t - t5t)

# print optimization results
println("Finished in : ", clkt, " (s)")
println("info: ", info)
println("end objective value: ", aep_wrapper(xopt)[1])

# extract final turbine locations
turbine_x = copy(xopt[1:nturbines])
turbine_y = copy(xopt[nturbines+1:end])

# add final turbine locations to plot
clf()
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C1", linestyle="--")) 
end

# add boundaries to plot
plt.gcf().gca().plot([boundary_vertices[1][:,1];boundary_vertices[1][1,1]],[boundary_vertices[1][:,2];boundary_vertices[1][1,2]], color="C2")
plt.gcf().gca().plot([boundary_vertices[2][:,1];boundary_vertices[2][1,1]],[boundary_vertices[2][:,2];boundary_vertices[2][1,2]], color="C2")
plt.gcf().gca().plot([boundary_hull_vertices[:,1];boundary_hull_vertices[1,1]],[boundary_hull_vertices[:,2];boundary_hull_vertices[1,2]], color="C3")

# set up plot window
axis("square")
xlim(-600, 600)
ylim(-400, 400)

# save current figure
savefig("../results/opt_plot4")

# write results to csv files
dataforcsv_xopt = DataFrame(xopt = xopt)
# dataforcsv_funceval = DataFrame(function_value = funcalls_AEP)
# CSV.write("ex2_functionvalue_log.csv", dataforcsv_funceval)
CSV.write("ex2_xopt.csv", dataforcsv_xopt)