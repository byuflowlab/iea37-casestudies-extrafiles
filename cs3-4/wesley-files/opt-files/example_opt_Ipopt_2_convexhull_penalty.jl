using Ipopt
using DelimitedFiles 
using PyPlot
using LazySets
import ForwardDiff

### uses a Julia interface to the Ipopt nonlinear solver

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
    return ff.ray_trace_boundary(boundary_hull_vertices, boundary_hull_normals, turbine_x, turbine_y)
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
    return ff.ray_trace_boundary(boundary_vertices, boundary_normals, turbine_x, turbine_y)
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
function aep_wrapper(x, params)
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

    # return the objective as an array
    return [AEP + penalty]
end

# objective function
function obj(x) 
    return -aep_wrapper(x)[1]
end

# constraint function
function con(x, g)
    # calculate spacing constraint value
    spacing_con = spacing_wrapper(x)

    # calculate non-discrete boundary constraint
    nondiscrete_boundary_con = nondiscrete_boundary_wrapper(x)

    # calculate discrete boudnary constraint
    # discrete_boundary_con = discrete_boundary_wrapper(x)

    # combine constraint values and jacobians into overall constaint value and jacobian arrays
    g[:] = [spacing_con; nondiscrete_boundary_con]#; discrete_boundary_con]
end

# objective gradient function
function obj_grad(x, grad_f)
    grad_f[:] = -ForwardDiff.jacobian(aep_wrapper,x)
end

# constraint gradients function
function con_grad(x, mode, rows, cols, values)
    if mode == :Structure
        # report the sparcity structure of the jacobian
        for i = 1:prob.m
            rows[(i-1)*prob.n+1:i*prob.n] = i*ones(Int,prob.n)
            cols[(i-1)*prob.n+1:i*prob.n] = 1:prob.n
        end
    else
        # calculate spacing constraint jacobian
        ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

        # calculate nondiscrete boundary constraint jacobian
        db1_dx = ForwardDiff.jacobian(nondiscrete_boundary_wrapper, x)

        # calculate discrete boundary constraint jacobian
        db2_dx = ForwardDiff.jacobian(discrete_boundary_wrapper, x)

        # combine constaint jacobians into overall constaint jacobian arrays
        for i = 1:prob.m
            for j = 1:prob.n
                values[(i-1)*prob.n+j] = [ds_dx; db1_dx; db2_dx][i, j]
            end
        end
    end
end

# import model set with wind farm and related details
include("./model_sets/model_set_6.jl")

# scale objective to be between 0 and 1
obj_scale = 1E-11

# set wind farm boundary parameters
# boundary_vertices = ([0 0; 1 0; 1 .75; .75 .75; .75 1; 0 1] .- .5).*500 # Utah-shape boundary
# boundary_normals = [0 1.0; -1 0; 0 -1; -1 0; 0 -1; 1 0]
boundary_vertices = ([0 1.5; .25 1.25; .75 1.25; .75 .875; 1.125 .5; .5 0; 1.75 0; 1.25 .25; 1.5 .5; 1 .75; 1.25 1; 1 1.5] .-.75).*350
boundary_normals = [1.0 1.0; 0 1; 1 0; 1 1; 4 -5; 0 1; -1 -2; -1 1; -1 -2; -1 1; -2 -1; 0 -1]
for i = 1:length(boundary_normals[:,1])
    boundary_normals[i,:] = boundary_normals[i,:] ./ sqrt(sum(boundary_normals[i,:].^2))
end
turbine_x .+= 100.0

# get the convex hull of wind farm boundary
all_boundary_vertices = boundary_vertices
nvertices = length(all_boundary_vertices[:,1])
v = fill(Float64[], nvertices)
for i = 1:nvertices
    v[i] = all_boundary_vertices[i,:]
end
boundary_hull_vertices_array = convex_hull(v)
boundary_hull_vertices = zeros(length(boundary_hull_vertices_array), 2)
for i = 1:length(boundary_hull_vertices_array)
    boundary_hull_vertices[i,:] = [boundary_hull_vertices_array[i][1] boundary_hull_vertices_array[i][2]]
end
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

global μ

params = params_struct2(model_set, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
    rotor_diameter, boundary_vertices, boundary_normals, boundary_hull_vertices, boundary_hull_normals, obj_scale, hub_height, turbine_yaw, 
    ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
    windresource, power_models)

# initialize design variable array
x = [copy(turbine_x);copy(turbine_y)]

# penalty parameters
μ = 0.0
ρ = 5.0

# report initial objective value
println("starting objective value: ", aep_wrapper(x, params)[1])

# add initial turbine location to plot
clf()
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C0"))
end

# get number of design variables
n_designvariables = length(x)

# get number of constraints
function numberofspacingconstraints(nturb)
    # calculates number of spacing constraints needed for given number of turbines
    ncon = 0
    for i = 1:nturb-1; ncon += i; end
    return ncon
end
n_spacingconstraints = numberofspacingconstraints(nturbines)
n_nondiscrete_boundaryconstraints = length(nondiscrete_boundary_wrapper(x, params))
n_discrete_boundaryconstraints = length(discrete_boundary_wrapper(x, params))
n_constraints = n_spacingconstraints + n_nondiscrete_boundaryconstraints# + n_discrete_boundaryconstraints

# set general lower and upper bounds for design variables
lb = ones(n_designvariables) * -Inf
ub = ones(n_designvariables) * Inf

# set lower and upper bounds for constraints
lb_g = ones(n_constraints) * -Inf
ub_g = zeros(n_constraints)

# create the problem
prob = createProblem(n_designvariables, lb, ub, n_constraints, lb_g, ub_g, n_designvariables*n_constraints, 0,
    obj, con, obj_grad, con_grad)
addOption(prob, "hessian_approximation", "limited-memory")
prob.x = x

# generate wrapper function surrogates
spacing_wrapper(x) = spacing_wrapper(x, params)
nondiscrete_boundary_wrapper(x) = nondiscrete_boundary_wrapper(x, params)
discrete_boundary_wrapper(x) = discrete_boundary_wrapper(x, params)
aep_wrapper(x) = aep_wrapper(x, params)

# run and time optimization
t1 = time()
global iter
iter = 1
while in(1,discrete_boundary_wrapper(prob.x) .> 1e-4) && iter < 20
    global μ
    global ρ
    global iter
    println(prob.x)
    status = solveProblem(prob)
    if μ == 0.0
        μ += 1
    else
        μ = μ*ρ
    end
    iter += 1
end

t2 = time()
clkt = t2-t1
xopt = prob.x
fopt = prob.obj_val
# info = Ipopt.ApplicationReturnStatus[status]

# print optimization results
println("Finished in : ", clkt, " (s)")
# println("info: ", info)
println("end objective value: ", aep_wrapper(xopt)[1])

# extract final turbine locations
turbine_x = copy(xopt[1:nturbines])
turbine_y = copy(xopt[nturbines+1:end])

# add final turbine locations to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C1", linestyle="--")) 
end

# add wind farm boundary to plot
plt.gcf().gca().plot([boundary_vertices[:,1];boundary_vertices[1,1]],[boundary_vertices[:,2];boundary_vertices[1,2]], color="C2")

# set up and show plot
axis("square")
xlim(minimum(boundary_vertices) - (maximum(boundary_vertices)-minimum(boundary_vertices))/5, maximum(boundary_vertices) + (maximum(boundary_vertices)-minimum(boundary_vertices))/5)
ylim(minimum(boundary_vertices) - (maximum(boundary_vertices)-minimum(boundary_vertices))/5, maximum(boundary_vertices) + (maximum(boundary_vertices)-minimum(boundary_vertices))/5)
savefig("opt_plot")
plt.show()
