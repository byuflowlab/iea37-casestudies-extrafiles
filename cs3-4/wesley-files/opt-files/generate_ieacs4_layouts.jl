function ieacs4_random_start(nturbines, boundary_vertices, boundary_normals, rotor_diameter; min_spacing=1.0)

    xrange = [minimum(boundary_vertices[:,1]); maximum(boundary_vertices[:,1])]
    yrange = [minimum(boundary_vertices[:,2]); maximum(boundary_vertices[:,2])]
    center = [(xrange[2]-xrange[1])/2 (yrange[2]-yrange[1])/2]

    min_spacing *= rotor_diameter
    locations = zeros(Float64,nturbines, 2)

    # generate random points within the wind farm boundary
    count = 0
    for i = 1:nturbines
        println(i)
        good_point = false
        while !good_point && count < 1e3
            count += 1
            # generate random point in the containing rectangle
            locations[i,:] = (rand(1,2).-0.5).*[xrange[2]-xrange[1] yrange[2]-yrange[1]] .+ center
            distances = ff.ray_trace_boundary(boundary_vertices, boundary_normals, locations[1:i,1], locations[1:i,2])
            # determine if the point is inside the wind farm boundary
            good_point = true
            for j = 1:length(distances)
                if distances[j] > 0.0
                    good_point = false
                end
            end
            # determine if the point is far enough away from other points
            n_bad_spacings = 0.0
            for turb = 1:nturbines
                if turb != i
                    spacing = sqrt((locations[turb,1]-locations[i,1])^2 + (locations[turb,2]-locations[i,2])^2)
                    if spacing < min_spacing
                        n_bad_spacings += 1
                    end
                end
                if n_bad_spacings > 0
                    good_point = false
                end
            end
        end
    end

    turbine_x = locations[:,1]
    turbine_y = locations[:,2]

    return turbine_x, turbine_y
end

import FlowFarm; const ff = FlowFarm
using PyPlot
include("boundary_normals_calculator.jl")
ieacs4_boundary_vertices_nondiscrete = [10363.8 6490.3; 9449.7 1602.2; 9387.0 1056.6; 9365.1 625.5; 9360.8 360.2; 9361.5 126.9; 9361.3 137.1; 7997.6 1457.9; 6098.3 3297.5; 
    5588.4 3791.3; 4670.7 4680.2; 4176.8 5158.6; 2047.8 7220.7; 1468.5 7781.7; 107.4 9100.0; 3267.1 10100.6; 4524.1 10498.7; 8953.7 11901.5; 7048.3 9531.5;
    6764.9 8399.7; 7274.9 7940.8; 7369.9 7896.2; 7455.1 7784.3; 7606.5 7713.0; 7638.9 7708.4; 8297.1 7398.9; 8450.3 6455.3; 8505.4 6422.3; 9133.0 6127.4; 
    9332.8 6072.6; 9544.2 6087.1; 9739.0 6171.2; 9894.9 6316.9; 10071.8 6552.5; 10106.9 6611.1]
ieacs4_boundary_normals_nondiscrete = boundary_normals_calculator(ieacs4_boundary_vertices_nondiscrete)
include("../../optimo-attempt-baker/Julia-files/baker_cs34_functions.jl")
strCase = "cs4"  # Which case study we're doing. 'cs3' or 'cs4'
bnry_file_name = "../../startup-files/iea37-boundary-" * strCase * ".yaml"
bndry_x, bndry_y = getBndryCs4YAML(bnry_file_name)
ieacs4_boundary_vertices_discrete = [[bndry_x[1][:] bndry_y[1][:]], [bndry_x[2][:] bndry_y[2][:]], [bndry_x[3][:] bndry_y[3][:]], [bndry_x[4][:] bndry_y[4][:]], [bndry_x[5][:] bndry_y[5][:]]]
ieacs4_boundary_normals_discrete = [boundary_normals_calculator(ieacs4_boundary_vertices_discrete[1]), boundary_normals_calculator(ieacs4_boundary_vertices_discrete[2]), boundary_normals_calculator(ieacs4_boundary_vertices_discrete[3]), boundary_normals_calculator(ieacs4_boundary_vertices_discrete[4]), boundary_normals_calculator(ieacs4_boundary_vertices_discrete[5])]

nturbines = 81
rotor_diameter = 198.0

start_layout_number = 101
end_layout_number = 120

for layout_number = start_layout_number:end_layout_number

    turbine_x, turbine_y = ieacs4_random_start(nturbines, ieacs4_boundary_vertices_nondiscrete, ieacs4_boundary_normals_nondiscrete, rotor_diameter; min_spacing=2.5)

    clf()
    for i = 1:length(turbine_x)
        plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C0", linestyle="--")) 
    end

    # add wind farm boundary to plot
    plt.gcf().gca().plot([ieacs4_boundary_vertices_discrete[1][:,1];ieacs4_boundary_vertices_discrete[1][1,1]],[ieacs4_boundary_vertices_discrete[1][:,2];ieacs4_boundary_vertices_discrete[1][1,2]], color="C2")
    plt.gcf().gca().plot([ieacs4_boundary_vertices_discrete[2][:,1];ieacs4_boundary_vertices_discrete[2][1,1]],[ieacs4_boundary_vertices_discrete[2][:,2];ieacs4_boundary_vertices_discrete[2][1,2]], color="C2")
    plt.gcf().gca().plot([ieacs4_boundary_vertices_discrete[3][:,1];ieacs4_boundary_vertices_discrete[3][1,1]],[ieacs4_boundary_vertices_discrete[3][:,2];ieacs4_boundary_vertices_discrete[3][1,2]], color="C2")
    plt.gcf().gca().plot([ieacs4_boundary_vertices_discrete[4][:,1];ieacs4_boundary_vertices_discrete[4][1,1]],[ieacs4_boundary_vertices_discrete[4][:,2];ieacs4_boundary_vertices_discrete[4][1,2]], color="C2")
    plt.gcf().gca().plot([ieacs4_boundary_vertices_discrete[5][:,1];ieacs4_boundary_vertices_discrete[5][1,1]],[ieacs4_boundary_vertices_discrete[5][:,2];ieacs4_boundary_vertices_discrete[5][1,2]], color="C2")

    # set up and show plot
    axis("square")
    xlim(0, 11000)
    ylim(-500, 13000)
    savefig("../initial-layouts/ieacs4_initial_layout-"*lpad(layout_number,3,"0")*".png")

    # save initial layout to YAML
    ff.write_turb_loc_YAML("../initial-layouts/ieacs4_initial_layout-"*lpad(layout_number,3,"0")*".yaml", turbine_x, turbine_y; title="IEA Case Study 4 Initial Wind Farm Layout $layout_number", titledescription="81 randomly generated turbine locations", 
    turbinefile="", locunits="m", wakemodelused="", windresourcefile="", aeptotal=[], 
    aepdirs=[], aepunits="MWh", baseyaml="default_cs4.yaml")
end
