# Test for LazySets package

using Plots, LazySets

include("../../optimo-attempt-baker/Julia-files/baker_cs34_functions.jl")
strCase = "cs4"  # Which case study we're doing. 'cs3' or 'cs4'
bnry_file_name = "../../startup-files/iea37-boundary-" * strCase * ".yaml"
bndry_x, bndry_y = getBndryCs4YAML(bnry_file_name)
bndry_x_clsd, bndry_y_clsd = closeBndryLists(bndry_x, bndry_y)
pyplot()
scatter(bndry_x, bndry_y)

all_boundaries_x_except_5 = [bndry_x_clsd[:][1]; bndry_x_clsd[:][2]; bndry_x_clsd[:][3]; bndry_x_clsd[:][4]]
all_boundaries_y_except_5 = [bndry_y_clsd[:][1]; bndry_y_clsd[:][2]; bndry_y_clsd[:][3]; bndry_y_clsd[:][4]]
nvertices = length(all_boundaries_x_except_5[:])
v1 = fill(Float64[], nvertices)
for i = 1:length(all_boundaries_x_except_5[:])
    v1[i] = [all_boundaries_x[i]; all_boundaries_y[i]]
end
hull1 = convex_hull(v1)
p = plot([Singleton(vi) for vi in v1])
plot!(p, VPolygon(hull1), alpha=0.2)

all_boundaries_x_45 = [bndry_x_clsd[:][3]; bndry_x_clsd[:][5]]
all_boundaries_y_45 = [bndry_y_clsd[:][3]; bndry_y_clsd[:][5]]
nvertices = length(all_boundaries_x_45[:])
v2 = fill(Float64[], nvertices)
for i = 1:length(all_boundaries_x_45[:])
    v2[i] = [all_boundaries_x_45[i]; all_boundaries_y_45[i]]
end
hull2 = convex_hull(v2)
plot!(p, VPolygon(hull2), alpha=0.2)



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

single_boundary = [10363.8 6490.3; 9449.7 1602.2; 9387.0 1056.6; 9365.1 625.5; 9360.8 360.2; 9361.5 126.9; 9361.3 137.1; 7997.6 1457.9; 6098.3 3297.5; 
    5588.4 3791.3; 4670.7 4680.2;
    4176.8 5158.6; 2047.8 7220.7;
    1468.5 7781.7; 107.4 9100.0; 3267.1 10100.6;
    4524.1 10498.7; 8953.7 11901.5; 7048.3 9531.5;
    6764.9 8399.7;
    7274.9 7940.8; 7369.9 7896.2; 7455.1 7784.3; 7606.5 7713.0; 7638.9 7708.4; 8297.1 7398.9;
    8450.3 6455.3; 8505.4 6422.3; 9133.0 6127.4; 9332.8 6072.6; 9544.2 6087.1; 9739.0 6171.2; 9894.9 6316.9; 10071.8 6552.5; 10106.9 6611.1]
plot!(single_boundary[:,1], single_boundary[:,2])

all_boundaries = [boundary_vertices_a; boundary_vertices_b; boundary_vertices_c; boundary_vertices_d; boundary_vertices_e]
nvertices = length(all_boundaries[:,1])
v = fill(Float64[], nvertices)
for i = 1:length(all_boundaries[:,1])
    v[i] = all_boundaries[i,:]
end
hull = convex_hull(v)

# points = N -> [randn(2) for i in 1:N]
# v = points(30)
# hull = convex_hull(v)

p = Plots.plot([Singleton(vi) for vi in v])
Plots.plot!(p, VPolygon(hull), alpha=0.2)