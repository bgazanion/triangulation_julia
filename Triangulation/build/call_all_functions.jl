using Base.Filesystem

include("../src/Triangulation.jl")

# create_nodes
nodes = Triangulation.create_nodes(10, [0., 1.], [0., 1.])

# in_circumcircle
Triangulation.in_circumcircle([1, 2, 3], 4, nodes)

# is_clockwise
Triangulation.is_clockwise([0., 0.], [1.,0.], [0., 1.])

# create_supertriangle
Triangulation.create_supertriangle(nodes)

# require_flip
AB = Triangulation.Edge([1, 2], [])
BC = Triangulation.Edge([2, 3], [])
CD = Triangulation.Edge([3, 4], [])
BD = Triangulation.Edge([2, 4], [])
AD = Triangulation.Edge([1, 4], [])
edges = [AB, BC, CD, BD, AD]
triangles = [Triangulation.Triangle([1, 4, 5]), Triangulation.Triangle([2, 3, 4])]
width = 1.
height = width/tan(pi/3)
nodes = [[0., -height], [width, 0.],  [0., height], [-width, 0.]]
Triangulation.require_flip(1, 2, triangles, edges, nodes) == true

# compute_angle
nodes_2 = [[0., 0.], [1., 0.], [0., tan(pi/3)]]
edge1 = Triangulation.Edge([1, 2], [])
edge2 = Triangulation.Edge([1, 3], [])
Triangulation.compute_angle(edge1, edge2, nodes_2)

# get_tri_nodes
edges = [
    Triangulation.Edge([1, 2], []),
    Triangulation.Edge([2, 3], []), 
    Triangulation.Edge([3, 1], [])
]
triangles = [Triangulation.Triangle([1, 2, 3])]
Triangulation.get_tri_node_ids(1, triangles, edges)

# inside_triangle
Triangulation.inside_triangle([1, 2, 3], 4, nodes)

# get_triangulation_nodes
Triangulation.get_triangulation_nodes(edges)

# flip_triangulation!
AB = Triangulation.Edge([1, 2], [1, 0])
BC = Triangulation.Edge([2, 3], [2, 0])
CD = Triangulation.Edge([3, 4], [2, 0])
BD = Triangulation.Edge([2, 4], [1, 2])
AD = Triangulation.Edge([1, 4], [1, 0])
edges = Dict{Int, Triangulation.Edge}(1=>AB, 2=>BC, 3=>CD, 4=>BD, 5=>AD)
triangles = Dict{Int, Triangulation.Triangle}(
    1 => Triangulation.Triangle([1, 4, 5]),
    2 => Triangulation.Triangle([2, 3, 4])
)
width = 1.
height = width/tan(pi/3)
nodes = [[0., -height], [width, 0.],  [0., height], [-width, 0.]]
Triangulation.flip_triangulation!(triangles, edges, nodes, true, "./_tmp_.txt", init_flip_step=0)

# plot_step
Triangulation.plot_step("./_tmp_flip_1.txt", false)

# triangulate
nodes = Triangulation.create_nodes(5, [0., 1.], [0., 1.])
Triangulation.triangulate(nodes, false, "./_tmp_.txt")

# parse_args
Triangulation.parse_args("-i 5 --long=3 --long2 3")

# print_help
Triangulation.print_help()

# get_arg_value
Triangulation.get_arg_value(["c", "d"], Dict("c" => "ok", "d" => "ok"))

# get_main_parameters
Triangulation.get_main_parameters("-c -o myOutput -h -p -s")
Triangulation.get_main_parameters("--create-nodes --output=myOutput --help --plot --save-steps")
Filesystem.touch("./_tmp_test.txt")
Triangulation.get_main_parameters("-f ./_tmp_test.txt")
Triangulation.get_main_parameters("-l ./_tmp_test*.txt")

# main
Triangulation.main(["--help"])
Triangulation.main(["-l ./results*.txt"])

# julia_main
Triangulation.julia_main()

# remove tmp files
Filesystem.rm("./_tmp_test.txt")
Filesystem.rm("./_tmp_flip_1.txt")
Filesystem.rm("./_tmp_flip_1.png")






