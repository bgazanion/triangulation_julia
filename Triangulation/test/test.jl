using Test

include("../src/lib.jl")

# Unit tests for functions related to delaunay triangulation
function test_delaunay()
    # create_nodes
    @testset "create_nodes" begin
        n = 100
        x_bounds = [-2, 3]
        y_bounds = [-4, 5]
        nodes = create_nodes(n, x_bounds, y_bounds)
        @test length(nodes) == (n)
        nodes_x = [a_node[1] for a_node in nodes]
        nodes_y = [a_node[2] for a_node in nodes]
        @test minimum(nodes_x) >= x_bounds[1]
        @test maximum(nodes_x) <= x_bounds[2]
        @test minimum(nodes_y) >= y_bounds[1]
        @test maximum(nodes_y) <= y_bounds[2]
    end

    # in_circumcircle
    @testset "in_circumcircle" begin
        nodes = [[0., 0.], [1., 0.], [0., 1.],   # tri nodes
                 [0.1, 0.1],               # node 4: in
                 [0.5, 2.]]                # node 5: out
        # tests are done with the triangle is counterclockwise and clockwise order
        @test in_circumcircle([1, 2, 3], 4, nodes) == true
        @test in_circumcircle([1, 3, 2], 4, nodes) == true
        @test in_circumcircle([1, 2, 3], 5, nodes) == false
        @test in_circumcircle([1, 3, 2], 5, nodes) == false
    end

    # is_clockwise(A, B, C)
    @testset "is_clockwise" begin
        A = [0., 0.]
        B = [1., 0.]
        C = [0., 1.]
        @test is_clockwise(A, B, C) == false
        @test is_clockwise(A, C, B) == true
    end

    # create_supertriangle
    @testset "create_supertriangle" begin
        n = 100
        x_bounds = [-2, 3]
        y_bounds = [-4, 5]
        nodes = create_nodes(n, x_bounds, y_bounds)
        triangle = create_supertriangle(nodes)
        @test length(triangle) == 3
        @test all([length(triangle_pt)==2 for triangle_pt in triangle])
        append!(nodes, triangle)
        triangle_ids = [n+1, n+2, n+3]
        @test all([inside_triangle(triangle_ids, node_id, nodes) for node_id=1:n])
    end

    # require_flip
    @testset "require_flip" begin
        AB = Edge([1, 2], [])
        BC = Edge([2, 3], [])
        CD = Edge([3, 4], [])
        BD = Edge([2, 4], [])
        AD = Edge([1, 4], [])
        edges = [AB, BC, CD, BD, AD]
        triangles = [Triangle([1, 4, 5]), Triangle([2, 3, 4])]
        # first try: both triangles have 2pi/3 opposite angles => flip required
        width = 1.
        height = width/tan(pi/3)
        nodes = [[0., -height], [width, 0.],  [0., height], [-width, 0.]]
        @test require_flip(1, 2, triangles, edges, nodes) == true
        # second try: both triangles have pi/3 opposite angles => flip not required
        height = width/tan(pi/6)
        nodes = [[0., -height], [width, 0.],  [0., height], [-width, 0.]]
        @test require_flip(1, 2, triangles, edges, nodes) == false
    end

    # compute_angle
    @testset "compute_angle" begin
        nodes = [[0., 0.], [1., 0.], [0., tan(pi/3)]]
        edge1 = Edge([1, 2], [])
        edge2 = Edge([1, 3], [])
        edge3 = Edge([2, 3], [])
        edge3_flip = Edge([3, 2], [])
        @test isapprox(compute_angle(edge1, edge2, nodes), pi/2)
        @test isapprox(compute_angle(edge1, edge3, nodes), 2*pi/3)
        @test isapprox(compute_angle(edge1, edge3_flip, nodes), pi/3)
        @test isapprox(compute_angle(edge2, edge3, nodes), pi/6)
    end

    # get_tri_node_ids
    @testset "get_tri_node_ids" begin
        edges = [Edge([1, 2], []), Edge([1, 3], []), Edge([2, 3], [])]
        triangles = [Triangle([1, 2, 3])]
        @test sort(get_tri_node_ids(1, triangles, edges)) == [1, 2, 3]
    end

    # inside_triangle
    @testset "inside_triangle" begin
        nodes = [[0., 0.], [1., 0.], [0., 1.],  # tri nodes
                    [0.1, 0.1],                 # node 4: in
                    [0., 0.5],                  # node 5: on an edge
                    [0., 1e-6],                 # node 6: slightly in
                    [0., -1e-6],                # node 7: slightly out
                ]
        # tests are done with the triangle is counterclockwise and clockwise order
        @test inside_triangle([1, 2, 3], 4, nodes) == true
        @test inside_triangle([1, 3, 2], 4, nodes) == true
        @test inside_triangle([1, 2, 3], 5, nodes) == false
        @test inside_triangle([1, 3, 2], 6, nodes) == true
        @test inside_triangle([1, 3, 2], 7, nodes) == false
    end

    # get_triangulation_nodes
    @testset "triangulation_nodes" begin
        edges = [
            Edge([1, 2], []),
            Edge([5, 1], []),
            Edge([3, 2], []),
            Edge([6, 3], []),
            Edge([7, 8], [])
        ]
        @test sort(get_triangulation_nodes(edges)) == [1, 2, 3, 5, 6, 7, 8]
    end
end


# Unit tests for functions related to the CLI application
function test_application()
    
    @testset "parse_args" begin
        args_string = "-f -g 5 -ij myFile --output-file=myOutput --use-dll 10"
        args = parse_args(args_string)
        @test "f" in keys(args)
        @test args["f"] === nothing
        @test "g" in keys(args)
        @test args["g"] == "5"
        @test "i" in keys(args)
        @test args["i"] === nothing
        @test "j" in keys(args)
        @test args["j"] === nothing
        @test "output-file" in keys(args)
        @test args["output-file"] == "myOutput"
        @test "use-dll" in keys(args)
        @test args["use-dll"] === nothing
        @test args["non_options"] == ["myFile", "10"]
    end

    @testset "print_help" begin
        original_stdout = stdout
        (read_end, ) = redirect_stdout()
        print_help()
        output = readline(read_end)
        @test output != ""
        redirect_stdout(original_stdout)
    end

    @testset "get_arg_value" begin
        args = Dict(
            "a" => nothing,
            "a_long" => nothing,
            "b" => nothing,
            "c_long" => nothing,
            "d" => "5",
            "e_long" => "10",
            "f" => "1",
            "f_long" => "2"
        )
        @test get_arg_value(["z", "z_long"], args) === false
        @test get_arg_value([], args) === false
        @test get_arg_value(["z", "a"], args) === nothing
        @test get_arg_value(["a", "a_long"], args) === nothing
        @test get_arg_value(["b", "b_long"], args) === nothing
        @test get_arg_value(["c", "c_long"], args) === nothing
        @test get_arg_value(["d", "d_long"], args) == "5"
        @test get_arg_value(["e", "e_long"], args) == "10"
        @test_throws ErrorException get_arg_value(["f", "f_long"], args)
    end

    @testset "get_main_parameters" begin
        # default values
        params = get_main_parameters("")
        @test params["create_nodes"] === false
        @test params["create_nodes_number"] == 100
        @test params["load_node_file"] === false
        @test params["load_result_files"] === false
        @test params["output"] == "./triangulation.json"
        @test params["plot_results"] === false
        @test params["save_steps"] === false
        @test params["help"] === true
        
        # set 1: short options
        args_string = "-c -o myOutput -h -p -s"
        params = get_main_parameters(args_string)
        @test params["create_nodes"] === true
        @test params["create_nodes_number"] == 100
        @test params["load_node_file"] === false
        @test params["load_result_files"] === false
        @test params["output"] == "myOutput"
        @test params["plot_results"] == "final"
        @test params["save_steps"] === true
        @test params["help"] === true

        # set 2: long options
        args_string = "--create-nodes --output=myOutput --help --plot --save-steps"
        params = get_main_parameters(args_string)
        @test params["create_nodes"] === true
        @test params["create_nodes_number"] == 100
        @test params["load_node_file"] === false
        @test params["load_result_files"] === false
        @test params["output"] == "myOutput"
        @test params["plot_results"] == "final"
        @test params["save_steps"] === true
        @test params["help"] === true        

        # -f and --from-file: options depending on the existence of the file
        args_string = "-f"
        @test_throws ErrorException get_main_parameters(args_string)

        args_string = "-f ./wrong_file__"
        @test_throws ErrorException get_main_parameters(args_string)
        
        (tmp_filename, tmp_io) = Filesystem.mktemp(".", cleanup=false)
        args_string = "-f $(tmp_filename)"
        params = get_main_parameters(args_string)
        @test params["load_node_file"] == tmp_filename
        close(tmp_io)
        Filesystem.rm(tmp_filename)

        # -l and --load-result: options depending on the existence of files matching the pattern
        args_string = "-l"
        @test_throws ErrorException get_main_parameters(args_string)

        args_string = "-l ./wrong_file__"
        @test_throws ErrorException get_main_parameters(args_string)
        
        tmp_filename_root = "./result_"*basename(Filesystem.tempname(cleanup=false))
        filenames = ["$(tmp_filename_root)_$(i)" for i in 1:3]
        Filesystem.touch.(filenames)
        file_pattern = "$(tmp_filename_root)_*"
        args_string = "-l $(file_pattern)"
        params = get_main_parameters(args_string)
        @test params["load_result_files"] == file_pattern
        Filesystem.rm.(filenames)

    end;

end


function test_all()
    test_delaunay()
    test_application()
end
;


# execution of the file from a terminal: run test_all function
if abspath(PROGRAM_FILE) == @__FILE__
    println("Running unit tests ...")
    test_all()
end
;