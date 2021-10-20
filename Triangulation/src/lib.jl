import Base.rand
import Base.Filesystem
import LinearAlgebra.det
import LinearAlgebra.dot
import LinearAlgebra.norm
import Plots
import JSON
using Logging
using LoggingExtras
import Glob.glob


#------------------------------------------------------------------------------
# Datastructures and functions related to delaunay triangulation
#------------------------------------------------------------------------------

struct Edge
    # Note: possible improvement using static arrays to nodes and in_triangles
    #       in_triangles::StaticArrays.SVector{2,Int}

    # indices of the nodes
    nodes::Array{Int, 1}
    # indices of the triangles containing this edge
    in_triangles::Array{Int, 1}
end

struct Triangle
    # Note: possible improvement using static arrays to nodes and in_triangles
    #       edges::StaticArrays.SVector{3,Int}

    # indices of the edges of the triangle
    edges::Array{Int, 1}
end

function create_nodes(number::Int, x_bounds, y_bounds)
    # number:   number of nodes created
    # x_bounds: X bounds of nodes - list of 2 reals
    # y_bounds: Y bounds of nodes - list of 2 reals
    nodes = [[rand(Float64), rand(Float64)] for i in 1:number]
    # apply X and Y bounds
    for i=1:number
        nodes[i][1] = nodes[i][1]*(x_bounds[2] - x_bounds[1]) + x_bounds[1]
        nodes[i][2] = nodes[i][2]*(y_bounds[2] - y_bounds[1]) + y_bounds[1]
    end
    return nodes
end


function in_circumcircle(triangle_ids, id, nodes)
    # Tell if node with if "id" is in the circumcircle of triangle
    # made of the nodes of "triangle_ids"
    # - triangle_ids : id of the nodes of the triangle - list of 3 int
    # - id : id of the node investigated - int
    # - nodes : coordinates of the nodes - list of [x, y] arrays of floats
    if (length(triangle_ids) != 3)
        error("Argument triangle_ids must be a list of 3 int")
    end
    # vertices of the triangle
    A = nodes[triangle_ids[1]]
    B = nodes[triangle_ids[2]]
    C = nodes[triangle_ids[3]]
    P = nodes[id]
    # switch vertices if the triangle is clockwise
    if (is_clockwise(A, B, C))
        B, C = C, B
    end
    # create the matrix (see https://fr.wikipedia.org/wiki/Triangulation_de_Delaunay)
    mat = zeros(Float64, 3, 3)
    mat[1, 1] = A[1] - P[1]
    mat[2, 1] = B[1] - P[1]
    mat[3, 1] = C[1] - P[1]
    mat[1, 2] = A[2] - P[2]
    mat[2, 2] = B[2] - P[2]
    mat[3, 2] = C[2] - P[2]
    mat[1, 3] = A[1]^2 - P[1]^2 + A[2]^2 - P[2]^2
    mat[2, 3] = B[1]^2 - P[1]^2 + B[2]^2 - P[2]^2
    mat[3, 3] = C[1]^2 - P[1]^2 + C[2]^2 - P[2]^2
    return (det(mat) > 0)
end

function is_clockwise(A, B, C)
    # A, B & C: list of 2 reals representing 2D points
    AB = B - A
    AC = C - A
    determinant = AB[1]*AC[2] - AB[2]*AC[1]
    return (determinant < 0)
end

function create_supertriangle(nodes)
    # Create a triangle containing all the nodes
    # - nodes : coordinates of the nodes - list of [x, y] arrays of floats

    # bounds of the nodes
    nodes_x = [a_node[1] for a_node in nodes]
    nodes_y = [a_node[2] for a_node in nodes]
    x_bounds = [minimum(nodes_x), maximum(nodes_x)]
    y_bounds = [minimum(nodes_y), maximum(nodes_y)]

    # take a margin (relative to the X & Y bounds)
    m = 0.05
    min_length = 1e-4
    lx = x_bounds[2] - x_bounds[1] + min_length
    x_bounds = [x_bounds[1]-m*lx, x_bounds[2]+m*lx]
    ly = y_bounds[2] - y_bounds[1] + min_length
    y_bounds = [y_bounds[1]-m*ly, y_bounds[2]+m*ly]

    # get the circle containing the nodes
    right_top = [x_bounds[2], y_bounds[2]]
    center = [0.5*(x_bounds[1]+x_bounds[2]), 0.5*(y_bounds[1]+y_bounds[2])]
    radius = sqrt( (right_top[1]-center[1])^2 + (right_top[2]-center[2])^2 )

    # build an equilateral triangle with edges adjacent to the circle
    # base's center = lowest point of the circle on the Y axis
    base_mid = [center[1], center[2]-radius]
    # radius projection on an edge other than the base
    edge_mid = [center[1]+radius*cos(pi/6), center[2]+radius*sin(pi/6)]
    edge_l = 4*abs(edge_mid[1] - base_mid[1])
    tri_lower_right = [base_mid[1]+edge_l/2, base_mid[2]]
    tri_lower_left = [base_mid[1]-edge_l/2, base_mid[2]]
    tri_top = [base_mid[1], base_mid[2]+edge_l*sin(pi/3)]

    return [tri_lower_right, tri_top, tri_lower_left]
end

function plot_start(nodes, tri)
    # plot the supertriangle and the nodes
    # - nodes : coordinates of the nodes - Array(n, 2) of floats
    # - tri : coordinates of the vertices of the triangle - array containing 3 doubles
    Plots.scatter(nodes[:,1], nodes[:,2], label="")
    Plots.plot!([tri[1][1], tri[2][1]], [tri[1][2], tri[2][2]], label="", color=:black)
    Plots.plot!([tri[2][1], tri[3][1]], [tri[2][2], tri[3][2]], label="", color=:black)
    Plots.plot!([tri[3][1], tri[1][1]], [tri[3][2], tri[1][2]], label="", color=:black)
    #Plots.plot!(supertriangle[2], supertriangle[3])
    #Plots.plot!(supertriangle[3], supertriangle[1])
end

function require_flip(tri1_id, tri2_id, triangles, edges, nodes)
    # Evaluate if the triangles must be flipped, that is to say if the opposite angles don't
    # match the delaunay criterion (sum <= 180°)
    # - tri1_id : id of the first triangle
    # - tri2_id : id of the second triangle
    # - triangles : array(n) containing all the triangles
    # - edges : array(n) containing all the edges
    # - nodes : coordinates of the all nodes - list of [x, y] arrays of floats

    triangles_pair = [triangles[tri1_id], triangles[tri2_id]]

    # common_edge
    common_edge = 0
    for triangle1_edge in triangles_pair[1].edges
        if triangle1_edge in triangles_pair[2].edges
            common_edge = triangle1_edge
            break
        end
    end
    if common_edge == 0
        error("Common edge not found")
    end

    # angle opposite to the common edge
    angles = [NaN, NaN]
    for i in [1, 2]
        angle_edges = [edges[e] for e in triangles_pair[i].edges if e != common_edge]
        angle_node = intersect(angle_edges[1].nodes, angle_edges[2].nodes)[1]

        # ensure that the first point of both edges is their common point
        if (angle_edges[1].nodes[1] != angle_node)
            angle_edges[1] = Edge(reverse(angle_edges[1].nodes), [])
        end
        if (angle_edges[2].nodes[1] != angle_node)
            angle_edges[2] = Edge(reverse(angle_edges[2].nodes), [])
        end

        # compute angle
        angles[i] = compute_angle(angle_edges[1], angle_edges[2], nodes)
    end
    return (abs(angles[1]) + abs(angles[2]) > pi)
end


function compute_angle(edge1, edge2, nodes)
    # - edge1 : first edge (type Edge)
    # - edge2 : second edge (type Edge)
    # - nodes : coordinates of the all nodes - list of [x, y] arrays of floats
    p1 = nodes[edge1.nodes[1]]
    p2 = nodes[edge1.nodes[2]]
    v1 = [p2[1]-p1[1], p2[2]-p1[2]]
    v1 = v1/norm(v1)
    
    p1 = nodes[edge2.nodes[1]]
    p2 = nodes[edge2.nodes[2]]
    v2 = [p2[1]-p1[1], p2[2]-p1[2]]
    v2 = v2/norm(v2)

    return acos(dot(v1, v2))
end

function get_tri_node_ids(triangle_id, triangles, edges)
    # - triangle_id : id of the triangle
    # - triangles : array(n) containing all the triangles
    # - edges : array(n) containing all the edges
    triangle_edges = [edges[edge_id] for edge_id in triangles[triangle_id].edges]
    all_nodes = [e.nodes for e in triangle_edges]
    return unique(Iterators.flatten(all_nodes))
end


function get_triangulation_nodes(edges)
    # - edges : array(n) containing all the edges
    node_ids = unique(Iterators.flatten([e.nodes for e in values(edges)]))
    return node_ids
end


function flip_triangulation!(triangles, edges, nodes, do_save_steps, step_file_pattern; init_flip_step=0)
    # Flip the triangles of the triangulation that don't passe delaunay's criterion
    # This function modifies the triangles and edges input arguments
    # - triangles : array(n) containing all the triangles
    # - edges : array(n) containing all the edges
    # - nodes : coordinates of the all nodes - list of [x, y] arrays of floats

    flip_step = init_flip_step
    nb_edges_flipped = 0
    for edge_id in keys(edges)
        common_edge = edges[edge_id]
        common_nodes = common_edge.nodes
        # find if neighbor triangles separated by this edge must be flipped
        neighbors_id = common_edge.in_triangles
        if 0 in neighbors_id
            # edge is on the border of the triangulation: no flip
            continue
        elseif require_flip(neighbors_id[1], neighbors_id[2], triangles, edges, nodes)
            flip_step = flip_step + 1
            nb_edges_flipped = nb_edges_flipped + 1
            # do flip
            # 1. analyse the neighbor triangles: opposite nodes and id of other edges
            tri1_nodes = get_tri_node_ids(neighbors_id[1], triangles, edges)
            tri2_nodes = get_tri_node_ids(neighbors_id[2], triangles, edges)
            @debug "triangle 1 : $(neighbor_tri_ids[1]) has nodes $(tri1_nodes)"
            @debug "triangle 2 : $(neighbor_tri_ids[2]) has nodes $(tri2_nodes)"
            opposite_nodes = [setdiff(tri1_nodes, common_nodes)[1], setdiff(tri2_nodes, common_nodes)[1]]
            #top_nodes = setdiff(union(triangles[neighbors_id[1]].nodes, triangles[neighbors_id[1]].nodes)
            #       common_edge.nodes)
            @debug "opposite_nodes : $(opposite_nodes)"
            if length(unique(opposite_nodes)) != 2
                @error "Flip: same opposite nodes"
                error("Flip: same opposite nodes")
            end
            other_edges_ids = setdiff(triangles[neighbors_id[1]].edges, edge_id)
            append!(other_edges_ids, setdiff(triangles[neighbors_id[2]].edges, edge_id))
            # TODO: rewrite with union
            # other_edges_ids = setdiff(union(triangles[neighbors_id[1]].edges, triangles[neighbors_id[2]].edges), common_edge)

            # 2. reconstruct the common edge and the neighbor triangles
            #    -> the common edge is now built with the opposite nodes
            #    -> the nodes of the former common edge are now the new opposite nodes
            
            # 2.1 common edge
            common_edge = Edge(opposite_nodes, neighbors_id)
            @debug "Common edge $(common_edge_id) updated: nodes = $(opposite_nodes)"
            edges[edge_id] = common_edge

            # 2.2 triangles        
            for i_tri in 1:2
                @debug "Rework triangle $(i_tri)"
                tri_id = neighbors_id[i_tri]
                # use a (former) common node and the new common edge
                top_node = common_nodes[i_tri]
                @debug "- Top node = $(top_node)"
                side_edges_ids = [e_id for e_id in other_edges_ids if top_node in edges[e_id].nodes]
                @debug "- Side edge 1 = $(side_edges_ids[1]): $(edges[side_edges_ids[1]])"
                @debug "- Side edge 2 = $(side_edges_ids[2]): $(edges[side_edges_ids[2]])"
                # - create new triangle (replace existing one with same id)
                triangles[tri_id] = Triangle([edge_id, side_edges_ids[1], side_edges_ids[2]])
                @debug "- Triangle $(tri_id) updated: $(triangles[tri_id].edges)"
                # - update the neighbor information (in_triangle) for the side edges
                for i_edge in side_edges_ids
                    # before the flip: side edge is in 
                    # - one of the triangles that will be flipped: an id in in neighbor_tri_ids
                    # - a triangle out of the parent triangle unchanged : outer_tri
                    #   -> outer_tri_id is 0 is this outer triangle does not exist
                    # after the flip: side edge is in
                    # - the triangle that has been flipped: id = tri_id
                    # - the outer triangle with id outer_tri_id
                    # => replace id of triangle-to-be-flipped by id of new triangle
                    in_triangles = copy(edges[i_edge].in_triangles)
                    if (in_triangles[1] in neighbors_id)
                        in_triangles[1] = tri_id
                    elseif (in_triangles[2] in neighbors_id)
                        in_triangles[2] = tri_id
                    else
                        @error "- Edge $(i_edge) is not in flipped triangle $(tri_id)"
                        error("- Edge ", i_edge, " is not in flipped triangle ", tri_id)
                    end
                    edge_nodes = edges[i_edge].nodes
                    #edges[i_edge] = Edge(edge_nodes, [tri_id, outer_tri_id])
                    edges[i_edge] = Edge(edge_nodes, in_triangles)
                    @debug "- Edge $(i_edge) now in triangles $(in_triangles)"
                end
            end

            # save result
            if do_save_steps
                new_edges_ids = union(other_edges_ids, edge_id)
                save_step(step_file_pattern, "flip_$(flip_step)", nodes, 0, edges, 
                new_edges_ids, [], [], triangles)
            end
        end
    end
    if nb_edges_flipped > 0
        @info "$(nb_edges_flipped) on $(length(edges)) edges have been flipped"
    end

    # check if there remain edges to flip after this first step
    nb_edges_to_flip = 0
    for edge_id in keys(edges)
        common_edge = edges[edge_id]
        # find if neighbor triangles separated by this edge must be flipped
        neighbors_id = common_edge.in_triangles
        if 0 in neighbors_id
            # edge is on the border of the triangulation: no flip
            continue
        elseif require_flip(neighbors_id[1], neighbors_id[2], triangles, edges, nodes)
            nb_edges_to_flip = nb_edges_to_flip + 1
        end
    end
    if (nb_edges_to_flip > 0)
        @info "Recursive call to flip_triangulation to fix $(nb_edges_to_flip) edges"
        flip_triangulation!(triangles, edges, nodes, do_save_steps, step_file_pattern, init_flip_step=flip_step)
    end
    
end


function save_step(file_pattern, file_suffix, nodes, last_node, edges, last_edges, 
                    supertriangle_node_ids, supertriangle_edge_ids, triangles)
    # Create the name of the step from the pattern by inserting a suffix
    # between the pattern and the file extension
    (base, ext) = splitext(file_pattern)
    file = base*file_suffix*ext
    output = Dict("nodes" => nodes,
                  "last_node" => last_node,
                  "edges" => edges,
                  "last_edges" => last_edges, 
                  "supertriangle_node_ids" => supertriangle_node_ids,
                  "supertriangle_edge_ids" => supertriangle_edge_ids,
                  "triangles" => triangles,
                  )
    open(file, "w") do io
        write(io, JSON.json(output))
    end
end

function plot_step(step_file, show_non_delaunay_triangles)
    # plot the nodes, edges and triangles for a step
    # warning: an existing file with the same name will be overwritten
    io = open(step_file, "r")
    string_content = read(io, String)
    close(io)
    json_content = JSON.parse(string_content)
    nodes = json_content["nodes"]
    last_node = json_content["last_node"]
    edges = Dict{Int, Edge}()
    for item in collect(json_content["edges"])
        key = parse(Int, item[1])
        edges[key] = Edge(item[2]["nodes"], item[2]["in_triangles"])
    end
    last_edges = json_content["last_edges"]
    supertriangle_edge_ids = json_content["supertriangle_edge_ids"]
    supertriangle_node_ids = json_content["supertriangle_node_ids"]
    triangles = Dict{Int, Triangle}()
    for item in collect(json_content["triangles"])
        key = parse(Int, item[1])
        triangles[key] = Triangle(item[2]["edges"])
    end

    # create a node to init figure - it will be reploted after
    Plots.scatter([nodes[1][1]], [nodes[1][2]], aspect_ratio=:equal, label="", color=:red)

    # non-delaunay triangles
    if show_non_delaunay_triangles
        for (tri_id, tri) in triangles
            tri_node_ids = get_tri_node_ids(tri_id, triangles, edges)
            # find if any node used in the triangulation is in this triangle's circumcircle
            for lookup_node_id in get_triangulation_nodes(edges)
                if !(lookup_node_id in tri_node_ids)
                    if in_circumcircle(tri_node_ids, lookup_node_id, nodes)
                        # the triangle is non-delaunay -> plot it
                        tri_points = [Tuple(nodes[i]) for i in tri_node_ids]
                        #append!(shapes, Plots.Shape(tri_points))
                        Plots.plot!([Plots.Shape(tri_points)], color=:red, label="")
                        # end nodes lookup
                        break
                    end
                end
            end
        end
    end

    # edges
    for e in edges
        if e[1] in last_edges
            # this edge will be plotted with the last edges - skip it
            continue
        elseif e[1] in supertriangle_edge_ids
            # this edge will be plotted with the supertriangle edges - skip it
            continue
        end
        edge_x = [nodes[n][1] for n in e[2].nodes]
        edge_y = [nodes[n][2] for n in e[2].nodes]
        Plots.plot!(edge_x, edge_y, label="", color=:black)
    end

    # supertriangle's edges (if any)
    for edge_id in supertriangle_edge_ids
        if edge_id in last_edges
            # this edge will be plotted with the last edges - skip it
            continue
        end
        edge_nodes = edges[edge_id].nodes
        edge_x = [nodes[n][1] for n in edge_nodes]
        edge_y = [nodes[n][2] for n in edge_nodes]
        Plots.plot!(edge_x, edge_y, label="", color=:chocolate3)
    end   

    # last edges (if any)
    for edge_id in last_edges
        edge_nodes = edges[edge_id].nodes
        edge_x = [nodes[n][1] for n in edge_nodes]
        edge_y = [nodes[n][2] for n in edge_nodes]
        Plots.plot!(edge_x, edge_y, label="", color=:green4)
    end

    # nodes
    nodes_x = [n[1] for n in nodes]
    nodes_y = [n[2] for n in nodes]
    Plots.scatter!(nodes_x, nodes_y, label="", color=:deepskyblue)

    # supertriangle's nodes
    nodes_x = [nodes[i][1] for i in supertriangle_node_ids]
    nodes_y = [nodes[i][2] for i in supertriangle_node_ids]
    Plots.scatter!(nodes_x, nodes_y, label="", color=:burlywood1)

    # last node inserted
    if last_node > 0
        Plots.scatter!([nodes[last_node][1]], [nodes[last_node][2]], label="", color=:red)
    end

    # title
    figure_title = basename(splitext(step_file)[1])
    Plots.title!(figure_title)

    # save fig
    figure_file = splitext(step_file)[1]*".png"
    #rm(figure_file, force=true)
    Plots.savefig(figure_file)
end


# TODO: remove -> use instead loop on globbed files + plot_step
function plot_several_steps(file_pattern, show_non_delaunay_triangles)
    # plot the step files matching the input pattern (eg "./steps/step_.txt")
    file_pattern_basename = splitext(basename(file_pattern))[1]
    file_pattern_end = splitext(basename(file_pattern))[2]
    for step_file in readdir(dirname(file_pattern), join=true)
        if startswith(basename(step_file), file_pattern_basename) && endswith(step_file, file_pattern_end)
            println("Plot step file: ", step_file)
            plot_step(step_file, show_non_delaunay_triangles)
        end
    end
end


function inside_triangle(triangle_ids, node_id, nodes)
    # Evaluate if a node is inside a triangle using barycentric coordinates method
    # https://en.wikipedia.org/wiki/Barycentric_coordinate_system
    #   See #Edge_approach for the calculation of the barycentric coordinates
    #   and #Determining_location_with_respect_to_a_triangle for the inside/outside criterion
    # - triangle_ids : id of the nodes of the triangle - list of 3 int
    # - node_id : id of the node investigated - int
    # - nodes : coordinates of the nodes - list of nodes coordinates [x, y]
    (x, y) = nodes[node_id]
    (x1, y1) = nodes[triangle_ids[1]]
    (x2, y2) = nodes[triangle_ids[2]]
    (x3, y3) = nodes[triangle_ids[3]]
    denominator = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3)
    lambda1 = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3))/denominator
    lambda2 = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3))/denominator
    lambda3 = 1 - lambda1 - lambda2
    inside = (0 < lambda1 < 1) && (0 < lambda2 < 1) && (0 < lambda3 < 1)
    return inside
end


function triangulate(input_nodes, must_save_steps=false, step_file_pattern="./steps/step_.txt")
    # - input_nodes : coordinates of the all nodes - list of [x, y] arrays of floats
    # - must_save_steps : save each step to a file (named from pattern)
    # - step_file_pattern : pattern of the step files - ex: "./steps/step_.txt"

    # nodes: work on a shallow copy of input nodes (we will add then remove the 3 supernodes)
    nodes = copy(input_nodes)
    number_of_nodes = length(nodes)
    if number_of_nodes < 3
        # TODO: return empty set instead of printing an error
        @error "The number of nodes must be 3 or more"
        error("The number of nodes must be 3 or more")
    end

    # remove existing step files
    if (must_save_steps)
        @info "Removing existing step files"
        step_file_pattern_basename = splitext(basename(step_file_pattern))[1]
        for step_file in readdir(dirname(step_file_pattern), join=true)
            if startswith(basename(step_file), step_file_pattern_basename)
                rm(step_file)
                @debug "Removed: "*step_file
            end
        end
    end
    
    # create dicts to store the triangles (1 triangle is described by the indices of its 3 nodes)
    # notice: dicts are used instead of arrays to avoid renumbering elements each time an element
    #         is removed 
    triangles = Dict{Int, Triangle}()
    current_triangle_index = 0
    edges = Dict{Int, Edge}()
    current_edge_index = 0
    # create the super triangle
    # -> its id is 1
    # -> its nodes have ids n+1, n+2 and n+3
    # -> its edges are numbered 1, 2 and 3
    supertriangle_nodes = create_supertriangle(nodes)
    supertriangle_node_ids = [number_of_nodes+1, number_of_nodes+2, number_of_nodes+3]
    append!(nodes, supertriangle_nodes)
    
    edge1 = Edge([supertriangle_node_ids[1], supertriangle_node_ids[2]], [1, 0])
    edges[1] = edge1
    edge2 = Edge([supertriangle_node_ids[2], supertriangle_node_ids[3]], [1, 0])
    edges[2] = edge2
    edge3 = Edge([supertriangle_node_ids[3], supertriangle_node_ids[1]], [1, 0])
    edges[3] = edge3
    current_edge_index = 3

    supertriangle_edge_ids = [1, 2, 3]
    supertriangle = Triangle(supertriangle_edge_ids)
    triangles[1] = supertriangle
    current_triangle_index = 1

    # iterate over the nodes to insert
    for node_id=1:number_of_nodes
        @info "Insert node $(node_id)"
        
        # find the parent triangle : triangle which circumcircle contains this node
        # --------------------------
        parent_tri_id = 0
        for tri_id in keys(triangles)
            # get nodes of the triangle
            tri_nodes = get_tri_node_ids(tri_id, triangles, edges)
            if inside_triangle(tri_nodes, node_id, nodes)
                parent_tri_id = tri_id
                break
            end
            # Warning: what to do if the node is ON an edge of the triangle ? 
        end
        if (parent_tri_id == 0)
            @error "parent triangle not found for node $(node_id)"
            error("parent triangle not found for node ", node_id)
        end
        # DEBUG: remove
        @debug "Id of the parent the triangle: $(parent_tri_id)"

        # update the triangulation: create the new triangles from the parent triangle and the inserted node
        # -------------------------
        # - edges of the parent triangle
        parent_edges = triangles[parent_tri_id].edges
        parent_nodes = get_tri_node_ids(parent_tri_id, triangles, edges)

        # - create the new edges (nodes of the parent tri + new node)
        new_edges = [Edge([parent_nodes[i], node_id], [0, 0]) for i in 1:3]
        new_edges_ids = [0, 0, 0]
        for i_new_edge in 1:3
            if length(unique(new_edges[i_new_edge].nodes)) != 2
                @error "triangulation: error edge with same node twice $(new_edges[i_new_edge])"
                error("triangulation: error edge with same node twice ", new_edges[i_new_edge])
            end
            current_edge_index = current_edge_index + 1
            new_edges_ids[i_new_edge] = current_edge_index
            edges[current_edge_index] = new_edges[i_new_edge]
        end

        # - create the new triangles with a bound edge and the two matching new edges
        new_triangles_indices = []
        for parent_edge_id in parent_edges
            @debug "Create new triangle from parent edge $(parent_edge_id) : $(edges[parent_edge_id])"
            # 
            ref_edge = edges[parent_edge_id]
            ref_nodes = ref_edge.nodes
            # find the edges of the new triangle
            tri_edges_ids = [parent_edge_id]
            for i_edge=1:3
                if (new_edges[i_edge].nodes[1] in ref_nodes) || (new_edges[i_edge].nodes[2] in ref_nodes)
                    append!(tri_edges_ids, new_edges_ids[i_edge])
                end
            end
            if length(tri_edges_ids) != 3
                @error "Newly created triangle has not 3 edges"
                error("Newly created triangle has not 3 edges")
            end

            # create the triangle and add it to the table of triangles
            current_triangle_index = current_triangle_index + 1
            new_triangle_id = current_triangle_index
            append!(new_triangles_indices, new_triangle_id)
            triangles[current_triangle_index] = Triangle(tri_edges_ids)
            @debug "Create triangle #$(new_triangle_id) with edges $(tri_edges_ids)"

            # - update neighbor info of the edge from the parent triangle
            ref_edge_triangles = ref_edge.in_triangles
            replace!(ref_edge_triangles, parent_tri_id => new_triangle_id)
            edges[parent_edge_id] = Edge(ref_nodes, ref_edge_triangles)
        end

        # - update the new edges: neighbor triangles
        for i_new_edge in new_edges_ids
            #tri_neighbors = [t_id for t_id in edges[i_new_edge].in_triangles 
            #                 if t_id != parent_tri_id &&  t_id != 0]
            
            # find the new triangle containing this edge
            tri_neighbors = []
            for tri_id in new_triangles_indices
                if i_new_edge in triangles[tri_id].edges
                    append!(tri_neighbors, tri_id)
                end
            end
            # update the edge
            @debug "Update edge #$(i_new_edge) : in_triangles = $(tri_neighbors)"
            edge_nodes = edges[i_new_edge].nodes
            edges[i_new_edge] = Edge(edge_nodes, tri_neighbors)
        end

        # - remove the parent triangle
        @debug "Delete parent triangle with id $(parent_tri_id)"
        pop!(triangles, parent_tri_id)

        # save step
        if must_save_steps
            save_step(step_file_pattern, string(node_id), nodes, node_id, edges, 
                new_edges_ids, supertriangle_node_ids, supertriangle_edge_ids, triangles)
        end

        # flip the new triangles, if needed
        # ---------------------------------
        # we iterate over the edges that made the former parent triangle
        # -> for each, we need to see if the triangles on each side respect the angle criterion

        flip_step_file_pattern = splitext(step_file_pattern)[1]*string(node_id)*"_"*splitext(step_file_pattern)[2]
        flip_triangulation!(triangles, edges, nodes, must_save_steps, flip_step_file_pattern)
    end

    #--------------------------------------------------------------------------
    # Remove the data related to the (now deleted) supertriangle: nodes, edges
    #--------------------------------------------------------------------------
    
    @info "Cleanup"

    # find the edges containg a supernode
    edges_to_delete_id = []
    triangles_to_delete_id = []
    for (edge_id, edge) in edges
        if intersect(edge.nodes, supertriangle_node_ids) != []
            @debug "Edge $(edge_id) contains a node from the supertriangle -> to delete"
            # mark the edge for deleting
            append!(edges_to_delete_id, edge_id)
            # mark the triangle(s) made with this edge for deletion
            for tri_id in edge.in_triangles
                if tri_id != 0
                    append!(triangles_to_delete_id, tri_id)
                end
            end
        end
    end
    triangles_to_delete_id = unique(triangles_to_delete_id)

    # delete the edges
    for edge_id in edges_to_delete_id
        pop!(edges, edge_id)
    end

    # update information for the edges in triangles about to be deleted
    for (edge_id, edge) in edges
        in_deleted_tri_ids = intersect(edge.in_triangles, triangles_to_delete_id)
        # 3 cases for the current edge:
        # - in 2 triangles marked for deletion -> it has *already* been deleted -> nothing to do
        # - in 1 triangle marked for deletion -> replace its id by 0
        # - in 0 triangles marked for deleting -> nothing to do
        if length(in_deleted_tri_ids) == 1
            in_triangles_updated = replace(edge.in_triangles, in_deleted_tri_ids[1] => 0)
            edges[edge_id] = Edge(edge.nodes, in_triangles_updated)
        end
    end

    # delete the triangles
    for tri_id in triangles_to_delete_id
        @debug "Delete triangle $(tri_id)"
        pop!(triangles, tri_id)
    end

    # delete the nodes of the supertriangle: 3 last elements of "nodes"
    nodes = nodes[1:number_of_nodes]

    # write the final state
    if must_save_steps
        save_step(step_file_pattern, "final", nodes, 0, edges, [], [], [], triangles)
    end

    # TODO: return results

    return edges, triangles
end


#------------------------------------------------------------------------------
# Functions related to the CLI application
# - parse_args: parses the options from a user's input string
# - print_help: displays the help of the CLI
# - get_arg_value: get the value of a couple of options (long/short) from 
#   parse_args's output
# - get_main_parameters: get the CLI parameters from the user's input string
# - main: perform the operation requested by the user
# - test_main: unit tests from the functions of this file
#------------------------------------------------------------------------------


function parse_args(args_string::String)
    # parse the input arguments
    output = Dict()
    output["non_options"] = String[]
    option_regex = r"^-[a-zA-Z]"
    long_option_regex = r"^--[a-zA-Z]"
    # note: elements of the result of split are of type SubString{String}
    #       -> convert to String to ease future handling
    args = [string(s) for s in split(args_string)]
    index = 1
    while index <= length(args)
        currentArg = string(args[index])
        
        if currentArg == "-"
            error("Invalid argument '-'")
        end

        if occursin(option_regex, currentArg)
            # options found -> extract each one are store it in the dict
            number_of_options = length(currentArg) - 1
            if number_of_options == 1
                # single option
                # -> the next argument can hold a value, or be another option
                option = string(currentArg[2])
                if index < length(args)
                    nextArg = args[index+1]
                    if occursin(option_regex, nextArg)
                        # option: skip it
                        value = nothing
                    else
                        # value: store it
                        value = nextArg
                        # increment index because arg[index+1] as been dealt with
                        index = index + 1
                    end
                else
                    # no following argument: empty value
                    value = nothing
                end
                # store option
                output[option] = value
            else
                # several options stacked -> they cannot have values
                for option in currentArg[2:end]
                    output[string(option)] = nothing
                end
            end
        elseif occursin(long_option_regex, currentArg)
            # long option
            if occursin("=", currentArg)
                # the long option is passed with a value
                option, value = split(currentArg[3:end], "=")
                value = string(value)   # convert from SubString{String} to String to ease handling
                if value == ""
                    error("No value for option $(option)")
                end
            else
                option = currentArg[3:end]
                value = nothing
            end
            output[option] = value
        else
            # value not related to an option
            push!(output["non_options"], currentArg)
        end

        # go to the next index in the arguments array
        index = index + 1
    end

    return output
end


function print_help()
    println("Help: Delaunay triangulation of a set of nodes")
    println("-----")
    println("")
    println("Syntax for options: -s short_value --long-option=long_value")
    println("")
    println("Options:")
    println("-c               Create a random set of nodes. The value is the number of nodes")
    println("                 (default: 100")
    println("--create-nodes   Same as -c")
    println("-f               Name of the text file containing the nodes (1 node per line, ")
    println("                 x et y coordinates separated by a comma)")
    println("                 [!] NOT YET IMPLEMENTED")
    println("--from-file      Same as -f")
    println("-h               Display the help")
    println("--help           Same as -h")
    println("-l               When used with -p (or --plot), do not generate a triangulation ")
    println("                 and plot from a result file. The value is the name of the ")
    println("                 result file. Give a regex to plot several files.")
    println("--load-result    Same as -l")
    println("-o               Name of the output file.")
    println("--output         Same as -o")
    println("-p               Plot the triangulation. The possible values are:")
    println("                   'final'   (default) only plot the final triangulation")
    println("                   'all'     plot the final triangulation and each step")
    println("--plot           Same as -p")
    println("-s               Save the steps. The file names are built from the output: ")
    println("                 <output>_step<stepId>_flip<flipId>.<output_ext>")
    println("--save-steps     Same as -s")
    
end


function get_arg_value(options, args_dict)
    # Get the value of an option which can have several names
    # - options: list of names corresponding to the same option.
    #            eg: ["c", "create"] to get short option "-c" and long option "--create"
    # - args_dict: dictionary containing the options and their values (nothing if no value)
    # Return
    #  - false if none of the options is in the arguments dict
    #  - the value (string or nothing) of the option (or of all options is they match)
    #  - error if several options are in the dict with different values
    options_in_args = [n for n in options if n in keys(args_dict)]
    if options_in_args == []
        # None of the options is in the arguments dict -> use the default value
        return false
    end

    if length(options_in_args) == 1
        return args_dict[options_in_args[1]]
    else
        # Redundant options: check that they have the same values
        if length(unique([args_dict[o] for o in options_in_args])) == 1
            return args_dict[options_in_args[1]]
        else
            # error
            error("Options passed with different value: $(options_in_args)")
        end
    end
end


function get_main_parameters(args_string::String)

    # initialize the parameters' dict with default values
    params = Dict(
        "create_nodes" => false,
        "create_nodes_number" => 100,
        "load_node_file" => false,
        "load_result_files" => false,
        "output" =>  "./triangulation.json",
        "plot_results" => false,
        "save_steps" => false,
        "help" => false
    )

    # parse arguments
    args = parse_args(args_string)

    # handle invalid options: signal, set help and return
    all_options = ["c", "create-nodes", "f", "from-file", "h", "help", "l", "load-result",
                     "o", "output", "p", "plot", "s", "save-steps", "non_options"]
    has_invalid_option = false
    for key in keys(args)
        if ! (key in all_options)
            pop!(args, key)
            if length(key) == 1
                option_name = "-$(key)"
            else
                option_name = "--$(key)"
            end
            println("Invalid option $(option_name)")
            has_invalid_option = true
        end
    end
    if has_invalid_option === true
        params["help"] = true
        return params
    end

    # if there is no option, set help and return
    if length(keys(args)) == 1 && "non_options" in keys(args)
        params["help"] = true
        return params
    end


    # get the parameters from input arguments & default value

    # * create nodes: -c and --create_nodes
    arg_create_nodes = get_arg_value(["c", "create-nodes"], args)
    if arg_create_nodes != false
        # either "c" or "create-nodes" is passed
        params["create_nodes"] = true
        if arg_create_nodes !== nothing
            # passed with a value: convert it to int
            params["create_nodes_number"] = parse(Int, arg_create_nodes)
        end
    end

    # * get nodes from a file: -f and --from-file
    arg_nodes_file = get_arg_value(["f", "from-file"], args)
    if arg_nodes_file === nothing
        error("Options -f and --from-file must be used with a value")
    elseif arg_nodes_file != false
        if Filesystem.isfile(arg_nodes_file)
            params["load_node_file"] = arg_nodes_file
        else
            error("Node file passed to option -f or --from-file does not exist")
        end
    end

    # * display help
    arg_help = get_arg_value(["h", "help"], args)
    if arg_help != false
        params["help"] = true
    end

    # * get the results to plot from a file: -l and --load-result
    arg_load_result = get_arg_value(["l", "load-result"], args)
    if arg_load_result === nothing
        error("Options -l and --load-result must be used with a value")
    elseif arg_load_result != false
        # check the existence of files matching the pattern
        if length(glob(arg_load_result)) > 0
            params["load_result_files"] = arg_load_result
        else
            error("No file match the pattern given for option -l or --load-result")
        end
    end

    # * output file: -o and --output
    arg_output = get_arg_value(["o", "output"], args)
    if typeof(arg_output) === String
        # use the value given
        params["output"] = arg_output
    end

    # * plot results: -p and --plot
    arg_plot_results = get_arg_value(["p", "plot"], args)
    if arg_plot_results === nothing
        params["plot_results"] = "final"
    elseif typeof(arg_plot_results) == String
        if arg_plot_results in ["final", "all"]
            # use the value given
            params["plot_results"] = arg_plot_results
        else
            error("Valid values for -p and --plot are \"final\" and \"all\"")
        end
    end

    # * save the triangulation steps: -s and --save
    arg_save_steps = get_arg_value(["s", "save-steps"], args)
    if arg_save_steps !== false
        params["save_steps"] = true
    end

    return params
end


function main(args)
    # get parameters from the user's arguments
    args_string = join(args, " ")
    params = Dict()     # need to create the variable at this level, before the try/catch
    try
        params = get_main_parameters(args_string)
    catch ErrorException
        # something is invalid with the input parameters
        # -> print help and exit
        print_help()
        return
    end

    # apply a scenario, depending on the inputs
    if params["help"] === true
        # scenario 1: display help
        print_help()
    elseif params["load_result_files"] !== false
        # scenario 2: plot the results of a previous triangulation from result files
        step_file_pattern = params["load_result_files"]
        for filename in glob(step_file_pattern)
            plot_step(filename, true)
        end
    else
        # scenario 3: calculate a triangulation
        # - nodes
        nodes = []
        if params["load_node_file"] !== false
            # load the nodes from a file
            error("Load the nodes from a file: not implemented")
        elseif params["create_nodes"] === true
            # create random nodes
            nodes = create_nodes(params["create_nodes_number"], [-1., 1.], [-1., 1.])
        else
            # no policy for the nodes: notify error and exit
            println("The nodes must be either created or loaded from a file")
            println("  - created: -c or --create-nodes | optional value: number of nodes")
            println("  - loaded: -f or --from-file | value: name of the file")
            return
        end
        # - triangulation
        output = params["output"]
        if params["save_steps"] === true
            do_save_steps = true
            # create the pattern for the step files
            step_file_pattern = splitext(output)[1]*"_step_"*splitext(output)[2]
        else
            do_save_steps = false
            step_file_pattern = ""
        end
        (edges, triangles) = triangulate(nodes, do_save_steps, step_file_pattern)
        save_step(output, "", nodes, 0, edges, [], [], [], triangles)
        # - plot
        if params["plot_results"] !== false
            if params["plot_results"] == "final"
                plot_step(params["output"], true)
            elseif params["plot_results"] == "all"
                output_glob_pattern = splitext(output)[1]*"*"*splitext(output)[2]
                for filename in glob(output_glob_pattern)
                    plot_step(filename, true)
                end
            end
        end
    end
end
;



function create_logger(log_file::String, log_min_level::Base.CoreLogging.LogLevel)
    # create a logger thats prints formatted output to console and file
    # see below for the combination of FormatLogger and other logged_test_triangulate
    # https://discourse.julialang.org/t/combination-of-formatlogger-and-filelogger-in-loggingextras/60654/2
    # - log_file: path of the log file
    function fmt(io, args)
        println(io, args._module, " | ", "[", args.level, "] ", args.message)
    end

    logger = TeeLogger(
        MinLevelLogger(FormatLogger(fmt, open(log_file, "w")), log_min_level),
        MinLevelLogger(FormatLogger(fmt, stdout), log_min_level)   # formatted but no color in REPL
        #global_logger()     # unformatted by color in REPL
    )
    return logger
end


# TODO: delete this function but keep the log redirection part for the main fct
function logged_test_triangulate(number_of_nodes, x_bounds, y_bounds, save_steps)
    # manual solution if the function outputs logs with prints
    #=
    open("test_triangulate.log", "w") do io
        redirect_stdout(io) do
            test_triangulate(number_of_nodes, x_bounds, y_bounds, save_steps)
        end
    end
    =#

    # clean solution
    log_file = "test_triangulate.log"
    logger = create_logger(log_file, Logging.Info)

    formatConsoleLogger = FormatLogger() do io, args
        println(io, args._module, " | ", "[", args.level, "] ", args.message)
    end

    with_logger(logger) do
        test_triangulate(number_of_nodes, x_bounds, y_bounds, save_steps)
    end
end


# CLI call from a terminal: call the main function with CLI args
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
;
