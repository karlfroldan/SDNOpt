struct SDNNode
    str_label :: String
    id :: Int
    loc_x :: Float64
    loc_y :: Float64
end

struct SDNEdge
    link_a :: Int
    link_b :: Int
    length :: Float64
end

load_dognet() = load_network("networks/dognet.dat")
load_coronet_conus() = load_network("networks/coronet_conus-l.dat")
load_cost266() = load_network("networks/net-cost266-l.dat")

# Codes - Internal integer representation for the graph
function codes_to_labels(g::MetaGraph, codes::Vector{Int})
    [label_for(g, c) for c in codes]
end

# Labels are the way we index the graph
function labels_to_codes(g::MetaGraph, node_labels)
    [code_for(g, l) for l in node_labels]
end

function load_network(filename :: AbstractString)
    contents = split(String(read(filename)), "\n")
    is_parsing_nodes = false

    mg = MetaGraph(
        SimpleWeightedGraph();
        label_type = Int,
        vertex_data_type = Tuple{Float64, Float64, String},
        edge_data_type = Float64,
        weight_function = identity,
        default_weight = Inf,
    )

    nodes = SDNNode[]
    edges = SDNEdge[]
    
    for line in contents
        if startswith(line, "param: sa_Nodes")
            # We are still parsing nodes
            is_parsing_nodes = true
            continue
        elseif startswith(line, "param: sa_Links")
            is_parsing_nodes = false
            continue
        elseif startswith(line, ";") || line == ""
            continue
        end

        lsplit = split(line)

        if is_parsing_nodes
            id = parse(Int, lsplit[1])
            loc_x = parse(Float64, lsplit[2])
            loc_y = parse(Float64, lsplit[3])
            str_label = lsplit[4]
            new_node = SDNNode(str_label, id, loc_x, loc_y)
            # @show new_node
            push!(nodes, new_node)
        else
            from_n = parse(Int, lsplit[2])
            to_n = parse(Int, lsplit[3])
            dist = parse(Float64, lsplit[4])
            new_edge = SDNEdge(from_n, to_n, dist)
            push!(edges, new_edge)
        end
    end

    for n in nodes
        mg[n.id] = (n.loc_x, n.loc_y, n.str_label)
    end
    
    for e in edges
        @printf("%d -> %d => %f\n", e.link_a, e.link_b, e.length)
        mg[e.link_a, e.link_b] = e.length
    end

    mg
end

"""
    attack_graph(g::MetaGraph, attacked_labels::AbstractArray{Int})

Given the `attacked_labels`, remove these nodes from the graph as if they are attacked.
"""
function attack_graph(g::MetaGraph, attacked_labels::AbstractArray{Int})
    attacked_indices = [code_for(g, l) for l in attacked_labels if haskey(g, l)]
    all_indices = vertices(g)
    surviving_indices = setdiff(all_indices, attacked_indices)
    new_g, _ = induced_subgraph(g, collect(surviving_indices))
    new_g
end

function plot_network(g::MetaGraph; title = nothing)
    locs = [g[label_for(g, v)] for v in vertices(g)]
    
    x_coords = [loc[1] for loc in locs]
    y_coords = [loc[2] for loc in locs]
    node_labels = [loc[3] for loc in locs]

    graphplot(
        g,
        x = x_coords,
        y = y_coords,
        names = node_labels,
        curves = false,
        nodeshape = :circle,
        nodesize = 0.30,
        fontsize = 8,
        linecolor = :gray,
        markercolor = :white,
        title = title,
        aspect_ratio = :equal
    )
end

function components(g::MetaGraph)
    component_indices = connected_components(g)
    [first(induced_subgraph(g, idx)) for idx in component_indices]
end

function distance_matrix(g::MetaGraph)
    n = nv(g)
    dist_matrix = fill(Inf, n, n)
    
    for i in 1:n
        ds = dijkstra_shortest_paths(g, i)
        dist_matrix[i, :] = ds.dists
    end
    
    return dist_matrix
end

function incident_edges(g :: MetaGraph, v :: Int)
    v_idx = code_for(g, v)
    [e for e ∈ edges(g) if e.src == v_idx || e.dst == v_idx]
end

function get_distance_matrix(g::AbstractGraph)
    V = nv(g)
    dists = fill(Inf, V, V)

    for i in 1:V
        algo_result = dijkstra_shortest_paths(g, i)
        dists[i, :] = algo_result.dists
    end

    return dists
end