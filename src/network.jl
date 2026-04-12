struct SDNNode
    str_label::String
    id::Int
    loc_x::Float64
    loc_y::Float64
end

struct SDNEdge
    link_a::Int
    link_b::Int
    length::Float64
end

load_dognet() = load_network("networks/dognet.dat")
load_coronet_conus() = load_network("networks/coronet-conus.dat")
load_cost266() = load_network("networks/cost266.dat")

# Codes - Internal integer representation for the graph
function codes_to_labels(g::MetaGraph, codes::Vector{Int})
    [label_for(g, c) for c in codes]
end

# Labels are the way we index the graph
function labels_to_codes(g::MetaGraph, node_labels)
    [code_for(g, l) for l in node_labels]
end

function load_network(filename::AbstractString)
    contents = split(String(read(filename)), "\n")
    is_parsing_nodes = false

    mg = MetaGraph(
        SimpleWeightedGraph();
        label_type = Int,
        vertex_data_type = Tuple{Float64,Float64,String},
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
            loc_x = parse(Float64, lsplit[3])
            loc_y = parse(Float64, lsplit[2])
            str_label = replace(lsplit[4], '_' => ' ')
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
        # The number of hops
        mg[e.link_a, e.link_b] = 1.0 # e.length
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
        aspect_ratio = :equal,
    )
end

"""
    components(g::MetaGraph)

Collect the connected components of a graph into a list of graphs.
"""
function components(g::MetaGraph)
    component_indices = connected_components(g)
    [first(induced_subgraph(g, idx)) for idx in component_indices]
end

"""
    incident_edges(g::MetaGraph, v::Int)

Given some node `v`, return the edges incident to think node.
"""
function incident_edges(g::MetaGraph, v::Int)
    v_idx = code_for(g, v)
    [e for e in edges(g) if e.src == v_idx || e.dst == v_idx]
end

"""
    get_distance_matrix(g::AbstractGraph)

Return the distance matrix of the graph using dijkstra's algorithm.
"""
function get_distance_matrix(g::AbstractGraph)
    V = nv(g)
    dists = fill(Inf, V, V)

    for i = 1:V
        algo_result = dijkstra_shortest_paths(g, i)
        dists[i, :] = algo_result.dists
    end

    return dists
end

function surviving_nodes(g::AbstractGraph, s::Vector{Int}, a::Attacks)
    # Return the surviving nodes (not components) given an attack 
    ag = attack_graph(g, a)

    Iterators.flatten([
        collect(labels(c)) for
        c in components(ag) if !isdisjoint([code_for(g, l) for l in labels(c)], s)
    ]) |>
    collect |>
    sort
end

function surviving_nodes(
    g::AbstractGraph,
    ps::Vector{Int}, # Primary controllers 
    bs::Vector{Int}, # Backup controllers 
    a::Attacks, # attacks
)
    ag = attack_graph(g, a)

    # Join the two controllers!
    cs = union(ps, bs)
    Iterators.flatten([
        collect(labels(c)) for
        c in components(ag) if !isdisjoint([code_for(g, l) for l in labels(c)], cs)
    ]) |>
    collect |>
    sort
end

function strings_to_indices(g::MetaGraph, str_labels::Vector{String})
    str_to_idx = Dict(g[label_for(g, v)][3] => v for v in vertices(g))
    return [str_to_idx[s] for s in str_labels if haskey(str_to_idx, s)]
end

function surviving_nodes(g::AbstractGraph, ps::Placements, a::Attacks)
    surviving_nodes(g, ps.pc, ps.bc, a)
end

function shortest_path_tree(g::MetaGraph, root::Int)
    root_idx = code_for(g, root)
    ds = dijkstra_shortest_paths(g, root_idx)

    tree = MetaGraph(
        SimpleWeightedGraph();
        label_type = Int,
        vertex_data_type = Tuple{Float64,Float64,String},
        edge_data_type = Float64,
        weight_function = identity,
        default_weight = Int,
    )

    for v_idx in vertices(g)
        v_label = label_for(g, v_idx)
        tree[v_label] = g[v_label]
    end

    for (v_idx, p_idx) in enumerate(ds.parents)
        # Skip source and unreachable nodes
        if p_idx != 0 && p_idx != v_idx
            w_label = label_for(g, p_idx)
            v_label = label_for(g, v_idx)
            # Retrieve the edge weight from the original MetaGraph 
            weight = g[w_label, v_label]
            tree[w_label, v_label] = weight
        end
    end
    tree
end


#### Experiments helpers

# Attack costs is based on betweenness centrality
function calculate_attack_costs(g::MetaGraph)
    costs = Dict{Int,Float64}()
    # Calculate the betweenness centrality for all vertices
    bc = betweenness_centrality(g)

    for v = 1:nv(g)
        costs[v] = Float64(bc[v])
    end
    return costs
end

# Attacker budget is the ratio of all attack costs 
function calculate_attacker_budget(costs::Dict{Int,Float64}, ρ::Float64)
    return ρ * sum(values(costs))
end

# Number of arrivals per graph
function calculate_arrivals(g::MetaGraph; multiplier::Float64 = 10.0)
    λ = Dict{Int,Float64}()
    for v = 1:nv(g)
        λ[v] = multiplier * degree(g, v)
    end
    return λ
end

# The value of the capacity for all controllers
function calculate_uniform_capacity(λ::Dict{Int,Float64}, P_prime::Int; κ = 1.0)
    Λ = sum(values(λ))
    return κ * Λ / P_prime
end

# Dictionary mapping for the controller capacity (uniform
# capacity using `calculate_uniform_capacity`)
function calculate_capacity_dict(g::MetaGraph, capacity_value::Float64)
    C_dict = Dict{Int,Float64}()
    for v = 1:nv(g)
        C_dict[v] = capacity_value
    end
    return C_dict
end
