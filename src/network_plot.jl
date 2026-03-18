
function get_active_nodes(g::MetaGraph, controllers::Vector{Int}, attacks::Vector{Int})
    attack_set = Set(attacks)
    controller_set = Set(controllers)
    active_nodes = Set{Int}()
    
    surviving_vertices = [v for v in vertices(g) if !(v in attack_set)]
    sub_g, vmap = induced_subgraph(g, surviving_vertices)
    
    for comp in connected_components(sub_g)
        orig_comp = [vmap[v] for v in comp]
        if !isdisjoint(orig_comp, controller_set)
            union!(active_nodes, orig_comp)
        end
    end
    
    return active_nodes
end
function get_active_components(g::MetaGraph, controllers::Vector{Int}, attacks::Vector{Int})
    attack_set = Set(attacks)
    controller_set = Set(controllers)
    active_components = Vector{Set{Int}}()
    
    surviving_vertices = [v for v in vertices(g) if !(v in attack_set)]
    sub_g, vmap = induced_subgraph(g, surviving_vertices)
    
    for comp in connected_components(sub_g)
        orig_comp = Set([vmap[v] for v in comp])
        if !isdisjoint(orig_comp, controller_set)
            push!(active_components, orig_comp)
        end
    end
    
    return active_components
end
function extract_coordinates(g::MetaGraph)
    locs = [g[label_for(g, v)] for v in vertices(g)]
    x_coords = [loc[1] for loc in locs]
    y_coords = [loc[2] for loc in locs]
    return x_coords, y_coords, locs
end

function determine_network_styles(g::MetaGraph, controllers::Vector{Int}, attacks::Vector{Int}, markercolor, markerstrokecolor, linecolor)
    inactive_color = plot_color("#E5E5E5")
    comp_palette = plot_color.([:black, :darkgreen, :darkorange, :purple, :teal, :saddlebrown, :magenta, :olive])
    
    n_colors = fill(plot_color(markercolor), nv(g))
    n_shapes = fill(:circle, nv(g))
    n_strokes = fill(plot_color(markerstrokecolor), nv(g))

    if isempty(controllers) && isempty(attacks)
        e_colors = (s, d, w) -> plot_color(linecolor)
        return n_colors, n_shapes, n_strokes, e_colors, inactive_color
    end

    active_components = get_active_components(g, controllers, attacks)
    attack_set = Set(attacks)
    controller_set = Set(controllers)

    node_to_comp = Dict{Int, Int}()
    for (i, comp) in enumerate(active_components)
        for v in comp
            node_to_comp[v] = i
        end
    end

    for v in vertices(g)
        if v in attack_set
            n_colors[v] = plot_color(v in controller_set ? :lightblue : inactive_color)
            n_shapes[v] = :rect
            n_strokes[v] = plot_color(:red)
        elseif haskey(node_to_comp, v)
            n_colors[v] = plot_color(v in controller_set ? :blue : markercolor)
            n_shapes[v] = :circle
            n_strokes[v] = plot_color(markerstrokecolor)
        else
            n_colors[v] = inactive_color
            n_shapes[v] = :circle
            n_strokes[v] = plot_color(:gray)
        end
    end

    e_color_dict = Dict{Tuple{Int, Int}, typeof(inactive_color)}()
    for e in edges(g)
        u, v = src(e), dst(e)
        if haskey(node_to_comp, u) && haskey(node_to_comp, v) && node_to_comp[u] == node_to_comp[v]
            c_idx = node_to_comp[u]
            c = comp_palette[mod1(c_idx, length(comp_palette))]
        else
            c = inactive_color
        end
        e_color_dict[(u, v)] = c
        e_color_dict[(v, u)] = c
    end
    
    e_colors = (s, d, w) -> get(e_color_dict, (s, d), inactive_color)

    return n_colors, n_shapes, n_strokes, e_colors, inactive_color
end

function add_network_legend!(p, markercolor, markerstrokecolor, inactive_color)
    scatter!(p, [NaN], [NaN], label="Controller", markercolor=:blue, markershape=:circle, markerstrokecolor=markerstrokecolor, framestyle=:none)
    scatter!(p, [NaN], [NaN], label="Switch", markercolor=markercolor, markershape=:circle, markerstrokecolor=markerstrokecolor, framestyle=:none)
    scatter!(p, [NaN], [NaN], label="Attack", markercolor=inactive_color, markershape=:rect, markerstrokecolor=:red, framestyle=:none)
    scatter!(p, [NaN], [NaN], label="Inactive Node", markercolor=inactive_color, markershape=:circle, markerstrokecolor=inactive_color, framestyle=:none)
    scatter!(p, [NaN], [NaN], label="Attacked Controller", markercolor=:lightblue, markershape=:rect, markerstrokecolor=:red, framestyle=:none)
end

function add_network_labels!(p, g::MetaGraph, locs, x_coords, y_coords, fontsize)
    y_spread = maximum(y_coords) - minimum(y_coords)
    offset = y_spread * 0.02

    for i in vertices(g)
        has_upper_edge = any(y_coords[j] > y_coords[i] for j in neighbors(g, i))
        
        if has_upper_edge
            lbl_y = y_coords[i] - offset
            valign = :top
        else
            lbl_y = y_coords[i] + offset
            valign = :bottom
        end
        
        annotate!(p, x_coords[i], lbl_y, text(locs[i][3], fontsize, :center, valign))
    end
end

function plot_network(
    g::MetaGraph; 
    title = "",
    show_labels::Bool = true,
    show_legend::Bool = false,
    nodesize = 2.30,
    fontsize = 4,
    legendfontsize = 6,
    markercolor = :white,
    markerstrokecolor = :black,
    linecolor = :gray,
    controllers::Vector{Int} = Int[],
    attacks::Vector{Int} = Int[],
    save_path::Union{String, Nothing} = nothing
)
    x_coords, y_coords, locs = extract_coordinates(g)
    
    n_colors, n_shapes, n_strokes, e_colors, inactive_color = determine_network_styles(
        g, controllers, attacks, markercolor, markerstrokecolor, linecolor
    )

    p = graphplot(
        g,
        x = x_coords,
        y = y_coords,
        curves = false,
        nodeshape = n_shapes,
        nodesize = nodesize,
        edgecolor = e_colors,
        markercolor = n_colors,
        markerstrokecolor = n_strokes,
        title = title,
        aspect_ratio = :equal,
        legend = show_legend ? :bottom : :none,
        legendfontsize = legendfontsize,
        legend_columns = 2,
        margin = -10Plots.mm,
        label = "" 
    )

    if show_legend
        add_network_legend!(p, markercolor, markerstrokecolor, inactive_color)
    end

    if show_labels
        add_network_labels!(p, g, locs, x_coords, y_coords, fontsize)
    end

    if !isnothing(save_path)
        savefig(p, save_path)
    end

    return p
end