
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
    
    # Custom Blue Gradient for Controller Load
    # lightblue = Low Load, darkblue = High Load
    blue_heatmap = cgrad([:lightblue, :darkblue]) 

    n_colors = fill(plot_color(markercolor), nv(g))
    n_shapes = fill(:circle, nv(g))
    n_strokes = fill(plot_color(markerstrokecolor), nv(g))

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
            if v in controller_set
                # Apply Blue Heatmap only to active controllers
                n_colors[v] = blue_heatmap[rand()] 
                n_shapes[v] = :hexagon
                n_strokes[v] = plot_color(:black)
            else
                # Regular switches remain as is
                n_colors[v] = plot_color(markercolor)
                n_shapes[v] = :circle
                n_strokes[v] = plot_color(markerstrokecolor)
            end
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
    # Heatmap specific legend entries
    scatter!(p, [NaN], [NaN], label="Low Load (Controller)", markercolor=:lightblue, markershape=:hexagon, markerstrokecolor=:black, framestyle=:none)
    scatter!(p, [NaN], [NaN], label="High Load (Controller)", markercolor=:darkblue, markershape=:hexagon, markerstrokecolor=:black, framestyle=:none)
    
    scatter!(p, [NaN], [NaN], label="Active Switch", markercolor=markercolor, markershape=:circle, markerstrokecolor=markerstrokecolor, framestyle=:none)
    scatter!(p, [NaN], [NaN], label="Attack Node", markercolor=inactive_color, markershape=:rect, markerstrokecolor=:red, framestyle=:none)
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

function calculate_tikz_label_position(g::MetaGraph, v::Int, x_coords::AbstractVector, y_coords::AbstractVector)
    neighbor_nodes = neighbors(g, v)
    
    # Default for isolated nodes
    if isempty(neighbor_nodes)
        return "above"
    end

    xv, yv = x_coords[v], y_coords[v]
    angles = Float64[]

    for u in neighbor_nodes
        xu, yu = x_coords[u], y_coords[u]
        push!(angles, atan(yu - yv, xu - xv))
    end

    target_angle = 0.0

    if length(angles) == 1
        # Place opposite to the single edge
        target_angle = mod2pi(angles[1] + π)
    else
        sort!(angles)
        max_gap = 0.0
        
        n_angles = length(angles)
        for i in 1:n_angles
            θ1 = angles[i]
            # Wrap around for the last angle
            θ2 = i == n_angles ? angles[1] + 2π : angles[i+1]
            gap = θ2 - θ1
            
            if gap > max_gap
                max_gap = gap
                target_angle = mod2pi(θ1 + gap / 2.0)
            end
        end
    end

    # Discretize angle to 8 compass directions
    # 0 is right, π/2 is above, π is left, 3π/2 is below
    octant = round(Int, target_angle / (π / 4)) % 8

    positions = Dict(
        0 => "right",
        1 => "above right",
        2 => "above",
        3 => "above left",
        4 => "left",
        5 => "below left",
        6 => "below",
        7 => "below right"
    )

    return positions[octant]
end


function generate_tikz(
    g::MetaGraph;
    controllers::Union{Placements, Nothing}=nothing, 
    attacks::Union{Vector{Int}, Nothing}=nothing, 
    show_labels::Bool=true,
    save_path::Union{String, Nothing}=nothing,
    rotate_deg::Real=0.0,
    flip_x::Bool=false,
    flip_y::Bool=false
)
    x_coords, y_coords, locs = extract_coordinates(g)
    
    # Transformations
    transformed_locs = Dict()
    rad = deg2rad(rotate_deg)
    
    for v in vertices(g)
        x, y = locs[v][1], locs[v][2]
        
        # Apply rotation
        if rotate_deg != 0.0
            x_rot = x * cos(rad) - y * sin(rad)
            y_rot = x * sin(rad) + y * cos(rad)
            x, y = x_rot, y_rot
        end
        
        # Apply flips
        if flip_x
            x = -x
        end
        if flip_y
            y = -y
        end

        transformed_locs[v] = tuple(x, y, locs[v][3:end]...)
    end
    
    locs = transformed_locs

    x_coords = [locs[v][1] for v in vertices(g)]
    y_coords = [locs[v][2] for v in vertices(g)]

    pc_set = Set{Int}()
    bc_set = Set{Int}()
    if !isnothing(controllers)
        union!(pc_set, controllers.pc)
        union!(bc_set, controllers.bc)
    end
    
    controller_set = union(pc_set, bc_set)
    ctrl_vec = collect(controller_set)
    atk_vec = isnothing(attacks) ? Int[] : collect(Int, attacks)
    
    active_components = get_active_components(g, ctrl_vec, atk_vec)
    attack_set = Set(atk_vec)

    node_to_comp = Dict{Int, Int}()
    for (i, comp) in enumerate(active_components)
        for v in comp
            node_to_comp[v] = i
        end
    end

    io = IOBuffer()
    println(io, "\\begin{tikzpicture}[scale=0.25]\n")
    
    println(io, "  % --- STYLES ---")
    println(io, "  \\tikzset{")
    println(io, "    every label/.append style={font=\\scriptsize, label distance=-1pt},")
    println(io, "    v/.style = {circle, fill=black, inner sep=1.2pt},")
    println(io, "    v_pc/.style = {circle, fill=cyan, inner sep=1.2pt},")
    println(io, "    v_bc/.style = {circle, fill=olive, inner sep=1,2pt},")
    println(io, "    v_dead/.style = {circle, fill=gray!50, inner sep=1.2pt},")
    println(io, "    attack_box/.style = {rectangle, draw=red, thick, inner sep=2.6pt, fill=none},")
    println(io, "    e/.style = {-, thin, black},")
    println(io, "    e_dead/.style = {-, thin, gray!30}")
    println(io, "  }\n")
    
    println(io, "  % --- NODES ---")
    
    sorted_vertices = sort(collect(vertices(g)))
    
    for v in sorted_vertices
        x = round(locs[v][1], digits=4)
        y = round(locs[v][2], digits=4)
        
        tikz_id = "node_$v"

        if v in pc_set
            style = "v_pc"
        elseif v in bc_set
            style = "v_bc"
        elseif v in attack_set
            style = "v_dead"
        elseif isempty(attack_set) || haskey(node_to_comp, v)
            style = "v"
        else
            style = "v_dead"
        end
        
        if show_labels
            label_pos = calculate_tikz_label_position(g, v, x_coords, y_coords)
            node_options = "$style, label={$label_pos:{$v}}"
        else
            node_options = "$style"
        end
        
        println(io, "  \\node[$node_options] ($tikz_id) at ($x, $y) {};")
        
        if v in attack_set
            println(io, "  \\node[attack_box] at ($x, $y) {};")
        end
    end
    
    println(io, "\n  % --- LINKS ---")
    for e in edges(g)
        u, v = src(e), dst(e)
        
        tikz_u = "node_$u"
        tikz_v = "node_$v"
        
        if isempty(attack_set) || (haskey(node_to_comp, u) && haskey(node_to_comp, v) && node_to_comp[u] == node_to_comp[v])
            edge_style = "e"
        else
            edge_style = "e_dead"
        end
        
        println(io, "  \\draw[$edge_style] ($tikz_u) -- ($tikz_v);")
    end

    println(io, "\n\\end{tikzpicture}")
    latex_str = String(take!(io))

    if !isnothing(save_path)
        open(save_path, "w") do f
            write(f, latex_str)
        end
    end

    return latex_str
end