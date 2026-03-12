const UNIQUE_NODE :: Int = 1

# Attacks
function add_uniqueness_constraint!(
    m::Model, 
    V::Int, 
    history::Vector{Attacks}, 
    diff::Int, 
    v1::AbstractArray{VariableRef}
)
    for T_nodes ∈ history
        not_T = setdiff(1:V, T_nodes)
        @constraint(m, 
            sum(v1[v] for v ∈ T_nodes) - sum(v1[v] for v ∈ not_T) ≤ length(T_nodes) - diff
        )
    end
end

# Placements with one variable 
function add_uniqueness_constraint!(
    m::Model, 
    V::Int, 
    history::Vector{Placements}, 
    diff::Int, 
    vp::AbstractArray{VariableRef}
)
    for sol ∈ history
        Tp = sol.pc
        not_Tp = setdiff(1:V, Tp)
        @constraint(m, 
            sum(vp[v] for v ∈ Tp) - sum(vp[v] for v ∈ not_Tp) ≤ length(Tp) - diff
        )
    end
end

# Placements with primary and backup controllers
function add_uniqueness_constraint!(
    m::Model, 
    V::Int, 
    history::Vector{Placements}, 
    diff::Int, 
    vp::AbstractArray{VariableRef}, 
    vb::AbstractArray{VariableRef}
)
    for sol ∈ history
        Tp = sol.pc
        Tb = sol.bc
        not_Tp = setdiff(1:V, Tp)
        not_Tb = setdiff(1:V, Tb)
        
        @constraint(m, 
            sum(vp[v] for v ∈ Tp) - sum(vp[v] for v ∈ not_Tp) +
            sum(vb[v] for v ∈ Tb) - sum(vb[v] for v ∈ not_Tb) ≤ length(Tp) + length(Tb) - diff
        )
    end
end

macro uniqueness_constraints(m, V, history, diff, vars...)
    return quote
        add_uniqueness_constraint!(
            $(esc(m)), 
            $(esc(V)), 
            $(esc(history)), 
            $(esc(diff)), 
            $(map(esc, vars)...)
        )
    end
end


function get_U(
    g :: MetaGraph, 
    BCC :: Float64,
    dists :: Matrix{Float64};
    tol = 1e-9
)
    V = nv(g)
    [
        (v, w) for (v, w) ∈ Base.Iterators.product(1:V, 1:V)
        if v ≤ w && dists[v, w] > BCC - tol
    ]
end

function get_Wv(
    g :: MetaGraph, 
    BSC :: Float64,
    dists :: Matrix{Float64};
    tol = 1e-9
)
    V = nv(g)
    W = Dict{Int, Vector{Int}}()

    for v ∈ 1:V
        vicinity = [w for w ∈ 1:V if dists[v, w] ≤ BSC + tol] # && v ≤ w]
        if !isempty(vicinity)
            W[v] = vicinity
        end
    end

    W
end

## Paper: Max-Min Optimization of Controller Placements
##        vs Min-Max Optimization of Attacks on Nodes in Service Networks

# Formulation 3.1: Max-min Controller Placement optimization problem (CPOP)
function cpop(
    g :: MetaGraph, 
    M :: Int, 
    # The nodes that we want to attack.
    attacks :: Vector{Vector{Int}};
    dists :: Union{Nothing, Matrix{Float64}} = nothing,
    BCC :: Union{Nothing, Float64}=nothing,
    BSC :: Union{Nothing, Float64}=nothing,
    optim = DEFAULT_OPTIM,
)
    # Find an M-node controller placement given the considered set of attacks
    model = Model(optim)
    set_silent(model)

    dists = isnothing(dists) ? get_distance_matrix(g) : dists

    V = nv(g)
    alen = length(attacks)

    # equals 1 iff controller is placed at loc v 
    @variable(model, s[1:V], Bin)

    # 1 if v survives attack a
    @variable(model, y[1:V, 1:alen], Bin)

    @variable(model, Y >= 0)

    # 3a
    @objective(model, Max, Y)

    # 3b - # of controllers should be M
    @constraint(model, sum(s) == M)

    # 3c - all attacked nodes are zeroed out
    for a ∈ 1:alen
        @constraint(model, y[attacks[a], a] .== 0)
    end

    # 3d - All nodes in components without controllers are zeroed out
    for a ∈ 1:alen
        ag = attack_graph(g, attacks[a])
        cs = components(ag)

        # The surviving nodes after the attack
        for comp ∈ cs
            # Vertices of comp are 1:size(comp)
            # Get their labels instead and then map those labels to codes.
            comp_labels = collect(labels(comp))
            V_c = length(comp_labels)

            comp_codes = labels_to_codes(g, comp_labels)
            @constraint(
                model,
                sum(y[comp_codes, a]) ≤ V_c * sum(s[comp_codes])
            )
        end
        # 3e - Bounding objective value
        @constraint(model, Y ≤ sum(y[:, a]))
    end

    # controllers must be within BCC distance of each other
    if !isnothing(BCC)
        U = get_U(g, BCC, dists)
        @constraint(model, [(v, w) ∈ U], s[v] + s[w] ≤ 1)
    end

    if !isnothing(BSC)
        W = get_Wv(g, BSC, dists)
        @constraint(model, [v ∈ 1:V], sum(s[W[v]]) ≥ 1)
    end

    time_taken = @solve_problem!(model)

    (
        controllers = to_indices(s),
        objective_value = objective_value(model),
        model = model,
        time = time_taken
    )
end

# Formulation 3.2: Min-Max Node Attack Optimization Problem (NAOP)
function naop(
    g :: MetaGraph, 
    K :: Int, 
    placements :: Vector{Vector{Int}};
    optim = DEFAULT_OPTIM,
)
    model = Model(optim)
    set_silent(model)

    V = nv(g)
    E = ne(g)
    S = length(placements)
    # Attack variable 
    @variable(model, a[1:V], Bin)
    @variable(model, t[1:E], Bin)
    @variable(model, z[1:V, 1:S], Bin)
    @variable(model, Z >= 0)

    @objective(model, Min, Z)

    # 5b - K-node attack 
    @constraint(model, sum(a) == K)

    # 5c - Link is down after attack a
    for v ∈ 1:V
        for (e, edge) ∈ enumerate(edges(g))
            if edge.src == v || edge.dst == v 
                @constraint(model, t[e] ≥ a[v])
            end
        end
    end

    # 5d - Link is down after attack a
    for (e, edge) ∈ enumerate(edges(g))
        α, β = edge.src, edge.dst
        @constraint(model, t[e] ≤ a[α] + a[β])
    end

    # 5e - Node v does not survive attack a if node v is attacked directly
    for s ∈ 1:S 
        @constraint(model, z[:, s] .≤ 1 .- a)
    end

    # 5f - Node v survives attack a when placement s is assumed if
    #      node v is not directly attacked and its location contains 
    #      a controller
    for (s, placement) ∈ enumerate(placements)
        @constraint(model, z[placement, s] .≥ 1 .- a[placement])
    end

    # 5g and 5h - Make sure that if link e is available after attack a, then
    #             its end nodes either simultaneaously suvive or both are out
    #             of service.
    for s ∈ 1:S
        for (e, edge) ∈ enumerate(edges(g))
            α, β = edge.src, edge.dst
            @constraint(model, z[α, s] ≥ z[β, s] - t[e])
            @constraint(model, z[β, s] ≥ z[α, s] - t[e])
        end
    end

    # 5i 
    for s ∈ 1:S 
        @constraint(model, Z ≥ sum(z[:, s]))
    end

    time_taken = @solve_problem!(model)
    
    (
        attack = to_indices(a),
        objective_value = objective_value(model),
        model = model,
        time = time_taken
    )
end

function maximum_sc_delay(
    g :: MetaGraph,
    P :: Union{Int, IntBound},
    dists :: Matrix{Float64},
    BSC :: Float64, 
    BCC :: Float64;
    tol = 1e-9,
    optim=DEFAULT_OPTIM,
)
    V = nv(g)

    @assert BSC > 0.0
    @assert BCC > 0.0

    P′, P″ = @unpack_bounds(P)

    m = Model(optim)
    set_silent(m)

    # A1i
    @variable(m, D)
    
    # A1h: Controller placements
    @variable(m, s[1:V], Bin)
    
    # A1h: Controller assignments
    Wv = get_Wv(g, BSC, dists; tol=tol)
    @variable(m, z[v in keys(Wv), w in Wv[v]], Bin)

    # A1a: Objective function
    @objective(m, Min, D)

    # A1b
    if P′ == P″
        @constraint(m, sum(s) == P′)
    else
        @constraint(m, P′ ≤ sum(s) ≤ P″)
    end

    # A1c
    U = get_U(g, BCC, dists; tol=tol)
    for (v, w) ∈ U 
        @constraint(m, s[v] + s[w] ≤ 1)
    end

    # A1d
    @constraint(m, [v ∈ keys(Wv)], sum(z[v, w] for w ∈ Wv[v]) == 1)

    # A1e
    @constraint(m, [v ∈ keys(Wv), w ∈ Wv[v]; w != v], z[v, w] ≤ s[w])

    # A1f
    @constraint(m, [v ∈ keys(Wv); v ∈ Wv[v]], z[v, v] == s[v])

    # A1g
    @constraint(m, [v ∈ keys(Wv)], D ≥ sum(dists[v, w] * z[v, w] for w ∈ Wv[v]))

    @solve_problem!(m)
    
    return objective_value(m)
end

### FEASIBILITY GENERATION 

# Primary controllers only, controllers capacity not considered
# max Z 
# s.t. sum {v in VERTICES} y[v] = M
#      {{v, w} in U} y[v] + y[w] ≤ 1
#      {v in VERTICES} sum {w ∈ W(v)} y[w] ≥ 1
#      l=1,2,…,L, sum {v ∈ T(l)} y[v] ≤ m 
#      {v in VERTICES} y[v] ∈ B
# where T(l) denotes the set of nodes where controllers of the generated
# placement number l∈1,2,…,L are located.

# For cost266 
# - P* = 3, tight CCD BCC = 1500, P* = 5 BCC = 2000 
function generate_controller_placement(
    g :: MetaGraph,
    P :: Union{Int, IntBound},
    B :: Union{Int, IntBound};
    C :: Union{Int, Nothing} = nothing,
    optim = DEFAULT_OPTIM,
    BCC :: Float64 = 0.0,
    BSC :: Float64 = 0.0,
    control_capacity :: Dict{Int, Float64} = Dict{Int, Float64}(),
    control_demand :: Dict{Int, Float64} = Dict{Int, Float64}(),
    placement_list :: Vector{Placements} = Placement[],
    placement_difference :: Int = 1,
    dists :: Union{Matrix{Float64}, Nothing} = nothing,
    history :: Vector{Placements} = Placements[]
)
    # Unfortunately, nvidia cuOpt does not support feasibility.
    m = Model(optim)
    set_silent(m)
    dists = isnothing(dists) ? get_distance_matrix(g) : dists
    V = nv(g)
    local x
    has_capacity = false

    if !isempty(control_capacity) || !isempty(control_demand)
        @assert (!isempty(control_capacity) && !isempty(control_demand)) "Both control capacity and control demand should not be empty"
        has_capacity = true
    end
    
    #  # if P isa Tuple

    P′, P″ = @unpack_bounds P 
    B′, B″ = @unpack_bounds B

    # Primary controller list
    @variable(m, y[1:V], Bin)

    # Backup controller list
    @variable(m, x[1:V], Bin)
    # Controller assignment list
    @variable(m, z[1:V, 1:V], Bin)

    if has_capacity
        @variable(m, Z)
        @objective(m, Min, Z)
    else
        @objective(m, FEASIBILITY_SENSE, 0)
    end
    
    # (2b) Total number of primary and backup controllers should equal M
    if !isnothing(C)
        @constraint(m, sum(x) + sum(y) == C)
    end 
    # (2c) Primary controllers constraint 
    @constraint(m, P′ ≤ sum(y) ≤ P″)

    # (2d) Backup controllers constraint
    @constraint(m, B′ ≤ sum(x) ≤ B″)

    if BCC > 0.0
        # println("(2e) activated")
        U = get_U(g, BCC, dists)

        # (2e) Controllers must be able to reach each other inside BCC delay.
        for (v, w) ∈ U
            @constraint(m, y[v] + y[w] ≤ 1)
        end
    end

    if BSC > 0.0
        # println("(2f) activated")
        W = get_Wv(g, BSC, dists)
        # @show W
        # (2f) Must have one controller in the vicinity of SCD.

        for (v, Wv) ∈ W
            # @show Wv
            @constraint(m, sum(y[Wv]) ≥ 1)
        end
    end

    # (2g) Can only be a primary or backup. Not both
    @constraint(m, y .+ x .≤ 1)

    # (2h) Controller assignment 
    for v ∈ 1:V
        @constraint(m, sum(z[v, :]) == 1)
    end

    # (2i)
    for v ∈ 1:V
        @constraint(m, z[v, :] .≤ y)
    end

    # (2j) Capacity constraint
    if has_capacity
        cd = [control_demand[v] for v ∈ 1:V]
        cc = [control_capacity[v] for v ∈ 1:V]
        for w ∈ 1:V
            @constraint(m, cd ⋅ z[:, w] ≤ cc[w] * y[w] + Z)
        end
    end

    # (2k) Uniqueness constraint
    @uniqueness_constraints(m, V, history, UNIQUE_NODE, y, x)

    time_taken = @solve_problem!(m)

    placement = Placements(to_indices(y), to_indices(x))
    
    (
        controllers = placement,
        time = time_taken,
        model = m,
        # objective_value = Z,
    )
end


### Paper:
### Finding Optimal Mixed-Strategies in a Matrix Game
### Between the Attacker and the Network Operator


## Master Problem 
function mixed_strategies_master(
    g :: MetaGraph,
    placementset :: Vector{Placements},
    attackset :: Vector{Attacks};
    optim = DEFAULT_OPTIM,
)
    m = Model(optim)
    set_silent(m)

    psize = length(placementset)
    asize = length(attackset)

    @variable(m, q[1:psize] ≥ 0)
    @variable(m, y)

    @objective(m, Max, y)
    
    # (22b) the probability should equal 1.
    @constraint(m, sum(q) == 1.0)

    # (22c)
    V_mat = zeros(Int, asize, psize)
    for i ∈ 1:asize
        for j ∈ 1:psize
            V_mat[i, j] = length(surviving_nodes(g, placementset[j], attackset[i]))
        end
    end

    @constraint(m, con[i=1:asize], y ≤ sum(V_mat[i, j] * q[j] for j ∈ 1:psize))
    time_taken = @solve_problem!(m)

    raw_duals = [dual(con[i]) for i ∈ 1:asize]

    (
        time = time_taken,
        model = m,
        objective = objective_value(m),
        q_star = safe_probs(value.(q)),
        p_star = safe_probs(abs.(raw_duals)), 
    )
end

function mixed_strategies_pricing_placement(
    g :: MetaGraph,
    P :: Union{Int, IntBound},
    B :: Union{Int, IntBound},
    attackset :: Vector{Attacks},
    p :: Vector{Float64};
    C :: Union{Int, Nothing} = nothing,
    optim = DEFAULT_OPTIM,
    BCC :: Float64 = 0.0,
    BSC :: Float64 = 0.0,
    dists :: Union{Nothing, Matrix{Float64}} = nothing,
    time_limit :: Union{Nothing, Float64} = nothing,
    history :: Vector{Placements} = Placements[],
)
    P′, P″ = @unpack_bounds P
    B′, B″ = @unpack_bounds B 

    V = nv(g)
    alen = length(attackset)

    dists = isnothing(dists) ? get_distance_matrix(g) : dists

    m = Model(optim) 
    set_silent(m)

    if !isnothing(time_limit)
        set_time_limit_sec(m, time_limit)
    end

    @variables(m, begin
        x[1:V], Bin # Backup controllers 
        y[1:V], Bin # Primary controllers 
        Y[1:alen]
    end)

    # Component survivability
    C_sets = Dict(a => components(attack_graph(g, attack)) for (a, attack) ∈ enumerate(attackset))

    # S=1 if component survives
    # @variable(m, S[1:alen, 1:length(C_sets)], Bin)
    @variable(m, S[a = 1:alen, c = 1:length(C_sets[a])], Bin)

    @objective(m, Max, p ⋅ Y)

    # Maximum number of controllers 
    if !isnothing(C)
        @constraint(m, sum(x) + sum(y) == C)
    end

    # Primary and backup controller bounding 
    @constraints(m, begin
        P′ ≤ sum(y) ≤ P″
        B′ ≤ sum(x) ≤ B″
    end)

    # controllers must be within BCC distance of each other
    if BCC > 0.0 
        U = get_U(g, BCC, dists)
        @constraint(m, [(v, w) ∈ U], y[v] + y[w] ≤ 1)
    end

    if BSC > 0.0 
        W = get_Wv(g, BSC, dists)
        @constraint(m, [v ∈ 1:V], sum(y[W[v]]) ≥ 1)
    end

    # controllers can only be one type.
    @constraint(m, y .+ x .≤ 1)
    
    # @constraint(m, [a in 1:alen], 
    # println(C_sets)
    for a ∈ 1:alen
        for (c, comp) ∈ enumerate(C_sets[a])
            vs = collect(labels(comp))
            @constraint(m, S[a, c] .≤ sum(x[vs] + y[vs]))
        end
    end

    @constraint(m,
        [a ∈ 1:alen],
        # length
        Y[a] == sum([
            let 
                vs = collect(labels(cs))
                length(vs) * S[a, c]
            end
            for (c, cs) ∈ enumerate(C_sets[a])
        ])
    )

    # @uniqueness_constraints(m, V, history, UNIQUE_NODE, y, x)

    time_taken = @solve_problem!(m)

    placements = B″ > 0 ? 
        Placements(to_indices(y), to_indices(x)) : 
        Placements(to_indices(y))

    (
        time = time_taken,
        objective = objective_value(m),
        model = m,
        s = placements,
    )
end

function mixed_strategies_pricing_attack(
    g :: MetaGraph,
    K :: Int,
    placementset :: Vector{Placements},
    q :: Vector{Float64};
    optim = DEFAULT_OPTIM,
    tol :: Float64 = 1e-9,
    time_limit :: Union{Float64} = nothing,
    history :: Vector{Attacks} = Attacks[],
)
    m = Model(optim)
    set_silent(m)

    if !isnothing(time_limit)
        set_time_limit_sec(m, time_limit)
    end

    V = nv(g)

    active_idxs = [i for (i, val) ∈ enumerate(q) if val > tol]
    S′ = placementset[active_idxs]
    q′ = q[active_idxs]

    slen = length(S′)

    @variable(m, F[1:slen] ≥ 0, Int)
    @variable(m, a[1:V], Bin)
    @variable(m, z[1:V, 1:slen], Bin)

    @objective(m, Min, q′ ⋅ F)

    # (25b)
    @constraint(m, sum(a) == K)

    # (25c)
    @constraint(m, [(s, ps) ∈ enumerate(S′)], 
        z[all_controllers(ps), s] .== 1 .- a[all_controllers(ps)])

    # (25d)
    for edge ∈ edges(g)
        α, β = edge.src, edge.dst
        @constraint(m, [s ∈ eachindex(S′)], 
            z[β, s] .≥ z[α, s] .- a[β])

        @constraint(m, [s ∈ eachindex(S′)],
            z[α, s] .≥ z[β, s] .- a[α])
    end

    # (25f)
    @constraint(m, [s ∈ eachindex(S′)], F[s] == sum(z[:, s]))

    # @uniqueness_constraints(m, V, history, UNIQUE_NODE, a)

    time_taken = @solve_problem!(m)

    (
        time = time_taken,
        objective = objective_value(m),
        model = m,
        a = to_indices(a),
    )
end
