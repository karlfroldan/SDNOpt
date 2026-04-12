const UNIQUE_NODE::Int = 1

# Attacks
function add_uniqueness_constraint!(
    m::Model,
    V::Int,
    history::Vector{Attacks},
    diff::Int,
    v1::AbstractArray{VariableRef},
)
    for T_nodes in history
        not_T = setdiff(1:V, T_nodes)
        @constraint(
            m,
            sum(v1[v] for v in T_nodes) - sum(v1[v] for v in not_T) ≤
            length(T_nodes) - diff
        )
    end
end

"""
    add_attack_constraint!(g, m, attack_config)

Add a constraint that bounds the number of attacks done on the network.
Introduces the variable `a[1:V]`.
"""
function add_attack_constraint!(g::MetaGraph, m::Model, attack_config::AttackConfig)
    V = nv(g)

    K = attack_config.K

    # Whether to attack a node or not.
    @variable(m, a[1:V], Bin)

    # The number of nodes to attack.
    if K isa Int
        @constraint(m, sum(a) == K)
    else
        K′, K″ = @unpack_bounds K
        @constraint(m, K′ ≤ sum(a) ≤ K″)
    end
end

"""
    add_edge_propagation_constraints!(g, m, attack_config, placements, q_star)

Adds constraints that will ensure that nodes that do not have access 
to controllers do not survive. This constraint also count the number 
of surviving nodes. Introduces the variables `F[1:|S|]`, `z[1:V, 1:|S|]`
"""
function add_edge_propagation_constraints!(
    g::MetaGraph,
    m::Model,
    attack_config::AttackConfig;
    tol::Float64 = 1e-9,
)
    V = nv(g)

    a = m[:a]

    # Reduce the problem size.
    q = attack_config.q
    active_idxs = [i for (i, val) in enumerate(q) if val > tol]
    S′ = attack_config.placements[active_idxs]

    slen = length(S′)

    @variable(m, F[1:slen] >= 0, Int)
    @variable(m, z[1:V, 1:slen], Bin)

    @constraint(
        m,
        [(s, ps) in enumerate(S′)],
        z[all_controllers(ps), s] .== 1 .- a[all_controllers(ps)]
    )

    for edge in edges(g)
        α, β = edge.src, edge.dst
        @constraint(m, [s in eachindex(S′)], z[β, s] .≥ z[α, s] .- a[β])
        @constraint(m, [s in eachindex(S′)], z[α, s] .≥ z[β, s] .- a[α])
    end

    @constraint(m, [s in eachindex(S′)], F[s] == sum(z[:, s]))
end

function add_attack_cost_constraints!(m::Model, cost_config::AttackCostConfig)
    R = cost_config.R
    attack_cost = cost_config.cost

    V = length(cost_config.cost)
    a = m[:a]

    @constraint(m, sum(attack_cost[v] * a[v] for v = 1:V) ≤ R)
end

# Placements with one variable 
function add_uniqueness_constraint!(
    m::Model,
    V::Int,
    history::Vector{Placements},
    diff::Int,
    vp::AbstractArray{VariableRef},
)
    for sol in history
        Tp = sol.pc
        not_Tp = setdiff(1:V, Tp)
        @constraint(
            m,
            sum(vp[v] for v in Tp) - sum(vp[v] for v in not_Tp) ≤ length(Tp) - diff
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
    vb::AbstractArray{VariableRef},
)
    for sol in history
        Tp = sol.pc
        Tb = sol.bc
        not_Tp = setdiff(1:V, Tp)
        not_Tb = setdiff(1:V, Tb)

        @constraint(
            m,
            sum(vp[v] for v in Tp) - sum(vp[v] for v in not_Tp) + sum(vb[v] for v in Tb) -
            sum(vb[v] for v in not_Tb) ≤ length(Tp) + length(Tb) - diff
        )
    end
end

function get_controller_placements(model::Model; with_backups = False)
    assert_is_solved_and_feasible(model)
    s = model[:s]
    if with_backups
        r = model[:r]
        return Placements(to_indices(s), to_indices(r))
    end

    return Placements(to_indices(s))
end

macro uniqueness_constraints(m, V, history, diff, vars...)
    return quote
        add_uniqueness_constraint!(
            $(esc(m)),
            $(esc(V)),
            $(esc(history)),
            $(esc(diff)),
            $(map(esc, vars)...),
        )
    end
end


function get_U(BCC::Float64, dists::Matrix{Float64}; tol = 1e-9)
    V, _ = size(dists)
    [
        (v, w) for
        (v, w) in Base.Iterators.product(1:V, 1:V) if v ≤ w && dists[v, w] > BCC - tol
    ]
end

function get_Wv(BSC::Float64, dists::Matrix{Float64}; tol = 1e-9)
    V, _ = size(dists)
    W = Dict{Int,Vector{Int}}()

    for v = 1:V
        vicinity = [w for w = 1:V if dists[v, w] ≤ BSC + tol] # && v ≤ w]
        if !isempty(vicinity)
            W[v] = vicinity
        end
    end

    W
end

function add_controllers!(
    g::MetaGraph,
    m::Model,
    placement_config::PlacementConfig;
    controller_constraints::Union{ControllerConstraints,Nothing} = nothing,
)
    V = nv(g)
    M = placement_config.M

    # Primary controller placement
    @variable(m, s[1:V], Bin)

    if isnothing(controller_constraints)
        # Total number of controllers must equal M.
        @constraint(m, sum(s) == M)
    else
        P′ = controller_constraints.P′
        P″ = controller_constraints.P″
        B′ = controller_constraints.B′
        B″ = controller_constraints.B″
        # There are backup controllers.
        @variable(m, r[1:V], Bin)

        # Total number of controllers 
        @constraint(m, sum(s) + sum(r) == M)

        # Bounds on the primary and backup controllers.
        @constraint(m, P′ ≤ sum(s) ≤ P″)
        @constraint(m, B′ ≤ sum(r) ≤ B″)

        # Ensure that a controller is of only one type.
        # * s[v] == 1 : primary controller
        # * r[v] == 1 : backup controller
        # * s[v] + r[v] == 0 : switch 
        @constraint(m, s .+ r .≤ 1)
    end
end

function add_delay_constraints!(
    g::MetaGraph,
    m::Model,
    delay_constraints::DelayConstraintsConfig;
    capacity_constraints::Union{CapacityConstraintsConfig,Nothing} = nothing,
)
    dists = delay_constraints.distance_matrix
    BCC = delay_constraints.BCC
    BSC = delay_constraints.BSC

    V = nv(g)

    s = m[:s]
    # Controller-to-controller delays
    U = get_U(BCC, dists)
    @constraint(m, [(v, w) in U], s[v] + s[w] ≤ 1)

    W = get_Wv(BSC, dists)
    if !isnothing(capacity_constraints)
        # We use the information about capacity constraints 
        h = m[:h] # h[i, j] = switch i is assigned to controller j
        @constraint(m, [v in 1:V], sum(h[v, W[v]]) == 1)
    else
        @constraint(m, [v in 1:V], sum(s[W[v]]) ≥ 1)
    end
end

function add_capacity_constraints!(
    g::MetaGraph,
    m::Model,
    capacity_constraints::CapacityConstraintsConfig;
    delay_constraints::Union{DelayConstraintsConfig,Nothing} = nothing,
)
    V = nv(g)
    s = m[:s]

    # h[i, j] == 1 means that switch i is assigned to controller j.
    @variable(m, h[1:V, 1:V], Bin)

    arrivals = capacity_constraints.arrivals
    controller_capacity = capacity_constraints.controller_capacity
    λ = [arrivals[v] for v = 1:V]
    cc = [controller_capacity[v] for v = 1:V]
    for w = 1:V
        @constraint(m, λ ⋅ h[:, w] ≤ cc[w] * s[w])
    end

    # But if there are no delay constraints, simply ensure that 
    # a switch is assigned to only one controller.
    if isnothing(delay_constraints)
        @constraint(m, [v in 1:V], sum(h[v, :]) == 1)
    end
end

function add_survivability_constraints!(
    g::MetaGraph,
    m::Model,
    placement_config::PlacementConfig;
    controller_constraints::Union{ControllerConstraints,Nothing} = nothing,
)
    attackset = placement_config.attacks
    alen = length(attackset)

    # The number of surviving nodes given attack a 
    @variable(m, Y[1:alen] >= 0, Int)

    C_sets = Dict(
        a => components(attack_graph(g, attack)) for (a, attack) in enumerate(attackset)
    )

    # Component Survivability variable 
    # S=1 if component survives
    @variable(m, S[a = 1:alen, c = 1:length(C_sets[a])], Bin)

    # The set of primary controller nodes 
    s = m[:s]

    # Component c survives if there is another controller in the 
    # given component.
    for a = 1:alen
        for (c, comp) in enumerate(C_sets[a])
            vs = collect(labels(comp))
            if isnothing(controller_constraints)
                @constraint(m, S[a, c] .≤ sum(s[vs]))
            else
                # The set of backup controller nodes 
                r = m[:r]
                @constraint(m, S[a, c] .≤ sum(s[vs] + r[vs]))
            end
        end
    end

    # Count the number of surviving nodes.
    @constraint(
        m,
        [a in 1:alen],
        # length
        Y[a] == sum([
            let
                vs = collect(labels(cs))
                length(vs) * S[a, c]
            end for (c, cs) in enumerate(C_sets[a])
        ])
    )
end

function maximum_sc_delay(
    g::MetaGraph,
    controller_config::ControllerConstraints,
    delay_constraints::DelayConstraintsConfig;
    tol = 1e-9,
    optim = DEFAULT_OPTIM,
)
    V = nv(g)

    BCC = delay_constraints.BCC
    BSC = delay_constraints.BSC

    @assert BSC > 0.0
    @assert BCC > 0.0

    P′, P″ = controller_config.P′, controller_config.P″
    dists = delay_constraints.distance_matrix

    m = Model(optim)
    set_silent(m)

    # A1i
    @variable(m, D)

    # A1h: Controller placements
    @variable(m, s[1:V], Bin)

    # A1h: Controller assignments
    Wv = get_Wv(BSC, dists; tol = tol)
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
    U = get_U(BCC, dists; tol = tol)
    for (v, w) in U
        @constraint(m, s[v] + s[w] ≤ 1)
    end

    # A1d
    @constraint(m, [v in keys(Wv)], sum(z[v, w] for w in Wv[v]) == 1)

    # A1e
    @constraint(m, [v in keys(Wv), w in Wv[v]; w != v], z[v, w] ≤ s[w])

    # A1f
    @constraint(m, [v in keys(Wv); v in Wv[v]], z[v, v] == s[v])

    # A1g
    @constraint(m, [v in keys(Wv)], D ≥ sum(dists[v, w] * z[v, w] for w in Wv[v]))

    @solve_problem!(m)

    return objective_value(m)
end

### FEASIBILITY GENERATION 

# Primary controllers only, controllers capacity not considered
# max Z 
# s.t. sum {v in VERTICES} y[v] = M
#      {{v, w} in U} y[v] + y[w] ≤ 1
#      {v in VERTICES} sum {w in W(v)} y[w] ≥ 1
#      l=1,2,…,L, sum {v in T(l)} y[v] ≤ m 
#      {v in VERTICES} y[v] in B
# where T(l) denotes the set of nodes where controllers of the generated
# placement number lin1,2,…,L are located.

# For cost266 
# - P* = 3, tight CCD BCC = 1500, P* = 5 BCC = 2000 
function generate_controller_placement(
    g::MetaGraph,
    placement_config::PlacementConfig;
    controller_constraints::Union{ControllerConstraints,Nothing} = nothing,
    delay_constraints::Union{DelayConstraintsConfig,Nothing} = nothing,
    capacity_constraints::Union{CapacityConstraintsConfig,Nothing} = nothing,
    optim = DEFAULT_OPTIM,
    time_limit::Union{Nothing,Float64} = nothing,
)
    m = Model(optim)
    set_silent(m)

    if !isnothing(time_limit)
        set_time_limit_sec(m, time_limit)
    end

    m = Model(optim)
    set_silent(m)

    if !isnothing(time_limit)
        set_time_limit_sec(m, time_limit)
    end

    add_controllers!(g, m, placement_config; controller_constraints)

    if !isnothing(capacity_constraints)
        add_capacity_constraints!(g, m, capacity_constraints)
    end

    if !isnothing(delay_constraints)
        add_delay_constraints!(g, m, delay_constraints; capacity_constraints)
    end

    @objective(m, FEASIBILITY_SENSE, 0)

    time_taken = @solve_problem!(m)

    placement =
        get_controller_placements(m; with_backups = !isnothing(controller_constraints))
    (controllers = placement, time = time_taken, model = m)
end


### Paper:
### Finding Optimal Mixed-Strategies in a Matrix Game
### Between the Attacker and the Network Operator


## Master Problem 
function mixed_strategies_master(
    g::MetaGraph,
    placementset::Vector{Placements},
    attackset::Vector{Attacks};
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
    for i = 1:asize
        for j = 1:psize
            V_mat[i, j] = length(surviving_nodes(g, placementset[j], attackset[i]))
        end
    end

    @constraint(m, con[i = 1:asize], y ≤ sum(V_mat[i, j] * q[j] for j = 1:psize))
    time_taken = @solve_problem!(m)

    raw_duals = [dual(con[i]) for i = 1:asize]

    (
        time = time_taken,
        model = m,
        objective = objective_value(m),
        q_star = safe_probs(value.(q)),
        p_star = safe_probs(abs.(raw_duals)),
    )
end

function mixed_strategies_pricing_placement(
    g::MetaGraph,
    placement_config::PlacementConfig;
    controller_constraints::Union{ControllerConstraints,Nothing} = nothing,
    delay_constraints::Union{DelayConstraintsConfig,Nothing} = nothing,
    capacity_constraints::Union{CapacityConstraintsConfig,Nothing} = nothing,
    optim = DEFAULT_OPTIM,
    time_limit::Union{Nothing,Float64} = nothing,
)
    @smart_assert length(placement_config.attacks) == length(placement_config.p)

    m = Model(optim)
    p = placement_config.p
    set_silent(m)

    if !isnothing(time_limit)
        set_time_limit_sec(m, time_limit)
    end

    add_controllers!(g, m, placement_config; controller_constraints)

    if !isnothing(capacity_constraints)
        add_capacity_constraints!(g, m, capacity_constraints)
    end

    if !isnothing(delay_constraints)
        add_delay_constraints!(g, m, delay_constraints; capacity_constraints)
    end



    add_survivability_constraints!(g, m, placement_config; controller_constraints)

    Y = m[:Y]
    # Maximize the number of surviving controllers (weighted by probability)
    @objective(m, Max, p ⋅ Y)

    time_taken = @solve_problem!(m)

    placements =
        get_controller_placements(m; with_backups = !isnothing(controller_constraints))

    (time = time_taken, objective = objective_value(m), model = m, s = placements)
end

function mixed_strategies_pricing_attack(
    g::MetaGraph,
    attack_config::AttackConfig;
    cost_config::Union{AttackCostConfig,Nothing} = nothing,
    optim = DEFAULT_OPTIM,
    tol::Float64 = 1e-9,
)
    @smart_assert length(attack_config.placements) == length(attack_config.q)

    placementset = attack_config.placements
    q = attack_config.q

    # Reduce the problem size.
    active_idxs = [i for (i, val) in enumerate(q) if val > tol]
    S′ = placementset[active_idxs]
    q′ = q[active_idxs]

    m = Model(optim)
    set_silent(m)

    # The number of attacks to perform 
    add_attack_constraint!(g, m, attack_config)

    # Edge propagation constraints 
    add_edge_propagation_constraints!(g, m, attack_config)

    # Attack cost constraints
    if !isnothing(cost_config)
        add_attack_cost_constraints!(m, cost_config)
    end

    F = m[:F]
    a = m[:a]
    @objective(m, Min, q′ ⋅ F)

    time_taken = @solve_problem!(m)
    (time = time_taken, objective = objective_value(m), model = m, a = to_indices(a))
end
