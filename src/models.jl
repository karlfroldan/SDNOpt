DEFAULT_OPTIM = SCIP.Optimizer

struct Placement
    pc :: Vector{Int}
    bc :: Vector{Int}

    Placement(pc :: Vector{Int}, bc :: Vector{Int}) = new(pc, bc)
    Placement(pc :: Vector{Int}) = new(pc, [])
end

function Base.show(io::IO, p::Placement)
    if !isempty(p.bc) && isempty(p.pc)
        error("Primary controller list is empty but backup is not empty")
    end

    if !isempty(p.pc) && isempty(p.bc)
        print(io, p.pc)
    else
        print(io, "primaries = $(p.pc), backups = $(p.bc)")
    end
end

macro elapsed_time(expr)
    return quote
        local start_time = time_ns()
        local result = $(esc(expr))
        local end_time = time_ns()
        local elapsed = (end_time - start_time) / 1e3
        
        if result === nothing
            elapsed
        else
            (elapsed, result)
        end
    end
end

function μs2s(ms)
    ms / 1e6
end

function randvec(k::Int, n::Int)
    if n > k
        error("n cannot be greater than k")
    end

    return sort(shuffle(1:k)[1:n])
end

function to_indices(controller_variable)
    findall(Int.(round.(value.(controller_variable))) .== 1)
end

## Paper: Max-Min Optimization of Controller Placements
##        vs Min-Max Optimization of Attacks on Nodes in Service Networks

# Formulation 3.1: Max-min Controller Placement optimization problem (CPOP)
function cpop(
    g :: MetaGraph, 
    M :: Int, 
    # The nodes that we want to attack.
    attacks :: Vector{Vector{Int}};
    optim = DEFAULT_OPTIM,
)
    # Find an M-node controller placement given the considered set of attacks
    model = Model(optim)
    set_silent(model)

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

    time_taken = @elapsed_time optimize!(model)
    if ! is_solved_and_feasible(model)
        return :infeasible 
    end
    
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

    time_taken = @elapsed_time optimize!(model)
    if ! is_solved_and_feasible(model)
        return :infeasible 
    end
    
    (
        attack = to_indices(a),
        objective_value = objective_value(model),
        model = model,
        time = time_taken
    )
end

# A1: Algorithm for controller placement optimization by means of
# attack generation 
function pure_controller_placement(
    g :: MetaGraph,
    M :: Int,
    K :: Int;
    optim = DEFAULT_OPTIM,
    tol = 1e-9
)
    V = nv(g)

    # Step 0: Generate a random M-Node controller placement s*
    s_star = randvec(V, M)
    attacks = Vector{Int}[]
    Y_star = Float64(V)

    cpop_time_ms = 0.0
    naop_time_ms = 0.0

    # Step 1: Solve NAOP to get the worst attack given the random placement
    #         This assures Z* surviving nodes.
    count = 0
    while true
        res = naop(g, K, [s_star]; optim=optim)
        @assert res != :infeasible "NAOP is infeasible"
        Z_star = res.objective_value

        naop_time_ms += res.time
        if Z_star ≥ Y_star - tol
            break
        end

        push!(attacks, res.attack)

        # Step 2: Solve CPOP to get better placement.
        res = cpop(g, M, attacks; optim=optim)
        @assert res != :infeasible "CPOP is infeasible"
        Y_star = res.objective_value
        s_star = res.controllers

        cpop_time_ms += res.time

        count += 1
    end

    (
        s_star = s_star,
        Y_star = Y_star,
        naop_time_ms = naop_time_ms,
        cpop_time_ms = cpop_time_ms,
        iterations = count,
        attacks = attacks,
    )
end


# A2: Pure Strategy Attack Generation by means of controller placement optimization
function pure_attack_generation(
    g :: MetaGraph,
    M :: Int,
    K :: Int;
    optim = DEFAULT_OPTIM,
    tol = 1e-9,
)
    V = nv(g)

    # Step 0
    a_star = randvec(V, K)
    placements = Vector{Int}[]
    Z_star = 0

    cpop_time_ms = 0.0
    naop_time_ms = 0.0

    # Step 1
    count = 0
    while true
        res = cpop(g, M, [a_star]; optim=optim)
        @assert res != :infeasible "CPOP is infeasible"
        Y_star = res.objective_value

        cpop_time_ms += res.time
        if Y_star ≤ Z_star + tol
            break
        end

        push!(placements, res.controllers)

        # Step 2
        res = naop(g, K, placements; optim=optim)
        @assert res != :infeasible "NAOP is infeasible"
        Z_star = res.objective_value
        a_star = res.attack

        naop_time_ms += res.time

        count += 1
    end

    (
        a_star = a_star,
        Z_star = Z_star,
        naop_time_ms = naop_time_ms,
        cpop_time_ms = cpop_time_ms,
        iterations = count,
        placements = placements,
    )
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
    M :: Int, 
    P :: Union{Int, Tuple{Int, Int}},
    B :: Union{Int, Tuple{Int, Int}};
    optim = DEFAULT_OPTIM,
    BCC :: Float64 = 0.0,
    BSC :: Float64 = 0.0,
    control_capacity :: Dict{Int, Float64} = Dict{Int, Float64}(),
    control_demand :: Dict{Int, Float64} = Dict{Int, Float64}(),
    placement_list :: Vector{Placement} = Placement[],
    placement_difference :: Int = 1,
    dists :: Union{Matrix{Float64}, Nothing} = nothing,
    tol = 1e-9
)
    m = Model(optim)
    dists = isnothing(dists) ? get_distance_matrix(g) : dists
    V = nv(g)
    local x
    has_capacity = false

    if !isempty(control_capacity) || !isempty(control_demand)
        @assert (!isempty(control_capacity) && !isempty(control_demand)) "Both control capacity and control demand should not be empty"
        has_capacity = true
    end
    
    #  # if P isa Tuple

    P′, P″ = P isa Tuple ? P : (P, P)
    B′, B″ = B isa Tuple ? B : (B, B)

    @assert P′ ≤ P″ "P′ should be less than or equal to P″"
    @assert B′ ≤ B″ "B′ should be less than or equal to B″"

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
    @constraint(m, sum(x) + sum(y) == M)

    # (2c) Primary controllers constraint 
    @constraint(m, P′ ≤ sum(y) ≤ P″)

    # (2d) Backup controllers constraint
    @constraint(m, B′ ≤ sum(x) ≤ B″)

    if BCC > 0.0
        # println("(2e) activated")
        U = get_U(g, BCC; dists=dists)

        # (2e) Controllers must be able to reach each other inside BCC delay.
        for (v, w) ∈ U
            @constraint(m, y[v] + y[w] ≤ 1)
        end
    end

    if BSC > 0.0
        # println("(2f) activated")
        W = get_Wv(g, BSC; dists=dists)
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
    for s ∈ placement_list 
        @constraint(m, sum(y[s.pc]) + sum(x[s.pc]) ≤ placement_difference)
    end

    time_taken = @elapsed_time optimize!(m)
    if ! is_solved_and_feasible(m)
        return :infeasible 
    end

    placement = Placement(to_indices(y), to_indices(x))
    
    (
        controllers = placement,
        time = time_taken,
        model = m,
        objective_value = Z,
    )
end

function get_U(
    g :: MetaGraph, 
    BCC :: Float64; 
    dists :: Union{Matrix{Float64}, Nothing},
    tol = 1e-9
)
    V = nv(g)
    dists = isnothing(dists) ? get_distance_matrix(g) : dists
    [
        (v, w) for (v, w) ∈ Base.Iterators.product(1:V, 1:V)
        if v ≤ w && dists[v, w] > BCC - tol
    ]
end

function get_Wv(
    g :: MetaGraph, 
    BSC :: Float64; 
    dists :: Union{Matrix{Float64}, Nothing},
    tol = 1e-9
)
    V = nv(g)
    W = Dict{Int, Vector{Int}}()
    dists = isnothing(dists) ? get_distance_matrix(g) : dists

    for v ∈ 1:V
        vicinity = [w for w ∈ 1:V if dists[v, w] ≤ BSC + tol] # && v ≤ w]
        if !isempty(vicinity)
            W[v] = vicinity
        end
    end

    W
end