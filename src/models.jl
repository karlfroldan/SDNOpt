DEFAULT_OPTIM = SCIP.Optimizer

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

function μ2s(ms)
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
        if Z_star ≥ Y_star
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
        if Y_star ≤ Z_star - tol
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