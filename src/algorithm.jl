const TIME_LIMIT = 2.0 * 60.0 * 60.0 # 2 Hours

function get_placement(placementlist::Vector{Vector{Int}}, mixed_strategy::Vector{Float64})
    @assert sum(mixed_strategy) ≈ 1.0 "Not a mixed Strategy!"
    @assert length(placementlist) == length(mixed_strategy) "Distribution does not equal the target set"

    dist = Categorical(mixed_strategy)
    idx = rand(dist)
    placementlist[idx]
end

# A1: Algorithm for controller placement optimization by means of
# attack generation 
function pure_controller_placement(
    g::MetaGraph,
    M::Int,
    K::Int;
    BCC::Union{Nothing,Float64} = nothing,
    BSC::Union{Nothing,Float64} = nothing,
    optim = DEFAULT_OPTIM,
    tol = 1e-9,
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
        res = naop(g, K, [s_star]; optim = optim)
        Z_star = res.objective_value

        naop_time_ms += res.time
        local_naop_time = res.time
        if Z_star ≥ Y_star - tol
            break
        end

        push!(attacks, res.attack)

        # Step 2: Solve CPOP to get better placement.
        res = cpop(g, M, attacks; optim = optim, BCC = BCC, BSC = BSC)
        Y_star = res.objective_value
        s_star = res.controllers

        cpop_time_ms += res.time
        local_cpop_time = res.time

        count += 1
        println(
            "Z_star: $Z_star, Y_star: $Y_star, naop_time: $local_naop_time, cpop_time: $local_cpop_time",
        )
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
    g::MetaGraph,
    M::Int,
    K::Int;
    BCC::Union{Nothing,Float64} = nothing,
    BSC::Union{Nothing,Float64} = nothing,
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
        res = cpop(g, M, [a_star]; optim = optim, BCC, BSC)
        Y_star = res.objective_value

        cpop_time_ms += res.time
        local_cpop_time = res.time
        if Y_star ≤ Z_star + tol
            break
        end

        push!(placements, res.controllers)

        # Step 2
        res = naop(g, K, placements; optim = optim)
        Z_star = res.objective_value
        a_star = res.attack

        naop_time_ms += res.time
        local_naop_time = res.time

        println(
            "Z_star: $Z_star, Y_star: $Y_star, naop_time: $local_naop_time, cpop_time: $local_cpop_time",
        )

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

function mixed_strategies_colgen(
    g::MetaGraph,
    placement_config::PlacementConfig,
    attack_config::AttackConfig;
    controller_constraints::Union{ControllerConstraints,Nothing} = nothing,
    delay_constraints::Union{DelayConstraintsConfig,Nothing} = nothing,
    capacity_constraints::Union{CapacityConstraintsConfig,Nothing} = nothing,
    cost_config::Union{AttackCostConfig,Nothing} = nothing,
    optim = DEFAULT_OPTIM,
    time_limit::Float64 = TIME_LIMIT,
)
    V = nv(g)
    _, K″ = @unpack_bounds attack_config.K

    res = generate_controller_placement(
        g,
        placement_config;
        controller_constraints,
        delay_constraints,
        capacity_constraints,
        optim,
    )

    placementset = Placements[res.controllers]
    attackset = Attacks[randvec(V, K″)]

    update_master() = begin
        res = mixed_strategies_master(g, placementset, attackset; optim)
        @assert res != :infeasible "Master problem is infeasible"
        (obj = res.objective, p_star = res.p_star, q_star = res.q_star)
    end

    has_changed = true
    obj, p_star, q_star = update_master()

    # Primal solution (controller's master problem)
    x_stars = Float64[]
    # Dual solution (attacker's master problem)
    y_stars = Float64[]

    # We also store the timestamps.
    master_times = Float64[]
    placement_times = Float64[]
    attack_times = Float64[]

    while has_changed
        has_changed = false
        # Step 1
        # Sovle the placement generation problem to get
        # placement s′
        placement_config = PlacementConfig(placement_config.M, attackset, p_star)
        p_res = mixed_strategies_pricing_placement(
            g,
            placement_config;
            controller_constraints,
            delay_constraints,
            capacity_constraints,
            optim,
        )

        s′ = p_res.s
        push!(placement_times, p_res.time)
        push!(x_stars, obj)

        if [length(surviving_nodes(g, s′, a)) for a in attackset] ⋅ p_star > obj
            if !(s′ in placementset)
                # We found a better placement.
                push!(placementset, s′)
                has_changed = true
            end
        end

        # Step 2, solve the attack generation problem 
        obj, p_star, q_star = update_master()

        # Generate a new attack.
        attack_config = AttackConfig(attack_config.K, placementset, q_star)
        a_res = mixed_strategies_pricing_attack(g, attack_config; cost_config)

        a′ = a_res.a # The new attack 
        push!(attack_times, a_res.time)
        push!(y_stars, obj)

        if obj > [length(surviving_nodes(g, s, a′)) for s in placementset] ⋅ q_star
            if !(a′ in attackset)
                push!(attackset, a′)
                has_changed = true
            end
        end

        obj, p_star, q_star = update_master()
    end

    (;
        attackset,
        placementset,
        p_star,
        q_star,
        placement_times,
        attack_times,
        objective = obj,

        # Objective value for statistics 
        x_stars = x_stars,
        y_stars = y_stars,
    )
end
