const TIME_LIMIT = 2.0 * 60.0 * 60.0 # 2 Hours

function get_placement(placementlist::Vector{Vector{Int}}, mixed_strategy::Vector{Float64})
    @assert sum(mixed_strategy) ≈ 1.0 "Not a mixed Strategy!"
    @assert length(placementlist) == length(mixed_strategy) "Distribution does not equal the target set"

    dist = Categorical(mixed_strategy)
    idx = rand(dist)
    placementlist[idx]
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
        mixed_strategies_master(g, placementset, attackset; optim)
    end

    has_changed = true

    # # Primal solution (controller's master problem)
    # x_stars = Float64[]
    # # Dual solution (attacker's master problem)
    # y_stars = Float64[]
    placement_results = SubResult{Placements}[]
    attack_results = SubResult{Attacks}[]
    master_results = SubResult{MasterResult}[]

    while has_changed
        has_changed = false
        master_res = update_master()
        p_star = master_res.results.p_star
        q_star = master_res.results.q_star
        obj = master_res.objective_value
        # Step 1
        # Sovle the placement generation problem to get
        # placement s′
        placement_config =
            PlacementConfig(placement_config.M, attackset, master_res.results.p_star)

        p_res = mixed_strategies_pricing_placement(
            g,
            placement_config;
            controller_constraints,
            delay_constraints,
            capacity_constraints,
            optim,
        )

        s′ = p_res.results
        if [length(surviving_nodes(g, s′, a)) for a in attackset] ⋅ p_star > obj
            if !(s′ in placementset)
                # We found a better placement.
                push!(placementset, s′)
                has_changed = true
            end
        end

        push!(placement_results, p_res)

        # Step 2, solve the attack generation problem 
        # obj, p_star, q_star, mt1 = update_master()
        master_res = update_master()
        push!(master_results, master_res)
        obj = master_res.objective_value
        p_star = master_res.results.p_star
        q_star = master_res.results.q_star

        # Generate a new attack.
        attack_config = AttackConfig(attack_config.K, placementset, q_star)
        a_res = mixed_strategies_pricing_attack(g, attack_config; cost_config)

        a′ = a_res.results # The new attack 
        if obj > [length(surviving_nodes(g, s, a′)) for s in placementset] ⋅ q_star
            if !(a′ in attackset)
                push!(attackset, a′)
                has_changed = true
            end
        end

        push!(attack_results, a_res)

        master_res = update_master()
        push!(master_results, master_res)
        p_star = master_res.results.p_star
        q_star = master_res.results.q_star
        obj = master_res.objective_value
    end

    MixedStrategyResult(
        master_results,
        placement_results,
        attack_results,
        placementset,
        attackset,
    )
end
