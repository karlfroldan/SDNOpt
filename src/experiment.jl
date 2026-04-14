# Method 1: Varying Primary Controllers (No B_list needed)
function run_experiments(
    g::MetaGraph,
    M::Int,
    P_list::Vector{Tuple{Int,Int}},
    K_list::Vector{Tuple{Int,Int}};
    dists = nothing,
    BCC::Float64 = 0.0,
    kwargs...,
)
    n_P = length(P_list)
    n_K = length(K_list)

    dists_matrix = isnothing(dists) ? SDNOpt.get_distance_matrix(g) : dists

    # Initialize 2D Arrays (P x K)
    objective_matrix = zeros(Float64, n_P, n_K)
    attackset_matrix = Matrix{Vector{Attacks}}(undef, n_P, n_K)
    placementset_matrix = Matrix{Vector{Placements}}(undef, n_P, n_K)
    p_star_matrix = Matrix{Vector{Float64}}(undef, n_P, n_K)
    q_star_matrix = Matrix{Vector{Float64}}(undef, n_P, n_K)
    
    master_times_matrix = Matrix{Vector{Float64}}(undef, n_P, n_K)
    placement_times_matrix = Matrix{Vector{Float64}}(undef, n_P, n_K)
    attack_times_matrix = Matrix{Vector{Float64}}(undef, n_P, n_K)
    
    x_stars_matrix = Matrix{Vector{Float64}}(undef, n_P, n_K)
    y_stars_matrix = Matrix{Vector{Float64}}(undef, n_P, n_K)

    master_results_matrix = Matrix{Vector{SubResult{MasterResult}}}(undef, n_P, n_K)
    placement_results_matrix = Matrix{Vector{SubResult{Placements}}}(undef, n_P, n_K)
    attack_results_matrix = Matrix{Vector{SubResult{Attacks}}}(undef, n_P, n_K)

    calculated_BSC_arr = zeros(Float64, n_P)
    placement_config = PlacementConfig(; M = M)

    for (i, P) in enumerate(P_list)
        # Recalculate minimum possible BSC bound for the current P
        println("Calculating minimum possible BSC bound for P=$(P), BCC=$(BCC)...")
        controller_constraints = ControllerConstraints(P[1], P[2], 0, 0)
        
        delay_config_init = DelayConstraintsConfig(;
            BCC = BCC,
            BSC = BCC^2,
            distance_matrix = dists_matrix,
        )

        calc_BSC = SDNOpt.maximum_sc_delay(g, controller_constraints, delay_config_init)
        calculated_BSC_arr[i] = calc_BSC
        println("Calculated BSC for P=$(P): $calc_BSC")

        delay_constraints = DelayConstraintsConfig(;
            BCC = BCC,
            BSC = calc_BSC,
            distance_matrix = dists_matrix,
        )

        for (j, K) in enumerate(K_list)
            println("Starting with P=$(P), B=(0, 0), K=$(K)")
            attack_config = AttackConfig(; K = K)

            res = mixed_strategies_colgen(
                g, placement_config, attack_config;
                controller_constraints = controller_constraints,
                delay_constraints = delay_constraints,
                kwargs...,
            )

            last_master = res.master[end]
            objective_matrix[i, j] = last_master.objective_value
            p_star_matrix[i, j] = last_master.results.p_star
            q_star_matrix[i, j] = last_master.results.q_star

            master_results_matrix[i, j] = res.master
            placement_results_matrix[i, j] = res.placement
            attack_results_matrix[i, j] = res.attack

            master_times_matrix[i, j] = solver_time(res, :master)
            placement_times_matrix[i, j] = solver_time(res, :placement)
            attack_times_matrix[i, j] = solver_time(res, :attack)

            attackset_matrix[i, j] = res.final_attacks
            placementset_matrix[i, j] = res.final_placements
            x_stars_matrix[i, j] = xstars(res)
            y_stars_matrix[i, j] = ystars(res)
        end
    end

    return (;
        objective = objective_matrix,
        attackset = attackset_matrix,
        placementset = placementset_matrix,
        p_star = p_star_matrix,
        q_star = q_star_matrix,
        master_times = master_times_matrix,
        placement_times = placement_times_matrix,
        attack_times = attack_times_matrix,
        x_stars = x_stars_matrix,
        y_stars = y_stars_matrix,
        master_results = master_results_matrix,
        placement_results = placement_results_matrix,
        attack_results = attack_results_matrix,
        calculated_BSC = calculated_BSC_arr,
    )
end


# Method 2: Varying Backup Controllers (Fixed P)
function run_experiments(
    g::MetaGraph,
    M::Int,
    P_fixed::Tuple{Int,Int},
    B_list::Vector{Tuple{Int,Int}},
    K_list::Vector{Tuple{Int,Int}};
    dists = nothing,
    BCC::Float64 = 0.0,
    kwargs...,
)
    n_B = length(B_list)
    n_K = length(K_list)

    dists_matrix = isnothing(dists) ? SDNOpt.get_distance_matrix(g) : dists

    println("Calculating minimum possible BSC bound for P=$(P_fixed), BCC=$(BCC)...")
    controller_config_init = ControllerConstraints(P_fixed[1], P_fixed[2], 0, 0)
    delay_config_init = DelayConstraintsConfig(;
        BCC = BCC,
        BSC = BCC^2,
        distance_matrix = dists_matrix,
    )

    calculated_BSC = SDNOpt.maximum_sc_delay(g, controller_config_init, delay_config_init)
    println("Calculated BSC: $calculated_BSC")

    delay_constraints = DelayConstraintsConfig(;
        BCC = BCC,
        BSC = calculated_BSC,
        distance_matrix = dists_matrix,
    )

    # Initialize 2D Arrays (B x K)
    objective_matrix = zeros(Float64, n_B, n_K)
    attackset_matrix = Matrix{Vector{Attacks}}(undef, n_B, n_K)
    placementset_matrix = Matrix{Vector{Placements}}(undef, n_B, n_K)
    p_star_matrix = Matrix{Vector{Float64}}(undef, n_B, n_K)
    q_star_matrix = Matrix{Vector{Float64}}(undef, n_B, n_K)
    
    master_times_matrix = Matrix{Vector{Float64}}(undef, n_B, n_K)
    placement_times_matrix = Matrix{Vector{Float64}}(undef, n_B, n_K)
    attack_times_matrix = Matrix{Vector{Float64}}(undef, n_B, n_K)
    
    x_stars_matrix = Matrix{Vector{Float64}}(undef, n_B, n_K)
    y_stars_matrix = Matrix{Vector{Float64}}(undef, n_B, n_K)

    master_results_matrix = Matrix{Vector{SubResult{MasterResult}}}(undef, n_B, n_K)
    placement_results_matrix = Matrix{Vector{SubResult{Placements}}}(undef, n_B, n_K)
    attack_results_matrix = Matrix{Vector{SubResult{Attacks}}}(undef, n_B, n_K)

    placement_config = PlacementConfig(; M = M)

    for (i, B) in enumerate(B_list)
        controller_constraints = ControllerConstraints(P_fixed[1], P_fixed[2], B[1], B[2])

        for (j, K) in enumerate(K_list)
            println("Starting with P=$(P_fixed), B=$(B), K=$(K)")
            attack_config = AttackConfig(; K = K)

            res = mixed_strategies_colgen(
                g, placement_config, attack_config;
                controller_constraints = controller_constraints,
                delay_constraints = delay_constraints,
                kwargs...,
            )

            last_master = res.master[end]
            objective_matrix[i, j] = last_master.objective_value
            p_star_matrix[i, j] = last_master.results.p_star
            q_star_matrix[i, j] = last_master.results.q_star

            master_results_matrix[i, j] = res.master
            placement_results_matrix[i, j] = res.placement
            attack_results_matrix[i, j] = res.attack

            master_times_matrix[i, j] = solver_time(res, :master)
            placement_times_matrix[i, j] = solver_time(res, :placement)
            attack_times_matrix[i, j] = solver_time(res, :attack)

            attackset_matrix[i, j] = res.final_attacks
            placementset_matrix[i, j] = res.final_placements
            x_stars_matrix[i, j] = xstars(res)
            y_stars_matrix[i, j] = ystars(res)
        end
    end

    return (;
        objective = objective_matrix,
        attackset = attackset_matrix,
        placementset = placementset_matrix,
        p_star = p_star_matrix,
        q_star = q_star_matrix,
        master_times = master_times_matrix,
        placement_times = placement_times_matrix,
        attack_times = attack_times_matrix,
        x_stars = x_stars_matrix,
        y_stars = y_stars_matrix,
        master_results = master_results_matrix,
        placement_results = placement_results_matrix,
        attack_results = attack_results_matrix,
        calculated_BSC = calculated_BSC,
    )
end