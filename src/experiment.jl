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

    # Calculate distance matrix if not provided
    dists_matrix = isnothing(dists) ? SDNOpt.get_distance_matrix(g) : dists

    # Calculate the minimum BSC
    println("Calculating mininum possible BSC bound for BCC = $BCC...")
    controller_config_init = ControllerConstraints(P_fixed[1], P_fixed[2], 0, 0)
    delay_config_init = DelayConstraintsConfig(;
        BCC = BCC,
        BSC = BCC^2,
        distance_matrix = dists_matrix,
    )

    calculated_BSC =
        SDNOpt.maximum_sc_delay(g, controller_config_init, delay_config_init)
    println("Calculated BSC: $calculated_BSC")

    # Define the final delay configuration with the calculated BSC
    delay_constraints = DelayConstraintsConfig(;
        BCC = BCC,
        BSC = calculated_BSC,
        distance_matrix = dists_matrix,
    )

    # The only matrix of scalars
    objective_matrix = zeros(Float64, n_B, n_K)

    # Matrices of Vectors
    attackset_matrix = Matrix{Vector{Attacks}}(undef, n_B, n_K)
    placementset_matrix = Matrix{Vector{Placements}}(undef, n_B, n_K)
    p_star_matrix = Matrix{Vector{Float64}}(undef, n_B, n_K)
    q_star_matrix = Matrix{Vector{Float64}}(undef, n_B, n_K)
    placement_times_matrix = Matrix{Vector{Float64}}(undef, n_B, n_K)
    attack_times_matrix = Matrix{Vector{Float64}}(undef, n_B, n_K)
    x_stars_matrix = Matrix{Vector{Float64}}(undef, n_B, n_K)
    y_stars_matrix = Matrix{Vector{Float64}}(undef, n_B, n_K)

    # Reusable placement configuration for the total number of placements M
    placement_config = PlacementConfig(; M = M)

    for (i, B) in enumerate(B_list)
        # Update controller constraints for the current backup budget B
        controller_constraints =
            ControllerConstraints(P_fixed[1], P_fixed[2], B[1], B[2])

        for (j, K) in enumerate(K_list)
            println("Starting with B=$(B) and K=$(K)")

            attack_config = AttackConfig(; K = K)

            res = mixed_strategies_colgen(
                g,
                placement_config,
                attack_config;
                controller_constraints = controller_constraints,
                delay_constraints = delay_constraints,
                kwargs...,
            )

            # Store the scalar objective
            objective_matrix[i, j] = res.objective

            # Store the vectors directly into the matrices
            attackset_matrix[i, j] = res.attackset
            placementset_matrix[i, j] = res.placementset
            p_star_matrix[i, j] = res.p_star
            q_star_matrix[i, j] = res.q_star
            placement_times_matrix[i, j] = res.placement_times
            attack_times_matrix[i, j] = res.attack_times
            x_stars_matrix[i, j] = res.x_stars
            y_stars_matrix[i, j] = res.y_stars
        end
    end

    return (;
        objective = objective_matrix,
        attackset = attackset_matrix,
        placementset = placementset_matrix,
        p_star = p_star_matrix,
        q_star = q_star_matrix,
        placement_times = placement_times_matrix,
        attack_times = attack_times_matrix,
        x_stars = x_stars_matrix,
        y_stars = y_stars_matrix,
        calculated_BSC = calculated_BSC,
    )
end
