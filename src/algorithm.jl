function get_placement(placementlist :: Vector{Vector{Int}}, mixed_strategy :: Vector{Float64})
    @assert sum(mixed_strategy) ≈ 1.0 "Not a mixed Strategy!"
    @assert length(placementlist) == length(mixed_strategy) "Distribution does not equal the target set"

    dist = Categorical(mixed_strategy)
    idx = rand(dist)
    placementlist[idx]
end

function mixed_strategies_colgen(
    g :: MetaGraph,
    P :: Union{Int, IntBound},
    B :: Union{Int, IntBound},
    K :: Int;
    C :: Union{Nothing, Int} = nothing,
    optim = DEFAULT_OPTIM,
    BCC :: Float64 = 0.0,
    BSC :: Float64 = 0.0,
    control_capacity :: Dict{Int, Float64} = Dict{Int, Float64}(),
    control_demand :: Dict{Int, Float64} = Dict{Int, Float64}(),
    placement_list :: Vector{Placements} = Placements[],
    placement_difference :: Int = 1,
    dists :: Union{Matrix{Float64}, Nothing} = nothing,
    time_limit :: Float64 = (2 * 60 * 60) # 2 hours
)
    V = nv(g)
    # Step 0 
    # Initialize list of placements and list of attacks
    res = generate_controller_placement(
        g, P, B; C, optim, BCC, BSC, control_capacity,
        control_demand, placement_list, placement_difference,
        dists
    )
    placementset = Placements[res.controllers]
    attackset = Attacks[randvec(V, K)]

    update_master() = begin
        res = mixed_strategies_master(g, placementset, attackset; optim=optim)
        @assert res != :infeasible "Master Problem is Infeasible"
        return res.objective, res.p_star, res.q_star
    end

    obj, p_star, q_star = update_master()

    has_changed = true

    xstars = Float64[]
    ystars = Float64[]

    placement_times = Float64[]
    attack_times = Float64[]

    while has_changed
        has_changed = false
        # Step 1
        # Solve the placement gen problem to get placement s′
        p_res = mixed_strategies_pricing_placement_backup(
            g, P, B, attackset, p_star;
            optim, C, BCC, BSC, dists,
            time_limit,
        )

        time_limit -= p_res.time

        # @assert p_res != :infeasible "Pricing Placement is infeasible"
        s′ = p_res.s

        push!(placement_times, p_res.time)
        push!(xstars, obj)

        if [length(surviving_nodes(g, s′, a)) for a ∈ attackset] ⋅ p_star > obj
            if !(s′ ∈ placementset)
                push!(placementset, s′)
                has_changed = true
            end
        end

        # Step 2 Solve the attack generation problem
        obj, p_star, q_star = update_master()

        a_res = mixed_strategies_pricing_attack(
            g, K, placementset, q_star; optim=optim,
            time_limit
        )

        time_limit -= p_res.time

        # @assert a_res != :infeasible "Pricing Attack is infeasible"
        a′ = a_res.a 

        push!(attack_times, a_res.time)
        push!(ystars, obj)

        if obj > [length(surviving_nodes(g, s, a′)) for s ∈ placementset] ⋅ q_star
            if !(a′ ∈ attackset)
                push!(attackset, a′)
                has_changed = true
            end
        end

        obj, p_star, q_star = update_master()
    end

    (
        ;
        attackset,
        placementset,
        p_star,
        q_star,
        placement_times,
        attack_times,
        objective = obj,

        # Objective value for statistics 
        x_stars = xstars,
        y_stars = ystars,
    )
end
