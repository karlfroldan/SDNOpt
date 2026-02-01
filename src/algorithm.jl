function mixed_strategies_colgen(
    g :: MetaGraph,
    M :: Int,
    K :: Int;
    optim = DEFAULT_OPTIM,
)
    V = nv(g)
    # Step 0 
    # Initialize list of placements and list of attacks
    placementset = Vector{Int}[randvec(V, M)]
    attackset = Vector{Int}[randvec(V, K)]
    res = mixed_strategies_master(g, placementset, attackset; optim=optim)

    @assert res != :infeasible "Initial solution is infeasible" 
    (obj, p_star, q_star) = res.objective, res.p_star, res.q_star
    placement_times = Float64[]
    attack_times = Float64[]

    has_changed = true

    xstars = Float64[]
    ystars = Float64[]

    while has_changed
        has_changed = false
        # Step 1
        # Solve the placement gen problem to get placement s′
        p_res = mixed_strategies_pricing_placement(
            g, M, attackset, p_star; optim=optim
        )

        # @show attackset
        # @show placementset
        @assert p_res != :infeasible "Pricing Placement is infeasible"
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
        res = mixed_strategies_master(
            g, placementset, attackset; optim=optim
        )

        @assert res != :infeasible "Master Problem is Infeasible"
        (obj, p_star, q_star) = res.objective, res.p_star, res.q_star

        a_res = mixed_strategies_pricing_attack(
            g, K, placementset, q_star; optim=optim,
        )

        @assert a_res != :infeasible "Pricing Attack is infeasible"
        a′ = a_res.a 

        push!(attack_times, a_res.time)
        push!(ystars, obj)

        if obj > [length(surviving_nodes(g, s, a′)) for s ∈ placementset] ⋅ q_star
            if !(a′ ∈ attackset)
                push!(attackset, a′)
                has_changed = true
            end
        end

        res = mixed_strategies_master(
            g, placementset, attackset; optim=optim
        )

        @assert res != :infeasible "Master Problem is Infeasible"
        (obj, p_star, q_star) = res.objective, res.p_star, res.q_star
        # println("")
    end

    (
        attackset = attackset,
        placementset = placementset,
        objective = obj,
        p_star = p_star,
        q_star = q_star,
        placement_times = placement_times,
        attack_times = attack_times,

        # Objective value for statistics 
        x_stars = xstars,
        y_stars = ystars,
    )
end