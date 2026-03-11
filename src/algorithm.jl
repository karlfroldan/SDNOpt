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
    g :: MetaGraph,
    M :: Int,
    K :: Int;
    BCC :: Union{Nothing, Float64} = nothing,
    BSC :: Union{Nothing, Float64} = nothing,
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
        Z_star = res.objective_value

        naop_time_ms += res.time
        if Z_star ≥ Y_star - tol
            break
        end

        push!(attacks, res.attack)

        # Step 2: Solve CPOP to get better placement.
        res = cpop(g, M, attacks; optim=optim, BCC=BCC, BSC=BSC)
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
    BCC :: Union{Nothing, Float64} = nothing,
    BSC :: Union{Nothing, Float64} = nothing,
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
        res = cpop(g, M, [a_star]; optim=optim, BCC, BSC)
        Y_star = res.objective_value

        cpop_time_ms += res.time
        if Y_star ≤ Z_star + tol
            break
        end

        push!(placements, res.controllers)

        # Step 2
        res = naop(g, K, placements; optim=optim)
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

function mixed_strategies_colgen(
    g::MetaGraph,
    P::Union{Int,IntBound},
    B::Union{Int,IntBound},
    K::Int;
    C::Union{Nothing,Int}=nothing,
    optim=DEFAULT_OPTIM,
    BCC::Float64=0.0,
    BSC::Float64=0.0,
    control_capacity::Dict{Int,Float64}=Dict{Int,Float64}(),
    control_demand::Dict{Int,Float64}=Dict{Int,Float64}(),
    placement_list::Vector{Placements}=Placements[],
    placement_difference::Int=1,
    dists::Union{Matrix{Float64},Nothing}=nothing,
    time_limit::Float64=TIME_LIMIT
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
        objective=obj,

        # Objective value for statistics 
        x_stars=xstars,
        y_stars=ystars,
    )
end
