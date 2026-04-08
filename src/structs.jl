struct Placements
    pc::Vector{Int}
    bc::Vector{Int}

    Placements(pc::Vector{Int}, bc::Vector{Int}) = new(pc, bc)
    Placements(pc::Vector{Int}) = new(pc, Int[])
end

function ==(a::Placements, b::Placements)
    return a.pc == b.pc && a.bc == b.bc
end

function hash(a::Placements, h::UInt)
    return hash(a.pc, hash(a.bc, h))
end

const IntBound = Tuple{Int,Int}
const Attacks = Vector{Int}

# # Base Config for the placement problem for mixed-strategies 
@with_kw struct PlacementConfig
    # Total number of placements
    M::Int
    # Current set of attacks 
    attacks::Vector{Attacks} = Attacks[]
    # The probability distribution of the attacks 
    p::Vector{Float64} = Float64[]
end

PlacementConfig(M::Int) = PlacementConfig(; M=M)

# Base config for the attack problem for mixed-strategies.
@with_kw struct AttackConfig
    # Total Number of attacks 
    K::Union{Int,IntBound}
    # Current set of placements so far.
    placements::Vector{Placements} = Placements[]
    # The probability distribution of the placements 
    q::Vector{Float64} = Float64[]

    # AttackConfig(K::Union{Int,IntBound}, placements::Vector{Placements}, q::Vector{Float64}) =
    #     AttackConfig(; K=K, placements=placements, q=q)
    # AttackConfig(K::Union{Int,IntBound}) = AttackConfig(K, Placements[], Float64[])
end

AttackConfig(K::Union{Int,IntBound}) = AttackConfig(; K=K)

# Backup controller config.
@with_kw struct ControllerConstraints
    # Bounds for the primary controller placement 
    # and backup controller placement 
    P′::Int
    P″::Int
    B′::Int
    B″::Int

    @assert P′ ≤ P″ "Primary controller lower bound should be less than or equal to the upper bound"
    @assert B′ ≤ B″ "Backup controller lower bound should be less tha or equal to the upper bound"
end

@kwdef struct DelayConstraintsConfig
    # Controller-Controller Delay bound 
    BCC::Float64
    # Switch-Controller Delay bound 
    BSC::Float64
    # Distance matrix 
    distance_matrix::Matrix{Float64}

    function DelayConstraintsConfig(g::AbstractGraph, BCC::Float64, BSC::Float64)
        distance_matrix = get_distance_matrix(g)
        return DelayConstraintsConfig(BCC, BSC, distance_matrix)
    end
end

@with_kw struct CapacityConstraintsConfig
    # Controller Capacity Mapping
    controller_capacity::Dict{Int,Float64}
    arrivals::Dict{Int,Float64}

    @assert !isempty(controller_capacity) && !isempty(arrivals)
end

@kwdef struct AttackCostConfig
    R :: Float64 # Total attacker's budget.
    cost::Dict{Int,Float64}
end