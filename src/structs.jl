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
@kwdef struct PlacementConfig
    # Total number of placements
    M::Int
    # Current set of attacks 
    attacks::Vector{Attacks}
    # The probability distribution of the attacks 
    p::Vector{Float64}
end

# Base config for the attack problem for mixed-strategies.
@kwdef struct AttackConfig
    # Total Number of attacks 
    K::Union{Int,IntBound}
    # Current set of placements so far.
    placements::Vector{Placements}
end

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

@with_kw struct DelayConstraintsConfig
    # Controller-Controller Delay bound 
    BCC::Float64
    # Switch-Controller Delay bound 
    BSC::Float64
    # Distance matrix 
    distance_matrix::Matrix{Float64}
end

@with_kw struct CapacityConstraintsConfig
    # Controller Capacity Mapping
    controller_capacity::Dict{Int,Float64}
    arrivals::Dict{Int,Float64}

    @assert !isempty(controller_capacity) && !isempty(arrivals)
end
