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

PlacementConfig(M::Int) = PlacementConfig(; M = M)

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

AttackConfig(K::Union{Int,IntBound}) = AttackConfig(; K = K)

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
end

function DelayConstraintsConfig(g::AbstractGraph, BCC::Float64, BSC::Float64)
    distance_matrix = get_distance_matrix(g)
    return DelayConstraintsConfig(BCC, BSC, distance_matrix)
end

@with_kw struct CapacityConstraintsConfig
    # Controller Capacity Mapping
    controller_capacity::Dict{Int,Float64}
    arrivals::Dict{Int,Float64}

    @assert !isempty(controller_capacity) && !isempty(arrivals)
end

@kwdef struct AttackCostConfig
    R::Float64 # Total attacker's budget.
    cost::Dict{Int,Float64}
end

struct MasterResult
    # Defender probability distribution
    q_star::Vector{Float64}
    # Attacker probability distribution
    p_star::Vector{Float64}
end

"""
The results of the pricing problem.
"""
mutable struct SubResult{T}
    time::Float64
    objective_value::Float64
    results::T
    # Whether this is added to the final list of generated controller
    # or attack placements.
end

function Base.show(io::IO, obj::SubResult{T}) where {T<:Any}
    res = obj.results
    obj_val = obj.objective_value
    time = obj.time

    return print(
        io,
        "Found: $res with objective value $obj_val with time $time",
    )
end

function Base.show(io::IO, obj::SubResult{MasterResult})
    obj_val = obj.objective_value
    time = obj.time

    return print(io, "Expected Payoff: $obj_val with time $time")
end

struct MixedStrategyResult
    master::Vector{SubResult{MasterResult}}
    placement::Vector{SubResult{Placements}}
    attack::Vector{SubResult{Attacks}}

    final_placements::Vector{Placements}
    final_attacks::Vector{Attacks}
end

function solver_time(res::MixedStrategyResult, problem::Symbol)
    return map(x -> x.time, getproperty(res, problem))
end

"""
Primal solution objective function. (Placement)
"""
function xstars(res::MixedStrategyResult)
    return [res.master[2i-1].objective_value for i in 1:length(res.placement)]
end

"""
Dual solution objective function. (Attacker)
"""
function ystars(res::MixedStrategyResult)
    return [res.master[2i].objective_value for i in 1:length(res.attack)]
end

function Base.show(io::IO, obj::MixedStrategyResult)
    n_placements = length(obj.final_placements)
    n_attacks = length(obj.final_attacks)

    master_time = sum(solver_time(obj, :master))
    placement_time = sum(solver_time(obj, :placement))
    attack_time = sum(solver_time(obj, :attack))

    println(io, "--- Solved and feasible ---")
    println(io, "n_placements: $n_placements")
    println(io, "n_attacks: $n_attacks")
    println(io, "time to solve master: $master_time")
    println(io, "time to generate placements: $placement_time")
    return println(io, "time to generate attacks: $attack_time")
end
