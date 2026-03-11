module SDNOpt

using Printf: @printf, @sprintf
using Graphs
using SimpleWeightedGraphs
using MetaGraphsNext
using Plots
using GraphRecipes
using Random: shuffle, rand
using LinearAlgebra
using JLD2
using Dates
using DataFrames

using Distributions

using Base.Iterators

using JuMP

# Optimizers
import SCIP
import HiGHS
import CPLEX

### Exports from network.jl
export attack_graph, plot_network, components, distance_matrix
export load_dognet, load_cost266, load_coronet_conus, load_network

export μs2s

### Exports from model.jl
# Paper 1: Max-min optimization of controller placements vs
#          min-max optimization of attacks on nodes in service networks.
export cpop, naop, pure_attack_generation, pure_controller_placement

export mixed_strategies_master, mixed_strategies_pricing_attack
export mixed_strategies_pricing_placement

### Exports from algorithm.jl 
export mixed_strategies_colgen

const IntBound = Tuple{Int, Int}


struct Placements
    pc :: Vector{Int}
    bc :: Vector{Int}

    Placements(pc :: Vector{Int}, bc :: Vector{Int}) = new(pc, bc)
    Placements(pc :: Vector{Int}) = new(pc, [])
end

const Attacks = Vector{Int}

DEFAULT_OPTIM = CPLEX.Optimizer

struct InfeasibleError <: Exception
    msg::String
end
Base.showerror(io::IO, e::InfeasibleError) = print(io, "InfeasibleError: ", e.msg)

struct TimeLimitError <: Exception
    msg::String
end
Base.showerror(io::IO, e::TimeLimitError) = print(io, "TimeLimitError: ", e.msg)

struct SolverError <: Exception
    msg::String
end
Base.showerror(io::IO, e::SolverError) = print(io, "SolverError: ", e.msg)

### Includes
include("macros.jl")
include("utils.jl")
include("network.jl")
include("models.jl")
include("algorithm.jl")
include("experimentation.jl")

end
