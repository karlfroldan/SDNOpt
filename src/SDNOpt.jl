module SDNOpt

using Base: @kwdef
using Printf: @printf, @sprintf
using Graphs
using SimpleWeightedGraphs
using MetaGraphsNext
using Plots
using Parameters: @with_kw
using GraphRecipes
using Random: shuffle, rand
using LinearAlgebra
using JLD2
using Dates
using DataFrames
using SmartAsserts: @smart_assert

using Distributions

using Base.Iterators

using JuMP

# Optimizers
import CPLEX

### Exports from network.jl
export attack_graph, plot_network, components, distance_matrix
export load_dognet, load_cost266, load_coronet_conus, load_network
export shortest_path_tree
export calculate_arrivals, calculate_attack_costs, calculate_attacker_budget
export calculate_capacity_dict

export μs2s

### Exports from model.jl
# Paper 1: Max-min optimization of controller placements vs
#          min-max optimization of attacks on nodes in service networks.
export cpop, naop, pure_attack_generation, pure_controller_placement

export mixed_strategies_master, mixed_strategies_pricing_attack
export mixed_strategies_pricing_placement
export mixed_strategies_pricing_attack_benders

### Exports from algorithm.jl 
export mixed_strategies_colgen

### Exports from network_plot.jl
export generate_tikz, generate_heatmap_tikz, calculate_node_probabilities

### Exports from structs.jl
export Placements, Attacks, PlacementConfig, AttackConfig, AttackCostConfig
export ControllerConstraints, DelayConstraintsConfig, CapacityConstraintsConfig

import Base: ==, hash

### Exports from experiment.jl
export run_experiments

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
include("structs.jl")
include("utils.jl")
include("network.jl")
include("network_plot.jl")
include("models.jl")
include("algorithm.jl")
include("experiment.jl")

end
