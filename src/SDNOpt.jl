module SDNOpt

using Printf: @printf, @sprintf
using Graphs
using SimpleWeightedGraphs
using MetaGraphsNext
using Plots
using GraphRecipes
using Random: shuffle
using LinearAlgebra


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

### Includes
include("network.jl")
include("models.jl")

g() = "hello"

end
