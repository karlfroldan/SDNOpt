module SDNOpt

using Printf: @printf, @sprintf
using Graphs
using SimpleWeightedGraphs
using MetaGraphsNext
using Plots
using GraphRecipes
using Random: shuffle

using JuMP

# Optimizers
import SCIP
import HiGHS
import CPLEX 

# Exports from network.jl
export attack_graph, plot_network, components, distance_matrix
export load_dognet, load_cost266, load_coronet_conus, load_network

export μ2s

include("network.jl")
include("models.jl")

f() = "hello1"

end
