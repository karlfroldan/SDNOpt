### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 9cb9a376-43c7-11f1-a9a7-593ccfdf5664
begin
	import Pkg;
	Pkg.activate(Base.current_project())
	Pkg.instantiate()

	cd(dirname(Base.current_project()))
	using SDNOpt
	using Plots
	using PlutoUI
	using Graphs
end

# ╔═╡ b940b140-ee8b-4b5f-b356-feba8dd658fa
md"""
This example uses the `cost266` network with $M=8$.

### Base Results
These results are when there are no additional constraints.
"""

# ╔═╡ a72ee2d6-bbdd-4b52-8403-e6455129197b
begin
	g = load_cost266()
	SDNOpt.plot_network(g; title="The cost266 network")
end

# ╔═╡ c46e9c70-22f7-4486-ab3a-d33cadf31dd8
@bind M Slider(3:20; default=8, show_value=true)

# ╔═╡ bc24e2f7-3cec-4522-a5f1-f8e734293044
begin
	base_placement_config = PlacementConfig(M)

	base_results = SDNOpt.MixedStrategyResult[]

	for k in [10, 11, 12]
		attack_config = AttackConfig(k)
		results = mixed_strategies_colgen(g, base_placement_config, attack_config)
		push!(base_results, results)
	end

	base_results
end

# ╔═╡ a16acee0-d807-4797-96f4-66134feba1fa


# ╔═╡ Cell order:
# ╠═9cb9a376-43c7-11f1-a9a7-593ccfdf5664
# ╠═b940b140-ee8b-4b5f-b356-feba8dd658fa
# ╠═a72ee2d6-bbdd-4b52-8403-e6455129197b
# ╠═c46e9c70-22f7-4486-ab3a-d33cadf31dd8
# ╠═bc24e2f7-3cec-4522-a5f1-f8e734293044
# ╠═a16acee0-d807-4797-96f4-66134feba1fa
