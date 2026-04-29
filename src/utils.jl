function Base.show(io::IO, p::Placements)
    if !isempty(p.bc) && isempty(p.pc)
        error("Primary controller list is empty but backup is not empty")
    end

    if !isempty(p.pc) && isempty(p.bc)
        print(io, p.pc)
    else
        print(io, "primaries = $(p.pc), backups = $(p.bc)")
    end
end

function safe_probs(p::AbstractArray; tol::Float64 = 1e-10)
    v = abs.(collect(p))
    v[v .< tol] .= 0.0

    # Normalize
    total_mass = sum(v)
    return v ./ total_mass
end

function μs2s(ms)
    return ms / 1e6
end

function randvec(k::Int, n::Int)
    if n > k
        error("n cannot be greater than k")
    end

    return sort(shuffle(1:k)[1:n])
end

function to_indices(controller_variable)
    return findall(Int.(round.(value.(controller_variable))) .== 1)
end

function all_controllers(p::Placements)
    return union(p.pc, p.bc)
end
