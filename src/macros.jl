macro elapsed_time(expr)
    return quote
        local start_time = time_ns()
        local result = $(esc(expr))
        local end_time = time_ns()
        local elapsed = (end_time - start_time) / 1e3
        
        if result === nothing
            elapsed
        else
            (elapsed, result)
        end
    end
end

macro unpack_bounds(var)
    # Generate the assertion message at parse time
    msg = "$(var)′ should be less than or equal to $(var)″"
    
    return quote
        let val = $(esc(var))
            tup = val isa Tuple ? val : (val, val)
            @assert tup[1] ≤ tup[2] $msg
            tup
        end
    end
end

macro solve_problem!(m)
    return quote
        local t = @elapsed optimize!($(esc(m)))
        if !is_solved_and_feasible($(esc(m)))
            error("Optimization failed: Model is infeasible or unsolved.")
        end
        t
    end
end