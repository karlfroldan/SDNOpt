experiment1() = run_experiment(
    load_cost266(),
    # Primary Controllers
    2, 10,
    # Backup Controllers 
    0, 2,
    # Attacks 
    2, 6,
    1500.0, # BCC 
    2000.0  # BSC
)

function experiment_pure_with_tight_delay_no_bc(g :: MetaGraph, BCC :: Float64; optim=DEFAULT_OPTIM)
    BSC_upper = BCC * BCC

    run_pure_strategies_experiment(
        g,
        # Primary controllers 
        3, 12,
        # Backup controllers 
        0, 0,
        # Attacks 
        2, 5,
        BCC,
        BSC_upper;
        optim=optim
    )
end

function experiment_with_tight_delay_no_bc(g :: MetaGraph, BCC :: Float64; optim=DEFAULT_OPTIM)
    BSC_upper :: Float64 = 3000.0

    run_experiment(
        g,
        # Primary Controllers 
        3, 12,
        # Backup controllers stay at 0
        0, 0,
        # Attacks 
        2, 5,
        BCC,
        BSC_upper;
        optim=optim
    )
end

function run_pure_strategies_experiment(
    g :: MetaGraph, 
    P_start :: Int,
    P_end :: Int, 
    B_start :: Int, 
    B_end :: Int, 
    K_start :: Int, 
    K_end :: Int,
    BCC :: Float64 = 0.0,
    BSC_upper :: Float64 = 0.0;
    optim=DEFAULT_OPTIM,
    tol = 1e-9
)
    dists = get_distance_matrix(g)
    
    # Create a unique timestamped directory for this batch of runs
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    outdir = "experiments/pure_results_$timestamp"
    mkpath(outdir)
    
    println("Starting pure strategies experiment...")
    println("Results will be saved to: $outdir/")

    # Keep a lightweight log of successes and failures
    run_log = []

    for P in P_start:P_end
        for B in B_start:B_end
            for K in K_start:K_end
                println("Running P = $P, B = $B, K = $K ...")
                
                try
                    C = P + B
                    # Tightest possible BSC 
                    BSC = maximum_sc_delay(g, P, dists, BSC_upper, BCC)

                    res_a1 = pure_controller_placement(
                        g, C, K;
                        BCC=BCC, BSC=BSC, optim=optim, tol=tol
                    )

                    res_a2 = pure_attack_generation(
                        g, C, K;
                        BCC=BCC, BSC=BSC, optim=optim, tol=tol
                    )

                    filename = joinpath(outdir, "run_P$(P)_B$(B)_K$(K).jld2")
                    
                    # Save metrics for both pure strategies
                    jldsave(filename; 
                        P=P, B=B, K=K, BCC=BCC, BSC=BSC, C=C,
                        # A1 outputs
                        a1_s_star = res_a1.s_star,
                        a1_Y_star = res_a1.Y_star,
                        a1_naop_time_ms = res_a1.naop_time_ms,
                        a1_cpop_time_ms = res_a1.cpop_time_ms,
                        a1_iterations = res_a1.iterations,
                        a1_attacks = res_a1.attacks,
                        # A2 outputs
                        a2_a_star = res_a2.a_star,
                        a2_Z_star = res_a2.Z_star,
                        a2_naop_time_ms = res_a2.naop_time_ms,
                        a2_cpop_time_ms = res_a2.cpop_time_ms,
                        a2_iterations = res_a2.iterations,
                        a2_placements = res_a2.placements
                    )
                    
                    push!(run_log, (P=P, B=B, K=K, C=C, BCC=BCC, BSC=BSC, 
                                    a1_Y_star=res_a1.Y_star, a2_Z_star=res_a2.Z_star, 
                                    status=:success, file=filename))
                    println("Saved to $filename")
                catch e
                    if e isa InfeasibleError
                        println("Infeasible P = $P, B = $B, K = $K")
                        push!(run_log, (P=P, B=B, K=K, C=C, BCC=BCC, BSC=BSC, 
                                        a1_Y_star=nothing, a2_Z_star=nothing, 
                                        status=:infeasible, error=string(e)))
                    elseif e isa TimeLimitError
                        println("Time Limit Exceeded P = $P, B = $B, K = $K")
                        push!(run_log, (P=P, B=B, K=K, C=C, BCC=BCC, BSC=BSC, 
                                        a1_Y_star=nothing, a2_Z_star=nothing, 
                                        status=:time_limit_exceeded, error=string(e)))
                    else
                        println("Other errors occurred")
                        push!(run_log, (P=P, B=B, K=K, C=C, BCC=BCC, BSC=BSC, 
                                        a1_Y_star=nothing, a2_Z_star=nothing, 
                                        status=:error, error=string(e)))
                    end
                    showerror(stdout, e)
                    println()
                end
            end
        end
    end
    
    # Save the summary log at the very end
    summary_file = joinpath(outdir, "experiment_summary.jld2")
    jldsave(summary_file; run_log=run_log)

    return run_log
end

function run_experiment(
    g :: MetaGraph, 
    P_start :: Int,
    P_end :: Int, 
    B_start :: Int, 
    B_end :: Int, 
    K_start :: Int, 
    K_end :: Int,
    BCC :: Float64 = 0.0,
    BSC_upper :: Float64 = 0.0;
    time_limit :: Float64 = 7200.0,
    optim=DEFAULT_OPTIM
)
    dists = get_distance_matrix(g)
    
    # Create a unique timestamped directory for this batch of runs
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    outdir = "experiments/results_$timestamp"
    mkpath(outdir)
    
    println("Starting overnight experiment...")
    println("Results will be saved to: $outdir/")

    # Keep a lightweight log of successes and failures
    run_log = []

    for P in P_start:P_end
        for B in B_start:B_end
            for K in K_start:K_end
                println("Running P = $P, B = $B, K = $K ...")
                
                try
                    C = P + B
                    # Tightest possible BSC 
                    BSC = maximum_sc_delay(g, P, dists, BSC_upper, BCC)
                    res = mixed_strategies_colgen(
                        g, P, B, K;
                        C=C, BCC=BCC, BSC=BSC, dists=dists,
                        optim=optim, time_limit=time_limit
                    )

                    filename = joinpath(outdir, "run_P$(P)_B$(B)_K$(K).jld2")
                    

                    jldsave(filename; 
                        P=P, B=B, K=K, BCC=BCC, BSC=BSC, C=C,
                        attackset = res.attackset,
                        placementset = res.placementset,
                        p_star = res.p_star,
                        q_star = res.q_star,
                        placement_times = res.placement_times,
                        attack_times = res.attack_times,
                        objective = res.objective,
                        x_stars = res.x_stars,
                        y_stars = res.y_stars
                    )

                    payoff=res.x_stars[end]
                    
                    push!(run_log, (P=P, B=B, K=K, C=C, BCC=BCC, BSC=BSC, payoff=payoff, status=:success, file=filename))
                    println("Saved to $filename")
                catch e
                    if e isa InfeasibleError
                        println("Infeasible P = $P, B = $B, K = $K")
                        push!(run_log, (P=P, B=B, K=K, C=C, BCC=BCC, BSC=BSC, payoff=nothing, status=:infeasible, error=string(e)))
                    elseif e isa TimeLimitError
                        println("Time Limit Exceeded P = $P, B = $B, K = $K")
                        push!(run_log, (P=P, B=B, K=K, C=C, BCC=BCC, BSC=BSC, payoff=nothing, status=:time_limit_exceeded, error=string(e)))
                    else
                        println("Other errors occurred")
                        push!(run_log, (P=P, B=B, K=K, C=C, BCC=BCC, BSC=BSC, payoff=nothing, status=:error, error=string(e)))
                    end
                    showerror(stdout, e)
                    println()
                end
            end
        end
    end
    
    # Save the summary log at the very end
    summary_file = joinpath(outdir, "experiment_summary.jld2")
    jldsave(summary_file; run_log=run_log)

    return run_log
end

function summary_to_dataframe(run_log)
    df = DataFrame()
    for entry in run_log
        push!(df, Dict(pairs(entry)), cols=:union)
    end
    
    return df
end