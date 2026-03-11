experiment1() = run_experiment(
    load_cost266(),
    # Primary Controllers
    2, 10,
    # Backup Controllers 
    0, 2,
    # Attacks 
    2, 5,
    1500.0, # BCC 
    2000.0  # BSC
)

function run_experiment(
    g :: MetaGraph, 
    P_start :: Int,
    P_end :: Int, 
    B_start :: Int, 
    B_end :: Int, 
    K_start :: Int, 
    K_end :: Int,
    BCC :: Float64 = 0.0,
    BSC :: Float64 = 0.0;
    time_limit :: Float64 = 7200.0,
    optim=DEFAULT_OPTIM
)
    dists = get_distance_matrix(g)
    
    # Create a unique timestamped directory for this batch of runs
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    outdir = "experiment_results_$timestamp"
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
                    
                    push!(run_log, (P=P, B=B, K=K, C=C, status=:success, file=filename))
                    println("Saved to $filename")
                catch e
                    if e isa InfeasibleError
                        println("Infeasible P = $P, B = $B, K = $K")
                        push!(run_log, (P=P, B=B, K=K, status=:infeasible, error=string(e)))
                    elseif e isa TimeLimitError
                        println("Time Limit Exceeded P = $P, B = $B, K = $K")
                        push!(run_log, (P=P, B=B, K=K, status=:time_limit_exceeded, error=string(e)))
                    else
                        println("Other errors occurred")
                        push!(run_log, (P=P, B=B, K=K, status=:error, error=string(e)))
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
    
    println("Experiment complete!")
    return run_log
end