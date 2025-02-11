# ADMM 
function ADMM!(results::Dict,ADMM::Dict,EOM::Dict,mdict::Dict,agents::Dict,scenario_overview_row::DataFrameRow,data::Dict,TO::TimerOutput, zones::Vector{String})
    convergence = 0
    iterations = ProgressBar(1:data["ADMM"]["max_iter"])

    for iter in iterations
        if convergence == 0
            # Multi-threaded version
            @sync for m in agents[:all] 
                # created subroutine to allow multi-treading to solve agents' decision problems
                @spawn ADMM_subroutine!(m,results,ADMM,EOM,mdict[m],agents,TO)
            end

            # Imbalances (for each zone)
            @timeit TO "Compute zonal imbalances" begin
                for z in zones
                    # push!(ADMM["Imbalances"]["EOM"][z], sum(results["g"][m][end] for m in agents[:zones][z]) - EOM["D"][z][:])
                    push!(ADMM["Imbalances"]["EOM"][z], sum(is_interconnector(m) ? sign_for_zone(m, z) * results["g"][m][end] : results["g"][m][end] for m in agents[:zones][z]))
                end                                
            end

            # Primal Residuals for each zone
            @timeit TO "Compute primal residuals" begin
                for z in zones
                    push!(ADMM["Residuals"]["Primal"]["EOM"][z], sqrt(sum(ADMM["Imbalances"]["EOM"][z][end].^2)))
                end
            end

            # Compute Dual Residuals # nAgents/n_zones?
            @timeit TO "Compute dual residuals" begin
                if iter > 1
                    for z in zones
                    push!(ADMM["Residuals"]["Dual"]["EOM"][z], sqrt(sum(sum((ADMM["ρ"]["EOM"][z][end]*((results["g"][m][end]-sum(results["g"][mstar][end] for mstar in agents[:eom])./(EOM["nAgents"]+1)) - (results["g"][m][end-1]-sum(results["g"][mstar][end-1] for mstar in agents[:eom])./(EOM["nAgents"]+1)))).^2 for m in agents[:eom]))))
                    end
                end
            end

            # Update prices for each zone
            @timeit TO "Update prices" begin
                for z in zones
                    push!(results["λ"]["EOM"][z], results["λ"]["EOM"][z][end] - ADMM["ρ"]["EOM"][z][end]/100*ADMM["Imbalances"]["EOM"][z][end])
                end
            end

            # Update ρ-values
            @timeit TO "Update ρ" begin
                 update_rho!(ADMM,iter)
            end

            # Progress bar for multi-zone
            @timeit TO "Progress bar" begin
                max_primal = maximum([ADMM["Residuals"]["Primal"]["EOM"][z][end] for z in zones])
                max_dual   = maximum([ADMM["Residuals"]["Dual"]["EOM"][z][end] for z in zones])
                set_description(iterations, @sprintf("Max Primal: %.3f, Max Dual: %.3f", max_primal, max_dual))
            end

            # Check convergence: primal and dual satisfy tolerance for each zone
            if all(ADMM["Residuals"]["Primal"]["EOM"][z][end] <= ADMM["Tolerance"]["EOM"] for z in zones) && all(ADMM["Residuals"]["Dual"]["EOM"][z][end] <= ADMM["Tolerance"]["EOM"] for z in zones)
                convergence = 1
            end

            # store number of iterations
            ADMM["n_iter"] = copy(iter)
        end
    end
end