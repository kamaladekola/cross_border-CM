# ADMM 
function ADMM!(results::Dict,ADMM::Dict,EOM::Dict,mdict::Dict,agents::Dict,scenario_overview_row::DataFrameRow,data::Dict,TO::TimerOutput)
    convergence = 0
    iterations = ProgressBar(1:data["ADMM"]["max_iter"])
    for iter in iterations
        if convergence == 0
            # Multi-threaded version
            @sync for m in agents[:all] 
                # created subroutine to allow multi-treading to solve agents' decision problems
                @spawn ADMM_subroutine!(m,results,ADMM,EOM,mdict[m],agents,TO)
            end

            # Imbalances 
            @timeit TO "Compute imbalances" begin
                push!(ADMM["Imbalances"]["EOM"], sum(results["g"][m][end] for m in agents[:eom]) - EOM["D"][:])
            end

            # Primal residuals 
            @timeit TO "Compute primal residuals" begin
                push!(ADMM["Residuals"]["Primal"]["EOM"], sqrt(sum(ADMM["Imbalances"]["EOM"][end].^2)))
            end

            # Dual residuals
            @timeit TO "Compute dual residuals" begin 
            if iter > 1
                push!(ADMM["Residuals"]["Dual"]["EOM"], sqrt(sum(sum((ADMM["ρ"]["EOM"][end]*((results["g"][m][end]-sum(results["g"][mstar][end] for mstar in agents[:eom])./(EOM["nAgents"]+1)) - (results["g"][m][end-1]-sum(results["g"][mstar][end-1] for mstar in agents[:eom])./(EOM["nAgents"]+1)))).^2 for m in agents[:eom]))))               
            end
            end

            # Price updates 
            @timeit TO "Update prices" begin
                push!(results[ "λ"]["EOM"], results[ "λ"]["EOM"][end] - ADMM["ρ"]["EOM"][end]/100*ADMM["Imbalances"]["EOM"][end])
            end

            # Update ρ-values
            @timeit TO "Update ρ" begin
                 update_rho!(ADMM,iter)
            end

            # Progress bar
            @timeit TO "Progress bar" begin
                set_description(iterations, string(@sprintf("Primal residual - EOM: %.3f -- Dual residual - EOM: %.3f",ADMM["Residuals"]["Primal"]["EOM"][end],ADMM["Residuals"]["Dual"]["EOM"][end])))
            end

            # Check convergence: primal and dual satisfy tolerance 
            if ADMM["Residuals"]["Primal"]["EOM"][end] <= ADMM["Tolerance"]["EOM"] && ADMM["Residuals"]["Dual"]["EOM"][end] <= ADMM["Tolerance"]["EOM"] 
                convergence = 1
            end

            # store number of iterations
            ADMM["n_iter"] = copy(iter)
        end
    end
end