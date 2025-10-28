# ADMM 
function ADMM!(results::Dict,ADMM::Dict,EOM::Dict,CM::Dict,mdict::Dict,agents::Dict,scenario_overview_row::DataFrameRow,data::Dict,TO::TimerOutput, zones::Vector{String})
    convergence = 0
    iterations = ProgressBar(1:data["ADMM"]["max_iter"])

    for iter in iterations
        if convergence == 0
            # Multi-threaded version
            @sync for m in agents[:all] 
                # created subroutine to allow multi-treading to solve agents' decision problems
                @spawn ADMM_subroutine!(m,results,ADMM,EOM,CM,mdict[m],agents,TO,zones,data)
            end

            # Imbalances
            @timeit TO "Compute zonal imbalances" begin
                for (zone_idx, z) in enumerate(zones)
                    push!(ADMM["Imbalances"]["EOM"][z],  sum(results["g"][m][end] for m in agents[:eom_Z][z]) + results["g"]["NetworkManager"][end][:, zone_idx])
                    ## Total generation capacity from generators in zone z (offered anywhere) - Capacity demand in zone z + capacity manager net position in zone Z
                    # push!(ADMM["Imbalances"]["CM"][z], sum(sum(results["cap_cm"][m][end]) for m in agents[:cm_Z][z] if m in agents[:Gen]) - sum(results["cap_cm"][m][end][zone_idx] for m in agents[:cm_Z][z] if m in agents[:Cons]) + results["cap_cm"]["CapacityManager"][end][zone_idx])

                   push!(ADMM["Imbalances"]["CM"][z], sum(results["cap_cm"][m][end][zone_idx]  for m in agents[:cm_Z][z] if m in agents[:Gen]) - sum(results["cap_cm"][m][end][zone_idx] for m in agents[:cm_Z][z] if m in agents[:Cons]) + results["cap_cm"]["CapacityManager"][end][zone_idx])                    
                end                                
            end

            # Primal Residuals for each zone
            @timeit TO "Compute primal residuals" begin
                for z in zones
                    push!(ADMM["Residuals"]["Primal"]["EOM"][z], sqrt(sum(ADMM["Imbalances"]["EOM"][z][end].^2)))
                    push!(ADMM["Residuals"]["Primal"]["CM"][z], sqrt(sum(ADMM["Imbalances"]["CM"][z][end].^2)))
                end
            end


            # Compute Dual Residuals for each zone
            @timeit TO "Compute dual residuals" begin                           
                if iter > 1
                    for z in zones
                        zone_idx = findfirst(isequal(z), zones)
                        NM_new = results["g"]["NetworkManager"][end][:, zone_idx]
                        NM_prev = results["g"]["NetworkManager"][end-1][:, zone_idx]
                        EOM_new = sum(results["g"][m][end] for m in agents[:eom_Z][z])
                        EOM_prev = sum(results["g"][m][end-1] for m in agents[:eom_Z][z])
                        
                         # why EOM["nAgents"]+1? 
                        push!(ADMM["Residuals"]["Dual"]["EOM"][z], sqrt(sum(sum((ADMM["ρ"]["EOM"][z][end]*((results["g"][m][end] - (EOM_new + NM_new)./(EOM["nAgents_z"][z]+1)) - (results["g"][m][end-1] - (EOM_prev + NM_prev)./(EOM["nAgents_z"][z]+1)))).^2 for m in agents[:eom_Z][z])) +
                                sum((ADMM["ρ"]["EOM"][z][end]*((NM_new .- (EOM_new + NM_new)./(EOM["nAgents_z"][z]+1)) - (NM_prev .- (EOM_prev + NM_prev)./(EOM["nAgents_z"][z]+1)))).^2)))
                        push!(ADMM["Residuals"]["Dual"]["CM"][z], sqrt(sum(sum((ADMM["ρ"]["CM"][z][end]*((results["cap_cm"][m][end]-sum(results["cap_cm"][mstar][end] for mstar in agents[:cm_Z][z])./(CM["nAgents_z"][z]+1)) - (results["cap_cm"][m][end-1]-sum(results["cap_cm"][mstar][end-1] for mstar in agents[:cm_Z][z])./(CM["nAgents_z"][z]+1)))).^2 for m in agents[:cm_Z][z]))))
                    end
                end
            end

            # Update prices for each zone #  ρ damper - 20 EOM; 50 CM; 
            @timeit TO "Update prices" begin
                for z in zones
                    push!(results["λ"]["EOM"][z], results["λ"]["EOM"][z][end] - ADMM["ρ"]["EOM"][z][end]/5*ADMM["Imbalances"]["EOM"][z][end])
                    push!(results["λ"]["CM"][z], results["λ"]["CM"][z][end] - ADMM["ρ"]["CM"][z][end]/10*ADMM["Imbalances"]["CM"][z][end])
                    # limit price update to investment cost of most expensive plant and no negative prices
                    # λ_CM_new = results["λ"]["CM"][z][end] - ADMM["ρ"]["CM"][z][end]/20*ADMM["Imbalances"]["CM"][z][end]
                    # λ_CM_new = max(λ_CM_new, data["CM"][z]["min_price"])
                    # λ_CM_new = min(λ_CM_new, data["CM"][z]["max_price"])
                    # push!(results["λ"]["CM"][z], λ_CM_new)
                end
            end

            # Update ρ-values
            @timeit TO "Update ρ" begin
                 update_rho!(ADMM,iter)
            end

            # Progress bar
            @timeit TO "Progress bar" begin
                max_primal_eom = maximum([ADMM["Residuals"]["Primal"]["EOM"][z][end] for z in zones])
                max_dual_eom   = maximum([ADMM["Residuals"]["Dual"]["EOM"][z][end] for z in zones])
                
                max_primal_cm = maximum([ADMM["Residuals"]["Primal"]["CM"][z][end] for z in zones])
                max_dual_cm   = maximum([ADMM["Residuals"]["Dual"]["CM"][z][end] for z in zones])
                
                set_description(iterations, @sprintf("Max EOM-Primal: %.3f, EOM-Dual: %.3f, CM-Primal: %.3f, CM-Dual: %.3f", max_primal_eom, max_dual_eom, max_primal_cm, max_dual_cm))
            end
            
            # Check convergence: primal and dual satisfy tolerance for each zone's EOM and CM
            if all(ADMM["Residuals"]["Primal"]["EOM"][z][end] ≤ ADMM["Tolerance"]["EOM"] for z in zones) &&
               all(ADMM["Residuals"]["Dual"]["EOM"][z][end] ≤ ADMM["Tolerance"]["EOM"] for z in zones) &&
               all(ADMM["Residuals"]["Primal"]["CM"][z][end] ≤ ADMM["Tolerance"]["CM"]  for z in zones) &&
               all(ADMM["Residuals"]["Dual"]["CM"][z][end] ≤ ADMM["Tolerance"]["CM"]  for z in zones)
                convergence = 1
            end
            

            # store number of iterations
            ADMM["n_iter"] = copy(iter)
        end
    end
end