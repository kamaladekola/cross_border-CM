function update_rho!(ADMM::Dict, iter::Int64)
    if mod(iter,1) == 0
        # ρ-updates following Boyd et al. (2011)
        if ADMM["Residuals"]["Primal"]["EOM"][end] > 2*ADMM["Residuals"]["Dual"]["EOM"][end]
            push!(ADMM["ρ"]["EOM"], minimum([1000,1.1*ADMM["ρ"]["EOM"][end]]))
        elseif ADMM["Residuals"]["Dual"]["EOM"][end] > 2*ADMM["Residuals"]["Primal"]["EOM"][end]
            push!(ADMM["ρ"]["EOM"], 1/1.1*ADMM["ρ"]["EOM"][end])
        end
    end
end