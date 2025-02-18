function update_rho!(ADMM::Dict, iter::Int64)
    # Update ρ for each zone
    if mod(iter, 1) == 0
        for zone in keys(ADMM["Residuals"]["Primal"]["EOM"])
            # ρ-updates following Boyd et al. (2011)
            if last(ADMM["Residuals"]["Primal"]["EOM"][zone]) > 2 * last(ADMM["Residuals"]["Dual"]["EOM"][zone])
                push!(ADMM["ρ"]["EOM"][zone], minimum([1000, 1.1 * last(ADMM["ρ"]["EOM"][zone])]))
            elseif last(ADMM["Residuals"]["Dual"]["EOM"][zone]) > 2 * last(ADMM["Residuals"]["Primal"]["EOM"][zone])
                push!(ADMM["ρ"]["EOM"][zone], 1/1.1 * last(ADMM["ρ"]["EOM"][zone]))
            end

            if last(ADMM["Residuals"]["Primal"]["CM"][zone]) > 2 * last(ADMM["Residuals"]["Dual"]["CM"][zone])
                push!(ADMM["ρ"]["CM"][zone], minimum([1000, 1.1 * last(ADMM["ρ"]["CM"][zone])]))
            elseif last(ADMM["Residuals"]["Dual"]["CM"][zone]) > 2 * last(ADMM["Residuals"]["Primal"]["CM"][zone])
                push!(ADMM["ρ"]["CM"][zone], 1/1.1 * last(ADMM["ρ"]["CM"][zone]))
            end
        end
    end
end