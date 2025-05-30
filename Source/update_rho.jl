# function update_rho!(ADMM::Dict, iter::Int64)
#     # Update ρ for each zone
#     if mod(iter, 1) == 0
#         for zone in keys(ADMM["Residuals"]["Primal"]["EOM"])
#             # ρ-updates following Boyd et al. (2011)
#             if last(ADMM["Residuals"]["Primal"]["EOM"][zone]) > 2 * last(ADMM["Residuals"]["Dual"]["EOM"][zone])
#                 push!(ADMM["ρ"]["EOM"][zone], minimum([1000, 1.1 * last(ADMM["ρ"]["EOM"][zone])]))
#             elseif last(ADMM["Residuals"]["Dual"]["EOM"][zone]) > 2 * last(ADMM["Residuals"]["Primal"]["EOM"][zone])
#                 push!(ADMM["ρ"]["EOM"][zone], 1/1.1 * last(ADMM["ρ"]["EOM"][zone]))
#             end

#             if last(ADMM["Residuals"]["Primal"]["CM"][zone]) > 2 * last(ADMM["Residuals"]["Dual"]["CM"][zone])
#                 push!(ADMM["ρ"]["CM"][zone], minimum([1000, 1.1 * last(ADMM["ρ"]["CM"][zone])]))
#             elseif last(ADMM["Residuals"]["Dual"]["CM"][zone]) > 2 * last(ADMM["Residuals"]["Primal"]["CM"][zone])
#                 push!(ADMM["ρ"]["CM"][zone], 1/1.1 * last(ADMM["ρ"]["CM"][zone]))
#             end
#         end
#     end
# end


function update_rho!(ADMM::Dict, iter::Int;
                     μ = 10.0,              # threshold ratio
                     τ = 1.1,               # scale factor
                     iter_skip = 1,         # change every k iterations
                     ρmax = 1_000.0)

    # only adjust every `iter_skip` iterations
    if iter_skip > 1 && (iter % iter_skip != 0)
        # still append unchanged ρ so vectors stay aligned
        for z in keys(ADMM["ρ"]["EOM"])
            push!(ADMM["ρ"]["EOM"][z],  last(ADMM["ρ"]["EOM"][z]))
            push!(ADMM["ρ"]["CM"][z],   last(ADMM["ρ"]["CM"][z]))
        end
        return
    end

    for z in keys(ADMM["Residuals"]["Primal"]["EOM"])

        # ----- EOM -----
        r_p = last(ADMM["Residuals"]["Primal"]["EOM"][z])
        r_d = last(ADMM["Residuals"]["Dual"]["EOM"][z])
        ρ_now = last(ADMM["ρ"]["EOM"][z])

        if r_p > μ * r_d
            ρ_new = min(ρmax, τ * ρ_now)
        elseif r_d > μ * r_p
            ρ_new = ρ_now / τ
        else
            ρ_new = ρ_now 
        end
        push!(ADMM["ρ"]["EOM"][z], ρ_new)

        # ----- CM -----
        r_p = last(ADMM["Residuals"]["Primal"]["CM"][z])
        r_d = last(ADMM["Residuals"]["Dual"]["CM"][z])
        ρ_now = last(ADMM["ρ"]["CM"][z])

        if r_p > μ * r_d
            ρ_new = min(ρmax, τ * ρ_now)
        elseif r_d > μ * r_p
            ρ_new = ρ_now / τ
        else
            ρ_new = ρ_now
        end
        push!(ADMM["ρ"]["CM"][z], ρ_new)
    end
end