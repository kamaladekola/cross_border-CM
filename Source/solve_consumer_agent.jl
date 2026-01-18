function solve_consumer_agent!(mod::Model, m::String, zones::Vector{String})
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JZ = mod.ext[:sets][:JZ]
    # Extract time series data
    D = mod.ext[:timeseries][:D]
    W = mod.ext[:parameters][:w] 

    # Extract parameters
    WTP = mod.ext[:parameters][:WTP]                  # value of lost load
    ela = mod.ext[:parameters][:ela]                  # fraction of demand that is elastic
    # CD = mod.ext[:parameters][:CD]                    # Capacity demand (administratively set?)
    # D_max = mod.ext[:parameters][:D_max]            # maximum demand (can be used in place of CD)
    WTP_CM = mod.ext[:parameters][:WTP_CM]          # Willingness to pay for capacity in the CM (price target)
    CD_margin = mod.ext[:parameters][:CD_margin]    # capacity demand margin
    σ_CM = mod.ext[:parameters][:σ_CM]                # 1 if capacity markets are active, 0 otherwise
    PM = mod.ext[:parameters][:participation_matrix]

    # ADMM parameters
    λ_EOM = mod.ext[:parameters][:λ_EOM]            # EOM prices
    g_bar = mod.ext[:parameters][:g_bar]            # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM]            # rho-value in ADMM related to EOM auctions
    λ_CM = mod.ext[:parameters][:λ_CM]              # CM prices
    ρ_CM = mod.ext[:parameters][:ρ_CM]              # rho-value in ADMM related to capacity markets
    cap_bar = mod.ext[:parameters][:cap_bar]        # element in ADMM penalty term related to capacity markets


    # Create variables
    g = mod.ext[:variables][:g]
    g_ela = mod.ext[:variables][:g_ela]
    cap_cm = mod.ext[:variables][:cap_cm]
    ens = mod.ext[:variables][:ens]                                                                              # negative capacity offered in capacity markets

    # Create affine expressions
    g_positive = mod.ext[:expressions][:g_positive] = @expression(mod, [jh=JH], -g[jh])  # consumption as positive value

    neg_utility = mod.ext[:expressions][:utility] = @expression(mod,                                                                                           
    sum(W[jh] * ((λ_EOM[jh] - WTP)*g_positive[jh] + (WTP/(2*ela*D[jh]))*(g_ela[jh])^2) for jh in JH)
    # + σ_CM * sum(λ_CM[jz] * cap_cm[jz] * PM[m][zones[jz]] for jz in JZ)
    # + σ_CM * sum((λ_CM[jz]) * cap_cm[jz] for jz in JZ)
    # + σ_CM * sum((λ_CM[jz] - WTP_CM) * cap_cm[jz] for jz in JZ)
    )

    # Objective => minimize negative utility (maximize utility)
    mod.ext[:objective] = @objective(mod, Min,
    neg_utility 
    + sum(W[jh] * ρ_EOM/2 * (g[jh] - g_bar[jh])^2 for jh in JH)
    # + σ_CM * sum(ρ_CM[jz]/2 * (cap_cm[jz] - cap_bar[jz])^2 for jz in JZ)
    )
    optimize!(mod);

    return mod
end

