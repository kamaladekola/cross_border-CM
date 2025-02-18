function solve_consumer_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    # Extract time series data
    D = mod.ext[:timeseries][:D] 

    # Extract parameters
    WTP = mod.ext[:parameters][:WTP]                  # value of lost load
    ela = mod.ext[:parameters][:ela]                  # fraction of demand that is elastic
    CD = mod.ext[:parameters][:CD]                    # Capacity demand (administratively set?)
    # D_max = mod.ext[:parameters][:D_max]            # maximum demand (can be used in place of CD)
    # WTP_CM = mod.ext[:parameters][:WTP_CM]          # Willingness to pay for capacity in the CM (price target)
    # CD_margin = mod.ext[:parameters][:CD_margin]    # capacity demand margin
    σ_CM = mod.ext[:parameters][:σ_CM]                # 1 if capacity markets are active, 0 otherwise

    # ADMM parameters
    λ_EOM = mod.ext[:parameters][:λ_EOM]            # EOM prices
    g_bar = mod.ext[:parameters][:g_bar]            # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM]            # rho-value in ADMM related to EOM auctions
    λ_CM = mod.ext[:parameters][:λ_CM]              # CM prices
    # ρ_CM = mod.ext[:parameters][:ρ_CM]              # rho-value in ADMM related to capacity markets
    # cap_bar = mod.ext[:parameters][:cap_bar]        # element in ADMM penalty term related to capacity markets




    # Create variables
    g = mod.ext[:variables][:g]
    g_VOLL = mod.ext[:variables][:g_VOLL]
    g_ela = mod.ext[:variables][:g_ela]
    # ens = mod.ext[:variables][:ens]
    cap_cm = mod.ext[:variables][:cap_cm]                                                                              # negative capacity offered in capacity markets

    # Create affine expressions
    g_positive = mod.ext[:expressions][:g_positive] = @expression(mod, [jh=JH], g_VOLL[jh] + g_ela[jh])

    utility = @expression(mod,                                                                                            # actually, negative utility
    sum((λ_EOM[jh] - WTP)*g_positive[jh] + (WTP/(2*ela*D[jh]))*(g_ela[jh])^2 for jh in JH)
    - σ_CM * λ_CM * cap_cm
    # - sum(WTP * ens[jh] for jh in JH)                                                                                     # cost of unserved energy (CUE)
    )

    # Objective => minimize negative utility (maximize utility)
    mod.ext[:objective] = @objective(mod, Min,
    utility 
    + sum(ρ_EOM/2 * (g[jh] - g_bar[jh])^2 for jh in JH)
    # + ρ_CM / 2 * (cap_cm - cap_bar)^2                                                                                  # may be used later
    )

    optimize!(mod);

    return mod
end

