function build_consumer_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]

    # Extract time series data
    D = mod.ext[:timeseries][:D] 

    # Extract parameters
    WTP = mod.ext[:parameters][:WTP]                  # value of lost load
    ela = mod.ext[:parameters][:ela]                  # fraction of demand that is elastic
    CD = mod.ext[:parameters][:CD]                    # if we want to use predefined Capacity demand
    CD_margin = mod.ext[:parameters][:CD_margin]      # capacity demand margin
    σ_CM = mod.ext[:parameters][:σ_CM]                # 1 if capacity markets are active, 0 otherwise
    D_max = mod.ext[:parameters][:D_max]            # maximum demand (can be used in place of CD)
    # WTP_CM = mod.ext[:parameters][:WTP_CM]          # Willingness to pay for capacity in the CM (price target)

    # ADMM parameters
    λ_EOM = mod.ext[:parameters][:λ_EOM]            # EOM prices
    g_bar = mod.ext[:parameters][:g_bar]            # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM]            # rho-value in ADMM related to EOM auctions
    λ_CM = mod.ext[:parameters][:λ_CM]              # CM prices
    ρ_CM = mod.ext[:parameters][:ρ_CM]              # rho-value in ADMM related to capacity markets
    cap_bar = mod.ext[:parameters][:cap_bar]        # element in ADMM penalty term related to capacity markets



    # Create variables
    g = mod.ext[:variables][:g] = @variable(mod, [jh=JH], base_name="generation")                                       # consumption as negative generation
    g_VOLL = mod.ext[:variables][:g_VOLL] = @variable(mod, [jh=JH], lower_bound = 0, base_name="inelastic demand")      # inelastic demand
    g_ela = mod.ext[:variables][:g_ela] = @variable(mod, [jh=JH], lower_bound = 0, base_name="elastic demand")          # elastic demand
    ens = mod.ext[:variables][:ens] = @variable(mod, [jh=JH], lower_bound = 0, base_name="unserved_energy")             # unserved energy
    cap_cm = mod.ext[:variables][:cap_cm] = @variable(mod, upper_bound = 0, base_name = "capacity offered")             # negative capacity offered in capacity markets
    # peak = mod.ext[:variables][:peak] = @variable(mod, base_name="peak_demand", lower_bound = 0 )                       # peak demand

    # Create affine expressions
    g_positive = mod.ext[:expressions][:g_positive] = @expression(mod, [jh=JH], g_VOLL[jh] + g_ela[jh])

    neg_utility =  mod.ext[:expressions][:utility] = @expression(mod,                                                                                            # actually, negative utility
    sum((λ_EOM[jh] - WTP)*g_positive[jh] + (WTP/(2*ela*D[jh]))*(g_ela[jh])^2 for jh in JH)
    + (- σ_CM * λ_CM * cap_cm)
    # + sum(WTP * ens[jh] for jh in JH)                                                                                       # penalty on energy not served (called cost of unserved energy) in Kaminski (PhD thesis)
    )

    # Objective => minimize negative utility (maximize utility)
    mod.ext[:objective] = @objective(mod, Min,
    neg_utility 
    + sum(ρ_EOM/2 * (g[jh] - g_bar[jh])^2 for jh in JH)
    + σ_CM * ρ_CM / 2 * (cap_cm - cap_bar)^2                # ADMM penalty term for capacity markets
    )

    # Constraints
    mod.ext[:constraints][:consumption] = @constraint(mod, [jh=JH], g[jh] == - g_positive[jh])                           # consumption as negative generation
    mod.ext[:constraints][:elastic_demand] = @constraint(mod, [jh=JH], g_ela[jh] <= ela * D[jh])                         # Elastic demand limit
    mod.ext[:constraints][:inelastic_demand] = @constraint(mod, [jh=JH], g_VOLL[jh] + ens[jh] == (1 - ela) * D[jh])      # Inelastic demand limit

    mod.ext[:constraints][:CM] = @constraint(mod, cap_cm <= - σ_CM * (1 + CD_margin) * D_max)

    # Battery / electrolyzer model 

    # find peak 
    # peak greater than g[jh]
    # peak greater than 0

    return mod
end

