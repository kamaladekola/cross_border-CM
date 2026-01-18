function build_consumer_agent!(mod::Model, m::String, zones::Vector{String})
    z, _ = parse_agent_name(m)

    # Extract sets
    JH = mod.ext[:sets][:JH]
    JZ = mod.ext[:sets][:JZ]
    JN = mod.ext[:sets][:JN]

    # Extract time series & parameters
    D = mod.ext[:timeseries][:D] 
    W = mod.ext[:parameters][:w]

    WTP = mod.ext[:parameters][:WTP]                  # value of lost load
    ela = mod.ext[:parameters][:ela]                  # fraction of demand that is elastic
    CD = mod.ext[:parameters][:CD]                    # if we want to use predefined Capacity demand
    CD_margin = mod.ext[:parameters][:CD_margin]      # capacity demand margin
    σ_CM = mod.ext[:parameters][:σ_CM]                # 1 if capacity markets are active, 0 otherwise
    # D_max = mod.ext[:parameters][:D_max]            # maximum demand (can be used in place of CD)
    WTP_CM = mod.ext[:parameters][:WTP_CM]            # Willingness to pay for capacity in the CM (price target)

    # ADMM parameters
    λ_EOM = mod.ext[:parameters][:λ_EOM]            # EOM prices
    g_bar = mod.ext[:parameters][:g_bar]            # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM]            # rho-value in ADMM related to EOM auctions
    λ_CM = mod.ext[:parameters][:λ_CM]              # CM prices
    ρ_CM = mod.ext[:parameters][:ρ_CM]              # rho-value in ADMM related to capacity markets
    cap_bar = mod.ext[:parameters][:cap_bar]        # element in ADMM penalty term related to capacity markets

    PM = mod.ext[:parameters][:participation_matrix]
    node_share = mod.ext[:parameters][:node_share_vec]
    nodes  = mod.ext[:parameters][:nodes]

    # Create variables
    g = mod.ext[:variables][:g] = @variable(mod, [jh=JH], base_name="generation")                                       # consumption as negative generation
    g_VOLL = mod.ext[:variables][:g_VOLL] = @variable(mod, [jh=JH], lower_bound = 0, base_name="inelastic demand")      # inelastic demand
    g_ela = mod.ext[:variables][:g_ela] = @variable(mod, [jh=JH], lower_bound = 0, base_name="elastic demand")          # elastic demand
    ens = mod.ext[:variables][:ens] = @variable(mod, [jh=JH], lower_bound = 0, base_name="unserved_energy")             # unserved energy
    cap_cm = mod.ext[:variables][:cap_cm] = @variable(mod, [jz=JZ], lower_bound = 0, base_name = "capacity offered")    # capacity offered in capacity markets

    # Create expressions
    g_positive = mod.ext[:expressions][:g_positive] = @expression(mod, [jh=JH], -g[jh])  # consumption as positive value
    neg_utility = mod.ext[:expressions][:utility] = @expression(mod,
    sum(W[jh] * ((λ_EOM[jh] - WTP) * g_positive[jh] + (WTP/(2 * ela * D[jh])) * (g_ela[jh])^2) for jh in JH)                    # Utility function for energy consumption
    # + σ_CM * sum((λ_CM[jz]) * cap_cm[jz] for jz in JZ)                                                    # payment for reliability through capacity markets
    # + σ_CM * sum((λ_CM[jz] - WTP_CM) * cap_cm[jz] for jz in JZ)                                                                 # cost of unserved energy
    )

    # Objective => minimize negative utility (maximize utility)
    mod.ext[:objective] = @objective(mod, Min,
    neg_utility 
    + sum(W[jh] * ρ_EOM/2 * (g[jh] - g_bar[jh])^2 for jh in JH)
    # + σ_CM * sum(ρ_CM[jz]/2 * (cap_cm[jz] - cap_bar[jz])^2 for jz in JZ)
    )

    # Constraints
    mod.ext[:constraints][:consumption] = @constraint(mod, [jh=JH], g[jh] == -1 * (g_VOLL[jh] + g_ela[jh]))                             # consumption as negative generation
    mod.ext[:constraints][:elastic_demand] = @constraint(mod, [jh=JH], g_ela[jh] <= ela * D[jh])                                        # Elastic demand limit
    mod.ext[:constraints][:inelastic_demand] = @constraint(mod, [jh=JH], g_VOLL[jh] + ens[jh] == (1 - ela) * D[jh])                     # Inelastic demand limit
    
    d_nodal = mod.ext[:expressions][:d_nodal] = @expression(mod, [jh in JH, jn in JN], node_share[jn] * g_positive[jh])

    for jz in JZ
        if zones[jz] == z
            # mod.ext[:constraints][Symbol("CD_$jz")] = @constraint(mod, cap_cm[jz] >= 0.0)
            mod.ext[:constraints][Symbol("CD_$jz")] = @constraint(mod, cap_cm[jz] >= σ_CM * PM[m][zones[jz]] * CD)
            mod.ext[:constraints][Symbol("CD_upper_$jz")] = @constraint(mod, cap_cm[jz] <= σ_CM * PM[m][zones[jz]] * (1 + CD_margin) * CD)
        else
            mod.ext[:constraints][Symbol("CD_$jz")] = @constraint(mod, cap_cm[jz] >= 0.0)
            mod.ext[:constraints][Symbol("CD_upper_$jz")] = @constraint(mod, cap_cm[jz] <= 0.0)
        end
    end

    return mod
end

