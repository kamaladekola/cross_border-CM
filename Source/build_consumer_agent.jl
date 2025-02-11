function build_consumer_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]

    # Extract time series data
    D = mod.ext[:timeseries][:D] 

    # Extract parameters
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EOM auctions
    WTP = mod.ext[:parameters][:WTP] # value of lost load
    ela = mod.ext[:parameters][:ela]   # fraction of demand that is elastic

    # Create variables
    g = mod.ext[:variables][:g] = @variable(mod, [jh=JH], base_name="generation")  # consumption as negative generation
    g_VOLL = mod.ext[:variables][:g_VOLL] = @variable(mod, [jh=JH], lower_bound = 0, base_name="inelastic_demand")
    g_ela = mod.ext[:variables][:g_ela] = @variable(mod, [jh=JH], lower_bound = 0, base_name="elastic_demand")
    ens = mod.ext[:variables][:ens] = @variable(mod, [jh=JH], lower_bound = 0, base_name="unserved_energy") # unserved energy
    # peak = mod.ext[:variables][:peak] = @variable(mod, base_name="peak_demand", lower_bound = 0 )  # peak demand

    # Create affine expressions
    g_positive = mod.ext[:expressions][:g_positive] = @expression(mod, [jh=JH], g_VOLL[jh] + g_ela[jh])

    # Constraints
    mod.ext[:constraints][:consumption] = @constraint(mod, [jh=JH], g[jh] == - g_positive[jh])  # consumption as negative generation
    mod.ext[:constraints][:elastic_demand] = @constraint(mod, [jh=JH], g_ela[jh] <= ela * D[jh]) # Elastic demand limit
    mod.ext[:constraints][:inelastic_demand] = @constraint(mod, [jh=JH], g_VOLL[jh] + ens[jh] == (1 - ela) * D[jh]) # Inelastic demand limit

    profit = @expression(mod,                                                                       # actually, negative profit
    sum((λ_EOM[jh] - WTP)*g_positive[jh] + (WTP/(2*ela*D[jh]))*(g_ela[jh])^2 for jh in JH)
    )

    # Objective - minimize negative profit (maximize profit)
    mod.ext[:objective] = @objective(mod, Min,
    profit + sum(ρ_EOM/2 * (g[jh] - g_bar[jh])^2 for jh in JH))

    # Battery model 
    # 

    # find peak 
    # peak greater than g[jh]
    # peak greater than 0

    return mod
end


    # Objective - without demand elasticity
    # mod.ext[:objective] = @objective(mod, Min,
    #     sum((λ_EOM[jh] - WTP)*g[jh] for jh in JH)
    #     + sum(ρ_EOM/2*(g[jh] - g_bar[jh])^2 for jh in JH)
    #     # + price_peak * peak 
    # )