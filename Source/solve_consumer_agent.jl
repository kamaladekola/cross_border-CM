function solve_consumer_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    # Extract time series data
    D = mod.ext[:timeseries][:D] 

    # Extract parameters
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EOM auctions
    WTP = mod.ext[:parameters][:WTP] # Willingness to pay
    ela = mod.ext[:parameters][:ela]   # fraction of demand that is elastic

    # Create variables
    g = mod.ext[:variables][:g]
    g_VOLL = mod.ext[:variables][:g_VOLL]
    g_ela = mod.ext[:variables][:g_ela]

    # Create affine expressions
    g_positive = mod.ext[:expressions][:g_positive] = @expression(mod, [jh=JH], g_VOLL[jh] + g_ela[jh])

    profit = @expression(mod, 
    sum((λ_EOM[jh] - WTP)*g_positive[jh] + (WTP/(2*ela*D[jh]))*(g_ela[jh])^2 for jh in JH)
    )
    # Objective - minimize negative profit (maximize profit)
    mod.ext[:objective] = @objective(mod, Min,
    profit + sum(ρ_EOM/2 * (g[jh] - g_bar[jh])^2 for jh in JH))

    optimize!(mod);

    return mod
end

