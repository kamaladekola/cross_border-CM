function build_consumer_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]

    # Extract time series data
    D = mod.ext[:timeseries][:D] 
    PV = mod.ext[:timeseries][:PV] 

    # Extract parameters
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EOM auctions

    # Create variables
    g = mod.ext[:variables][:g] = @variable(mod, [jh=JH], base_name="generation")  # positive if consumer feeds power to the rest of the system 
    #peak = mod.ext[:variables][:peak] = @variable(mod, base_name="peak_demand", lower_bound = 0 )  # peak demand

    # Create affine expressions 

    # Objective 
    mod.ext[:objective] = @objective(mod, Min,
        - sum(λ_EOM[jh]*g[jh] for jh in JH)
        + sum(ρ_EOM/2*(g[jh] - g_bar[jh])^2 for jh in JH)
        # + price_peak * peak 
    )

    # Demand balance
    mod.ext[:constraints][:energy_balance] = @constraint(mod, [jh=JH],
        g[jh] <= PV[jh] - D[jh] # - charge +discharge
    )

    # Battery model 
    # 

    # find peak 
    # peak greater than g[jh]
    # peak greather than 0

    return mod
end