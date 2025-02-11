function build_generator_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]

    # Extract time series data
    AC = mod.ext[:timeseries][:AC]

    # Extract parameters
    A = mod.ext[:parameters][:A] 
    B = mod.ext[:parameters][:B]  
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EOM auctions
    
    # Investment parameters
    # I = mod.ext[:parameters][:I]
    # max_cap  = mod.ext[:parameters][:max_cap]
    # C_init   = mod.ext[:parameters][:C]   # existing capacity

    # Create variables
    g = mod.ext[:variables][:g] = @variable(mod, [jh=JH], lower_bound=0, base_name="generation") # generation
    # cp = mod.ext[:variables][:cp] = @variable(mod, lower_bound=C_init, upper_bound=max_cap, base_name="InstalledCapacity") # installed capacity

    # Objective
    mod.ext[:objective] = @objective(mod, Min,
        + sum(A/2*g[jh]^2 for jh in JH) # cost function for generation
        + sum(B*g[jh] for jh in JH)
        - sum(λ_EOM[jh]*g[jh] for jh in JH) # revenue from EOM
        # + I * cp # investment cost
        - 0 # include capacity market revenue here (placeholder)
        + sum(ρ_EOM/2*(g[jh] - g_bar[jh])^2 for jh in JH) # ADMM penalty term
        
    )

    # capacity constraints
    mod.ext[:constraints][:cap_limit] = @constraint(mod, [jh=JH], 
        g[jh] <=  AC[jh] # capacity limit
    )

    # offered capacity in capacity markets <= installed capacity

    return mod
end