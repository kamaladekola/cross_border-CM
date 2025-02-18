function build_generator_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]

    # Extract time series data
    # AC = mod.ext[:timeseries][:AC] # deprecated
    AF = mod.ext[:timeseries][:AF]

    # Extract parameters
    A = mod.ext[:parameters][:A] 
    B = mod.ext[:parameters][:B]  
    λ_EOM = mod.ext[:parameters][:λ_EOM]                # EOM prices
    g_bar = mod.ext[:parameters][:g_bar]                # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM]                # rho-value in ADMM related to EOM auctions
    cap_bar = mod.ext[:parameters][:cap_bar]            # element in ADMM penalty term related to capacity markets
    ρ_CM = mod.ext[:parameters][:ρ_CM]
    λ_CM = mod.ext[:parameters][:λ_CM]
    σ_CM = mod.ext[:parameters][:σ_CM]                  # 1 if capacity markets are active, 0 otherwise


    
    # Investment parameters
    I = mod.ext[:parameters][:I]
    max_cap  = mod.ext[:parameters][:max_cap]
    y_init   = mod.ext[:parameters][:C]   # existing capacity

    # Create variables
    g = mod.ext[:variables][:g] = @variable(mod, [jh=JH], lower_bound=0, base_name="generation")                                # generation
    y = mod.ext[:variables][:y] = @variable(mod, lower_bound=y_init, upper_bound=max_cap, base_name="InstalledCapacity")       # installed capacity
    cap_cm = mod.ext[:variables][:cap_cm] = @variable(mod, lower_bound = 0, base_name = "capacity offered")                     # capacity offered in capacity markets

    # Objective => minimize GenCo costs
    mod.ext[:objective] = @objective(mod, Min,
        + sum(A/2*g[jh]^2 for jh in JH)                       # cost function for generation
        + sum(B*g[jh] for jh in JH)
        - sum(λ_EOM[jh]*g[jh] for jh in JH)                     # revenue from EOM
        + I * (y - y_init)                                      # investment cost
        - σ_CM * λ_CM * cap_cm                                  # revenue from capacity markets
        + sum(ρ_EOM/2*(g[jh] - g_bar[jh])^2 for jh in JH)       # ADMM penalty term for EOM clearing
        + σ_CM * ρ_CM / 2 * (cap_cm - cap_bar)^2                # ADMM penalty term for capacity markets
    )

    # generation <= installed capacity
    mod.ext[:constraints][:cap_limit] = @constraint(mod, [jh=JH], 
        # g[jh] <=  AC[jh] # capacity limit
        g[jh] <= y * AF[jh] # capacity limit
    )

    # offered capacity in capacity markets
    if σ_CM == 1
        mod.ext[:constraints][:CM] = @constraint(mod, 
            cap_cm <= y
        )
    else
        mod.ext[:constraints][:CM] = @constraint(mod, 
            cap_cm == 0
        )
    end

    return mod
end