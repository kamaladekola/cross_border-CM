function solve_generator_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]

    # Extract parameters
    A = mod.ext[:parameters][:A] 
    B = mod.ext[:parameters][:B]  
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EOM auctions

    # Investment parameters
    # I = mod.ext[:parameters][:I]

    # Create variables
    g = mod.ext[:variables][:g]  
    # cp = mod.ext[:variables][:cp]

    # Objective 
    mod.ext[:objective] = @objective(mod, Min,
        + sum(A/2*g[jh]^2 for jh in JH) # cost function for generation
        + sum(B*g[jh] for jh in JH)
        - sum(λ_EOM[jh]*g[jh] for jh in JH) # revenue from EOM
        # + I * cp # investment cost
        - 0 # include capacity market revenue here (placeholder)
        + sum(ρ_EOM/2*(g[jh] - g_bar[jh])^2 for jh in JH) # ADMM penalty term
        
    )


    optimize!(mod);

    return mod
end

