function solve_interconnector_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]

    # Extract parameters
    λ1 = mod.ext[:parameters][:λ1]
    λ2 = mod.ext[:parameters][:λ2]

    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EOM auctions

    # Create variables
    g = mod.ext[:variables][:g]

    # Objective
    mod.ext[:objective] = @objective(mod, Min,
    - sum((λ2[jh] - λ1[jh]) * g[jh] for jh in JH)
    + sum(ρ_EOM/2 * (g[jh] - g_bar[jh])^2 for jh in JH)
)

    optimize!(mod);

    return mod
end
