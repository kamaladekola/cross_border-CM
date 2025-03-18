function solve_generator_agent!(mod::Model, m::String, zones::Vector{String})
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JZ = mod.ext[:sets][:JZ]
    # Extract parameters
    A = mod.ext[:parameters][:A] 
    B = mod.ext[:parameters][:B]  
    λ_EOM = mod.ext[:parameters][:λ_EOM]        # EOM prices
    g_bar = mod.ext[:parameters][:g_bar]        # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM]        # rho-value in ADMM related to EOM auctions
    cap_bar = mod.ext[:parameters][:cap_bar]    # element in ADMM penalty term related to capacity markets
    ρ_CM = mod.ext[:parameters][:ρ_CM]          # rho-value in ADMM related to capacity markets
    λ_CM = mod.ext[:parameters][:λ_CM]          # CM prices
    I = mod.ext[:parameters][:I]                # investment cost
    y_init   = mod.ext[:parameters][:C]         # existing capacity
    σ_CM = mod.ext[:parameters][:σ_CM]          # 1 if capacity markets are active, 0 otherwise

    PM = mod.ext[:parameters][:participation_matrix]

    # Create variables
    g = mod.ext[:variables][:g]  
    y = mod.ext[:variables][:y]
    cap_cm = mod.ext[:variables][:cap_cm]


    # Objective => minimize GenCo costs
    mod.ext[:objective] = @objective(mod, Min,
        # + sum(A/2*g[jh]^2 for jh in JH)                       # cost function for generation
        + sum(B*g[jh] for jh in JH)
        - sum(λ_EOM[jh]*g[jh] for jh in JH)                     # revenue from EOM
        + I * (y - y_init)                                      # investment cost
        - σ_CM * sum(λ_CM[jz] * cap_cm[jz] * PM[m][zones[jz]] for jz in JZ)        # revenue from capacity markets
        + sum(ρ_EOM/2*(g[jh] - g_bar[jh])^2 for jh in JH)       # ADMM penalty term for EOM clearing
        + σ_CM * sum(ρ_CM[jz]/2 * PM[m][zones[jz]] * (cap_cm[jz] - cap_bar[jz])^2 for jz in JZ) # ADMM penalty term for capacity markets
    )


    optimize!(mod);

    return mod
end

